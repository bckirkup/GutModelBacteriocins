/* -----------------------------------------------------------------------
   GutIBM – Chemical field implementation
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "domain.h"
#include "error.h"
#include "species_names.h"
#include "tridiagonal_factorization.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {

namespace {

void solve_tridiagonal_with_diagonal(
    std::vector<Real>& values, const std::vector<Real>& diagonal,
    Real alpha, Real first_source = 0.0) {
  if (values.empty()) return;
  std::vector<Real> diag = diagonal;
  std::vector upper(diag.size() > 1 ? diag.size() - 1 : 0, -alpha);
  values.front() += first_source;
  for (size_t i = 1; i < diag.size(); ++i) {
    const Real multiplier = -alpha / diag[i - 1];
    diag[i] -= multiplier * upper[i - 1];
    values[i] -= multiplier * values[i - 1];
  }
  values.back() /= diag.back();
  for (size_t i = values.size() - 1; i > 0; --i) {
    values[i - 1] =
        (values[i - 1] - upper[i - 1] * values[i]) / diag[i - 1];
  }
}

void solve_periodic_with_sink(
    std::vector<Real>& values, const std::vector<Real>& sink, Real alpha) {
  const size_t n = values.size();
  if (n == 0) return;
  if (n == 1) {
    values[0] /= 1.0 + sink[0];
    return;
  }
  if (n == 2) {
    const Real d0 = 1.0 + 2.0 * alpha + sink[0];
    const Real d1 = 1.0 + 2.0 * alpha + sink[1];
    const Real off = 2.0 * alpha;
    const Real determinant = d0 * d1 - off * off;
    const Real first = (d1 * values[0] + off * values[1]) / determinant;
    const Real second = (off * values[0] + d0 * values[1]) / determinant;
    values[0] = first;
    values[1] = second;
    return;
  }
  const Real corner = -alpha;
  const Real gamma = -(1.0 + 2.0 * alpha + sink[0]);
  std::vector diagonal(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    diagonal[i] = 1.0 + 2.0 * alpha + sink[i];
  }
  diagonal[0] -= gamma;
  diagonal[n - 1] -= corner * corner / gamma;
  TridiagonalFactorization factorization;
  std::vector lower(n - 1, -alpha);
  std::vector upper(n - 1, -alpha);
  factorization.factorize(lower, diagonal, upper);
  std::vector correction(n, 0.0);
  correction[0] = gamma;
  correction[n - 1] = corner;
  factorization.solve_in_place(correction);
  const Real denominator = 1.0 + correction[0]
      + corner * correction[n - 1] / gamma;
  factorization.solve_in_place(values);
  const Real adjustment =
      (values[0] + corner * values[n - 1] / gamma) / denominator;
  for (size_t i = 0; i < n; ++i) values[i] -= adjustment * correction[i];
}

class PeriodicLineSolver {
 public:
  PeriodicLineSolver(Int size, Real alpha) : size_(size), alpha_(alpha) {
    if (size_ < 3) return;

    const Real diagonal_value = 1.0 + 2.0 * alpha_;
    const Real gamma = -diagonal_value;
    const Real corner = -alpha_;
    std::vector lower(static_cast<size_t>(size_ - 1), -alpha_);
    std::vector upper(static_cast<size_t>(size_ - 1), -alpha_);
    std::vector diagonal(static_cast<size_t>(size_), diagonal_value);
    diagonal.front() -= gamma;
    diagonal.back() -= corner * corner / gamma;
    factorization_.factorize(lower, diagonal, upper);

    correction_.assign(static_cast<size_t>(size_), 0.0);
    correction_.front() = gamma;
    correction_.back() = corner;
    factorization_.solve_in_place(correction_);
    denominator_ = 1.0 + correction_.front()
        + corner * correction_.back() / gamma;
  }

  void solve(std::vector<Real>& values) const {
    if (size_ <= 1) return;
    if (size_ == 2) {
      const Real diagonal = 1.0 + 2.0 * alpha_;
      const Real off_diagonal = 2.0 * alpha_;
      const Real determinant = 1.0 + 4.0 * alpha_;
      const Real first =
          (diagonal * values[0] + off_diagonal * values[1]) / determinant;
      const Real second =
          (off_diagonal * values[0] + diagonal * values[1]) / determinant;
      values[0] = first;
      values[1] = second;
      return;
    }

    const Real diagonal_value = 1.0 + 2.0 * alpha_;
    const Real gamma = -diagonal_value;
    const Real corner = -alpha_;
    factorization_.solve_in_place(values);
    const Real numerator = values.front() + corner * values.back() / gamma;
    const Real adjustment = numerator / denominator_;
    for (size_t i = 0; i < values.size(); ++i) {
      values[i] -= adjustment * correction_[i];
    }
  }

 private:
  Int size_ = 0;
  Real alpha_ = 0.0;
  Real denominator_ = 1.0;
  TridiagonalFactorization factorization_;
  std::vector<Real> correction_;
};

class NeumannTopLineSolver {
 public:
  NeumannTopLineSolver(Int size, Real alpha) : size_(size), alpha_(alpha) {
    if (size_ <= 0) return;
    std::vector lower(static_cast<size_t>(std::max(size_ - 1, 0)), -alpha_);
    std::vector upper(static_cast<size_t>(std::max(size_ - 1, 0)), -alpha_);
    std::vector diagonal(static_cast<size_t>(size_), 1.0 + 2.0 * alpha_);
    diagonal.back() = 1.0 + alpha_;
    factorization_.factorize(lower, diagonal, upper);
  }

  void solve(std::vector<Real>& values, Real boundary_conc) const {
    if (size_ <= 0) return;
    values.front() += alpha_ * boundary_conc;
    factorization_.solve_in_place(values);
  }

 private:
  Int size_ = 0;
  Real alpha_ = 0.0;
  TridiagonalFactorization factorization_;
};

void diffuse_periodic_x(std::vector<Real>& concentration,
                        const Domain& domain,
                        Real alpha) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const PeriodicLineSolver solver(nx, alpha);

  #ifdef GUTIBM_OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<Real> line(static_cast<size_t>(nx));
    #ifdef GUTIBM_OPENMP
    #pragma omp for collapse(2) schedule(static)
    #endif
    for (Int iz = 0; iz < nz; ++iz) {
      for (Int iy = 0; iy < ny; ++iy) {
        for (Int ix = 0; ix < nx; ++ix) {
          line[static_cast<size_t>(ix)] =
              concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))];
        }
        solver.solve(line);
        for (Int ix = 0; ix < nx; ++ix) {
          concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
              line[static_cast<size_t>(ix)];
        }
      }
    }
  }
}

void diffuse_periodic_y(std::vector<Real>& concentration,
                        const Domain& domain,
                        Real alpha) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const PeriodicLineSolver solver(ny, alpha);

  #ifdef GUTIBM_OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<Real> line(static_cast<size_t>(ny));
    #ifdef GUTIBM_OPENMP
    #pragma omp for collapse(2) schedule(static)
    #endif
    for (Int iz = 0; iz < nz; ++iz) {
      for (Int ix = 0; ix < nx; ++ix) {
        for (Int iy = 0; iy < ny; ++iy) {
          line[static_cast<size_t>(iy)] =
              concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))];
        }
        solver.solve(line);
        for (Int iy = 0; iy < ny; ++iy) {
          concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
              line[static_cast<size_t>(iy)];
        }
      }
    }
  }
}

void diffuse_periodic_x_delivery(
    std::vector<Real>& concentration, const std::vector<Real>& sink_rate,
    std::vector<Real>& realized, const Domain& domain, Real alpha,
    Real sink_dt, Real cell_volume) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int iy = 0; iy < ny; ++iy) {
      std::vector<Real> line(static_cast<size_t>(nx));
      std::vector<Real> sink(static_cast<size_t>(nx));
      for (Int ix = 0; ix < nx; ++ix) {
        const Int cell = domain.cell_index(ix, iy, iz);
        line[static_cast<size_t>(ix)] = concentration[static_cast<size_t>(cell)];
        sink[static_cast<size_t>(ix)] =
            sink_rate[static_cast<size_t>(cell)] * sink_dt;
      }
      solve_periodic_with_sink(line, sink, alpha);
      for (Int ix = 0; ix < nx; ++ix) {
        const Int cell = domain.cell_index(ix, iy, iz);
        concentration[static_cast<size_t>(cell)] = line[static_cast<size_t>(ix)];
        realized[static_cast<size_t>(cell)] +=
            sink[static_cast<size_t>(ix)] * line[static_cast<size_t>(ix)]
            * cell_volume;
      }
    }
  }
}

void diffuse_periodic_y_delivery(
    std::vector<Real>& concentration, const std::vector<Real>& sink_rate,
    std::vector<Real>& realized, const Domain& domain, Real alpha,
    Real sink_dt, Real cell_volume) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int ix = 0; ix < nx; ++ix) {
      std::vector<Real> line(static_cast<size_t>(ny));
      std::vector<Real> sink(static_cast<size_t>(ny));
      for (Int iy = 0; iy < ny; ++iy) {
        const Int cell = domain.cell_index(ix, iy, iz);
        line[static_cast<size_t>(iy)] = concentration[static_cast<size_t>(cell)];
        sink[static_cast<size_t>(iy)] =
            sink_rate[static_cast<size_t>(cell)] * sink_dt;
      }
      solve_periodic_with_sink(line, sink, alpha);
      for (Int iy = 0; iy < ny; ++iy) {
        const Int cell = domain.cell_index(ix, iy, iz);
        concentration[static_cast<size_t>(cell)] = line[static_cast<size_t>(iy)];
        realized[static_cast<size_t>(cell)] +=
            sink[static_cast<size_t>(iy)] * line[static_cast<size_t>(iy)]
            * cell_volume;
      }
    }
  }
}

Real diffuse_bounded_z(std::vector<Real>& concentration,
                       const Domain& domain,
                       Real alpha,
                       Real boundary_conc,
                       Real cell_volume) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  if (nz <= 1) return 0.0;

  const NeumannTopLineSolver solver(nz - 1, alpha);
  Real face_exchange = 0.0;
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<Real> line(static_cast<size_t>(nz - 1));
    #ifdef GUTIBM_OPENMP
    #pragma omp for collapse(2) schedule(static)
    #endif
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        for (Int iz = 1; iz < nz; ++iz) {
          line[static_cast<size_t>(iz - 1)] =
              concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))];
        }
        solver.solve(line, boundary_conc);
        #ifdef GUTIBM_OPENMP
        #pragma omp atomic
        #endif
        face_exchange += alpha * (boundary_conc - line.front()) * cell_volume;
        for (Int iz = 1; iz < nz; ++iz) {
          concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
              line[static_cast<size_t>(iz - 1)];
        }
      }
    }
  }
  return face_exchange;
}

struct DeliveryBoundaryParameters {
  Real boundary_conc = 0.0;
  Real beta = 0.0;
  Real flux_source = 0.0;
  EpithelialBoundaryMode mode = EpithelialBoundaryMode::Flux;
  Real cell_volume = 0.0;
};

Real diffuse_bounded_z_delivery(
    std::vector<Real>& concentration, const std::vector<Real>& sink_rate,
    std::vector<Real>& realized, const Domain& domain, Real alpha,
    const DeliveryBoundaryParameters& params, Real sink_dt) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  if (nz <= 1) return 0.0;
  for (Int iy = 0; iy < ny; ++iy) {
    for (Int ix = 0; ix < nx; ++ix) {
      std::vector<Real> line(static_cast<size_t>(nz - 1));
      std::vector<Real> sink(static_cast<size_t>(nz - 1));
      for (Int iz = 1; iz < nz; ++iz) {
        const Int cell = domain.cell_index(ix, iy, iz);
        line[static_cast<size_t>(iz - 1)] =
            concentration[static_cast<size_t>(cell)];
        sink[static_cast<size_t>(iz - 1)] =
            sink_rate[static_cast<size_t>(cell)] * sink_dt;
      }
      std::vector<Real> diagonal(static_cast<size_t>(nz - 1));
      for (Int iz = 1; iz < nz; ++iz) {
        diagonal[static_cast<size_t>(iz - 1)] =
            1.0 + 2.0 * alpha + sink[static_cast<size_t>(iz - 1)];
      }
      diagonal.back() = 1.0 + alpha + sink.back();
      solve_tridiagonal_with_diagonal(line, diagonal, alpha,
                                      alpha * params.boundary_conc);
      for (Int iz = 1; iz < nz; ++iz) {
        const Int cell = domain.cell_index(ix, iy, iz);
        concentration[static_cast<size_t>(cell)] =
            line[static_cast<size_t>(iz - 1)];
        realized[static_cast<size_t>(cell)] +=
            sink[static_cast<size_t>(iz - 1)]
            * line[static_cast<size_t>(iz - 1)] * params.cell_volume;
      }
    }
  }
  return 0.0;
}

class DeliveryBottomLineSolver {
 public:
  DeliveryBottomLineSolver(Int size, Real alpha,
                           EpithelialBoundaryMode mode, Real beta)
      : size_(size), alpha_(alpha), beta_(beta), mode_(mode) {
    if (size_ <= 0) return;
    std::vector lower(static_cast<size_t>(std::max(size_ - 1, 0)), -alpha_);
    std::vector upper(static_cast<size_t>(std::max(size_ - 1, 0)), -alpha_);
    std::vector diagonal(static_cast<size_t>(size_), 1.0 + 2.0 * alpha_);
    diagonal.front() = mode_ == EpithelialBoundaryMode::Robin
        ? 1.0 + alpha_ + beta_ : 1.0 + alpha_;
    diagonal.back() = 1.0 + alpha_;
    factorization_.factorize(lower, diagonal, upper);
  }

  void solve(std::vector<Real>& values, Real source) const {
    if (size_ <= 0) return;
    values.front() += source;
    factorization_.solve_in_place(values);
  }

 private:
  Int size_ = 0;
  Real alpha_ = 0.0;
  Real beta_ = 0.0;
  EpithelialBoundaryMode mode_ = EpithelialBoundaryMode::Flux;
  TridiagonalFactorization factorization_;
};

struct DeliveryGridParameters {
  Int nx = 0;
  Int ny = 0;
  Int nz = 0;
  Real alpha = 0.0;
};

template <typename LoadLine, typename StoreLine, typename LoadSink,
          typename AddRealized>
struct DeliveryLineOperations {
  LoadLine load_line;
  StoreLine store_line;
  LoadSink load_sink;
  AddRealized add_realized;
};

template <typename LoadLine, typename StoreLine, typename LoadSink,
          typename AddRealized>
void solve_delivery_z_line(
    Int ix, Int iy, const DeliveryGridParameters& grid,
    const DeliveryBoundaryParameters& params, Real sink_dt,
    const DeliveryLineOperations<LoadLine, StoreLine, LoadSink, AddRealized>&
        operations,
    Real& face_exchange) {
  std::vector<Real> line(static_cast<size_t>(grid.nz));
  std::vector<Real> sink(static_cast<size_t>(grid.nz));
  std::vector diagonal(static_cast<size_t>(grid.nz), 0.0);
  operations.load_line(ix, iy, line);
  operations.load_sink(ix, iy, sink);
  for (Int iz = 0; iz < grid.nz; ++iz) {
    diagonal[static_cast<size_t>(iz)] =
        1.0 + 2.0 * grid.alpha
        + sink[static_cast<size_t>(iz)] * sink_dt;
  }
  diagonal.front() += params.mode == EpithelialBoundaryMode::Robin
      ? -grid.alpha + params.beta : -grid.alpha;
  diagonal.back() -= grid.alpha;
  const Real source = params.mode == EpithelialBoundaryMode::Robin
      ? params.beta * params.boundary_conc : params.flux_source;
  solve_tridiagonal_with_diagonal(line, diagonal, grid.alpha, source);
  const Real boundary_realized =
      params.mode == EpithelialBoundaryMode::Robin
          ? params.beta * (params.boundary_conc - line.front())
              * params.cell_volume
          : params.flux_source * params.cell_volume;
  #ifdef GUTIBM_OPENMP
  #pragma omp atomic
  #endif
  face_exchange += boundary_realized;
  for (Int iz = 0; iz < grid.nz; ++iz) {
    const auto index = static_cast<size_t>(iz);
    operations.store_line(ix, iy, iz, line[index]);
    operations.add_realized(
        ix, iy, iz, sink[index] * sink_dt * line[index] * params.cell_volume);
  }
}

template <typename LoadLine, typename StoreLine>
Real diffuse_bounded_z_delivery_impl(
    Int nx, Int ny, Int nz, Real alpha,
    const DeliveryBoundaryParameters& params,
    LoadLine load_line, StoreLine store_line) {
  if (nz <= 0) return 0.0;

  const DeliveryBottomLineSolver solver(nz, alpha, params.mode, params.beta);
  Real face_exchange = 0.0;
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<Real> line(static_cast<size_t>(nz));
    #ifdef GUTIBM_OPENMP
    #pragma omp for collapse(2) schedule(static)
    #endif
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        load_line(ix, iy, line);
        const Real source = params.mode == EpithelialBoundaryMode::Robin
            ? params.beta * params.boundary_conc : params.flux_source;
        solver.solve(line, source);
        const Real realized = params.mode == EpithelialBoundaryMode::Robin
            ? params.beta * (params.boundary_conc - line.front())
                * params.cell_volume
            : params.flux_source * params.cell_volume;
        #ifdef GUTIBM_OPENMP
        #pragma omp atomic
        #endif
        face_exchange += realized;
        store_line(ix, iy, line);
      }
    }
  }
  return face_exchange;
}

Real diffuse_bounded_z_delivery(
    std::vector<Real>& concentration, const Domain& domain, Real alpha,
    const DeliveryBoundaryParameters& params) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const auto load_line = [&concentration, &domain, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] =
          concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))];
    }
  };
  const auto store_line = [&concentration, &domain, nz](
                              Int ix, Int iy,
                              const std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
          line[static_cast<size_t>(iz)];
    }
  };
  return diffuse_bounded_z_delivery_impl(
      nx, ny, nz, alpha, params, load_line, store_line);
}

template <typename LoadLine, typename StoreLine, typename LoadSink,
          typename AddRealized>
Real diffuse_bounded_z_delivery_with_sink_impl(
    const DeliveryGridParameters& grid,
    const DeliveryBoundaryParameters& params, Real sink_dt,
    LoadLine load_line, StoreLine store_line, LoadSink load_sink,
    AddRealized add_realized) {
  const Int nx = grid.nx;
  const Int ny = grid.ny;
  const Int nz = grid.nz;
  if (nz <= 0) return 0.0;

  Real face_exchange = 0.0;
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel
  #endif
  {
    const DeliveryLineOperations<LoadLine, StoreLine, LoadSink, AddRealized>
        operations{load_line, store_line, load_sink, add_realized};
    #ifdef GUTIBM_OPENMP
    #pragma omp for collapse(2) schedule(static)
    #endif
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        solve_delivery_z_line(
            ix, iy, grid, params, sink_dt, operations, face_exchange);
      }
    }
  }
  return face_exchange;
}

Real diffuse_bounded_z_delivery_with_sink(
    std::vector<Real>& concentration, const std::vector<Real>& sink_rate,
    std::vector<Real>& realized, const Domain& domain, Real alpha,
    const DeliveryBoundaryParameters& params, Real sink_dt) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const auto load_line = [&concentration, &domain, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] =
          concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))];
    }
  };
  const auto store_line = [&concentration, &domain](
                              Int ix, Int iy, Int iz, Real value) {
    concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
        value;
  };
  const auto load_sink = [&sink_rate, &domain, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] =
          sink_rate[static_cast<size_t>(domain.cell_index(ix, iy, iz))];
    }
  };
  const auto add_realized = [&realized, &domain](
                                Int ix, Int iy, Int iz, Real amount) {
    realized[static_cast<size_t>(domain.cell_index(ix, iy, iz))] += amount;
  };
  return diffuse_bounded_z_delivery_with_sink_impl(
      {nx, ny, nz, alpha}, params, sink_dt, load_line, store_line, load_sink,
      add_realized);
}

Real set_epithelial_boundary(std::vector<Real>& concentration,
                             const Domain& domain,
                             Real boundary_conc,
                             Real cell_volume) {
  Real amount = 0.0;
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.nx(); ++ix) {
      const auto index = static_cast<size_t>(
          domain.cell_index(ix, iy, 0));
      amount += (boundary_conc - concentration[index]) * cell_volume;
      concentration[index] = boundary_conc;
    }
  }
  return amount;
}

void set_luminal_neumann_boundary(std::vector<Real>& concentration,
                                  const Domain& domain) {
  if (domain.nz() < 2) return;
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.nx(); ++ix) {
      const auto top = static_cast<size_t>(
          domain.cell_index(ix, iy, domain.nz() - 1));
      const auto below = static_cast<size_t>(
          domain.cell_index(ix, iy, domain.nz() - 2));
      concentration[top] = concentration[below];
    }
  }
}

Real z_gradient_reference(const ChemicalSpec& spec,
                          const Domain& domain,
                          Int iz) {
  if (iz == 0) return spec.boundary_conc;
  const Int profile_iz = (domain.nz() >= 2 && iz == domain.nz() - 1)
      ? domain.nz() - 2 : iz;
    const Real z_rel = (profile_iz + 0.5) * domain.dx_z();
  return spec.initial_conc * std::exp(-z_rel / spec.z_gradient_lambda);
}

void shift_z_gradient(std::vector<Real>& concentration,
                      const ChemicalSpec& spec,
                      const Domain& domain,
                      Real scale) {
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    const Real shift = scale * z_gradient_reference(spec, domain, iz);
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.nx(); ++ix) {
        concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] += shift;
      }
    }
  }
}

void clamp_nonnegative(std::vector<Real>& concentration) {
  const auto size = static_cast<Int>(concentration.size());
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Int cell = 0; cell < size; ++cell) {
    concentration[static_cast<size_t>(cell)] =
        std::max(concentration[static_cast<size_t>(cell)], 0.0);
  }
}

void validate_chemical_decomposition(const Domain& domain,
                                    std::string_view decomposition) {
  if (decomposition != "replicated" && decomposition != "slab") {
    throw ConfigError("invalid chemistry decomposition mode");
  }
  if (decomposition != "slab") return;

  const auto required_halo = static_cast<Int>(
      std::ceil(domain.ghost_width() / domain.dx_x()));
  if (domain.grid_halo_width() < required_halo) {
    throw ConfigError(
        "slab chemistry requires grid_halo_width >= ceil(ghost_width / dx)");
  }
  if (domain.local_grid_nx() < domain.grid_halo_width()) {
    throw ConfigError(
        "slab chemistry requires each owned x-slab to be at least "
        "grid_halo_width cells wide");
  }
}

void validate_epithelial_delivery(const std::vector<ChemicalSpec>& specs) {
  for (const auto& spec : specs) {
    if (spec.epithelial_transfer_coeff < 0.0
        || spec.epithelial_flux < 0.0) {
      throw ConfigError(
          "epithelial delivery coefficients and fluxes must be nonnegative");
    }
    if (spec.epithelial_boundary_mode != EpithelialBoundaryMode::Dirichlet
        && spec.z_gradient_enabled) {
      throw ConfigError(
          "z_gradient_enabled cannot be combined with Robin or flux "
          "epithelial boundary modes");
    }
  }
}

}  // namespace

void ChemicalField::init(const Domain& domain,
                          const std::vector<ChemicalSpec>& specs,
                          std::string_view decomposition) {
  validate_chemical_decomposition(domain, decomposition);
  validate_epithelial_delivery(specs);

  domain_ = &domain;
  mode_ = decomposition == "slab"
      ? DecompositionMode::Slab : DecompositionMode::Replicated;
  global_nx_ = domain.nx();
  global_ny_ = domain.ny();
  global_nz_ = domain.nz();
  global_ncells_ = domain.ncells();
  owned_x_begin_ = decomposition == "slab"
      ? domain.local_grid_x_begin() : 0;
  owned_x_end_ = decomposition == "slab"
      ? domain.local_grid_x_end() : domain.nx();
  halo_width_ = decomposition == "slab" ? domain.grid_halo_width() : 0;
  storage_nx_ = decomposition == "slab"
      ? domain.local_grid_storage_nx() : domain.nx();
  specs_  = specs;
  nspec_  = static_cast<Int>(specs.size());
  ncells_ = storage_nx_ * global_ny_ * global_nz_;
  flux_accounting_.init(specs.size());

  conc_.resize(nspec_);
  reac_.resize(nspec_);
  sink_rate_.assign(static_cast<size_t>(ncells_), 0.0);
  sink_realized_.assign(static_cast<size_t>(ncells_), 0.0);
  for (Int s = 0; s < nspec_; ++s) {
    conc_[s].assign(ncells_, specs_[s].initial_conc);
    reac_[s].assign(ncells_, 0.0);

    if (specs_[s].z_gradient_enabled) {
      for (Int storage_cell = 0; storage_cell < ncells_; ++storage_cell) {
        const Int global_cell = storage_to_global_cell(storage_cell);
        if (global_cell < 0) continue;
        const Int iz = global_cell / (global_nx_ * global_ny_);
        const Real z_rel = (iz + 0.5) * domain.dx_z();
        conc_[s][static_cast<size_t>(storage_cell)] =
            specs_[s].initial_conc
            * std::exp(-z_rel / specs_[s].z_gradient_lambda);
      }
    }
  }
}

Int ChemicalField::global_to_storage_cell(Int global_cell) const {
  assert(global_cell >= 0 && global_cell < global_ncells_);
  if (mode_ == DecompositionMode::Replicated) return global_cell;

  const Int ix = global_cell % global_nx_;
  const Int yz = global_cell / global_nx_;
  const Int local_ix = domain_->global_to_local_grid_x(ix);
  if (local_ix < 0) return -1;
  return yz * storage_nx_ + local_ix;
}

Int ChemicalField::storage_to_global_cell(Int storage_cell) const {
  assert(storage_cell >= 0 && storage_cell < ncells_);
  if (mode_ == DecompositionMode::Replicated) return storage_cell;

  const Int local_ix = storage_cell % storage_nx_;
  const Int yz = storage_cell / storage_nx_;
  const Int global_ix = domain_->local_to_global_grid_x(local_ix);
  if (global_ix < 0) return -1;
  return yz * global_nx_ + global_ix;
}

bool ChemicalField::owns_global_cell(Int global_cell) const {
  if (global_cell < 0 || global_cell >= global_ncells_) return false;
  const Int ix = global_cell % global_nx_;
  return mode_ != DecompositionMode::Slab
      || (ix >= owned_x_begin_ && ix < owned_x_end_);
}

bool ChemicalField::global_cell_in_halo(Int global_cell) const {
  if (global_cell < 0 || global_cell >= global_ncells_) return false;
  return mode_ == DecompositionMode::Slab
      && !owns_global_cell(global_cell)
      && global_to_storage_cell(global_cell) >= 0;
}

void ChemicalField::exchange_concentration_halos() {
  if (mode_ != DecompositionMode::Slab) return;

  const Int plane_cells = global_ny_ * global_nz_;
  const Int halo_cells = halo_width_ * plane_cells;
  const Int local_nx = owned_x_end_ - owned_x_begin_;
  if (halo_cells <= 0 || local_nx <= 0) return;

  if (domain_->nprocs() <= 1) {
    for (Int s = 0; s < nspec_; ++s) {
      for (Int iz = 0; iz < global_nz_; ++iz) {
        for (Int iy = 0; iy < global_ny_; ++iy) {
          for (Int h = 0; h < halo_width_; ++h) {
            const Int lo_storage = iz * storage_nx_ * global_ny_
                + iy * storage_nx_ + h;
            const Int hi_storage = iz * storage_nx_ * global_ny_
                + iy * storage_nx_ + halo_width_ + local_nx + h;
            const Int lo_source = iz * storage_nx_ * global_ny_
                + iy * storage_nx_ + halo_width_ + local_nx
                - halo_width_ + h;
            const Int hi_source = iz * storage_nx_ * global_ny_
                + iy * storage_nx_ + halo_width_ + h;
            if (!domain_->config().periodic[0]) continue;
            conc_[static_cast<size_t>(s)][static_cast<size_t>(lo_storage)] =
                conc_[static_cast<size_t>(s)][static_cast<size_t>(lo_source)];
            conc_[static_cast<size_t>(s)][static_cast<size_t>(hi_storage)] =
                conc_[static_cast<size_t>(s)][static_cast<size_t>(hi_source)];
          }
        }
      }
    }
    return;
  }

#ifdef GUTIBM_MPI
  const int rank_lo = domain_->rank_lo() >= 0
      ? domain_->rank_lo() : MPI_PROC_NULL;
  const int rank_hi = domain_->rank_hi() >= 0
      ? domain_->rank_hi() : MPI_PROC_NULL;
  for (Int s = 0; s < nspec_; ++s) {
    std::vector<Real> send_lo(static_cast<size_t>(halo_cells));
    std::vector<Real> send_hi(static_cast<size_t>(halo_cells));
    std::vector<Real> recv_lo(static_cast<size_t>(halo_cells));
    std::vector<Real> recv_hi(static_cast<size_t>(halo_cells));
    for (Int iz = 0; iz < global_nz_; ++iz) {
      for (Int iy = 0; iy < global_ny_; ++iy) {
        for (Int h = 0; h < halo_width_; ++h) {
          const Int lo_storage = iz * storage_nx_ * global_ny_
              + iy * storage_nx_ + halo_width_ + h;
          const Int hi_storage = iz * storage_nx_ * global_ny_
              + iy * storage_nx_ + halo_width_ + local_nx
              - halo_width_ + h;
          const Int offset = (iz * global_ny_ + iy) * halo_width_ + h;
          send_lo[static_cast<size_t>(offset)] =
              conc_[static_cast<size_t>(s)][static_cast<size_t>(lo_storage)];
          send_hi[static_cast<size_t>(offset)] =
              conc_[static_cast<size_t>(s)][static_cast<size_t>(hi_storage)];
        }
      }
    }

    const int tag = 700 + 2 * static_cast<int>(s);
    MPI_Sendrecv(send_lo.data(), halo_cells, MPI_DOUBLE, rank_lo, tag,
                 recv_hi.data(), halo_cells, MPI_DOUBLE, rank_hi, tag,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(send_hi.data(), halo_cells, MPI_DOUBLE, rank_hi, tag + 1,
                 recv_lo.data(), halo_cells, MPI_DOUBLE, rank_lo, tag + 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (Int iz = 0; iz < global_nz_; ++iz) {
      for (Int iy = 0; iy < global_ny_; ++iy) {
        for (Int h = 0; h < halo_width_; ++h) {
          const Int lo_storage = iz * storage_nx_ * global_ny_
              + iy * storage_nx_ + h;
          const Int hi_storage = iz * storage_nx_ * global_ny_
              + iy * storage_nx_ + halo_width_ + local_nx + h;
          const Int offset = (iz * global_ny_ + iy) * halo_width_ + h;
          conc_[static_cast<size_t>(s)][static_cast<size_t>(lo_storage)] =
              recv_lo[static_cast<size_t>(offset)];
          conc_[static_cast<size_t>(s)][static_cast<size_t>(hi_storage)] =
              recv_hi[static_cast<size_t>(offset)];
        }
      }
    }
  }
#endif
}

void ChemicalField::zero_reactions() {
  for (Int s = 0; s < nspec_; ++s) {
    std::ranges::fill(reac_[s], 0.0);
  }
  std::ranges::fill(sink_rate_, 0.0);
  std::ranges::fill(sink_realized_, 0.0);
}

void ChemicalField::add_sink_rate_global(Int cell, Real rate) {
  const Int storage_cell = global_to_storage_cell(cell);
  if (storage_cell < 0 || rate <= 0.0) return;
  #ifdef GUTIBM_OPENMP
  #pragma omp atomic
  #endif
  sink_rate_[static_cast<size_t>(storage_cell)] += rate;
}

Real ChemicalField::sink_realized_global(Int cell) const {
  const Int storage_cell = global_to_storage_cell(cell);
  return storage_cell >= 0 ? sink_realized_[static_cast<size_t>(storage_cell)]
                           : 0.0;
}

bool ChemicalField::has_sink_rate() const {
  const bool local = std::ranges::any_of(
      sink_rate_, [](Real value) { return value > 0.0; });
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized && !finalized) {
    int local_value = local ? 1 : 0;
    int global_value = 0;
    MPI_Allreduce(&local_value, &global_value, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    return global_value != 0;
  }
#endif
  return local;
}

void ChemicalField::sum_reactions_across_ranks() {
#ifdef GUTIBM_MPI
  if (mode_ == DecompositionMode::Slab) return;
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return;

  int ranks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  if (ranks <= 1) return;
  for (auto& reaction : reac_) {
    MPI_Allreduce(MPI_IN_PLACE, reaction.data(), ncells_, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
  }
  MPI_Allreduce(MPI_IN_PLACE, sink_rate_.data(), ncells_, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#endif
}

void ChemicalField::sum_agent_uptake_across_ranks() {
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return;
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.agent_uptake_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.maintenance_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.maintenance_shortfall_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.maintenance_limited_agents_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.uptake_demand_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.uptake_shortfall_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.uptake_limited_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}

void ChemicalField::sum_accounting_across_ranks() {
#ifdef GUTIBM_MPI
  if (mode_ != DecompositionMode::Slab) return;
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return;
  const int count = static_cast<int>(nspec_);
  auto reduce = [count](std::vector<Real>& values) {
    MPI_Allreduce(MPI_IN_PLACE, values.data(), count, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
  };
  reduce(flux_accounting_.boundary_step);
  reduce(flux_accounting_.reaction_clip_step);
#endif
}

namespace {

Int slab_storage_index(
    Int local_ix, Int iy, Int iz, Int storage_nx, Int ny) {
  return iz * storage_nx * ny + iy * storage_nx + local_ix;
}

struct SlabTransportContext {
  std::vector<Real>& concentration;
  const Domain& domain;
  Int storage_nx;
  Int halo_width;
  Real alpha;
  const std::vector<Real>* sink_rate = nullptr;
  std::vector<Real>* realized = nullptr;
  Real sink_dt = 0.0;
  Real cell_volume = 0.0;
};

void diffuse_periodic_x_slab_single(
    const SlabTransportContext& context, const PeriodicLineSolver& solver) {
  auto& concentration = context.concentration;
  const auto& domain = context.domain;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  const Real alpha = context.alpha;
  const auto* sink_rate = context.sink_rate;
  auto* realized = context.realized;
  const Real sink_dt = context.sink_dt;
  const Real cell_volume = context.cell_volume;
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int line_count = ny * domain.nz();
  for (Int line_id = 0; line_id < line_count; ++line_id) {
    const Int iy = line_id % ny;
    const Int iz = line_id / ny;
    std::vector<Real> line(static_cast<size_t>(nx));
    std::vector sink(static_cast<size_t>(nx), 0.0);
    for (Int ix = 0; ix < nx; ++ix) {
      const Int cell = slab_storage_index(
          halo_width + ix, iy, iz, storage_nx, ny);
      line[static_cast<size_t>(ix)] =
          concentration[static_cast<size_t>(cell)];
      if (sink_rate != nullptr) {
        sink[static_cast<size_t>(ix)] =
            (*sink_rate)[static_cast<size_t>(cell)] * sink_dt;
      }
    }
    if (sink_rate != nullptr) {
      solve_periodic_with_sink(line, sink, alpha);
    } else {
      solver.solve(line);
    }
    for (Int ix = 0; ix < nx; ++ix) {
      const Int cell = slab_storage_index(
          halo_width + ix, iy, iz, storage_nx, ny);
      concentration[static_cast<size_t>(cell)] = line[static_cast<size_t>(ix)];
      if (realized != nullptr) {
        (*realized)[static_cast<size_t>(cell)] +=
            sink[static_cast<size_t>(ix)] * line[static_cast<size_t>(ix)]
            * cell_volume;
      }
    }
  }
}

#ifdef GUTIBM_MPI
struct SlabPeriodicXLayout {
  Int nx = 0;
  Int ny = 0;
  Int local_nx = 0;
  Int line_count = 0;
  Int process_count = 0;
  Int local_rank = 0;
  std::vector<Int> line_counts;
  std::vector<Int> line_displacements;
  std::vector<Int> x_counts;
  std::vector<Int> x_displacements;
  std::vector<Int> send_counts;
  std::vector<Int> send_displacements;
  std::vector<Int> recv_counts;
  std::vector<Int> recv_displacements;
};

struct SlabPeriodicXBuffers {
  std::vector<Real> send;
  std::vector<Real> recv;
  std::vector<Real> solved;
  std::vector<Real> sink_send;
  std::vector<Real> sink_recv;
};

SlabPeriodicXLayout make_slab_periodic_x_layout(
    const Domain& domain, Int nx) {
  SlabPeriodicXLayout layout;
  layout.nx = nx;
  layout.ny = domain.ny();
  layout.local_nx = domain.local_grid_nx();
  layout.line_count = layout.ny * domain.nz();
  layout.process_count = domain.nprocs();
  layout.local_rank = domain.rank();
  const auto count = static_cast<size_t>(layout.process_count);
  layout.line_counts.assign(count, 0);
  layout.line_displacements.assign(count, 0);
  layout.x_counts.assign(count, 0);
  layout.x_displacements.assign(count, 0);
  for (Int rank = 0; rank < layout.process_count; ++rank) {
    const auto [begin, end] =
        Domain::grid_x_range_for_rank(nx, layout.process_count, rank);
    layout.x_counts[static_cast<size_t>(rank)] = end - begin;
    layout.x_displacements[static_cast<size_t>(rank)] = begin;
    layout.line_counts[static_cast<size_t>(rank)] =
        (layout.line_count + layout.process_count - rank - 1)
        / layout.process_count;
    if (rank > 0) {
      layout.line_displacements[static_cast<size_t>(rank)] =
          layout.line_displacements[static_cast<size_t>(rank - 1)]
          + layout.line_counts[static_cast<size_t>(rank - 1)];
    }
  }
  layout.send_counts.assign(count, 0);
  layout.send_displacements.assign(count, 0);
  layout.recv_counts.assign(count, 0);
  layout.recv_displacements.assign(count, 0);
  for (Int rank = 0; rank < layout.process_count; ++rank) {
    layout.send_counts[static_cast<size_t>(rank)] =
        layout.line_counts[static_cast<size_t>(rank)] * layout.local_nx;
    layout.recv_counts[static_cast<size_t>(rank)] =
        layout.line_counts[static_cast<size_t>(layout.local_rank)]
        * layout.x_counts[static_cast<size_t>(rank)];
    if (rank > 0) {
      layout.send_displacements[static_cast<size_t>(rank)] =
          layout.send_displacements[static_cast<size_t>(rank - 1)]
          + layout.send_counts[static_cast<size_t>(rank - 1)];
      layout.recv_displacements[static_cast<size_t>(rank)] =
          layout.recv_displacements[static_cast<size_t>(rank - 1)]
          + layout.recv_counts[static_cast<size_t>(rank - 1)];
    }
  }
  return layout;
}

void pack_slab_periodic_x_line(
    const SlabTransportContext& context,
    const SlabPeriodicXLayout& layout, Int iy, Int iz, Int offset,
    SlabPeriodicXBuffers& buffers) {
  const auto& concentration = context.concentration;
  const auto* sink_rate = context.sink_rate;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  for (Int ix = 0; ix < layout.local_nx; ++ix) {
    const Int cell = slab_storage_index(
        halo_width + ix, iy, iz, storage_nx, layout.ny);
    buffers.send[static_cast<size_t>(offset + ix)] =
        concentration[static_cast<size_t>(cell)];
    if (sink_rate != nullptr) {
      buffers.sink_send[static_cast<size_t>(offset + ix)] =
          (*sink_rate)[static_cast<size_t>(cell)] * context.sink_dt;
    }
  }
}

void pack_slab_periodic_x(
    const SlabTransportContext& context,
    const SlabPeriodicXLayout& layout, SlabPeriodicXBuffers& buffers) {
  for (Int destination = 0; destination < layout.process_count;
       ++destination) {
    Int offset = layout.send_displacements[static_cast<size_t>(destination)];
    for (Int line_id = destination; line_id < layout.line_count;
         line_id += layout.process_count) {
      pack_slab_periodic_x_line(
          context, layout, line_id % layout.ny, line_id / layout.ny, offset,
          buffers);
      offset += layout.local_nx;
    }
  }
}

void solve_slab_periodic_x_lines(
    const SlabTransportContext& context,
    const SlabPeriodicXLayout& layout, const PeriodicLineSolver& solver,
    SlabPeriodicXBuffers& buffers) {
  const auto* sink_rate = context.sink_rate;
  std::vector sink(static_cast<size_t>(layout.nx), 0.0);
  const auto gathered_displacements = layout.recv_displacements;
  for (Int line_index = 0;
       line_index < layout.line_counts[static_cast<size_t>(layout.local_rank)];
       ++line_index) {
    std::vector line(static_cast<size_t>(layout.nx), 0.0);
    for (Int source = 0; source < layout.process_count; ++source) {
      const Int segment_offset =
          gathered_displacements[static_cast<size_t>(source)]
          + line_index * layout.x_counts[static_cast<size_t>(source)];
      std::copy_n(
          buffers.recv.begin() + segment_offset,
          layout.x_counts[static_cast<size_t>(source)],
          line.begin() + layout.x_displacements[static_cast<size_t>(source)]);
      if (sink_rate != nullptr) {
        std::copy_n(
            buffers.sink_recv.begin() + segment_offset,
            layout.x_counts[static_cast<size_t>(source)],
            sink.begin() + layout.x_displacements[static_cast<size_t>(source)]);
      }
    }
    if (sink_rate != nullptr) {
      solve_periodic_with_sink(line, sink, context.alpha);
    } else {
      solver.solve(line);
    }
    for (Int destination = 0; destination < layout.process_count;
         ++destination) {
      const Int segment_offset =
          gathered_displacements[static_cast<size_t>(destination)]
          + line_index * layout.x_counts[static_cast<size_t>(destination)];
      std::copy_n(
          line.begin() + layout.x_displacements[static_cast<size_t>(destination)],
          layout.x_counts[static_cast<size_t>(destination)],
          buffers.solved.begin() + segment_offset);
    }
  }
}

void store_slab_periodic_x_line(
    const SlabTransportContext& context,
    const SlabPeriodicXLayout& layout, Int iy, Int iz, Int offset,
    const SlabPeriodicXBuffers& buffers) {
  auto& concentration = context.concentration;
  auto* realized = context.realized;
  const auto* sink_rate = context.sink_rate;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  for (Int ix = 0; ix < layout.local_nx; ++ix) {
    const Int cell = slab_storage_index(
        halo_width + ix, iy, iz, storage_nx, layout.ny);
    concentration[static_cast<size_t>(cell)] =
        buffers.recv[static_cast<size_t>(offset + ix)];
    if (realized != nullptr) {
      (*realized)[static_cast<size_t>(cell)] +=
          (*sink_rate)[static_cast<size_t>(cell)] * context.sink_dt
          * concentration[static_cast<size_t>(cell)] * context.cell_volume;
    }
  }
}

void diffuse_periodic_x_slab_mpi(
    const SlabTransportContext& context, const PeriodicLineSolver& solver) {
  const auto& domain = context.domain;
  const SlabPeriodicXLayout layout =
      make_slab_periodic_x_layout(domain, domain.nx());
  SlabPeriodicXBuffers buffers;
  const auto solved_total = std::accumulate(
      layout.recv_counts.begin(), layout.recv_counts.end(), 0);
  const auto output_total = std::accumulate(
      layout.send_counts.begin(), layout.send_counts.end(), 0);
  buffers.send.resize(static_cast<size_t>(output_total));
  buffers.recv.resize(static_cast<size_t>(solved_total));
  buffers.solved.resize(static_cast<size_t>(solved_total));
  if (context.sink_rate != nullptr) {
    buffers.sink_send.resize(static_cast<size_t>(output_total));
    buffers.sink_recv.resize(static_cast<size_t>(solved_total));
  }
  pack_slab_periodic_x(context, layout, buffers);
  MPI_Alltoallv(
      buffers.send.data(), layout.send_counts.data(),
      layout.send_displacements.data(), MPI_DOUBLE, buffers.recv.data(),
      layout.recv_counts.data(), layout.recv_displacements.data(), MPI_DOUBLE,
      MPI_COMM_WORLD);
  if (context.sink_rate != nullptr) {
    MPI_Alltoallv(
        buffers.sink_send.data(), layout.send_counts.data(),
        layout.send_displacements.data(), MPI_DOUBLE, buffers.sink_recv.data(),
        layout.recv_counts.data(), layout.recv_displacements.data(), MPI_DOUBLE,
        MPI_COMM_WORLD);
  }
  solve_slab_periodic_x_lines(context, layout, solver, buffers);
  const auto output_counts = layout.recv_counts;
  const auto output_displacements = layout.recv_displacements;
  const auto solved_counts = layout.send_counts;
  const auto solved_displacements = layout.send_displacements;
  const auto gathered_displacements = layout.recv_displacements;
  buffers.send.assign(
      static_cast<size_t>(std::accumulate(
          output_counts.begin(), output_counts.end(), 0)), 0.0);
  buffers.recv.assign(
      static_cast<size_t>(std::accumulate(
          solved_counts.begin(), solved_counts.end(), 0)), 0.0);
  for (Int destination = 0; destination < layout.process_count;
       ++destination) {
    const Int segment_length =
        layout.x_counts[static_cast<size_t>(destination)];
    Int offset = output_displacements[static_cast<size_t>(destination)];
    for (Int line_index = 0;
         line_index < layout.line_counts[static_cast<size_t>(layout.local_rank)];
         ++line_index) {
      const Int source_offset =
          gathered_displacements[static_cast<size_t>(destination)]
          + line_index * segment_length;
      std::copy_n(
          buffers.solved.begin() + source_offset, segment_length,
          buffers.send.begin() + offset);
      offset += segment_length;
    }
  }
  MPI_Alltoallv(
      buffers.send.data(), output_counts.data(), output_displacements.data(),
      MPI_DOUBLE, buffers.recv.data(), solved_counts.data(),
      solved_displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  for (Int source = 0; source < layout.process_count; ++source) {
    Int offset = layout.send_displacements[static_cast<size_t>(source)];
    for (Int line_id = source; line_id < layout.line_count;
         line_id += layout.process_count) {
      store_slab_periodic_x_line(
          context, layout, line_id % layout.ny, line_id / layout.ny, offset,
          buffers);
      offset += layout.local_nx;
    }
  }
}
#endif

void diffuse_periodic_x_slab(
    const SlabTransportContext& context) {
  const PeriodicLineSolver solver(
      context.domain.nx(), context.alpha);
  if (context.domain.nprocs() == 1) {
    diffuse_periodic_x_slab_single(context, solver);
    return;
  }
#ifdef GUTIBM_MPI
  diffuse_periodic_x_slab_mpi(context, solver);
#else
  (void)solver;
#endif
}

void diffuse_periodic_y_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width, Real alpha) {
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const Int local_nx = domain.local_grid_nx();
  const PeriodicLineSolver solver(ny, alpha);
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int ix = 0; ix < local_nx; ++ix) {
      std::vector<Real> line(static_cast<size_t>(ny));
      for (Int iy = 0; iy < ny; ++iy) {
        line[static_cast<size_t>(iy)] = concentration[
            static_cast<size_t>(slab_storage_index(
                halo_width + ix, iy, iz, storage_nx, ny))];
      }
      solver.solve(line);
      for (Int iy = 0; iy < ny; ++iy) {
        concentration[static_cast<size_t>(slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny))] =
            line[static_cast<size_t>(iy)];
      }
    }
  }
}

void diffuse_periodic_y_slab_delivery(
    const SlabTransportContext& context) {
  auto& concentration = context.concentration;
  const auto& domain = context.domain;
  const auto& sink_rate = *context.sink_rate;
  auto& realized = *context.realized;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  const Real alpha = context.alpha;
  const Real sink_dt = context.sink_dt;
  const Real cell_volume = context.cell_volume;
  const Int nx = domain.local_grid_nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int ix = 0; ix < nx; ++ix) {
      std::vector<Real> line(static_cast<size_t>(ny));
      std::vector<Real> sink(static_cast<size_t>(ny));
      for (Int iy = 0; iy < ny; ++iy) {
        const Int cell = slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny);
        line[static_cast<size_t>(iy)] = concentration[static_cast<size_t>(cell)];
        sink[static_cast<size_t>(iy)] =
            sink_rate[static_cast<size_t>(cell)] * sink_dt;
      }
      solve_periodic_with_sink(line, sink, alpha);
      for (Int iy = 0; iy < ny; ++iy) {
        const Int cell = slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny);
        concentration[static_cast<size_t>(cell)] = line[static_cast<size_t>(iy)];
        realized[static_cast<size_t>(cell)] +=
            sink[static_cast<size_t>(iy)] * line[static_cast<size_t>(iy)]
            * cell_volume;
      }
    }
  }
}

Real diffuse_bounded_z_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width, Real alpha, Real boundary_conc,
    Real cell_volume) {
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const Int local_nx = domain.local_grid_nx();
  if (nz <= 1) return 0.0;
  const NeumannTopLineSolver solver(nz - 1, alpha);
  Real face_exchange = 0.0;
  for (Int iy = 0; iy < ny; ++iy) {
    for (Int ix = 0; ix < local_nx; ++ix) {
      std::vector<Real> line(static_cast<size_t>(nz - 1));
      for (Int iz = 1; iz < nz; ++iz) {
        line[static_cast<size_t>(iz - 1)] = concentration[
            static_cast<size_t>(slab_storage_index(
                halo_width + ix, iy, iz, storage_nx, ny))];
      }
      solver.solve(line, boundary_conc);
      face_exchange += alpha * (boundary_conc - line.front()) * cell_volume;
      for (Int iz = 1; iz < nz; ++iz) {
        concentration[static_cast<size_t>(slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny))] =
            line[static_cast<size_t>(iz - 1)];
      }
    }
  }
  return face_exchange;
}

Real diffuse_bounded_z_slab_delivery(
    const SlabTransportContext& context, Real boundary_conc) {
  auto& concentration = context.concentration;
  const auto& domain = context.domain;
  const auto& sink_rate = *context.sink_rate;
  auto& realized = *context.realized;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  const Real alpha = context.alpha;
  const Real cell_volume = context.cell_volume;
  const Real sink_dt = context.sink_dt;
  const Int nx = domain.local_grid_nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  if (nz <= 1) return 0.0;
  for (Int iy = 0; iy < ny; ++iy) {
    for (Int ix = 0; ix < nx; ++ix) {
      std::vector<Real> line(static_cast<size_t>(nz - 1));
      std::vector<Real> sink(static_cast<size_t>(nz - 1));
      for (Int iz = 1; iz < nz; ++iz) {
        const Int cell = slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny);
        line[static_cast<size_t>(iz - 1)] = concentration[
            static_cast<size_t>(cell)];
        sink[static_cast<size_t>(iz - 1)] =
            sink_rate[static_cast<size_t>(cell)] * sink_dt;
      }
      std::vector<Real> diagonal(static_cast<size_t>(nz - 1));
      for (Int iz = 1; iz < nz; ++iz) {
        diagonal[static_cast<size_t>(iz - 1)] =
            1.0 + 2.0 * alpha + sink[static_cast<size_t>(iz - 1)];
      }
      diagonal.back() = 1.0 + alpha + sink.back();
      solve_tridiagonal_with_diagonal(line, diagonal, alpha,
                                      alpha * boundary_conc);
      for (Int iz = 1; iz < nz; ++iz) {
        const Int cell = slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny);
        concentration[static_cast<size_t>(cell)] =
            line[static_cast<size_t>(iz - 1)];
        realized[static_cast<size_t>(cell)] +=
            sink[static_cast<size_t>(iz - 1)]
            * line[static_cast<size_t>(iz - 1)] * cell_volume;
      }
    }
  }
  return 0.0;
}

Real diffuse_bounded_z_delivery_slab(
    std::vector<Real>& concentration, const Domain& domain, Int storage_nx,
    Int halo_width, Real alpha,
    const DeliveryBoundaryParameters& params) {
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const Int local_nx = domain.local_grid_nx();
  const auto load_line = [&concentration, storage_nx, halo_width, ny, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] = concentration[
          static_cast<size_t>(slab_storage_index(
              halo_width + ix, iy, iz, storage_nx, ny))];
    }
  };
  const auto store_line = [&concentration, storage_nx, halo_width, ny, nz](
                              Int ix, Int iy,
                              const std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      concentration[static_cast<size_t>(slab_storage_index(
          halo_width + ix, iy, iz, storage_nx, ny))] =
          line[static_cast<size_t>(iz)];
    }
  };
  return diffuse_bounded_z_delivery_impl(
      local_nx, ny, nz, alpha, params, load_line, store_line);
}

Real diffuse_bounded_z_delivery_with_sink_slab(
    const SlabTransportContext& context,
    const DeliveryBoundaryParameters& params) {
  auto& concentration = context.concentration;
  const auto& sink_rate = *context.sink_rate;
  auto& realized = *context.realized;
  const auto& domain = context.domain;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  const Real alpha = context.alpha;
  const Real sink_dt = context.sink_dt;
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const Int local_nx = domain.local_grid_nx();
  const auto load_line = [&concentration, storage_nx, halo_width, ny, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] = concentration[
          static_cast<size_t>(slab_storage_index(
              halo_width + ix, iy, iz, storage_nx, ny))];
    }
  };
  const auto store_line = [&concentration, storage_nx, halo_width, ny](
                              Int ix, Int iy, Int iz, Real value) {
    concentration[static_cast<size_t>(slab_storage_index(
        halo_width + ix, iy, iz, storage_nx, ny))] = value;
  };
  const auto load_sink = [&sink_rate, storage_nx, halo_width, ny, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] = sink_rate[
          static_cast<size_t>(slab_storage_index(
              halo_width + ix, iy, iz, storage_nx, ny))];
    }
  };
  const auto add_realized = [&realized, storage_nx, halo_width, ny](
                                Int ix, Int iy, Int iz, Real amount) {
    realized[static_cast<size_t>(slab_storage_index(
        halo_width + ix, iy, iz, storage_nx, ny))] += amount;
  };
  return diffuse_bounded_z_delivery_with_sink_impl(
      {local_nx, ny, nz, alpha}, params, sink_dt, load_line, store_line,
      load_sink, add_realized);
}

Real set_epithelial_boundary_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width, Real boundary_conc, Real cell_volume) {
  Real amount = 0.0;
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.local_grid_nx(); ++ix) {
      const Int index = slab_storage_index(
          halo_width + ix, iy, 0, storage_nx, domain.ny());
      amount += (boundary_conc - concentration[static_cast<size_t>(index)])
          * cell_volume;
      concentration[static_cast<size_t>(index)] = boundary_conc;
    }
  }
  return amount;
}

void set_luminal_neumann_boundary_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width) {
  if (domain.nz() < 2) return;
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.local_grid_nx(); ++ix) {
      const Int top = slab_storage_index(
          halo_width + ix, iy, domain.nz() - 1, storage_nx, domain.ny());
      const Int below = slab_storage_index(
          halo_width + ix, iy, domain.nz() - 2, storage_nx, domain.ny());
      concentration[static_cast<size_t>(top)] =
          concentration[static_cast<size_t>(below)];
    }
  }
}

void shift_z_gradient_slab(
    std::vector<Real>& concentration, const ChemicalSpec& spec,
    const Domain& domain, Int storage_nx, Int halo_width, Real scale) {
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    const Int profile_iz = domain.nz() >= 2 && iz == domain.nz() - 1
        ? domain.nz() - 2 : iz;
  const Real z_rel = (profile_iz + 0.5) * domain.dx_z();
    const Real shift = scale * spec.initial_conc
        * std::exp(-z_rel / spec.z_gradient_lambda);
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.local_grid_nx(); ++ix) {
        concentration[static_cast<size_t>(slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, domain.ny()))] += shift;
      }
    }
  }
}

void clamp_nonnegative_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width) {
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.local_grid_nx(); ++ix) {
        const Int index = slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, domain.ny());
        concentration[static_cast<size_t>(index)] =
            std::max(concentration[static_cast<size_t>(index)], 0.0);
      }
    }
  }
}

struct ReplicatedDiffusionContext {
  std::vector<Real>& concentration;
  const Domain& domain;
  const ChemicalSpec& chemical;
  NutrientFluxAccounting& flux;
  std::vector<Real>& sink_rate;
  std::vector<Real>& sink_realized;
  Int spec;
  Real dt;
  Real alpha_x;
  Real alpha_y;
  Real alpha_z;
  Real cell_volume;
  Real diffusion_boundary;
  bool preserve_gradient;
  bool delivery;
};

void prepare_replicated_diffusion(ReplicatedDiffusionContext& context) {
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    context.flux.add_boundary(
        context.spec, set_epithelial_boundary(
            context.concentration, context.domain,
            context.chemical.boundary_conc, context.cell_volume));
  }
  if (!context.preserve_gradient) return;
  set_luminal_neumann_boundary(context.concentration, context.domain);
  shift_z_gradient(
      context.concentration, context.chemical, context.domain, -1.0);
  context.diffusion_boundary = 0.0;
}

void transport_replicated_periodic(
    const ReplicatedDiffusionContext& context) {
  if (context.delivery) {
    const Real sink_dt = context.dt / 3.0;
    diffuse_periodic_x_delivery(
        context.concentration, context.sink_rate, context.sink_realized,
        context.domain, context.alpha_x, sink_dt, context.cell_volume);
    diffuse_periodic_y_delivery(
        context.concentration, context.sink_rate, context.sink_realized,
        context.domain, context.alpha_y, sink_dt, context.cell_volume);
  } else {
    diffuse_periodic_x(
        context.concentration, context.domain, context.alpha_x);
    diffuse_periodic_y(
        context.concentration, context.domain, context.alpha_y);
  }
}

Real transport_replicated_z(
    const ReplicatedDiffusionContext& context) {
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    return context.delivery
        ? diffuse_bounded_z_delivery(
              context.concentration, context.sink_rate,
              context.sink_realized, context.domain, context.alpha_z,
              {context.diffusion_boundary, 0.0, 0.0,
               EpithelialBoundaryMode::Dirichlet, context.cell_volume},
              context.dt / 3.0)
        : diffuse_bounded_z(
              context.concentration, context.domain, context.alpha_z,
              context.diffusion_boundary, context.cell_volume);
  }
  const Real beta = context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Robin
      ? context.chemical.epithelial_transfer_coeff * context.dt
          / context.domain.dx_z() : 0.0;
  const Real flux_source = context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Flux
      ? context.chemical.epithelial_flux * context.dt
          / context.domain.dx_z() : 0.0;
  const Real boundary = context.delivery
      ? diffuse_bounded_z_delivery_with_sink(
            context.concentration, context.sink_rate,
            context.sink_realized, context.domain, context.alpha_z,
            {context.chemical.boundary_conc, beta, flux_source,
             context.chemical.epithelial_boundary_mode, context.cell_volume},
            context.dt / 3.0)
      : diffuse_bounded_z_delivery(
            context.concentration, context.domain, context.alpha_z,
            {context.chemical.boundary_conc, beta, flux_source,
             context.chemical.epithelial_boundary_mode, context.cell_volume});
  return boundary;
}

void transport_replicated_diffusion(
    ReplicatedDiffusionContext& context) {
  transport_replicated_periodic(context);
  context.flux.add_boundary(
      context.spec, transport_replicated_z(context));
}

void finish_replicated_diffusion(ReplicatedDiffusionContext& context) {
  if (context.preserve_gradient) {
    shift_z_gradient(
        context.concentration, context.chemical, context.domain, 1.0);
    set_luminal_neumann_boundary(context.concentration, context.domain);
  }
  clamp_nonnegative(context.concentration);
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    context.flux.add_boundary(
        context.spec, set_epithelial_boundary(
            context.concentration, context.domain,
            context.chemical.boundary_conc, context.cell_volume));
  }
}

struct SlabDiffusionContext {
  std::vector<Real>& concentration;
  const Domain& domain;
  const ChemicalSpec& chemical;
  NutrientFluxAccounting& flux;
  std::vector<Real>& sink_rate;
  std::vector<Real>& sink_realized;
  Int storage_nx;
  Int halo_width;
  Int spec;
  Real dt;
  Real alpha_x;
  Real alpha_y;
  Real alpha_z;
  Real cell_volume;
  Real diffusion_boundary;
  bool preserve_gradient;
  bool delivery;
};

void prepare_slab_diffusion(SlabDiffusionContext& context) {
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    context.flux.add_boundary(
        context.spec, set_epithelial_boundary_slab(
            context.concentration, context.domain, context.storage_nx,
            context.halo_width, context.chemical.boundary_conc,
            context.cell_volume));
  }
  if (!context.preserve_gradient) return;
  set_luminal_neumann_boundary_slab(
      context.concentration, context.domain, context.storage_nx,
      context.halo_width);
  shift_z_gradient_slab(
      context.concentration, context.chemical, context.domain,
      context.storage_nx, context.halo_width, -1.0);
  context.diffusion_boundary = 0.0;
}

void transport_slab_periodic(const SlabDiffusionContext& context) {
  const SlabTransportContext transport{
      context.concentration, context.domain, context.storage_nx,
      context.halo_width, context.alpha_x,
      context.delivery ? &context.sink_rate : nullptr,
      context.delivery ? &context.sink_realized : nullptr,
      context.dt / 3.0, context.cell_volume};
  diffuse_periodic_x_slab(transport);
  if (context.delivery) {
    diffuse_periodic_y_slab_delivery(
        {context.concentration, context.domain, context.storage_nx,
         context.halo_width, context.alpha_y, &context.sink_rate,
         &context.sink_realized, context.dt / 3.0, context.cell_volume});
  } else {
    diffuse_periodic_y_slab(
        context.concentration, context.domain, context.storage_nx,
        context.halo_width, context.alpha_y);
  }
}

Real transport_slab_z(const SlabDiffusionContext& context) {
  const SlabTransportContext z_transport{
      context.concentration, context.domain, context.storage_nx,
      context.halo_width, context.alpha_z,
      context.delivery ? &context.sink_rate : nullptr,
      context.delivery ? &context.sink_realized : nullptr,
      context.dt / 3.0, context.cell_volume};
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    const Real boundary = context.delivery
        ? diffuse_bounded_z_slab_delivery(
              z_transport, context.diffusion_boundary)
        : diffuse_bounded_z_slab(
              context.concentration, context.domain, context.storage_nx,
              context.halo_width, context.alpha_z,
              context.diffusion_boundary, context.cell_volume);
    return boundary;
  }
  const Real beta = context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Robin
      ? context.chemical.epithelial_transfer_coeff * context.dt
          / context.domain.dx_z() : 0.0;
  const Real flux_source = context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Flux
      ? context.chemical.epithelial_flux * context.dt
          / context.domain.dx_z() : 0.0;
  const DeliveryBoundaryParameters boundary_params{
      context.chemical.boundary_conc, beta, flux_source,
      context.chemical.epithelial_boundary_mode, context.cell_volume};
  const Real boundary = context.delivery
      ? diffuse_bounded_z_delivery_with_sink_slab(
            z_transport, boundary_params)
      : diffuse_bounded_z_delivery_slab(
            context.concentration, context.domain, context.storage_nx,
            context.halo_width, context.alpha_z, boundary_params);
  return boundary;
}

void transport_slab_diffusion(SlabDiffusionContext& context) {
  transport_slab_periodic(context);
  context.flux.add_boundary(context.spec, transport_slab_z(context));
}

void finish_slab_diffusion(SlabDiffusionContext& context) {
  if (context.preserve_gradient) {
    shift_z_gradient_slab(
        context.concentration, context.chemical, context.domain,
        context.storage_nx, context.halo_width, 1.0);
    set_luminal_neumann_boundary_slab(
        context.concentration, context.domain, context.storage_nx,
        context.halo_width);
  }
  clamp_nonnegative_slab(
      context.concentration, context.domain, context.storage_nx,
      context.halo_width);
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    context.flux.add_boundary(
        context.spec, set_epithelial_boundary_slab(
            context.concentration, context.domain, context.storage_nx,
            context.halo_width, context.chemical.boundary_conc,
            context.cell_volume));
  }
}

}  // namespace

void ChemicalField::apply_diffusion(const Domain& domain, Real dt) {
  if (dt <= 0.0 || domain.dx_x() <= 0.0 || domain.dx_y() <= 0.0
      || domain.dx_z() <= 0.0) return;
  if (mode_ == DecompositionMode::Slab) {
    apply_diffusion_slab(domain, dt);
    return;
  }
  for (Int s = 0; s < nspec_; ++s) {
    apply_diffusion_species(domain, dt, s);
  }
}

void ChemicalField::apply_diffusion_species(
    const Domain& domain, Real dt, Int s) {
  const ChemicalSpec& chemical = specs_[s];
  if (!chemical.diffusion_enabled || chemical.diff_coeff <= 0.0
      || chemical.retardation <= 0.0) {
    return;
  }
  const Real effective_diffusion = chemical.diff_coeff / chemical.retardation;
  const Real alpha_x = effective_diffusion * dt
      / (domain.dx_x() * domain.dx_x());
  const Real alpha_y = effective_diffusion * dt
      / (domain.dx_y() * domain.dx_y());
  const Real alpha_z = effective_diffusion * dt
      / (domain.dx_z() * domain.dx_z());
  ReplicatedDiffusionContext context{
      conc_[s], domain, chemical, flux_accounting_, sink_rate_,
      sink_realized_, s, dt, alpha_x, alpha_y, alpha_z,
      domain.cell_volume(), chemical.boundary_conc,
      chemical.z_gradient_enabled && chemical.z_gradient_lambda > 0.0,
      s == find(species::CARBON) && has_sink_rate()};
  prepare_replicated_diffusion(context);
  transport_replicated_diffusion(context);
  finish_replicated_diffusion(context);
}

void ChemicalField::apply_periodic_x_diffusion(const Domain& domain, Real dt) {
  if (dt <= 0.0 || domain.dx_x() <= 0.0) return;
  for (Int s = 0; s < nspec_; ++s) {
    apply_periodic_x_diffusion(domain, dt, s);
  }
}

void ChemicalField::apply_periodic_x_diffusion(const Domain& domain, Real dt,
                                               Int spec) {
  if (dt <= 0.0 || domain.dx_x() <= 0.0 || spec < 0 || spec >= nspec_) return;
  const ChemicalSpec& chemical = specs_[static_cast<size_t>(spec)];
  if (!chemical.diffusion_enabled || chemical.diff_coeff <= 0.0
      || chemical.retardation <= 0.0) {
    return;
  }
  const Real alpha = (chemical.diff_coeff / chemical.retardation) * dt
      / (domain.dx_x() * domain.dx_x());
  if (mode_ == DecompositionMode::Slab) {
    diffuse_periodic_x_slab(
        {conc_[static_cast<size_t>(spec)], domain, storage_nx_, halo_width_,
         alpha});
  } else {
    diffuse_periodic_x(conc_[static_cast<size_t>(spec)], domain, alpha);
  }
}

void ChemicalField::apply_periodic_y_diffusion(const Domain& domain, Real dt,
                                               Int spec) {
  if (dt <= 0.0 || domain.dx_y() <= 0.0 || spec < 0 || spec >= nspec_) return;
  const ChemicalSpec& chemical = specs_[static_cast<size_t>(spec)];
  if (!chemical.diffusion_enabled || chemical.diff_coeff <= 0.0
      || chemical.retardation <= 0.0) {
    return;
  }
  const Real alpha = (chemical.diff_coeff / chemical.retardation) * dt
      / (domain.dx_y() * domain.dx_y());
  diffuse_periodic_y(conc_[static_cast<size_t>(spec)], domain, alpha);
}

void ChemicalField::apply_bounded_z_diffusion(const Domain& domain, Real dt,
                                              Int spec) {
  if (dt <= 0.0 || domain.dx_z() <= 0.0 || spec < 0 || spec >= nspec_) return;
  const ChemicalSpec& chemical = specs_[static_cast<size_t>(spec)];
  if (!chemical.diffusion_enabled || chemical.diff_coeff <= 0.0
      || chemical.retardation <= 0.0) {
    return;
  }
  const Real effective_diffusion =
      chemical.diff_coeff / chemical.retardation;
  const Real alpha = effective_diffusion * dt
      / (domain.dx_z() * domain.dx_z());
  if (chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    flux_accounting_.add_boundary(
        spec, diffuse_bounded_z(
                  conc_[static_cast<size_t>(spec)], domain, alpha,
                  chemical.boundary_conc, domain.cell_volume()));
  } else {
    const Real beta = chemical.epithelial_boundary_mode
        == EpithelialBoundaryMode::Robin
        ? chemical.epithelial_transfer_coeff * dt / domain.dx_z() : 0.0;
    const Real flux_source = chemical.epithelial_boundary_mode
        == EpithelialBoundaryMode::Flux
        ? chemical.epithelial_flux * dt / domain.dx_z() : 0.0;
    flux_accounting_.add_boundary(
        spec, diffuse_bounded_z_delivery(
                  conc_[static_cast<size_t>(spec)], domain, alpha,
                  {chemical.boundary_conc, beta, flux_source,
                   chemical.epithelial_boundary_mode, domain.cell_volume()}));
  }
}

void ChemicalField::apply_diffusion_slab(const Domain& domain, Real dt) {
  for (Int s = 0; s < nspec_; ++s) {
    apply_diffusion_slab_species(domain, dt, s);
  }
}

void ChemicalField::apply_diffusion_slab_species(
    const Domain& domain, Real dt, Int s) {
  const ChemicalSpec& chemical = specs_[s];
  if (!chemical.diffusion_enabled || chemical.diff_coeff <= 0.0
      || chemical.retardation <= 0.0) {
    return;
  }
  const Real effective_diffusion = chemical.diff_coeff / chemical.retardation;
  const Real alpha_x = effective_diffusion * dt
      / (domain.dx_x() * domain.dx_x());
  const Real alpha_y = effective_diffusion * dt
      / (domain.dx_y() * domain.dx_y());
  const Real alpha_z = effective_diffusion * dt
      / (domain.dx_z() * domain.dx_z());
  SlabDiffusionContext context{
      conc_[static_cast<size_t>(s)], domain, chemical, flux_accounting_,
      sink_rate_, sink_realized_, storage_nx_, halo_width_, s, dt,
      alpha_x, alpha_y, alpha_z, domain.cell_volume(),
      chemical.boundary_conc,
      chemical.z_gradient_enabled && chemical.z_gradient_lambda > 0.0,
      s == find(species::CARBON) && has_sink_rate()};
  prepare_slab_diffusion(context);
  transport_slab_diffusion(context);
  finish_slab_diffusion(context);
}

namespace {

void apply_epithelial_boundary_layer(
    std::vector<std::vector<Real>>& concentration,
    const Domain& domain, Int species_index, Real boundary_conc,
    NutrientFluxAccounting& flux_accounting) {
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.nx(); ++ix) {
      const Int idx = domain.cell_index(ix, iy, 0);
      const Real old = concentration[static_cast<size_t>(species_index)]
          [static_cast<size_t>(idx)];
      concentration[static_cast<size_t>(species_index)]
          [static_cast<size_t>(idx)] = boundary_conc;
      flux_accounting.add_boundary(
          species_index,
          (boundary_conc - old) * domain.cell_volume());
    }
  }
}

void mirror_non_diffusing_top_layer(std::vector<std::vector<Real>>& concentration,
                                    const Domain& domain, Int species_index) {
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.nx(); ++ix) {
      const Int top = domain.cell_index(ix, iy, domain.nz() - 1);
      const Int below = domain.cell_index(ix, iy, domain.nz() - 2);
      concentration[static_cast<size_t>(species_index)]
          [static_cast<size_t>(top)] =
              concentration[static_cast<size_t>(species_index)]
                  [static_cast<size_t>(below)];
    }
  }
}

}  // namespace

void ChemicalField::apply_boundaries(const Domain& domain) {
  if (mode_ == DecompositionMode::Slab) {
    apply_boundaries_slab(domain);
    return;
  }
  const Int nz = domain.nz();

  for (Int s = 0; s < nspec_; ++s) {
    const Real bc = specs_[s].boundary_conc;

    // z=0 (epithelial surface): Dirichlet for nutrients. When a z-gradient is
    // configured, this is the peak concentration at the epithelium.
    if (specs_[s].epithelial_boundary_mode
        == EpithelialBoundaryMode::Dirichlet) {
      apply_epithelial_boundary_layer(
          conc_, domain, s, bc, flux_accounting_);
    }

    // The implicit z solve enforces the luminal zero-flux condition directly.
    // Non-diffusing fields retain the legacy mirrored top layer.
    if (!specs_[s].diffusion_enabled && nz >= 2) {
      mirror_non_diffusing_top_layer(conc_, domain, s);
    }
  }
}

void ChemicalField::apply_boundaries_slab(const Domain& domain) {
  const Real cell_volume = domain.cell_volume();
  for (Int s = 0; s < nspec_; ++s) {
    auto& concentration = conc_[static_cast<size_t>(s)];
    if (specs_[s].epithelial_boundary_mode
        == EpithelialBoundaryMode::Dirichlet) {
      flux_accounting_.add_boundary(
          s, set_epithelial_boundary_slab(
                 concentration, domain, storage_nx_, halo_width_,
                 specs_[s].boundary_conc, cell_volume));
    }
    if (!specs_[s].diffusion_enabled && domain.nz() >= 2) {
      set_luminal_neumann_boundary_slab(
          concentration, domain, storage_nx_, halo_width_);
    }
  }
}

Int ChemicalField::find(std::string_view name) const {
  for (Int i = 0; i < nspec_; ++i) {
    if (specs_[i].name == name) return i;
  }
  return -1;
}

}  // namespace gutibm
