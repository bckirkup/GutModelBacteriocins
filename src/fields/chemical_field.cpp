/* -----------------------------------------------------------------------
   GutIBM – Chemical field implementation
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "domain.h"
#include <algorithm>
#include <cmath>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {

namespace {

class TridiagonalFactorization {
 public:
  void factorize(const std::vector<Real>& lower,
                 const std::vector<Real>& diagonal,
                 const std::vector<Real>& upper) {
    diagonal_ = diagonal;
    upper_ = upper;
    multipliers_.assign(lower.size(), 0.0);
    for (size_t i = 1; i < diagonal_.size(); ++i) {
      const Real multiplier = lower[i - 1] / diagonal_[i - 1];
      multipliers_[i - 1] = multiplier;
      diagonal_[i] -= multiplier * upper_[i - 1];
    }
  }

  void solve_in_place(std::vector<Real>& values) const {
    for (size_t i = 1; i < values.size(); ++i) {
      values[i] -= multipliers_[i - 1] * values[i - 1];
    }
    values.back() /= diagonal_.back();
    for (size_t i = values.size() - 1; i > 0; --i) {
      values[i - 1] =
          (values[i - 1] - upper_[i - 1] * values[i]) / diagonal_[i - 1];
    }
  }

 private:
  std::vector<Real> diagonal_;
  std::vector<Real> upper_;
  std::vector<Real> multipliers_;
};

class PeriodicLineSolver {
 public:
  PeriodicLineSolver(Int size, Real alpha) : size_(size), alpha_(alpha) {
    if (size_ < 3) return;

    const Real diagonal_value = 1.0 + 2.0 * alpha_;
    const Real gamma = -diagonal_value;
    const Real corner = -alpha_;
    std::vector<Real> lower(static_cast<size_t>(size_ - 1), -alpha_);
    std::vector<Real> upper(static_cast<size_t>(size_ - 1), -alpha_);
    std::vector<Real> diagonal(static_cast<size_t>(size_), diagonal_value);
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
    std::vector<Real> lower(static_cast<size_t>(std::max(size_ - 1, 0)), -alpha_);
    std::vector<Real> upper(static_cast<size_t>(std::max(size_ - 1, 0)), -alpha_);
    std::vector<Real> diagonal(static_cast<size_t>(size_), 1.0 + 2.0 * alpha_);
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

void diffuse_bounded_z(std::vector<Real>& concentration,
                       const Domain& domain,
                       Real alpha,
                       Real boundary_conc) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  if (nz <= 1) return;

  const NeumannTopLineSolver solver(nz - 1, alpha);
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
        for (Int iz = 1; iz < nz; ++iz) {
          concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
              line[static_cast<size_t>(iz - 1)];
        }
      }
    }
  }
}

void set_epithelial_boundary(std::vector<Real>& concentration,
                             const Domain& domain,
                             Real boundary_conc) {
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.nx(); ++ix) {
      concentration[static_cast<size_t>(domain.cell_index(ix, iy, 0))] = boundary_conc;
    }
  }
}

void set_luminal_neumann_boundary(std::vector<Real>& concentration,
                                  const Domain& domain) {
  if (domain.nz() < 2) return;
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.nx(); ++ix) {
      const size_t top = static_cast<size_t>(
          domain.cell_index(ix, iy, domain.nz() - 1));
      const size_t below = static_cast<size_t>(
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
  const Real z_rel = (profile_iz + 0.5) * domain.dx();
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
  const Int size = static_cast<Int>(concentration.size());
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Int cell = 0; cell < size; ++cell) {
    concentration[static_cast<size_t>(cell)] =
        std::max(concentration[static_cast<size_t>(cell)], 0.0);
  }
}

void apply_z_gradient(std::vector<Real>& conc_row,
                      const ChemicalSpec& spec,
                      const Domain& domain) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  for (Int iz = 0; iz < nz; ++iz) {
    Real z_rel = (iz + 0.5) * domain.dx();
    Real factor = std::exp(-z_rel / spec.z_gradient_lambda);
    Real conc = spec.initial_conc * factor;
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        conc_row[domain.cell_index(ix, iy, iz)] = conc;
      }
    }
  }
}

}  // namespace

void ChemicalField::init(const Domain& domain,
                          const std::vector<ChemicalSpec>& specs) {
  specs_  = specs;
  nspec_  = static_cast<Int>(specs.size());
  ncells_ = domain.ncells();

  conc_.resize(nspec_);
  reac_.resize(nspec_);
  for (Int s = 0; s < nspec_; ++s) {
    conc_[s].assign(ncells_, specs_[s].initial_conc);
    reac_[s].assign(ncells_, 0.0);

    if (specs_[s].z_gradient_enabled) {
      apply_z_gradient(conc_[s], specs_[s], domain);
    }
  }
}

void ChemicalField::zero_reactions() {
  for (Int s = 0; s < nspec_; ++s) {
    std::ranges::fill(reac_[s], 0.0);
  }
}

void ChemicalField::sum_reactions_across_ranks() {
#ifdef GUTIBM_MPI
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
#endif
}

void ChemicalField::apply_diffusion(const Domain& domain, Real dt) {
  if (dt <= 0.0 || domain.dx() <= 0.0) return;

  const Real dx2 = domain.dx() * domain.dx();
  for (Int s = 0; s < nspec_; ++s) {
    const ChemicalSpec& chemical = specs_[s];
    if (!chemical.diffusion_enabled || chemical.diff_coeff <= 0.0
        || chemical.retardation <= 0.0) {
      continue;
    }

    // Explicit diffusion would be catastrophically unstable at the biological
    // timestep (alpha is 5,040 for O2 at dx=5 um, dt=60 s). Backward-Euler
    // directional splitting is L-stable, positivity-preserving, and O(N).
    const Real effective_diffusion = chemical.diff_coeff / chemical.retardation;
    const Real alpha = effective_diffusion * dt / dx2;
    auto& concentration = conc_[s];
    const bool preserve_gradient =
        chemical.z_gradient_enabled && chemical.z_gradient_lambda > 0.0;
    Real diffusion_boundary = chemical.boundary_conc;

    set_epithelial_boundary(concentration, domain, chemical.boundary_conc);
    if (preserve_gradient) {
      // The configured z-gradient is an environmental background profile.
      // Diffuse reaction-driven departures from it rather than erasing it.
      set_luminal_neumann_boundary(concentration, domain);
      shift_z_gradient(concentration, chemical, domain, -1.0);
      diffusion_boundary = 0.0;
    }

    diffuse_periodic_x(concentration, domain, alpha);
    diffuse_periodic_y(concentration, domain, alpha);
    diffuse_bounded_z(concentration, domain, alpha, diffusion_boundary);

    if (preserve_gradient) {
      shift_z_gradient(concentration, chemical, domain, 1.0);
      set_luminal_neumann_boundary(concentration, domain);
    }
    clamp_nonnegative(concentration);
    set_epithelial_boundary(concentration, domain, chemical.boundary_conc);
  }
}

void ChemicalField::apply_boundaries(const Domain& domain) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();

  for (Int s = 0; s < nspec_; ++s) {
    const Real bc = specs_[s].boundary_conc;

    // z=0 (epithelial surface): Dirichlet for nutrients. When a z-gradient is
    // configured, this is the peak concentration at the epithelium.
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        const Int idx = domain.cell_index(ix, iy, 0);
        conc_[s][idx] = bc;
      }
    }

    // The implicit z solve enforces the luminal zero-flux condition directly.
    // Non-diffusing fields retain the legacy mirrored top layer.
    if (!specs_[s].diffusion_enabled && nz >= 2) {
      for (Int iy = 0; iy < ny; ++iy) {
        for (Int ix = 0; ix < nx; ++ix) {
          const Int top = domain.cell_index(ix, iy, nz - 1);
          const Int below = domain.cell_index(ix, iy, nz - 2);
          conc_[s][top] = conc_[s][below];
        }
      }
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
