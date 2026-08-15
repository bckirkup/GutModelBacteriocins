/* -----------------------------------------------------------------------
   GutIBM – Chemical field implementation
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "domain.h"
#include "error.h"
#include <algorithm>
#include <cmath>
#include <numeric>
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
  const auto size = static_cast<Int>(concentration.size());
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Int cell = 0; cell < size; ++cell) {
    concentration[static_cast<size_t>(cell)] =
        std::max(concentration[static_cast<size_t>(cell)], 0.0);
  }
}

}  // namespace

void ChemicalField::init(const Domain& domain,
                          const std::vector<ChemicalSpec>& specs,
                          std::string_view decomposition) {
  if (decomposition != "replicated" && decomposition != "slab") {
    throw ConfigError("invalid chemistry decomposition mode");
  }
  if (decomposition == "slab") {
    const auto required_halo = static_cast<Int>(
        std::ceil(domain.ghost_width() / domain.dx()));
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
  for (Int s = 0; s < nspec_; ++s) {
    conc_[s].assign(ncells_, specs_[s].initial_conc);
    reac_[s].assign(ncells_, 0.0);

    if (specs_[s].z_gradient_enabled) {
      for (Int storage_cell = 0; storage_cell < ncells_; ++storage_cell) {
        const Int global_cell = storage_to_global_cell(storage_cell);
        if (global_cell < 0) continue;
        const Int iz = global_cell / (global_nx_ * global_ny_);
        const Real z_rel = (iz + 0.5) * domain.dx();
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

void diffuse_periodic_x_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width, Real alpha) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const Int local_nx = domain.local_grid_nx();
  const PeriodicLineSolver solver(nx, alpha);
  const Int line_count = ny * nz;
  const Int process_count = domain.nprocs();
  std::vector<Int> line_counts(static_cast<size_t>(process_count), 0);
  std::vector<Int> line_displacements(static_cast<size_t>(process_count), 0);
  std::vector<Int> x_counts(static_cast<size_t>(process_count), 0);
  std::vector<Int> x_displacements(static_cast<size_t>(process_count), 0);
  for (Int rank = 0; rank < process_count; ++rank) {
    const auto [begin, end] =
        Domain::grid_x_range_for_rank(nx, process_count, rank);
    x_counts[static_cast<size_t>(rank)] = end - begin;
    x_displacements[static_cast<size_t>(rank)] = begin;
    line_counts[static_cast<size_t>(rank)] =
        (line_count + process_count - rank - 1) / process_count;
    line_displacements[static_cast<size_t>(rank)] =
        rank == 0 ? 0 : line_displacements[static_cast<size_t>(rank - 1)]
            + line_counts[static_cast<size_t>(rank - 1)];
  }

  if (process_count == 1) {
    for (Int line_id = 0; line_id < line_count; ++line_id) {
      const Int iy = line_id % ny;
      const Int iz = line_id / ny;
      std::vector<Real> line(static_cast<size_t>(nx));
      for (Int ix = 0; ix < local_nx; ++ix) {
        line[static_cast<size_t>(ix)] = concentration[
            static_cast<size_t>(slab_storage_index(
                halo_width + ix, iy, iz, storage_nx, ny))];
      }
      solver.solve(line);
      for (Int ix = 0; ix < local_nx; ++ix) {
        concentration[static_cast<size_t>(slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny))] =
            line[static_cast<size_t>(ix)];
      }
    }
    return;
  }

#ifdef GUTIBM_MPI
  std::vector<Int> send_counts(static_cast<size_t>(process_count), 0);
  std::vector<Int> send_displacements(
      static_cast<size_t>(process_count), 0);
  std::vector<Int> recv_counts(static_cast<size_t>(process_count), 0);
  std::vector<Int> recv_displacements(
      static_cast<size_t>(process_count), 0);
  for (Int rank = 0; rank < process_count; ++rank) {
    send_counts[static_cast<size_t>(rank)] =
        line_counts[static_cast<size_t>(rank)] * local_nx;
    send_displacements[static_cast<size_t>(rank)] =
        rank == 0 ? 0 : send_displacements[static_cast<size_t>(rank - 1)]
            + send_counts[static_cast<size_t>(rank - 1)];
    recv_counts[static_cast<size_t>(rank)] =
        line_counts[static_cast<size_t>(domain.rank())]
        * x_counts[static_cast<size_t>(rank)];
    recv_displacements[static_cast<size_t>(rank)] =
        rank == 0 ? 0 : recv_displacements[static_cast<size_t>(rank - 1)]
            + recv_counts[static_cast<size_t>(rank - 1)];
  }
  const std::vector<Int> gathered_displacements = recv_displacements;

  const Int send_total = std::accumulate(send_counts.begin(),
                                         send_counts.end(), 0);
  const Int recv_total = std::accumulate(recv_counts.begin(),
                                         recv_counts.end(), 0);
  std::vector<Real> send_buffer(static_cast<size_t>(send_total));
  std::vector<Real> recv_buffer(static_cast<size_t>(recv_total));
  for (Int destination = 0; destination < process_count; ++destination) {
    Int offset = send_displacements[static_cast<size_t>(destination)];
    for (Int line_id = destination; line_id < line_count;
         line_id += process_count) {
      const Int iy = line_id % ny;
      const Int iz = line_id / ny;
      for (Int ix = 0; ix < local_nx; ++ix) {
        send_buffer[static_cast<size_t>(offset + ix)] = concentration[
            static_cast<size_t>(slab_storage_index(
                halo_width + ix, iy, iz, storage_nx, ny))];
      }
      offset += local_nx;
    }
  }
  MPI_Alltoallv(send_buffer.data(), send_counts.data(),
                send_displacements.data(), MPI_DOUBLE,
                recv_buffer.data(), recv_counts.data(),
                recv_displacements.data(), MPI_DOUBLE,
                MPI_COMM_WORLD);

  std::vector<Real> solved_buffer(static_cast<size_t>(recv_total));
  const Int local_rank = domain.rank();
  for (Int line_index = 0; line_index < line_counts[static_cast<size_t>(local_rank)];
       ++line_index) {
    std::vector<Real> line(static_cast<size_t>(nx));
    for (Int source = 0; source < process_count; ++source) {
      const Int segment_offset =
          gathered_displacements[static_cast<size_t>(source)]
          + line_index * x_counts[static_cast<size_t>(source)];
      std::copy_n(recv_buffer.begin() + segment_offset,
                  x_counts[static_cast<size_t>(source)],
                  line.begin() + x_displacements[static_cast<size_t>(source)]);
    }
    solver.solve(line);
    for (Int destination = 0; destination < process_count; ++destination) {
      const Int segment_offset =
          gathered_displacements[static_cast<size_t>(destination)]
          + line_index * x_counts[static_cast<size_t>(destination)];
      std::copy_n(line.begin() + x_displacements[static_cast<size_t>(destination)],
                  x_counts[static_cast<size_t>(destination)],
                  solved_buffer.begin() + segment_offset);
    }
  }

  std::fill(send_counts.begin(), send_counts.end(), 0);
  std::fill(send_displacements.begin(), send_displacements.end(), 0);
  std::fill(recv_counts.begin(), recv_counts.end(), 0);
  std::fill(recv_displacements.begin(), recv_displacements.end(), 0);
  for (Int rank = 0; rank < process_count; ++rank) {
    send_counts[static_cast<size_t>(rank)] =
        line_counts[static_cast<size_t>(local_rank)]
        * x_counts[static_cast<size_t>(rank)];
    recv_counts[static_cast<size_t>(rank)] =
        line_counts[static_cast<size_t>(rank)] * local_nx;
    if (rank > 0) {
      send_displacements[static_cast<size_t>(rank)] =
          send_displacements[static_cast<size_t>(rank - 1)]
          + send_counts[static_cast<size_t>(rank - 1)];
      recv_displacements[static_cast<size_t>(rank)] =
          recv_displacements[static_cast<size_t>(rank - 1)]
          + recv_counts[static_cast<size_t>(rank - 1)];
    }
  }
  const Int solved_total = std::accumulate(send_counts.begin(),
                                           send_counts.end(), 0);
  const Int output_total = std::accumulate(recv_counts.begin(),
                                           recv_counts.end(), 0);
  send_buffer.assign(static_cast<size_t>(solved_total), 0.0);
  recv_buffer.assign(static_cast<size_t>(output_total), 0.0);
  for (Int destination = 0; destination < process_count; ++destination) {
    const Int segment_length = x_counts[static_cast<size_t>(destination)];
    Int offset = send_displacements[static_cast<size_t>(destination)];
    for (Int line_index = 0;
         line_index < line_counts[static_cast<size_t>(local_rank)];
         ++line_index) {
      const Int source_offset =
          gathered_displacements[static_cast<size_t>(destination)]
          + line_index * segment_length;
      std::copy_n(solved_buffer.begin() + source_offset, segment_length,
                  send_buffer.begin() + offset);
      offset += segment_length;
    }
  }
  MPI_Alltoallv(send_buffer.data(), send_counts.data(),
                send_displacements.data(), MPI_DOUBLE,
                recv_buffer.data(), recv_counts.data(),
                recv_displacements.data(), MPI_DOUBLE,
                MPI_COMM_WORLD);
  for (Int source = 0; source < process_count; ++source) {
    Int offset = recv_displacements[static_cast<size_t>(source)];
    for (Int line_id = source; line_id < line_count;
         line_id += process_count) {
      const Int iy = line_id % ny;
      const Int iz = line_id / ny;
      for (Int ix = 0; ix < local_nx; ++ix) {
        concentration[static_cast<size_t>(slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny))] =
            recv_buffer[static_cast<size_t>(offset + ix)];
      }
      offset += local_nx;
    }
  }
#else
  (void)concentration;
  (void)domain;
  (void)storage_nx;
  (void)halo_width;
  (void)alpha;
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
    const Real z_rel = (profile_iz + 0.5) * domain.dx();
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

}  // namespace

void ChemicalField::apply_diffusion(const Domain& domain, Real dt) {
  if (dt <= 0.0 || domain.dx() <= 0.0) return;
  if (mode_ == DecompositionMode::Slab) {
    apply_diffusion_slab(domain, dt);
    return;
  }

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

    const Real cell_volume = domain.dx() * domain.dx() * domain.dx();
    flux_accounting_.add_boundary(
        s, set_epithelial_boundary(concentration, domain,
                                   chemical.boundary_conc, cell_volume));
    if (preserve_gradient) {
      // The configured z-gradient is an environmental background profile.
      // Diffuse reaction-driven departures from it rather than erasing it.
      // Its reintroduction is a prescribed bulk profile, not epithelial
      // Dirichlet injection, so it is intentionally excluded from the
      // boundary flux counter.
      set_luminal_neumann_boundary(concentration, domain);
      shift_z_gradient(concentration, chemical, domain, -1.0);
      diffusion_boundary = 0.0;
    }

    diffuse_periodic_x(concentration, domain, alpha);
    diffuse_periodic_y(concentration, domain, alpha);
    flux_accounting_.add_boundary(
        s, diffuse_bounded_z(concentration, domain, alpha,
                             diffusion_boundary, cell_volume));

    if (preserve_gradient) {
      shift_z_gradient(concentration, chemical, domain, 1.0);
      set_luminal_neumann_boundary(concentration, domain);
    }
    clamp_nonnegative(concentration);
    flux_accounting_.add_boundary(
        s, set_epithelial_boundary(concentration, domain,
                                   chemical.boundary_conc, cell_volume));
  }
}

void ChemicalField::apply_periodic_x_diffusion(const Domain& domain, Real dt) {
  if (dt <= 0.0 || domain.dx() <= 0.0) return;
  for (Int s = 0; s < nspec_; ++s) {
    apply_periodic_x_diffusion(domain, dt, s);
  }
}

void ChemicalField::apply_periodic_x_diffusion(const Domain& domain, Real dt,
                                               Int spec) {
  if (dt <= 0.0 || domain.dx() <= 0.0 || spec < 0 || spec >= nspec_) return;
  const ChemicalSpec& chemical = specs_[static_cast<size_t>(spec)];
  if (!chemical.diffusion_enabled || chemical.diff_coeff <= 0.0
      || chemical.retardation <= 0.0) {
    return;
  }
  const Real dx2 = domain.dx() * domain.dx();
  const Real alpha = (chemical.diff_coeff / chemical.retardation) * dt / dx2;
  if (mode_ == DecompositionMode::Slab) {
    diffuse_periodic_x_slab(conc_[static_cast<size_t>(spec)], domain,
                            storage_nx_, halo_width_, alpha);
  } else {
    diffuse_periodic_x(conc_[static_cast<size_t>(spec)], domain, alpha);
  }
}

void ChemicalField::apply_diffusion_slab(const Domain& domain, Real dt) {
  const Real dx2 = domain.dx() * domain.dx();
  for (Int s = 0; s < nspec_; ++s) {
    const ChemicalSpec& chemical = specs_[s];
    if (!chemical.diffusion_enabled || chemical.diff_coeff <= 0.0
        || chemical.retardation <= 0.0) {
      continue;
    }
    const Real alpha =
        (chemical.diff_coeff / chemical.retardation) * dt / dx2;
    auto& concentration = conc_[static_cast<size_t>(s)];
    const bool preserve_gradient =
        chemical.z_gradient_enabled && chemical.z_gradient_lambda > 0.0;
    Real diffusion_boundary = chemical.boundary_conc;
    const Real cell_volume = dx2 * domain.dx();
    flux_accounting_.add_boundary(
        s, set_epithelial_boundary_slab(
               concentration, domain, storage_nx_, halo_width_,
               chemical.boundary_conc, cell_volume));
    if (preserve_gradient) {
      set_luminal_neumann_boundary_slab(
          concentration, domain, storage_nx_, halo_width_);
      shift_z_gradient_slab(
          concentration, chemical, domain, storage_nx_, halo_width_, -1.0);
      diffusion_boundary = 0.0;
    }
    diffuse_periodic_x_slab(
        concentration, domain, storage_nx_, halo_width_, alpha);
    diffuse_periodic_y_slab(
        concentration, domain, storage_nx_, halo_width_, alpha);
    flux_accounting_.add_boundary(
        s, diffuse_bounded_z_slab(
               concentration, domain, storage_nx_, halo_width_, alpha,
               diffusion_boundary, cell_volume));
    if (preserve_gradient) {
      shift_z_gradient_slab(
          concentration, chemical, domain, storage_nx_, halo_width_, 1.0);
      set_luminal_neumann_boundary_slab(
          concentration, domain, storage_nx_, halo_width_);
    }
    clamp_nonnegative_slab(
        concentration, domain, storage_nx_, halo_width_);
    flux_accounting_.add_boundary(
        s, set_epithelial_boundary_slab(
               concentration, domain, storage_nx_, halo_width_,
               chemical.boundary_conc, cell_volume));
  }
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
          (boundary_conc - old) * domain.dx() * domain.dx() * domain.dx());
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
    apply_epithelial_boundary_layer(conc_, domain, s, bc, flux_accounting_);

    // The implicit z solve enforces the luminal zero-flux condition directly.
    // Non-diffusing fields retain the legacy mirrored top layer.
    if (!specs_[s].diffusion_enabled && nz >= 2) {
      mirror_non_diffusing_top_layer(conc_, domain, s);
    }
  }
}

void ChemicalField::apply_boundaries_slab(const Domain& domain) {
  const Real cell_volume = domain.dx() * domain.dx() * domain.dx();
  for (Int s = 0; s < nspec_; ++s) {
    auto& concentration = conc_[static_cast<size_t>(s)];
    flux_accounting_.add_boundary(
        s, set_epithelial_boundary_slab(
               concentration, domain, storage_nx_, halo_width_,
               specs_[s].boundary_conc, cell_volume));
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
