/* -----------------------------------------------------------------------
   GutIBM – Domain implementation
   ----------------------------------------------------------------------- */

#include "domain.h"
#include "error.h"
#include <cmath>
#include <algorithm>

namespace gutibm {

void Domain::init(const DomainConfig& cfg) {
  cfg_      = cfg;
  lo_       = cfg.lo;
  hi_       = cfg.hi;
  periodic_ = cfg.periodic;
  dx_       = cfg.grid_dx;

  Vec3 sz = size();
  nx_ = std::max(1, static_cast<Int>(std::round(sz[0] / dx_)));
  ny_ = std::max(1, static_cast<Int>(std::round(sz[1] / dx_)));
  nz_ = std::max(1, static_cast<Int>(std::round(sz[2] / dx_)));
  grid_halo_width_ = std::max(0, cfg.grid_halo_width);

#ifdef GUTIBM_MPI
  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);
  if (mpi_initialized) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_);
  }
#endif

  decompose();

  hash_.init(lo_, hi_, cfg.hash_cell_size);
}

void Domain::pos_to_grid(const Vec3& pos, Int& ix, Int& iy, Int& iz) const {
  ix = static_cast<Int>(std::floor((pos[0] - lo_[0]) / dx_));
  iy = static_cast<Int>(std::floor((pos[1] - lo_[1]) / dx_));
  iz = static_cast<Int>(std::floor((pos[2] - lo_[2]) / dx_));
  ix = std::clamp(ix, 0, nx_ - 1);
  iy = std::clamp(iy, 0, ny_ - 1);
  iz = std::clamp(iz, 0, nz_ - 1);
}

Vec3 Domain::cell_center(Int ix, Int iy, Int iz) const {
  return {
    lo_[0] + (ix + 0.5) * dx_,
    lo_[1] + (iy + 0.5) * dx_,
    lo_[2] + (iz + 0.5) * dx_
  };
}

void Domain::apply_pbc(Vec3& pos) const {
  Vec3 sz = size();
  for (int d = 0; d < 3; ++d) {
    if (!periodic_[d]) {
      pos[d] = std::clamp(pos[d], lo_[d], hi_[d]);
      continue;
    }
    while (pos[d] < lo_[d]) pos[d] += sz[d];
    while (pos[d] >= hi_[d]) pos[d] -= sz[d];
  }
}

Vec3 Domain::min_image_delta(const Vec3& pos_i, const Vec3& pos_j) const {
  Vec3 delta;
  Vec3 sz = size();
  for (int d = 0; d < 3; ++d) {
    delta[d] = pos_j[d] - pos_i[d];
    if (periodic_[d]) {
      if (delta[d] >  0.5 * sz[d]) delta[d] -= sz[d];
      if (delta[d] < -0.5 * sz[d]) delta[d] += sz[d];
    }
  }
  return delta;
}

Real Domain::min_image_dist_sq(const Vec3& a, const Vec3& b) const {
  Vec3 d = min_image_delta(a, b);
  return d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
}

void Domain::decompose() {
  Int axis = cfg_.mpi_decomp_axis;  // 0 = x
  Real global_lo = lo_[axis];
  Real global_hi = hi_[axis];
  const auto grid_range = grid_x_range_for_rank(nx_, nprocs_, rank_);
  local_lo_x_ = global_lo + grid_range.first * dx_;
  local_hi_x_ = global_lo + grid_range.second * dx_;
  if (rank_ == nprocs_ - 1) local_hi_x_ = global_hi;

  ghost_width_ = cfg_.ghost_width;

  // Determine neighbor ranks
  bool axis_periodic = periodic_[axis];
  rank_lo_ = rank_ - 1;
  rank_hi_ = rank_ + 1;

  if (rank_lo_ < 0) {
    rank_lo_ = axis_periodic ? (nprocs_ - 1) : -1;
  }
  if (rank_hi_ >= nprocs_) {
    rank_hi_ = axis_periodic ? 0 : -1;
  }

  local_grid_x_begin_ = grid_range.first;
  local_grid_x_end_ = grid_range.second;
}

std::pair<Int, Int> Domain::grid_x_range_for_rank(
    Int global_nx, Int nprocs, Int rank) {
  if (global_nx < 0) {
    throw ConfigError("grid partition requires a non-negative cell count");
  }
  if (nprocs <= 0) {
    throw ConfigError("grid partition requires at least one rank");
  }
  if (rank < 0 || rank >= nprocs) {
    throw ConfigError("grid partition rank is outside the process range");
  }
  const Int begin = (global_nx * rank) / nprocs;
  const Int end = (global_nx * (rank + 1)) / nprocs;
  return {begin, end};
}

Int Domain::grid_owner_rank_for_cell(
    Int global_nx, Int nprocs, Int global_ix) {
  if (global_nx <= 0) {
    throw ConfigError("grid ownership requires at least one cell");
  }
  if (nprocs <= 0) {
    throw ConfigError("grid ownership requires at least one rank");
  }
  if (global_ix < 0 || global_ix >= global_nx) {
    throw ConfigError("grid ownership cell is outside the global range");
  }
  const Int rank = ((global_ix + 1) * nprocs - 1) / global_nx;
  return std::clamp(rank, 0, nprocs - 1);
}

Int Domain::global_to_local_grid_x(Int global_ix) const {
  if (nx_ <= 0) return -1;

  Int candidate = global_ix;
  if (!periodic_[0] && (candidate < 0 || candidate >= nx_)) return -1;
  if (periodic_[0]) {
    const Int lower = local_grid_x_begin_ - grid_halo_width_;
    const Int upper = local_grid_x_end_ + grid_halo_width_;
    if (candidate >= lower && candidate < upper) {
      return candidate - local_grid_x_begin_ + grid_halo_width_;
    }
    candidate %= nx_;
    if (candidate < 0) candidate += nx_;
    const Int periods = static_cast<Int>(std::floor(
        static_cast<Real>(local_grid_x_begin_ - candidate) / nx_));
    candidate += periods * nx_;
    if (candidate < lower) candidate += nx_;
    if (candidate >= upper) candidate -= nx_;
  }

  if (candidate < local_grid_x_begin_ - grid_halo_width_
      || candidate >= local_grid_x_end_ + grid_halo_width_) {
    return -1;
  }
  return candidate - local_grid_x_begin_ + grid_halo_width_;
}

Int Domain::local_to_global_grid_x(Int local_ix) const {
  const Int storage_nx = local_grid_storage_nx();
  if (local_ix < 0 || local_ix >= storage_nx) return -1;

  Int global_ix = local_grid_x_begin_ + local_ix - grid_halo_width_;
  if (periodic_[0]) {
    global_ix %= nx_;
    if (global_ix < 0) global_ix += nx_;
    return global_ix;
  }
  if (global_ix < 0 || global_ix >= nx_) return -1;
  return global_ix;
}

bool Domain::is_local(const Vec3& pos) const {
  if (nprocs_ <= 1) return true;
  Int axis = cfg_.mpi_decomp_axis;
  return pos[axis] >= local_lo_x_ && pos[axis] < local_hi_x_;
}

Int Domain::owner_rank(const Vec3& pos) const {
  if (nprocs_ <= 1) return 0;
  Int axis = cfg_.mpi_decomp_axis;
  Real global_lo = lo_[axis];
  auto ix = static_cast<Int>(std::floor((pos[axis] - global_lo) / dx_));
  ix = std::clamp(ix, 0, nx_ - 1);
  return grid_owner_rank_for_cell(nx_, nprocs_, ix);
}

}  // namespace gutibm
