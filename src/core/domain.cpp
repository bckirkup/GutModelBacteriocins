/* -----------------------------------------------------------------------
   GutIBM – Domain implementation
   ----------------------------------------------------------------------- */

#include "domain.h"
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
  Real total_len = global_hi - global_lo;
  Real slab_width = total_len / nprocs_;

  local_lo_x_ = global_lo + rank_ * slab_width;
  local_hi_x_ = global_lo + (rank_ + 1) * slab_width;

  // Last rank absorbs any floating-point remainder
  if (rank_ == nprocs_ - 1) {
    local_hi_x_ = global_hi;
  }

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
  Real total_len = hi_[axis] - global_lo;
  Real slab_width = total_len / nprocs_;

  auto r = static_cast<Int>(std::floor((pos[axis] - global_lo) / slab_width));
  return std::clamp(r, 0, nprocs_ - 1);
}

}  // namespace gutibm
