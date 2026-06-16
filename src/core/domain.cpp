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

  hash_.init(lo_, hi_, cfg.hash_cell_size);

#ifdef GUTIBM_MPI
  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);
  if (mpi_initialized) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_);
  }
#endif
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

}  // namespace gutibm
