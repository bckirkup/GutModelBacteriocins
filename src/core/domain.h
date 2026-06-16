/* -----------------------------------------------------------------------
   GutIBM – Simulation domain with MPI decomposition
   Models a 3D section of the colonic mucus layer.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_DOMAIN_H
#define GUTIBM_DOMAIN_H

#include "types.h"
#include "spatial_hash.h"

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {

struct DomainConfig {
  // Physical domain (meters) — colonic mucus slab
  Vec3 lo = {0.0, 0.0, 0.0};
  Vec3 hi = {1.0e-3, 1.0e-3, 100.0e-6};  // 1mm x 1mm x 100um mucus

  // Periodicity: x,y periodic (lateral), z non-periodic (epithelium→lumen)
  std::array<bool, 3> periodic = {true, true, false};

  // Grid resolution for chemical fields
  Real grid_dx = 2.0e-6;  // 2 um

  // Spatial hash cell size (should be >= max interaction range)
  Real hash_cell_size = 10.0e-6;  // 10 um
};

class Domain {
 public:
  Domain() = default;

  void init(const DomainConfig& cfg);

  // Global domain bounds
  const Vec3& lo() const { return lo_; }
  const Vec3& hi() const { return hi_; }
  Vec3 size() const {
    return {hi_[0] - lo_[0], hi_[1] - lo_[1], hi_[2] - lo_[2]};
  }

  // Grid dimensions
  Int nx() const { return nx_; }
  Int ny() const { return ny_; }
  Int nz() const { return nz_; }
  Int ncells() const { return nx_ * ny_ * nz_; }
  Real dx() const { return dx_; }

  // Cell index from grid coordinates
  Int cell_index(Int ix, Int iy, Int iz) const {
    return iz * (nx_ * ny_) + iy * nx_ + ix;
  }

  // Grid coords from position
  void pos_to_grid(const Vec3& pos, Int& ix, Int& iy, Int& iz) const;

  // Cell center position
  Vec3 cell_center(Int ix, Int iy, Int iz) const;

  // Apply periodic boundary conditions
  void apply_pbc(Vec3& pos) const;

  // Minimum image displacement (pos_j - pos_i) with PBC
  Vec3 min_image_delta(const Vec3& pos_i, const Vec3& pos_j) const;

  Real min_image_dist_sq(const Vec3& a, const Vec3& b) const;

  // MPI decomposition
  Int rank() const { return rank_; }
  Int nprocs() const { return nprocs_; }

  // Spatial hash
  SpatialHash& spatial_hash() { return hash_; }
  const SpatialHash& spatial_hash() const { return hash_; }

  const DomainConfig& config() const { return cfg_; }

 private:
  DomainConfig cfg_;
  Vec3 lo_{}, hi_{};
  Int nx_ = 0, ny_ = 0, nz_ = 0;
  Real dx_ = 2.0e-6;

  std::array<bool, 3> periodic_{};

  SpatialHash hash_;

  Int rank_   = 0;
  Int nprocs_ = 1;
};

}  // namespace gutibm

#endif  // GUTIBM_DOMAIN_H
