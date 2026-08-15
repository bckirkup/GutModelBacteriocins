/* -----------------------------------------------------------------------
   GutIBM – Simulation domain with MPI decomposition
   Models a 3D section of the colonic mucus layer.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_DOMAIN_H
#define GUTIBM_DOMAIN_H

#include "types.h"
#include "spatial_hash.h"

#include <utility>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {

inline bool is_first_periodic_offset(Int offset, Int count, Int span,
                                     bool periodic) {
  return !periodic || offset - count < -span;
}

struct DomainConfig {
  // Physical domain (meters) — colonic mucus slab
  Vec3 lo = {0.0, 0.0, 0.0};
  Vec3 hi = {1.0e-3, 1.0e-3, 100.0e-6};  // 1mm x 1mm x 100um mucus

  // Periodicity: x,y periodic (lateral), z non-periodic (epithelium→lumen)
  std::array<bool, 3> periodic = {true, true, false};

  // Grid resolution for chemical fields
  Real grid_dx = 2.0e-6;  // 2 um
  // Integer coarsening of the chemistry grid in x, y, and z.  Agent
  // mechanics and hashing remain on physical-length scales.
  std::array<Int, 3> chemistry_stride = {1, 1, 1};

  // Spatial hash cell size (should be >= max interaction range)
  Real hash_cell_size = 10.0e-6;  // 10 um

  // MPI decomposition axis: 0 = x (distal flow direction)
  Int mpi_decomp_axis = 0;

  // Ghost layer thickness (should be >= hash_cell_size for correct neighbor queries)
  Real ghost_width = 10.0e-6;  // 10 um default

  // Chemical grid halo width in cells. Stage 2a validates this metadata for
  // slab chemistry; storage remains replicated until the local-storage stage.
  Int grid_halo_width = 1;
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
  Real dx_x() const { return dx_[0]; }
  Real dx_y() const { return dx_[1]; }
  Real dx_z() const { return dx_[2]; }
  Real cell_volume() const { return dx_[0] * dx_[1] * dx_[2]; }

  // Rank-local x-grid ownership, half-open in global cell coordinates.
  Int local_grid_x_begin() const { return local_grid_x_begin_; }
  Int local_grid_x_end() const { return local_grid_x_end_; }
  Int local_grid_nx() const { return local_grid_x_end_ - local_grid_x_begin_; }
  Int grid_halo_width() const { return grid_halo_width_; }
  Int local_grid_storage_nx() const {
    return local_grid_nx() + 2 * grid_halo_width_;
  }

  // Convert between global x indices and this rank's owned-plus-halo storage.
  // Global indices wrap in periodic x; -1 means the global index is outside
  // this rank's configured halo (or an invalid non-periodic halo).
  Int global_to_local_grid_x(Int global_ix) const;
  Int local_to_global_grid_x(Int local_ix) const;

  // Integer partition helper shared by decomposition metadata and tests.
  static std::pair<Int, Int> grid_x_range_for_rank(
      Int global_nx, Int nprocs, Int rank);
  static Int grid_owner_rank_for_cell(
      Int global_nx, Int nprocs, Int global_ix);

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

  // Local slab bounds (this rank's partition along decomp axis)
  Real local_lo_x() const { return local_lo_x_; }
  Real local_hi_x() const { return local_hi_x_; }
  Real ghost_width() const { return ghost_width_; }

  // Check if a position falls within this rank's local domain
  bool is_local(const Vec3& pos) const;

  // Determine which rank owns a position
  Int owner_rank(const Vec3& pos) const;

  // Neighbor ranks (-1 if no neighbor; wraps for periodic axis)
  Int rank_lo() const { return rank_lo_; }
  Int rank_hi() const { return rank_hi_; }

  // True when periodic wrap maps both slab faces to the same neighbor rank
  // (e.g. 2-rank decomposition along a periodic axis).
  bool neighbors_collapsed() const {
    return nprocs_ > 1 && rank_lo_ >= 0 && rank_lo_ == rank_hi_;
  }

  // Spatial hash
  SpatialHash& spatial_hash() { return hash_; }
  const SpatialHash& spatial_hash() const { return hash_; }

  const DomainConfig& config() const { return cfg_; }

 private:
  void decompose();

  DomainConfig cfg_;
  Vec3 lo_{};
  Vec3 hi_{};
  Int nx_ = 0;
  Int ny_ = 0;
  Int nz_ = 0;
  Vec3 dx_ = {2.0e-6, 2.0e-6, 2.0e-6};

  std::array<bool, 3> periodic_{};

  SpatialHash hash_;

  Int rank_   = 0;
  Int nprocs_ = 1;

  // Slab decomposition along x-axis
  Real local_lo_x_ = 0.0;
  Real local_hi_x_ = 0.0;
  Real ghost_width_ = 0.0;
  Int rank_lo_ = -1;  // neighbor rank in -x direction
  Int rank_hi_ = -1;  // neighbor rank in +x direction

  // Rank-local chemical grid metadata along x.
  Int local_grid_x_begin_ = 0;
  Int local_grid_x_end_ = 0;
  Int grid_halo_width_ = 1;
};

}  // namespace gutibm

#endif  // GUTIBM_DOMAIN_H
