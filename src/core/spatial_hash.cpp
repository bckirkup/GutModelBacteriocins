/* -----------------------------------------------------------------------
   GutIBM – Spatial hashing implementation
   ----------------------------------------------------------------------- */

#include "spatial_hash.h"
#include <cmath>
#include <algorithm>

namespace gutibm {

void SpatialHash::init(Vec3 domain_lo, Vec3 domain_hi, Real cell_size) {
  lo_        = domain_lo;
  hi_        = domain_hi;
  cell_size_ = cell_size;
  cells_.clear();
}

void SpatialHash::clear() {
  cells_.clear();
}

void SpatialHash::insert(Int agent_idx, const Vec3& pos) {
  CellKey key = pos_to_cell(pos);
  cells_[key].push_back(agent_idx);
}

SpatialHash::CellKey SpatialHash::pos_to_cell(const Vec3& pos) const {
  return {
    static_cast<Int>(std::floor((pos[0] - lo_[0]) / cell_size_)),
    static_cast<Int>(std::floor((pos[1] - lo_[1]) / cell_size_)),
    static_cast<Int>(std::floor((pos[2] - lo_[2]) / cell_size_))
  };
}

std::vector<Int> SpatialHash::query_radius(const Vec3& pos, Real radius) const {
  std::vector<Int> result;

  auto cells_span = static_cast<Int>(std::ceil(radius / cell_size_));
  CellKey center = pos_to_cell(pos);

  for (Int dz = -cells_span; dz <= cells_span; ++dz) {
    for (Int dy = -cells_span; dy <= cells_span; ++dy) {
      for (Int dx = -cells_span; dx <= cells_span; ++dx) {
        CellKey key{center.ix + dx, center.iy + dy, center.iz + dz};
        auto it = cells_.find(key);
        if (it == cells_.end()) continue;
        for (Int idx : it->second) {
          result.push_back(idx);
        }
      }
    }
  }
  return result;
}

std::vector<Int> SpatialHash::query_neighbors(const Vec3& pos) const {
  std::vector<Int> result;
  CellKey center = pos_to_cell(pos);

  for (Int dz = -1; dz <= 1; ++dz) {
    for (Int dy = -1; dy <= 1; ++dy) {
      for (Int dx = -1; dx <= 1; ++dx) {
        CellKey key{center.ix + dx, center.iy + dy, center.iz + dz};
        auto it = cells_.find(key);
        if (it == cells_.end()) continue;
        for (Int idx : it->second) {
          result.push_back(idx);
        }
      }
    }
  }
  return result;
}

}  // namespace gutibm
