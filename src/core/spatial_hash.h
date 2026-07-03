/* -----------------------------------------------------------------------
   GutIBM – Spatial hashing for O(N) neighbor lookups
   Replaces naive O(N^2) pairwise search with cell-list hashing.
   Critical for scaling to 10^6–10^7 agents.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_SPATIAL_HASH_H
#define GUTIBM_SPATIAL_HASH_H

#include "types.h"
#include <cstddef>
#include <vector>
#include <unordered_map>

namespace gutibm {

class SpatialHash {
 public:
  SpatialHash() = default;

  void init(Vec3 domain_lo, Vec3 domain_hi, Real cell_size);
  void clear();
  void insert(Int agent_idx, const Vec3& pos);

  // Return indices of agents within `radius` of `pos`
  std::vector<Int> query_radius(const Vec3& pos, Real radius) const;

  // Return indices of agents in the same cell and 26 neighbors
  std::vector<Int> query_neighbors(const Vec3& pos) const;

  Real cell_size() const { return cell_size_; }

 private:
  struct CellKey {
    Int ix;
    Int iy;
    Int iz;
    bool operator==(const CellKey& o) const = default;
  };

  struct CellKeyHash {
    size_t operator()(const CellKey& k) const {
      // FNV-1a inspired combine
      size_t h = 14695981039346656037ULL;
      h ^= static_cast<size_t>(k.ix); h *= 1099511628211ULL;
      h ^= static_cast<size_t>(k.iy); h *= 1099511628211ULL;
      h ^= static_cast<size_t>(k.iz); h *= 1099511628211ULL;
      return h;
    }
  };

  CellKey pos_to_cell(const Vec3& pos) const;
  void append_cell_agents(CellKey key, std::vector<Int>& result) const;

  Real cell_size_ = 10.0e-6;  // default 10 um
  Vec3 lo_{};
  Vec3 hi_{};
  std::unordered_map<CellKey, std::vector<Int>, CellKeyHash> cells_;
};

}  // namespace gutibm

#endif  // GUTIBM_SPATIAL_HASH_H
