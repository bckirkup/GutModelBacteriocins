/* -----------------------------------------------------------------------
   GutIBM – Barnes-Hut Octree for O(N log N) far-field evaluation

   Accelerates QSSA Green's function superposition by aggregating
   distant source contributions via monopole approximation.  Near-field
   sources (within the existing cutoff radius) are still evaluated
   exactly with the advection-diffusion kernel.

   Opening-angle parameter theta controls the accuracy/speed trade-off:
     theta → 0 : exact (all interactions direct)
     theta ≈ 0.5: good balance for bacteriocin fields
     theta → 1 : fast but lower accuracy
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_OCTREE_H
#define GUTIBM_OCTREE_H

#include "types.h"
#include "greens_function.h"
#include <array>
#include <vector>

namespace gutibm {

class Domain;
class AdvectionField;

struct OctreeNode {
  Vec3 center;                     // geometric center of this cell
  Real half_size;                  // half the side length
  Real total_source_strength;      // sum of source rates in subtree
  Vec3 center_of_source;           // source-weighted centroid
  std::array<int, 8> children;     // child indices (-1 = no child)
  std::vector<int> sources;        // leaf source indices
  bool is_leaf;
};

class Octree {
 public:
  Octree() = default;

  // Build octree from source positions and strengths.
  // The domain is used only to determine the bounding box.
  void build(const std::vector<Vec3>& positions,
             const std::vector<Real>& strengths,
             const Domain& domain);

  // Evaluate aggregate field at a target point using Barnes-Hut
  // traversal.  Sources closer than near_cutoff are skipped (they
  // will be evaluated exactly by the caller).
  // gf + adv + avg_params are used for the far-field monopole kernel.
  Real evaluate_far_field(const Vec3& target,
                          Real theta,
                          Real near_cutoff,
                          const GreensFunction& gf,
                          const GreensFunctionParams& avg_params) const;

  // Evaluate the full field (near + far) at a single target.
  // Near-field uses per-source params; far-field uses monopole.
  Real evaluate_field(const Vec3& target,
                      Real theta,
                      Real near_cutoff,
                      const GreensFunction& gf,
                      const std::vector<Vec3>& positions,
                      const std::vector<GreensFunctionParams>& params,
                      const GreensFunctionParams& avg_params) const;

  int num_nodes() const { return static_cast<int>(nodes_.size()); }
  bool empty() const { return nodes_.empty(); }
  const OctreeNode& node(int i) const { return nodes_[i]; }

  static constexpr int MAX_LEAF_SOURCES = 8;

 private:
  // Recursive build helper — returns node index
  int build_recursive(const std::vector<Vec3>& positions,
                      const std::vector<Real>& strengths,
                      const std::vector<int>& indices,
                      const Vec3& center, Real half_size);

  // Recursive far-field traversal
  void traverse_far(int node_idx,
                    const Vec3& target,
                    Real theta,
                    Real near_cutoff2,
                    const GreensFunction& gf,
                    const GreensFunctionParams& avg_params,
                    Real& accumulator) const;

  // Which octant does a position fall into?
  static int octant(const Vec3& pos, const Vec3& center);

  std::vector<OctreeNode> nodes_;
  const Domain* domain_ = nullptr;
};

}  // namespace gutibm

#endif  // GUTIBM_OCTREE_H
