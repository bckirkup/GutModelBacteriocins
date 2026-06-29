/* -----------------------------------------------------------------------
   GutIBM – Barnes-Hut Octree implementation

   Provides O(N log N) far-field aggregation for QSSA Green's function
   superposition.  Distant source clusters are replaced by a single
   monopole at the source-weighted centroid.
   ----------------------------------------------------------------------- */

#include "octree.h"
#include "domain.h"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace gutibm {

// ── Octant assignment ────────────────────────────────────────────────────

int Octree::octant(const Vec3& pos, const Vec3& center) {
  int oct = 0;
  if (pos[0] >= center[0]) oct |= 1;
  if (pos[1] >= center[1]) oct |= 2;
  if (pos[2] >= center[2]) oct |= 4;
  return oct;
}

// ── Build ────────────────────────────────────────────────────────────────

void Octree::build(const std::vector<Vec3>& positions,
                   const std::vector<Real>& strengths,
                   const Domain& domain) {
  nodes_.clear();
  domain_ = &domain;

  if (positions.empty()) return;

  // Bounding box from domain (not positions — keeps tree aligned)
  Vec3 lo = domain.lo();
  Vec3 hi = domain.hi();

  Vec3 center;
  Real half_size = 0.0;
  for (int d = 0; d < 3; ++d) {
    center[d] = 0.5 * (lo[d] + hi[d]);
    Real hs = 0.5 * (hi[d] - lo[d]);
    if (hs > half_size) half_size = hs;
  }
  // Pad slightly so all sources are strictly inside
  half_size *= 1.001;

  std::vector<int> all_indices(positions.size());
  std::iota(all_indices.begin(), all_indices.end(), 0);

  build_recursive(positions, strengths, all_indices, center, half_size);
}

int Octree::build_recursive(const std::vector<Vec3>& positions,
                            const std::vector<Real>& strengths,
                            const std::vector<int>& indices,
                            const Vec3& center, Real half_size) {
  auto node_idx = static_cast<int>(nodes_.size());
  nodes_.emplace_back();
  OctreeNode& node = nodes_[node_idx];

  node.center = center;
  node.half_size = half_size;
  for (int i = 0; i < 8; ++i) node.children[i] = -1;

  // Compute monopole: total strength and weighted centroid
  Real total_s = 0.0;
  Vec3 weighted_pos = {0.0, 0.0, 0.0};
  for (int idx : indices) {
    Real s = strengths[idx];
    total_s += s;
    for (int d = 0; d < 3; ++d)
      weighted_pos[d] += s * positions[idx][d];
  }
  node.total_source_strength = total_s;
  if (total_s > 0.0) {
    for (int d = 0; d < 3; ++d)
      node.center_of_source[d] = weighted_pos[d] / total_s;
  } else {
    node.center_of_source = center;
  }

  // Leaf condition
  if (static_cast<int>(indices.size()) <= MAX_LEAF_SOURCES) {
    node.is_leaf = true;
    node.sources = indices;
    return node_idx;
  }

  node.is_leaf = false;

  // Partition into 8 octants
  std::array<std::vector<int>, 8> child_indices;
  for (int idx : indices) {
    int oct = octant(positions[idx], center);
    child_indices[oct].push_back(idx);
  }

  Real child_half = half_size * 0.5;
  for (int oct = 0; oct < 8; ++oct) {
    if (child_indices[oct].empty()) continue;

    Vec3 child_center;
    child_center[0] = center[0] + ((oct & 1) ? child_half : -child_half);
    child_center[1] = center[1] + ((oct & 2) ? child_half : -child_half);
    child_center[2] = center[2] + ((oct & 4) ? child_half : -child_half);

    // Recursive call may invalidate node reference (vector reallocation)
    int child_idx = build_recursive(positions, strengths,
                                    child_indices[oct],
                                    child_center, child_half);
    nodes_[node_idx].children[oct] = child_idx;
  }

  return node_idx;
}

// ── Far-field traversal ──────────────────────────────────────────────────

void Octree::traverse_far(int node_idx,
                          const Vec3& target,
                          Real theta,
                          Real near_cutoff2,
                          const GreensFunction& gf,
                          const GreensFunctionParams& avg_params,
                          Real& accumulator) const {
  const OctreeNode& node = nodes_[node_idx];

  if (node.total_source_strength <= 0.0) return;

  // Distance from target to node's source centroid
  Vec3 delta = domain_->min_image_delta(node.center_of_source, target);
  Real dist2 = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];

  // If this entire node is within near_cutoff, skip it —
  // those sources will be evaluated exactly by the caller.
  // Check: if centroid + diagonal is within cutoff, all sources are near.
  Real diag = node.half_size * std::sqrt(3.0);
  Real dist = std::sqrt(dist2);
  if (dist + diag < std::sqrt(near_cutoff2)) return;

  // Barnes-Hut opening criterion: cell_size / distance < theta
  Real cell_size = 2.0 * node.half_size;
  if (!node.is_leaf && (cell_size / std::max(dist, 1e-30)) >= theta) {
    // Open this node — recurse into children
    for (int c = 0; c < 8; ++c) {
      if (node.children[c] >= 0) {
        traverse_far(node.children[c], target, theta, near_cutoff2,
                     gf, avg_params, accumulator);
      }
    }
    return;
  }

  // Accept monopole approximation: treat the cluster as a single source
  // at center_of_source with total strength
  if (dist2 <= near_cutoff2) {
    // This leaf/cell overlaps with the near-field region.
    // For leaf nodes, skip individual sources that are within cutoff
    // (the caller handles them exactly).
    if (node.is_leaf) {
      // Nothing to do — leaf sources within cutoff handled by caller
      return;
    }
    // For internal nodes that pass opening criterion but overlap
    // near-field: recurse deeper to separate near/far sources
    for (int c = 0; c < 8; ++c) {
      if (node.children[c] >= 0) {
        traverse_far(node.children[c], target, theta, near_cutoff2,
                     gf, avg_params, accumulator);
      }
    }
    return;
  }

  // Far-field monopole: use Green's function with aggregate source
  GreensFunctionParams mono = avg_params;
  mono.source_rate = node.total_source_strength;
  accumulator += gf.concentration_bounded(node.center_of_source, target, mono);
}

Real Octree::evaluate_far_field(const Vec3& target,
                                Real theta,
                                Real near_cutoff,
                                const GreensFunction& gf,
                                const GreensFunctionParams& avg_params) const {
  if (nodes_.empty()) return 0.0;

  Real result = 0.0;
  Real near_cutoff2 = near_cutoff * near_cutoff;
  traverse_far(0, target, theta, near_cutoff2, gf, avg_params, result);
  return std::max(result, 0.0);
}

Real Octree::evaluate_field(const Vec3& target,
                            Real theta,
                            Real near_cutoff,
                            const GreensFunction& gf,
                            const std::vector<Vec3>& positions,
                            const std::vector<GreensFunctionParams>& params,
                            const GreensFunctionParams& avg_params) const {
  if (nodes_.empty()) return 0.0;

  Real near_cutoff2 = near_cutoff * near_cutoff;

  // Near-field: exact evaluation for sources within cutoff
  Real near = 0.0;
  for (size_t s = 0; s < positions.size(); ++s) {
    Real d2 = domain_->min_image_dist_sq(positions[s], target);
    if (d2 <= near_cutoff2) {
      near += gf.concentration_bounded(positions[s], target, params[s]);
    }
  }

  // Far-field: Barnes-Hut monopole approximation
  Real far = evaluate_far_field(target, theta, near_cutoff, gf, avg_params);

  return std::max(near + far, 0.0);
}

}  // namespace gutibm
