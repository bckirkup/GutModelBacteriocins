/* -----------------------------------------------------------------------
   GutIBM – Kernel-independent FMM implementation
   ----------------------------------------------------------------------- */

#include "fmm.h"
#include "fmm_kernel.h"
#include "domain.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace gutibm {

int FMM::num_coefficients(int order) {
  return fmm_detail::num_coefficients(order);
}

int FMM::octant(const Vec3& pos, const Vec3& center) {
  int oct = 0;
  if (pos[0] >= center[0]) oct |= 1;
  if (pos[1] >= center[1]) oct |= 2;
  if (pos[2] >= center[2]) oct |= 4;
  return oct;
}

bool FMM::well_separated(const FMMNode& source,
                         const FMMNode& target,
                         Real theta) {
  Vec3 delta = {source.center_of_source[0] - target.center[0],
                source.center_of_source[1] - target.center[1],
                source.center_of_source[2] - target.center[2]};
  Real dist = std::sqrt(delta[0] * delta[0] + delta[1] * delta[1]
                        + delta[2] * delta[2]);
  Real src_size = 2.0 * source.half_size;
  Real tgt_size = 2.0 * target.half_size;
  return dist > (src_size + tgt_size) / std::max(theta, 1e-30);
}

void FMM::build(const std::vector<Vec3>& positions,
                const std::vector<Real>& strengths,
                const Domain& domain,
                int expansion_order) {
  nodes_.clear();
  domain_ = &domain;
  locals_ready_ = false;
  expansion_order_ = std::max(1, std::min(expansion_order, MAX_EXPANSION_ORDER));

  if (positions.empty()) return;

  Vec3 lo = domain.lo();
  Vec3 hi = domain.hi();

  Vec3 center;
  Real half_size = 0.0;
  for (int d = 0; d < 3; ++d) {
    center[d] = 0.5 * (lo[d] + hi[d]);
    Real hs = 0.5 * (hi[d] - lo[d]);
    if (hs > half_size) half_size = hs;
  }
  half_size *= 1.001;

  std::vector<int> all_indices(positions.size());
  std::iota(all_indices.begin(), all_indices.end(), 0);

  build_recursive(positions, strengths, all_indices, center, half_size);
  upward_pass(positions, strengths);
}

int FMM::build_recursive(const std::vector<Vec3>& positions,
                           const std::vector<Real>& strengths,
                           const std::vector<int>& indices,
                           const Vec3& center, Real half_size) {
  int node_idx = static_cast<int>(nodes_.size());
  nodes_.emplace_back();
  FMMNode& node = nodes_[node_idx];

  node.center = center;
  node.half_size = half_size;
  for (int i = 0; i < 8; ++i) node.children[i] = -1;

  const int ncoeff = num_coefficients(expansion_order_);
  node.multipole.assign(ncoeff, 0.0);
  node.local.assign(ncoeff, 0.0);

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

  if (static_cast<int>(indices.size()) <= MAX_LEAF_SOURCES) {
    node.is_leaf = true;
    node.sources = indices;
    return node_idx;
  }

  node.is_leaf = false;

  std::vector<int> child_indices[8];
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

    int child_idx = build_recursive(positions, strengths,
                                    child_indices[oct],
                                    child_center, child_half);
    nodes_[node_idx].children[oct] = child_idx;
  }

  return node_idx;
}

void FMM::upward_pass(const std::vector<Vec3>& positions,
                      const std::vector<Real>& strengths) {
  // Post-order traversal: leaves first, then internal nodes.
  for (int i = static_cast<int>(nodes_.size()) - 1; i >= 0; --i) {
    FMMNode& node = nodes_[i];
    std::fill(node.multipole.begin(), node.multipole.end(), 0.0);

    if (node.is_leaf) {
      for (int idx : node.sources) {
        fmm_detail::add_particle(node.multipole, expansion_order_,
                                 strengths[idx], positions[idx],
                                 node.center_of_source);
      }
      continue;
    }

    for (int c = 0; c < 8; ++c) {
      if (node.children[c] < 0) continue;
      const FMMNode& child = nodes_[node.children[c]];
      fmm_detail::add_shifted_moments(
          node.multipole, child.multipole, expansion_order_,
          child.center_of_source, node.center_of_source);
    }
  }
}

void FMM::m2l_at_node(int node_idx,
                      Real theta,
                      const GreensFunction& gf,
                      const GreensFunctionParams& avg_params) {
  FMMNode& target = nodes_[node_idx];

  for (int src_idx = 0; src_idx < static_cast<int>(nodes_.size()); ++src_idx) {
    if (src_idx == node_idx) continue;
    const FMMNode& source = nodes_[src_idx];
    if (source.total_source_strength <= 0.0) continue;
    if (!well_separated(source, target, theta)) continue;

    std::vector<Real> contrib = multipole_to_local(
        source.multipole, expansion_order_,
        source.center_of_source, target.center,
        gf, avg_params);

    for (size_t k = 0; k < target.local.size(); ++k)
      target.local[k] += contrib[k];
  }
}

void FMM::l2l_downward(int node_idx) {
  FMMNode& node = nodes_[node_idx];

  for (int c = 0; c < 8; ++c) {
    if (node.children[c] < 0) continue;
    FMMNode& child = nodes_[node.children[c]];

    std::vector<Real> shifted = shift_local_expansion(
        node.local, expansion_order_, node.center, child.center);
    for (size_t k = 0; k < child.local.size(); ++k)
      child.local[k] += shifted[k];

    l2l_downward(node.children[c]);
  }
}

void FMM::compute_local_expansions(Real theta,
                                   const GreensFunction& gf,
                                   const GreensFunctionParams& avg_params) {
  if (nodes_.empty()) return;

  for (auto& node : nodes_)
    std::fill(node.local.begin(), node.local.end(), 0.0);

  // M2L at every node, then L2L from root downward.
  for (int i = 0; i < static_cast<int>(nodes_.size()); ++i)
    m2l_at_node(i, theta, gf, avg_params);

  l2l_downward(0);
  locals_ready_ = true;
}

int FMM::find_containing_node(const Vec3& target) const {
  int node_idx = 0;
  while (!nodes_[node_idx].is_leaf) {
    const FMMNode& node = nodes_[node_idx];
    int child_oct = octant(target, node.center);
    if (node.children[child_oct] < 0) break;
    node_idx = node.children[child_oct];
  }
  return node_idx;
}

Real FMM::evaluate_local_at(const Vec3& target,
                            const GreensFunction& gf,
                            const GreensFunctionParams& avg_params) const {
  if (!locals_ready_ || nodes_.empty()) return 0.0;

  int leaf = find_containing_node(target);
  return evaluate_local(nodes_[leaf].local, expansion_order_,
                        nodes_[leaf].center, target);
}

void FMM::traverse_far(int node_idx,
                       const Vec3& target,
                       Real theta,
                       Real near_cutoff2,
                       const GreensFunction& gf,
                       const GreensFunctionParams& avg_params,
                       Real& accumulator) const {
  const FMMNode& node = nodes_[node_idx];

  if (node.total_source_strength <= 0.0) return;

  Vec3 delta = domain_->min_image_delta(node.center_of_source, target);
  Real dist2 = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];

  Real diag = node.half_size * std::sqrt(3.0);
  Real dist = std::sqrt(dist2);
  if (dist + diag < std::sqrt(near_cutoff2)) return;

  Real cell_size = 2.0 * node.half_size;
  if (!node.is_leaf && (cell_size / std::max(dist, 1e-30)) >= theta) {
    for (int c = 0; c < 8; ++c) {
      if (node.children[c] >= 0) {
        traverse_far(node.children[c], target, theta, near_cutoff2,
                     gf, avg_params, accumulator);
      }
    }
    return;
  }

  if (dist2 <= near_cutoff2) {
    if (node.is_leaf) return;
    for (int c = 0; c < 8; ++c) {
      if (node.children[c] >= 0) {
        traverse_far(node.children[c], target, theta, near_cutoff2,
                     gf, avg_params, accumulator);
      }
    }
    return;
  }

  // Higher-order multipole evaluation at target.
  KernelTaylorCoeffs kernel = kernel_taylor_at_source(
      gf, node.center_of_source, target, avg_params, expansion_order_);
  accumulator += evaluate_multipole(node.multipole, expansion_order_, kernel);
}

Real FMM::evaluate_total_field(const Vec3& target,
                               Real theta,
                               const GreensFunction& gf,
                               const GreensFunctionParams& avg_params) const {
  if (nodes_.empty()) return 0.0;

  if (locals_ready_)
    return std::max(evaluate_local_at(target, gf, avg_params), 0.0);

  Real result = 0.0;
  const Real huge_cutoff2 = 1.0e20;
  traverse_far(0, target, theta, huge_cutoff2, gf, avg_params, result);
  return std::max(result, 0.0);
}

Real FMM::evaluate_far_field(const Vec3& target,
                             Real theta,
                             Real near_cutoff,
                             const GreensFunction& gf,
                             const GreensFunctionParams& avg_params) const {
  if (nodes_.empty()) return 0.0;

  Real far = 0.0;
  Real near_cutoff2 = near_cutoff * near_cutoff;
  traverse_far(0, target, theta, near_cutoff2, gf, avg_params, far);
  return std::max(far, 0.0);
}

Real FMM::evaluate_field(const Vec3& target,
                         Real theta,
                         Real near_cutoff,
                         const GreensFunction& gf,
                         const std::vector<Vec3>& positions,
                         const std::vector<GreensFunctionParams>& params,
                         const GreensFunctionParams& avg_params) const {
  if (nodes_.empty()) return 0.0;

  Real near_cutoff2 = near_cutoff * near_cutoff;

  Real near = 0.0;
  for (size_t s = 0; s < positions.size(); ++s) {
    Real d2 = domain_->min_image_dist_sq(positions[s], target);
    if (d2 <= near_cutoff2) {
      near += gf.concentration_bounded(positions[s], target, params[s]);
    }
  }

  Real far = evaluate_far_field(target, theta, near_cutoff, gf, avg_params);
  return std::max(near + far, 0.0);
}

}  // namespace gutibm
