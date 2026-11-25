/* -----------------------------------------------------------------------
   GutIBM – Kernel-independent FMM implementation
   ----------------------------------------------------------------------- */

#include "fmm.h"
#include "fmm_kernel.h"
#include "domain.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#ifdef GUTIBM_OPENMP
#include <omp.h>
#endif

namespace gutibm {

namespace {

constexpr Real kSqrt3 = 1.732050807568877293527446341505872367;

Real vec3_norm(const Vec3& v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

}  // namespace

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

bool FMM::has_empty_child(const FMMNode& node) {
  if (node.is_leaf) return false;
  for (int child : node.children) {
    if (child < 0) return true;
  }
  return false;
}

Real FMM::stencil_margin(Real dist) const {
  const Real L = std::max(dist, 1.0e-9);
  if (expansion_order_ >= 3) return 3.0 * fmm_detail::fd_step(3, L);
  if (expansion_order_ == 2) return 2.0 * fmm_detail::fd_step(2, L);
  return fmm_detail::fd_step(1, L);
}

bool FMM::crosses_periodic_branch_cut(const Vec3& delta, Real pad) const {
  const auto& periodic = domain_->config().periodic;
  const Vec3 lo = domain_->lo();
  const Vec3 hi = domain_->hi();
  for (int d = 0; d < 3; ++d) {
    if (!periodic[static_cast<size_t>(d)]) continue;
    const Real span = hi[static_cast<size_t>(d)] - lo[static_cast<size_t>(d)];
    if (std::abs(delta[static_cast<size_t>(d)]) + pad >= 0.5 * span) {
      return true;
    }
  }
  return false;
}

bool FMM::is_source_in_near_region(const Vec3& src,
                                   const Vec3& tgt,
                                   Real near_cutoff,
                                   FmmNearCutoffShape cutoff_shape) const {
  if (near_cutoff <= 0.0) return false;
  if (cutoff_shape == FmmNearCutoffShape::Sphere) {
    return domain_->min_image_dist_sq(src, tgt) <= near_cutoff * near_cutoff;
  }
  Int sx = 0;
  Int sy = 0;
  Int sz = 0;
  Int tx = 0;
  Int ty = 0;
  Int tz = 0;
  domain_->pos_to_grid(src, sx, sy, sz);
  domain_->pos_to_grid(tgt, tx, ty, tz);
  const Int span_x = static_cast<Int>(std::ceil(near_cutoff / domain_->dx_x()));
  const Int span_y = static_cast<Int>(std::ceil(near_cutoff / domain_->dx_y()));
  const Int span_z = static_cast<Int>(std::ceil(near_cutoff / domain_->dx_z()));
  const auto& periodic = domain_->config().periodic;
  auto axis_dist = [](Int a, Int b, Int n, bool per) {
    const Int d = std::abs(a - b);
    return per ? std::min(d, n - d) : d;
  };
  return axis_dist(tx, sx, domain_->nx(), periodic[0]) <= span_x
      && axis_dist(ty, sy, domain_->ny(), periodic[1]) <= span_y
      && axis_dist(tz, sz, domain_->nz(), periodic[2]) <= span_z;
}

bool FMM::pair_completely_inside_near(const Vec3& box_delta,
                                      Real combined_half,
                                      Real near_cutoff,
                                      FmmNearCutoffShape cutoff_shape) const {
  if (near_cutoff <= 0.0) return false;
  if (cutoff_shape == FmmNearCutoffShape::Sphere) {
    return vec3_norm(box_delta) + combined_half * kSqrt3 <= near_cutoff;
  }
  const std::array<Real, 3> dx = {domain_->dx_x(), domain_->dx_y(), domain_->dx_z()};
  const std::array<Int, 3> n = {domain_->nx(), domain_->ny(), domain_->nz()};
  const auto& periodic = domain_->config().periodic;
  for (size_t d = 0; d < 3; ++d) {
    const Int span_d = static_cast<Int>(std::ceil(near_cutoff / dx[d]));
    if (periodic[d] && (2 * span_d + 1 >= n[d])) continue;
    if (std::abs(box_delta[d]) + combined_half + dx[d]
        > static_cast<Real>(span_d) * dx[d]) {
      return false;
    }
  }
  return true;
}

bool FMM::pair_completely_outside_near(const Vec3& box_delta,
                                       Real combined_half,
                                       Real near_cutoff,
                                       FmmNearCutoffShape cutoff_shape) const {
  if (near_cutoff <= 0.0) return true;
  if (cutoff_shape == FmmNearCutoffShape::Sphere) {
    return vec3_norm(box_delta) - combined_half * kSqrt3 > near_cutoff;
  }
  const std::array<Real, 3> dx = {domain_->dx_x(), domain_->dx_y(), domain_->dx_z()};
  const std::array<Int, 3> n = {domain_->nx(), domain_->ny(), domain_->nz()};
  const auto& periodic = domain_->config().periodic;
  for (size_t d = 0; d < 3; ++d) {
    const Int span_d = static_cast<Int>(std::ceil(near_cutoff / dx[d]));
    if (periodic[d] && (2 * span_d + 1 >= n[d])) continue;
    if (std::abs(box_delta[d]) - combined_half
        > static_cast<Real>(span_d + 1) * dx[d]) {
      return true;
    }
  }
  return false;
}

bool FMM::well_separated(const FMMNode& source,
                         const FMMNode& target,
                         Real theta,
                         Real near_cutoff,
                         FmmNearCutoffShape cutoff_shape) const {
  const Vec3 cos_delta = domain_->min_image_delta(source.center_of_source, target.center);
  const Real dist = vec3_norm(cos_delta);
  // Standard dual-tree MAC: (r_src + r_tgt) / dist < theta, where r = half_size.
  // For equal-sized boxes of side L = 2r at theta = 0.5, this requires dist > 2L
  // (classic Greengard-Rokhlin one-box separation), matching traverse_far's
  // 2*half_size / dist < theta single-box criterion.
  const Real combined_half = source.half_size + target.half_size;
  if (dist <= combined_half / std::max(theta, 1e-30)) {
    return false;
  }

  const Vec3 box_delta = domain_->min_image_delta(source.center, target.center);
  // Reject M2L pairs whose central-difference stencil around cos_delta spans a
  // periodic half-domain branch cut where min_image_delta switches sign.
  const Real fd_pad = 2.0 * stencil_margin(dist);
  if (crosses_periodic_branch_cut(cos_delta, fd_pad)) {
    return false;
  }

  if (near_cutoff > 0.0
      && !pair_completely_outside_near(
          box_delta, combined_half, near_cutoff, cutoff_shape)) {
    return false;
  }
  return true;
}

bool FMM::locals_match_query(Real theta,
                             Real near_cutoff,
                             FmmNearCutoffShape cutoff_shape) const {
  if (!locals_ready_) return false;
  if (cutoff_shape != locals_cutoff_shape_) return false;
  const Real theta_tol = 1e-12 * std::max<Real>(1.0, std::abs(theta));
  const Real cutoff_tol = 1e-12 * std::max<Real>(1.0, std::abs(near_cutoff));
  return std::abs(theta - locals_theta_) <= theta_tol
      && std::abs(near_cutoff - locals_near_cutoff_) <= cutoff_tol;
}

void FMM::build(const std::vector<Vec3>& positions,
                const std::vector<Real>& strengths,
                const Domain& domain,
                int expansion_order) {
  nodes_.clear();
  source_positions_ = positions;
  source_strengths_ = strengths;
  domain_ = &domain;
  locals_ready_ = false;
  locals_theta_ = 0.0;
  locals_near_cutoff_ = 0.0;
  locals_cutoff_shape_ = FmmNearCutoffShape::Sphere;
  m2l_interaction_count_ = 0;
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

  build_recursive(positions, strengths, all_indices, center, half_size, 0);
  upward_pass(positions, strengths);
}

int FMM::build_recursive(const std::vector<Vec3>& positions,
                           const std::vector<Real>& strengths,
                           const std::vector<int>& indices,
                           const Vec3& center, Real half_size, int depth) {
  auto node_idx = static_cast<int>(nodes_.size());
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

  std::array<std::vector<int>, 8> child_indices;
  for (int idx : indices) {
    int oct = octant(positions[idx], center);
    child_indices[oct].push_back(idx);
  }

  bool subdivision_progress = false;
  for (const auto& child : child_indices) {
    if (child.size() < indices.size()) {
      subdivision_progress = true;
      break;
    }
  }
  if (!subdivision_progress || depth >= MAX_TREE_DEPTH) {
    // Preserve all sources in an oversized leaf when subdivision cannot
    // separate them or has reached the recursion backstop.
    node.is_leaf = true;
    node.sources = indices;
    return node_idx;
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
                                    child_center, child_half, depth + 1);
    nodes_[node_idx].children[oct] = child_idx;
  }

  return node_idx;
}

void FMM::upward_pass(const std::vector<Vec3>& positions,
                      const std::vector<Real>& strengths) {
  // Post-order traversal: leaves first, then internal nodes.
  for (auto i = static_cast<int>(nodes_.size()) - 1; i >= 0; --i) {
    FMMNode& node = nodes_[i];
    std::ranges::fill(node.multipole, 0.0);

    if (node.is_leaf) {
      for (int idx : node.sources) {
        fmm_detail::add_particle(node.multipole, expansion_order_,
                                 strengths[idx], positions[idx],
                                 node.center_of_source);
      }
    } else {
      for (int c = 0; c < 8; ++c) {
        if (node.children[c] < 0) continue;
        const FMMNode& child = nodes_[node.children[c]];
        fmm_detail::add_shifted_moments(
            node.multipole, child.multipole, expansion_order_,
            child.center_of_source, node.center_of_source);
      }
    }

    // Multipole moments are centered at the weighted source centroid
    // (center_of_source), so the dipole moments vanish identically.
    if (expansion_order_ >= 1 && node.total_source_strength > 0.0) {
      node.multipole[1] = 0.0;
      node.multipole[2] = 0.0;
      node.multipole[3] = 0.0;
    }
  }
}

void FMM::collect_dual_tree_interactions(
    int target_idx,
    int src_idx,
    Real theta,
    Real near_cutoff,
    FmmNearCutoffShape cutoff_shape,
    std::vector<std::vector<int>>& m2l_sources) {
  const FMMNode& source = nodes_[src_idx];
  FMMNode& target = nodes_[target_idx];
  if (source.total_source_strength <= 0.0) return;

  if (near_cutoff > 0.0) {
    const Vec3 box_delta = domain_->min_image_delta(source.center, target.center);
    if (pair_completely_inside_near(
            box_delta, source.half_size + target.half_size,
            near_cutoff, cutoff_shape)) {
      return;
    }
  }

  if (src_idx != target_idx
      && well_separated(source, target, theta, near_cutoff, cutoff_shape)) {
    m2l_sources[static_cast<size_t>(target_idx)].push_back(src_idx);
    ++m2l_interaction_count_;
    return;
  }

  if (target.is_leaf && source.is_leaf) {
    target.near_list.push_back(src_idx);
    return;
  }

  if (target.is_leaf || (!source.is_leaf && source.half_size > target.half_size)) {
    for (int child : source.children) {
      if (child >= 0) {
        collect_dual_tree_interactions(
            target_idx, child, theta, near_cutoff, cutoff_shape, m2l_sources);
      }
    }
    return;
  }

  if (has_empty_child(target)) {
    target.near_list.push_back(src_idx);
  }
  for (int child : target.children) {
    if (child >= 0) {
      collect_dual_tree_interactions(
          child, src_idx, theta, near_cutoff, cutoff_shape, m2l_sources);
    }
  }
}

void FMM::l2l_downward(int node_idx) {
  FMMNode& node = nodes_[node_idx];

  for (int child_idx : node.children) {
    if (child_idx < 0) continue;
    FMMNode& child = nodes_[child_idx];

    const std::vector<Real> shifted = shift_local_expansion(
        node.local, expansion_order_, node.center, child.center);
    for (size_t k = 0; k < child.local.size(); ++k) {
      child.local[k] += shifted[k];
    }

    l2l_downward(child_idx);
  }
}

void FMM::compute_local_expansions(Real theta,
                                   const GreensFunction& gf,
                                   const GreensFunctionParams& avg_params,
                                   Real near_cutoff,
                                   FmmNearCutoffShape cutoff_shape) {
  if (nodes_.empty()) return;

  for (auto& node : nodes_) {
    std::ranges::fill(node.local, 0.0);
    node.near_list.clear();
  }
  m2l_interaction_count_ = 0;

  // Sub-quadratic O(N_nodes) dual-tree interaction walk from (root, root).
  std::vector<std::vector<int>> m2l_sources(nodes_.size());
  collect_dual_tree_interactions(
      0, 0, theta, near_cutoff, cutoff_shape, m2l_sources);

  const auto n_nodes = static_cast<int>(nodes_.size());
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (int tgt_idx = 0; tgt_idx < n_nodes; ++tgt_idx) {
    const auto& sources_for_tgt = m2l_sources[static_cast<size_t>(tgt_idx)];
    if (sources_for_tgt.empty()) continue;
    FMMNode& target = nodes_[tgt_idx];
    for (int src_idx : sources_for_tgt) {
      const FMMNode& source = nodes_[src_idx];
      const Real sep = vec3_norm(
          domain_->min_image_delta(source.center_of_source, target.center));
      // Keep physical source_center so spatially varying advection/peristalsis
      // reads the true producer phase; concentration_bounded wraps displacements
      // via min_image_delta internally.
      const std::vector<Real> contrib = multipole_to_local(
          source.multipole, expansion_order_,
          source.center_of_source, target.center,
          gf, avg_params, sep);
      for (size_t k = 0; k < target.local.size(); ++k) {
        target.local[k] += contrib[k];
      }
    }
  }

  l2l_downward(0);
  locals_theta_ = theta;
  locals_near_cutoff_ = near_cutoff;
  locals_cutoff_shape_ = cutoff_shape;
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

Real FMM::evaluate_local_with_residual(
    const Vec3& target,
    Real theta,
    Real near_cutoff,
    FmmNearCutoffShape cutoff_shape,
    const GreensFunction& gf,
    const GreensFunctionParams& avg_params) const {
  if (!locals_ready_ || nodes_.empty()) return 0.0;

  const int node_idx = find_containing_node(target);
  const FMMNode& node = nodes_[node_idx];
  Real acc = evaluate_local(node.local, expansion_order_, node.center, target);
  for (int src_idx : node.near_list) {
    traverse_far(src_idx, target, theta, near_cutoff, cutoff_shape,
                 gf, avg_params, acc);
  }
  return acc;
}

void FMM::traverse_far(int node_idx,
                       const Vec3& target,
                       Real theta,
                       Real near_cutoff,
                       FmmNearCutoffShape cutoff_shape,
                       const GreensFunction& gf,
                       const GreensFunctionParams& avg_params,
                       Real& accumulator) const {
  const FMMNode& node = nodes_[node_idx];
  if (node.total_source_strength <= 0.0) return;

  const Vec3 box_delta = domain_->min_image_delta(node.center, target);
  if (pair_completely_inside_near(
          box_delta, node.half_size, near_cutoff, cutoff_shape)) {
    return;
  }

  const Vec3 cos_delta = domain_->min_image_delta(node.center_of_source, target);
  const Real dist = vec3_norm(cos_delta);
  const bool straddles_near = near_cutoff > 0.0
      && !pair_completely_outside_near(
          box_delta, node.half_size, near_cutoff, cutoff_shape);
  const Real pad = stencil_margin(dist);
  const bool straddles_branch = crosses_periodic_branch_cut(cos_delta, pad);

  const Real cell_size = 2.0 * node.half_size;
  if (!node.is_leaf) {
    const bool fails_mac = (cell_size / std::max(dist, 1e-30)) >= theta;
    if (fails_mac || straddles_near || straddles_branch) {
      for (int child : node.children) {
        if (child >= 0) {
          traverse_far(child, target, theta, near_cutoff, cutoff_shape,
                       gf, avg_params, accumulator);
        }
      }
      return;
    }
  } else if (straddles_near || straddles_branch) {
    // Leaf straddles the near cutoff or a periodic half-domain branch cut:
    // evaluate each individual far source directly to avoid cutoff
    // misclassification or finite-difference differentiation across a jump.
    GreensFunctionParams unit = avg_params;
    for (int src_id : node.sources) {
      const Vec3& src_pos = source_positions_[static_cast<size_t>(src_id)];
      if (is_source_in_near_region(src_pos, target, near_cutoff, cutoff_shape)) {
        continue;
      }
      unit.source_rate = source_strengths_[static_cast<size_t>(src_id)];
      if (unit.source_rate <= 0.0) continue;
      accumulator += gf.concentration_bounded(src_pos, target, unit);
    }
    return;
  }

  // Strictly far, well-separated, non-branch-cut-straddling node: evaluate
  // multipole expansion at physical center_of_source with periodic separation dist.
  accumulator += evaluate_multipole_at_source(
      gf, node.multipole, expansion_order_,
      node.center_of_source, target, avg_params, dist);
}

Real FMM::evaluate_total_field(const Vec3& target,
                               Real theta,
                               const GreensFunction& gf,
                               const GreensFunctionParams& avg_params) const {
  if (nodes_.empty()) return 0.0;

  if (locals_match_query(theta, 0.0, FmmNearCutoffShape::Sphere)) {
    return std::max(
        evaluate_local_with_residual(
            target, theta, 0.0, FmmNearCutoffShape::Sphere, gf, avg_params),
        0.0);
  }

  Real result = 0.0;
  traverse_far(0, target, theta, 0.0, FmmNearCutoffShape::Sphere,
               gf, avg_params, result);
  return std::max(result, 0.0);
}

Real FMM::evaluate_far_field(const Vec3& target,
                             Real theta,
                             Real near_cutoff,
                             const GreensFunction& gf,
                             const GreensFunctionParams& avg_params,
                             FmmNearCutoffShape cutoff_shape) const {
  if (nodes_.empty()) return 0.0;

  if (locals_match_query(theta, near_cutoff, cutoff_shape)) {
    return std::max(
        evaluate_local_with_residual(
            target, theta, near_cutoff, cutoff_shape, gf, avg_params),
        0.0);
  }

  Real far = 0.0;
  traverse_far(0, target, theta, near_cutoff, cutoff_shape,
               gf, avg_params, far);
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

  Real far = evaluate_far_field(
      target, theta, near_cutoff, gf, avg_params, FmmNearCutoffShape::Sphere);
  return std::max(near + far, 0.0);
}

}  // namespace gutibm
