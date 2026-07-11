/* -----------------------------------------------------------------------
   GutIBM – Kernel-independent Fast Multipole Method (FMM)

   Extends the Barnes-Hut octree (#9) with higher-order Cartesian
   multipole expansions and local-to-local translations for O(N+M)
   far-field evaluation of the advection-diffusion Green's function.

   Expansion order p controls accuracy: error ~ theta^p for well-separated
   clusters.  Order 1 recovers monopole Barnes-Hut; order 2 adds dipole
   and quadrupole moments; order 3 adds octupole terms.

   The kernel is treated as a black box via source-position Taylor
   coefficients of concentration_bounded(), so the method applies to the
   exponentially-weighted advection-diffusion Green's function without
   assuming a 1/r Laplace kernel.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FMM_H
#define GUTIBM_FMM_H

#include "types.h"
#include "greens_function.h"
#include <array>
#include <vector>

namespace gutibm {

class Domain;

struct FMMNode {
  Vec3 center;
  Real half_size;
  Real total_source_strength;
  Vec3 center_of_source;
  std::array<int, 8> children;
  std::vector<int> sources;
  bool is_leaf;

  // Multipole (P2M + M2M) and local (M2L + L2L) Taylor coefficients.
  // Indexed by multi-index (ix, iy, iz) with ix+iy+iz <= expansion_order.
  std::vector<Real> multipole;
  std::vector<Real> local;
};

class FMM {
 public:
  FMM() = default;

  // Build octree and run upward P2M/M2M pass.
  void build(const std::vector<Vec3>& positions,
             const std::vector<Real>& strengths,
             const Domain& domain,
             int expansion_order);

  // Tree-walk M2L at every node, then L2L downward to populate local expansions.
  // Call once per solve before batch far-field evaluation.
  void compute_local_expansions(Real theta,
                                const GreensFunction& gf,
                                const GreensFunctionParams& avg_params);

  // Total field from all sources (near + far) via multipole expansions.
  Real evaluate_total_field(const Vec3& target,
                            Real theta,
                            const GreensFunction& gf,
                            const GreensFunctionParams& avg_params) const;

  // Far-field only: total minus exact near-source contributions within cutoff.
  Real evaluate_far_field(const Vec3& target,
                          Real theta,
                          Real near_cutoff,
                          const GreensFunction& gf,
                          const GreensFunctionParams& avg_params) const;

  // Near + far field at a single target (for unit tests).
  Real evaluate_field(const Vec3& target,
                      Real theta,
                      Real near_cutoff,
                      const GreensFunction& gf,
                      const std::vector<Vec3>& positions,
                      const std::vector<GreensFunctionParams>& params,
                      const GreensFunctionParams& avg_params) const;

  int num_nodes() const { return static_cast<int>(nodes_.size()); }
  bool locals_ready() const { return locals_ready_; }
  bool empty() const { return nodes_.empty(); }
  const FMMNode& node(int i) const { return nodes_[i]; }
  int expansion_order() const { return expansion_order_; }

  static constexpr int MAX_LEAF_SOURCES = 8;
  static constexpr int MAX_EXPANSION_ORDER = 3;

  // Number of Cartesian Taylor coefficients for total degree <= order.
  static int num_coefficients(int order);

 private:
  int build_recursive(const std::vector<Vec3>& positions,
                      const std::vector<Real>& strengths,
                      const std::vector<int>& indices,
                      const Vec3& center, Real half_size);

  void upward_pass(const std::vector<Vec3>& positions,
                   const std::vector<Real>& strengths);

  void m2l_traverse_into_target(int target_idx,
                                int src_idx,
                                Real theta,
                                const GreensFunction& gf,
                                const GreensFunctionParams& avg_params);

  void m2l_visit_children(int target_idx,
                          int src_idx,
                          Real theta,
                          const GreensFunction& gf,
                          const GreensFunctionParams& avg_params);

  void l2l_downward(int node_idx);

  void traverse_far(int node_idx,
                    const Vec3& target,
                    Real theta,
                    Real near_cutoff2,
                    const GreensFunction& gf,
                    const GreensFunctionParams& avg_params,
                    Real& accumulator) const;

  Real evaluate_local_at(const Vec3& target,
                         const GreensFunction& gf,
                         const GreensFunctionParams& avg_params) const;

  int find_containing_node(const Vec3& target) const;

  static int octant(const Vec3& pos, const Vec3& center);
  static bool well_separated(const FMMNode& source,
                             const FMMNode& target,
                             Real theta);

  std::vector<FMMNode> nodes_;
  const Domain* domain_ = nullptr;
  int expansion_order_ = 1;
  bool locals_ready_ = false;
};

}  // namespace gutibm

#endif  // GUTIBM_FMM_H
