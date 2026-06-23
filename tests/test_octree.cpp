/* -----------------------------------------------------------------------
   GutIBM – Octree / Barnes-Hut tests
   Verify tree construction and field evaluation accuracy vs exact method.
   ----------------------------------------------------------------------- */

#include "octree.h"
#include "greens_function.h"
#include "qssa_solver.h"
#include "domain.h"
#include "advection.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <random>

using namespace gutibm;

// Shared test domain setup
static void make_test_domain(Domain& domain, AdvectionField& adv,
                             bool with_flow = false) {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {1e-3, 1e-3, 100e-6};
  dcfg.grid_dx = 5e-6;
  domain.init(dcfg);

  AdvectionConfig acfg;
  if (with_flow) {
    acfg.radial_turnover = 5400.0;
    acfg.distal_transit_time = 43200.0;
  } else {
    acfg.radial_turnover = 1e20;
    acfg.distal_transit_time = 1e20;
  }
  acfg.mucus_thickness = 100e-6;
  acfg.distal_length = 1e-3;
  adv.init(acfg, domain);
}

void test_empty_octree() {
  Domain domain;
  AdvectionField adv;
  make_test_domain(domain, adv);

  Octree tree;
  std::vector<Vec3> positions;
  std::vector<Real> strengths;
  tree.build(positions, strengths, domain);

  assert(tree.empty());
  assert(tree.num_nodes() == 0);

  GreensFunction gf;
  gf.init(domain, adv);

  GreensFunctionParams avg;
  avg.diff_coeff = 4e-11;
  avg.source_rate = 0;
  avg.pI = 7.0;
  avg.retardation = 1.0;

  Vec3 target = {500e-6, 500e-6, 50e-6};
  Real val = tree.evaluate_far_field(target, 0.5, 200e-6, gf, avg);
  assert(val == 0.0);

  std::cout << "  test_empty_octree: PASSED\n";
}

void test_single_source() {
  Domain domain;
  AdvectionField adv;
  make_test_domain(domain, adv);

  GreensFunction gf;
  gf.init(domain, adv);

  Vec3 src = {500e-6, 500e-6, 50e-6};
  Real strength = 1e-18;
  std::vector<Vec3> positions = {src};
  std::vector<Real> strengths = {strength};

  Octree tree;
  tree.build(positions, strengths, domain);

  assert(!tree.empty());
  assert(tree.num_nodes() >= 1);

  // The single source should be in the root leaf
  const OctreeNode& root = tree.node(0);
  assert(root.is_leaf);
  assert(root.sources.size() == 1);
  assert(std::abs(root.total_source_strength - strength) < 1e-30);

  std::cout << "  test_single_source: PASSED\n";
}

void test_octree_construction() {
  Domain domain;
  AdvectionField adv;
  make_test_domain(domain, adv);

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> dist_xy(0.0, 1e-3);
  std::uniform_real_distribution<double> dist_z(0.0, 100e-6);

  int N = 100;
  std::vector<Vec3> positions(N);
  std::vector<Real> strengths(N, 1e-18);
  for (int i = 0; i < N; ++i) {
    positions[i] = {dist_xy(rng), dist_xy(rng), dist_z(rng)};
  }

  Octree tree;
  tree.build(positions, strengths, domain);

  assert(!tree.empty());

  // Root node should contain total strength
  const OctreeNode& root = tree.node(0);
  Real expected_total = N * 1e-18;
  Real actual_total = root.total_source_strength;
  Real err = std::abs(actual_total - expected_total) / expected_total;
  assert(err < 1e-10);

  std::cout << "  test_octree_construction: PASSED (nodes="
            << tree.num_nodes() << ")\n";
}

void test_accuracy_vs_exact() {
  // Compare Barnes-Hut field evaluation against exact brute-force
  // at several target points.  With theta=0.3, error should be small.
  Domain domain;
  AdvectionField adv;
  make_test_domain(domain, adv);

  GreensFunction gf;
  gf.init(domain, adv);

  std::mt19937 rng(123);
  std::uniform_real_distribution<double> dist_xy(100e-6, 900e-6);
  std::uniform_real_distribution<double> dist_z(10e-6, 90e-6);

  int N = 50;
  std::vector<Vec3> positions(N);
  std::vector<GreensFunctionParams> params(N);
  std::vector<Real> strengths(N);

  for (int i = 0; i < N; ++i) {
    positions[i] = {dist_xy(rng), dist_xy(rng), dist_z(rng)};
    params[i].diff_coeff  = 4e-11;
    params[i].source_rate = 1e-18;
    params[i].pI          = 7.0;
    params[i].retardation = 5.0;
    strengths[i] = params[i].source_rate;
  }

  GreensFunctionParams avg_params;
  avg_params.diff_coeff  = 4e-11;
  avg_params.source_rate = 0.0;
  avg_params.pI          = 7.0;
  avg_params.retardation = 5.0;

  Octree tree;
  tree.build(positions, strengths, domain);

  // Test at several target points away from sources
  Real cutoff = 200e-6;
  Real theta = 0.3;  // tight opening angle for accuracy

  std::vector<Vec3> targets = {
    {200e-6, 200e-6, 50e-6},
    {800e-6, 800e-6, 50e-6},
    {500e-6, 500e-6, 25e-6},
  };

  for (const Vec3& tgt : targets) {
    // Exact: sum all sources
    Real exact = 0.0;
    for (int i = 0; i < N; ++i) {
      exact += gf.concentration_bounded(positions[i], tgt, params[i]);
    }

    // Barnes-Hut: octree evaluate_field (near + far)
    Real approx = tree.evaluate_field(
        tgt, theta, cutoff, gf, positions, params, avg_params);

    Real rel_err = (exact > 1e-30)
        ? std::abs(approx - exact) / exact
        : std::abs(approx - exact);

    // With theta=0.3, relative error should be well under 50%
    // (Barnes-Hut monopole gives ~theta^2 relative error for smooth fields)
    std::cout << "    target=(" << tgt[0]*1e6 << "," << tgt[1]*1e6
              << "," << tgt[2]*1e6 << ") exact=" << exact
              << " approx=" << approx << " rel_err=" << rel_err << "\n";
    assert(rel_err < 0.5);
  }

  std::cout << "  test_accuracy_vs_exact: PASSED\n";
}

void test_far_field_only() {
  // With a very small cutoff, most sources are far-field.
  // The octree far-field evaluation should be non-negative.
  Domain domain;
  AdvectionField adv;
  make_test_domain(domain, adv);

  GreensFunction gf;
  gf.init(domain, adv);

  std::mt19937 rng(456);
  std::uniform_real_distribution<double> dist_xy(100e-6, 900e-6);
  std::uniform_real_distribution<double> dist_z(10e-6, 90e-6);

  int N = 30;
  std::vector<Vec3> positions(N);
  std::vector<Real> strengths(N, 1e-18);
  for (int i = 0; i < N; ++i) {
    positions[i] = {dist_xy(rng), dist_xy(rng), dist_z(rng)};
  }

  GreensFunctionParams avg;
  avg.diff_coeff  = 4e-11;
  avg.source_rate = 0.0;
  avg.pI          = 7.0;
  avg.retardation = 5.0;

  Octree tree;
  tree.build(positions, strengths, domain);

  Vec3 target = {50e-6, 50e-6, 50e-6};
  Real far_val = tree.evaluate_far_field(target, 0.5, 10e-6, gf, avg);
  assert(far_val >= 0.0);

  std::cout << "  test_far_field_only: PASSED (far=" << far_val << ")\n";
}

void test_qssa_config_fmm_flag() {
  // Verify QSSAConfig default has use_fmm=false and theta=0.5
  QSSAConfig cfg;
  assert(cfg.use_fmm == false);
  assert(std::abs(cfg.fmm_theta - 0.5) < 1e-10);

  cfg.use_fmm = true;
  cfg.fmm_theta = 0.3;
  assert(cfg.use_fmm == true);
  assert(std::abs(cfg.fmm_theta - 0.3) < 1e-10);

  std::cout << "  test_qssa_config_fmm_flag: PASSED\n";
}

int main() {
  std::cout << "=== Octree / Barnes-Hut Tests ===\n";
  test_empty_octree();
  test_single_source();
  test_octree_construction();
  test_accuracy_vs_exact();
  test_far_field_only();
  test_qssa_config_fmm_flag();
  std::cout << "All octree tests passed.\n";
  return 0;
}
