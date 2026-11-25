/* -----------------------------------------------------------------------
   GutIBM – FMM tests
   Verify higher-order multipole expansions, M2M/M2L/L2L passes, and
   accuracy improvement over monopole Barnes-Hut.
   ----------------------------------------------------------------------- */

#include "fmm.h"
#include "fmm_kernel.h"
#include "octree.h"
#include "greens_function.h"
#include "qssa_solver.h"
#include "domain.h"
#include "advection.h"
#include "random.h"
#include "test_shard_dispatch.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

using namespace gutibm;

static void make_test_domain(Domain& domain, AdvectionField& adv) {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {1e-3, 1e-3, 100e-6};
  dcfg.grid_dx = 5e-6;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 1e20;
  acfg.distal_transit_time = 1e20;
  acfg.mucus_thickness = 100e-6;
  acfg.distal_length = 1e-3;
  adv.init(acfg, domain);
}

static void make_random_sources(int N, RNG& rng,
                                std::vector<Vec3>& positions,
                                std::vector<GreensFunctionParams>& params,
                                std::vector<Real>& strengths) {
  positions.resize(N);
  params.resize(N);
  strengths.resize(N);

  for (int i = 0; i < N; ++i) {
    positions[i] = {rng.uniform(100e-6, 900e-6), rng.uniform(100e-6, 900e-6),
                    rng.uniform(10e-6, 90e-6)};
    params[i].diff_coeff  = 4e-11;
    params[i].source_rate = 1e-18;
    params[i].pI          = 7.0;
    params[i].retardation = 5.0;
    strengths[i] = params[i].source_rate;
  }
}

void test_fmm_coefficient_count() {
  assert(FMM::num_coefficients(0) == 1);
  assert(FMM::num_coefficients(1) == 4);
  assert(FMM::num_coefficients(2) == 10);
  assert(FMM::num_coefficients(3) == 20);
  std::cout << "  test_fmm_coefficient_count: PASSED\n";
}

void test_fmm_build_and_moments() {
  Domain domain;
  AdvectionField adv;
  make_test_domain(domain, adv);

  RNG rng(7);
  std::vector<Vec3> positions;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strengths;
  make_random_sources(40, rng, positions, params, strengths);

  FMM fmm;
  fmm.build(positions, strengths, domain, 2);
  assert(!fmm.empty());

  const FMMNode& root = fmm.node(0);
  Real expected_total = 40 * 1e-18;
  assert(std::abs(root.total_source_strength - expected_total) / expected_total < 1e-10);
  assert(static_cast<int>(root.multipole.size()) == FMM::num_coefficients(2));
  assert(root.multipole[0] > 0.0);  // monopole

  std::cout << "  test_fmm_build_and_moments: PASSED\n";
}

void test_fmm_degenerate_clusters_match_direct_field() {
  Domain domain;
  AdvectionField adv;
  make_test_domain(domain, adv);

  GreensFunction gf;
  gf.init(domain, adv);

  GreensFunctionParams params;
  params.diff_coeff = 4e-11;
  params.source_rate = 1e-18;
  params.pI = 7.0;
  params.retardation = 5.0;
  const Vec3 target = {300e-6, 300e-6, 50e-6};

  for (const bool tight_cluster : {false, true}) {
    std::vector<Vec3> positions(9);
    std::vector strengths(9, params.source_rate);
    for (int i = 0; i < 9; ++i) {
      const Real offset = tight_cluster
                            ? static_cast<Real>(i) * 1e-12
                            : 0.0;
      positions[static_cast<size_t>(i)] = {
          500e-6 + offset, 500e-6, 50e-6};
    }

    FMM fmm;
    fmm.build(positions, strengths, domain, 2);

    Real exact = 0.0;
    for (const Vec3& source : positions)
      exact += gf.concentration_bounded(source, target, params);

    const GreensFunctionParams avg = params;
    const Real approx = fmm.evaluate_field(
        target, 0.5, 10e-6, gf, positions,
        std::vector<GreensFunctionParams>(positions.size(), params), avg);
    const Real relative_error = std::abs(approx - exact) / exact;
    assert(relative_error < 1e-6);
  }

  std::cout << "  test_fmm_degenerate_clusters_match_direct_field: PASSED\n";
}

void test_fmm_accuracy_order2_vs_exact() {
  Domain domain;
  AdvectionField adv;
  make_test_domain(domain, adv);

  GreensFunction gf;
  gf.init(domain, adv);

  RNG rng(321);
  std::vector<Vec3> positions;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strengths;
  make_random_sources(50, rng, positions, params, strengths);

  GreensFunctionParams avg;
  avg.diff_coeff  = 4e-11;
  avg.pI          = 7.0;
  avg.retardation = 5.0;
  avg.source_rate = 0.0;

  Real cutoff = 200e-6;
  Real theta = 0.5;

  FMM fmm;
  fmm.build(positions, strengths, domain, 2);
  fmm.compute_local_expansions(theta, gf, avg);

  std::vector<Vec3> targets = {
    {200e-6, 200e-6, 50e-6},
    {800e-6, 800e-6, 50e-6},
    {500e-6, 500e-6, 25e-6},
  };

  for (const Vec3& tgt : targets) {
    Real exact = 0.0;
    for (size_t i = 0; i < positions.size(); ++i)
      exact += gf.concentration_bounded(positions[i], tgt, params[i]);

    Real approx = fmm.evaluate_field(tgt, theta, cutoff, gf, positions, params, avg);
    Real rel_err = (exact > 1e-30) ? std::abs(approx - exact) / exact
                                   : std::abs(approx - exact);

    std::cout << "    order=2 target rel_err=" << rel_err << "\n";
    assert(rel_err < 0.5);
  }

  std::cout << "  test_fmm_accuracy_order2_vs_exact: PASSED\n";
}

void test_fmm_higher_order_more_accurate_than_monopole() {
  Domain domain;
  AdvectionField adv;
  make_test_domain(domain, adv);

  GreensFunction gf;
  gf.init(domain, adv);

  RNG rng(999);
  std::vector<Vec3> positions;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strengths;
  make_random_sources(60, rng, positions, params, strengths);

  GreensFunctionParams avg;
  avg.diff_coeff  = 4e-11;
  avg.pI          = 7.0;
  avg.retardation = 5.0;
  avg.source_rate = 0.0;

  Real cutoff = 200e-6;
  Real theta = 0.6;
  Vec3 target = {300e-6, 700e-6, 40e-6};

  Real exact = 0.0;
  for (size_t i = 0; i < positions.size(); ++i)
    exact += gf.concentration_bounded(positions[i], target, params[i]);

  // Order 1: Barnes-Hut-style tree walk (no local expansion preprocessing)
  FMM fmm1;
  fmm1.build(positions, strengths, domain, 1);
  Real err1 = std::abs(fmm1.evaluate_field(target, theta, cutoff, gf,
                                           positions, params, avg) - exact);
  if (exact > 1e-30) err1 /= exact;

  // Order 2: full FMM with local expansions
  FMM fmm2;
  fmm2.build(positions, strengths, domain, 2);
  fmm2.compute_local_expansions(theta, gf, avg);
  Real err2 = std::abs(fmm2.evaluate_field(target, theta, cutoff, gf,
                                           positions, params, avg) - exact);
  if (exact > 1e-30) err2 /= exact;

  std::cout << "    monopole rel_err=" << err1 << " order2 rel_err=" << err2 << "\n";

  // Higher order should not be worse; typically tighter for clustered sources.
  assert(err2 <= err1 * 1.5 + 1e-12);

  std::cout << "  test_fmm_higher_order_more_accurate_than_monopole: PASSED\n";
}

void test_fmm_local_expansion_nonnegative() {
  Domain domain;
  AdvectionField adv;
  make_test_domain(domain, adv);

  GreensFunction gf;
  gf.init(domain, adv);

  RNG rng(55);
  std::vector<Vec3> positions;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strengths;
  make_random_sources(30, rng, positions, params, strengths);

  GreensFunctionParams avg;
  avg.diff_coeff  = 4e-11;
  avg.pI          = 7.0;
  avg.retardation = 5.0;
  avg.source_rate = 0.0;

  FMM fmm;
  fmm.build(positions, strengths, domain, 2);
  fmm.compute_local_expansions(0.5, gf, avg);

  Vec3 target = {50e-6, 50e-6, 50e-6};
  Real far = fmm.evaluate_far_field(target, 0.5, 10e-6, gf, avg);
  assert(far >= 0.0);

  std::cout << "  test_fmm_local_expansion_nonnegative: PASSED\n";
}

// Spec 0 §3: the finite-difference derivative machinery must be numerically
// stable up to order 3. Validate against the analytic 1/r kernel, whose
// single-axis derivatives have closed form. A tiny fixed step (the former
// 1e-9) makes the order-3 stencil (÷2h^3) lose all precision; the order- and
// length-scaled step keeps every order accurate.
void test_fmm_kernel_derivative_accuracy() {
  const Vec3 x0 = {0.0, 0.0, 0.0};
  const auto phi = [&x0](const Vec3& p) {
    const Real dx = p[0] - x0[0];
    const Real dy = p[1] - x0[1];
    const Real dz = p[2] - x0[2];
    return 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz);
  };

  const Vec3 c = {3.0e-5, 1.0e-5, 2.0e-5};
  const Real dx = c[0] - x0[0];
  const Real r2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
  const Real r = std::sqrt(r2);

  // Analytic single-axis (x) derivatives of 1/r.
  const Real d1 = -dx / (r2 * r);                                  // φ_x
  const Real d2 = (3.0 * dx * dx - r2) / (r2 * r2 * r);            // φ_xx
  const Real d3 = 9.0 * dx / (r2 * r2 * r)                         // φ_xxx
                  - 15.0 * dx * dx * dx / (r2 * r2 * r2 * r);

  const Real fd1 = fd_axis_derivative(phi, c, 0, 1, r);
  const Real fd2 = fd_axis_derivative(phi, c, 0, 2, r);
  const Real fd3 = fd_axis_derivative(phi, c, 0, 3, r);

  const Real e1 = std::abs(fd1 - d1) / std::abs(d1);
  const Real e2 = std::abs(fd2 - d2) / std::abs(d2);
  const Real e3 = std::abs(fd3 - d3) / std::abs(d3);
  std::cout << "    1/r derivative rel_err: order1=" << e1
            << " order2=" << e2 << " order3=" << e3 << std::endl;

  // Order 1/2 reach ~1e-8 or better; order 3 (÷h^3) is near the double-
  // precision floor for a central difference of 1/r (~1e-5) but is now
  // accurate, vs O(1e-3) with the former fixed 1e-9 step at this micron scale.
  assert(e1 < 1e-6);
  assert(e2 < 1e-6);
  assert(e3 < 1e-4);

  // Demonstrate the fix: the former fixed 1e-9 step is far less accurate for
  // the order-3 derivative at this micron length scale (measured ~200x worse).
  const Real h_old = 1.0e-9;
  const auto at = [&phi, &c](Real off) { Vec3 p = c; p[0] += off; return phi(p); };
  const Real fd3_old = (at(2 * h_old) - 2 * at(h_old) + 2 * at(-h_old) - at(-2 * h_old))
                       / (2 * h_old * h_old * h_old);
  const Real e3_old = std::abs(fd3_old - d3) / std::abs(d3);
  std::cout << "    order3 rel_err with old fixed h=1e-9: " << e3_old << std::endl;
  assert(e3_old > e3 * 10.0);  // length-scaled step is >10x more accurate

  std::cout << "  test_fmm_kernel_derivative_accuracy: PASSED\n";
}

void test_fmm_large_tree_locals_ready() {
  Domain domain;
  AdvectionField adv;
  make_test_domain(domain, adv);

  GreensFunction gf;
  gf.init(domain, adv);

  RNG rng(2024);
  std::vector<Vec3> positions;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strengths;
  make_random_sources(1000, rng, positions, params, strengths);

  GreensFunctionParams avg;
  avg.diff_coeff  = 4e-11;
  avg.pI          = 7.0;
  avg.retardation = 5.0;
  avg.source_rate = 0.0;

  Real cutoff = 200e-6;
  Real theta = 0.5;
  Vec3 target = {400e-6, 600e-6, 35e-6};

  FMM fmm;
  fmm.build(positions, strengths, domain, 2);
  assert(fmm.num_nodes() > 256);
  fmm.compute_local_expansions(theta, gf, avg);
  assert(fmm.locals_ready());

  Real exact = 0.0;
  for (size_t i = 0; i < positions.size(); ++i) {
    exact += gf.concentration_bounded(positions[i], target, params[i]);
  }

  Real approx = fmm.evaluate_field(target, theta, cutoff, gf,
                                 positions, params, avg);
  Real rel_err = (exact > 1e-30) ? std::abs(approx - exact) / exact
                                 : std::abs(approx - exact);
  std::cout << "    large_tree nodes=" << fmm.num_nodes()
            << " rel_err=" << rel_err << "\n";
  assert(rel_err < 0.5);

  std::cout << "  test_fmm_large_tree_locals_ready: PASSED\n";
}

void test_fmm_config_defaults() {
  QSSAConfig cfg;
  assert(cfg.fmm_expansion_order == 2);
  std::cout << "  test_fmm_config_defaults: PASSED\n";
}

void test_fmm_periodic_cross_seam_accuracy() {
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {200.0e-6, 200.0e-6, 100.0e-6};
  dcfg.periodic = {true, true, false};
  dcfg.grid_dx = 5.0e-6;
  Domain domain;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 1e20;
  acfg.distal_transit_time = 1e20;
  acfg.mucus_thickness = 100.0e-6;
  acfg.distal_length = 200.0e-6;
  AdvectionField adv;
  adv.init(acfg, domain);

  GreensFunction gf;
  gf.init(domain, adv);

  RNG rng(15901);
  std::vector<Vec3> positions;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strengths;
  auto add_cluster = [&](int count, Real x0, Real x1, Real y0, Real y1,
                         Real z0, Real z1) {
    for (int i = 0; i < count; ++i) {
      positions.push_back({rng.uniform(x0, x1), rng.uniform(y0, y1),
                           rng.uniform(z0, z1)});
      GreensFunctionParams p;
      p.diff_coeff = 4.0e-11;
      p.source_rate = 1.0e-18;
      p.pI = 7.0;
      p.retardation = 5.0;
      p.decay_rate = 1.0e-3;
      params.push_back(p);
      strengths.push_back(p.source_rate);
    }
  };

  // Producers within 10 um of x = lo[0] (cross-seam neighbors of x ~ hi[0]),
  // plus far clusters around domain center and x ~ hi[0] so targets near hi[0]
  // exercise both cross-seam near exclusion and periodic M2L far-field.
  add_cluster(24, 2.0e-6, 9.0e-6, 85.0e-6, 115.0e-6, 35.0e-6, 65.0e-6);
  add_cluster(16, 2.0e-6, 8.0e-6, 2.0e-6, 9.0e-6, 40.0e-6, 60.0e-6);
  add_cluster(24, 92.0e-6, 108.0e-6, 92.0e-6, 108.0e-6, 35.0e-6, 65.0e-6);
  add_cluster(16, 191.0e-6, 198.0e-6, 140.0e-6, 160.0e-6, 40.0e-6, 60.0e-6);

  GreensFunctionParams avg = params[0];
  avg.source_rate = 0.0;

  // Targets near x = hi[0] spanning cross-seam near-field and empty/populated octants.
  const std::vector<Vec3> targets = {
      {195.0e-6, 100.0e-6, 50.0e-6},
      {192.0e-6, 95.0e-6, 45.0e-6},
      {197.0e-6, 105.0e-6, 55.0e-6},
      {196.0e-6, 196.0e-6, 50.0e-6},
      {194.0e-6, 5.0e-6, 50.0e-6},
      {193.0e-6, 150.0e-6, 50.0e-6},
      {196.0e-6, 148.0e-6, 48.0e-6},
      {190.0e-6, 60.0e-6, 40.0e-6},
      {195.0e-6, 180.0e-6, 60.0e-6},
      {198.0e-6, 100.0e-6, 35.0e-6},
  };

  std::vector<Real> exact(targets.size(), 0.0);
  Real exact_l2_sq = 0.0;
  for (size_t t = 0; t < targets.size(); ++t) {
    Real sum = 0.0;
    for (size_t s = 0; s < positions.size(); ++s) {
      sum += gf.concentration_bounded(positions[s], targets[t], params[s]);
    }
    exact[t] = sum;
    exact_l2_sq += sum * sum;
  }
  const Real exact_l2 = std::sqrt(exact_l2_sq);
  assert(exact_l2 > 0.0);

  constexpr Real kCutoff = 35.0e-6;
  constexpr Real kTheta = 0.5;
  for (int order : {2, 3}) {
    FMM fmm_tree;
    fmm_tree.build(positions, strengths, domain, order);
    assert(!fmm_tree.locals_ready());

    FMM fmm_local;
    fmm_local.build(positions, strengths, domain, order);
    fmm_local.compute_local_expansions(kTheta, gf, avg, kCutoff);
    assert(fmm_local.locals_ready());

    Real err_tree_sq = 0.0;
    Real err_local_sq = 0.0;
    for (size_t t = 0; t < targets.size(); ++t) {
      const Real val_tree = fmm_tree.evaluate_field(
          targets[t], kTheta, kCutoff, gf, positions, params, avg);
      const Real val_local = fmm_local.evaluate_field(
          targets[t], kTheta, kCutoff, gf, positions, params, avg);
      assert(std::isfinite(val_tree) && val_tree >= 0.0);
      assert(std::isfinite(val_local) && val_local >= 0.0);
      const Real dt = val_tree - exact[t];
      const Real dl = val_local - exact[t];
      err_tree_sq += dt * dt;
      err_local_sq += dl * dl;
    }

    const Real rel_l2_tree = std::sqrt(err_tree_sq) / exact_l2;
    const Real rel_l2_local = std::sqrt(err_local_sq) / exact_l2;
    std::cout << "    periodic cross-seam order=" << order
              << " rel_l2_tree=" << rel_l2_tree
              << " rel_l2_local=" << rel_l2_local
              << " m2l_count=" << fmm_local.m2l_interaction_count() << "\n";
    assert(rel_l2_tree <= 2.0e-2);
    assert(rel_l2_local <= 2.0e-2);
  }

  // Single-source VoxelBox near+far disjoint partition check:
  // Diagonal voxels inside the [-span..span] near box but outside the radial
  // cutoff sphere must not be double-counted when combining superpose_to_grid
  // with VoxelBox FMM far-field.
  {
    const std::vector<Vec3> single_pos = {{100.0e-6, 100.0e-6, 50.0e-6}};
    const std::vector<GreensFunctionParams> single_params = {params[0]};
    const std::vector<Real> single_str = {params[0].source_rate};
    constexpr Real kBoxCutoff = 40.0e-6;
    std::vector<Real> near_grid;
    gf.superpose_to_grid(single_pos, single_params, near_grid, kBoxCutoff);

    FMM fmm_box;
    fmm_box.build(single_pos, single_str, domain, 2);
    fmm_box.compute_local_expansions(
        kTheta, gf, avg, kBoxCutoff, FmmNearCutoffShape::VoxelBox);

    // Check diagonal inside box (dx=30um, dy=30um -> r=42.4um > 40um) and far cell.
    for (const std::array<Int, 3>& cell :
         std::vector<std::array<Int, 3>>{{13, 13, 5}, {16, 10, 5}, {2, 10, 5}}) {
      const Vec3 tgt = domain.cell_center(cell[0], cell[1], cell[2]);
      const Int flat = cell[2] * domain.nx() * domain.ny()
                     + cell[1] * domain.nx() + cell[0];
      const Real near = near_grid[static_cast<size_t>(flat)];
      const Real far = fmm_box.evaluate_far_field(
          tgt, kTheta, kBoxCutoff, gf, avg, FmmNearCutoffShape::VoxelBox);
      const Real ref = gf.concentration_bounded(single_pos[0], tgt, single_params[0]);
      assert(ref > 0.0);
      assert(std::abs((near + far) - ref) / ref <= 1.0e-12);
    }
  }

  // Nonzero distal advection across periodic half-domain branch cut (|dx| ~ L_x/2).
  {
    AdvectionConfig flow_cfg = acfg;
    flow_cfg.distal_transit_time = 200.0;  // v_distal_max = 1 um/s along +x
    flow_cfg.radial_turnover = 1e20;
    AdvectionField flow_adv;
    flow_adv.init(flow_cfg, domain);
    GreensFunction flow_gf;
    flow_gf.init(domain, flow_adv);

    FMM fmm_flow;
    fmm_flow.build(positions, strengths, domain, 2);
    fmm_flow.compute_local_expansions(kTheta, flow_gf, avg, kCutoff);

    const std::vector<Vec3> cut_targets = {
        {102.0e-6, 100.0e-6, 50.0e-6},
        {105.0e-6, 100.0e-6, 50.0e-6},
        {108.0e-6, 100.0e-6, 50.0e-6},
        {195.0e-6, 100.0e-6, 50.0e-6},
    };
    Real flow_err_sq = 0.0;
    Real flow_ref_sq = 0.0;
    for (const Vec3& tgt : cut_targets) {
      Real ref = 0.0;
      for (size_t s = 0; s < positions.size(); ++s) {
        ref += flow_gf.concentration_bounded(positions[s], tgt, params[s]);
      }
      const Real est = fmm_flow.evaluate_field(
          tgt, kTheta, kCutoff, flow_gf, positions, params, avg);
      assert(std::isfinite(est) && est >= 0.0);
      flow_err_sq += (est - ref) * (est - ref);
      flow_ref_sq += ref * ref;
    }
    const Real flow_rel_l2 = std::sqrt(flow_err_sq / flow_ref_sq);
    std::cout << "    nonzero-flow branch-cut rel_l2=" << flow_rel_l2 << "\n";
    assert(flow_rel_l2 <= 2.0e-2);
  }

  std::cout << "  test_fmm_periodic_cross_seam_accuracy: PASSED\n";
}

void test_fmm_subquadratic_m2l_scaling() {
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {400.0e-6, 400.0e-6, 400.0e-6};
  dcfg.periodic = {true, true, false};
  dcfg.grid_dx = 10.0e-6;
  Domain domain;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 1e20;
  acfg.distal_transit_time = 1e20;
  acfg.mucus_thickness = 400.0e-6;
  acfg.distal_length = 400.0e-6;
  AdvectionField adv;
  adv.init(acfg, domain);

  GreensFunction gf;
  gf.init(domain, adv);

  GreensFunctionParams avg;
  avg.diff_coeff = 4.0e-11;
  avg.pI = 7.0;
  avg.retardation = 5.0;
  avg.source_rate = 0.0;
  avg.decay_rate = 1.0e-3;

  const std::array<int, 3> counts = {1000, 4000, 16000};
  std::array<int, 3> node_counts = {};
  std::array<uint64_t, 3> m2l_counts = {};
  std::array<double, 3> elapsed_ms = {};

  for (size_t k = 0; k < counts.size(); ++k) {
    RNG rng(15900 + static_cast<uint64_t>(k));
    std::vector<Vec3> positions(static_cast<size_t>(counts[k]));
    std::vector<Real> strengths(static_cast<size_t>(counts[k]), 1.0e-18);
    for (int i = 0; i < counts[k]; ++i) {
      positions[static_cast<size_t>(i)] = {
          rng.uniform(0.0, 400.0e-6),
          rng.uniform(0.0, 400.0e-6),
          rng.uniform(20.0e-6, 380.0e-6)};
    }

    FMM fmm;
    fmm.build(positions, strengths, domain, 1);
    assert(fmm.num_nodes() > 256);

    const auto t0 = std::chrono::steady_clock::now();
    fmm.compute_local_expansions(0.5, gf, avg, 0.0);
    const auto t1 = std::chrono::steady_clock::now();
    assert(fmm.locals_ready());

    node_counts[k] = fmm.num_nodes();
    m2l_counts[k] = fmm.m2l_interaction_count();
    elapsed_ms[k] = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double per_node =
        static_cast<double>(m2l_counts[k]) / static_cast<double>(node_counts[k]);
    std::cout << "    N_sources=" << counts[k]
              << " N_nodes=" << node_counts[k]
              << " M2L=" << m2l_counts[k]
              << " M2L/node=" << per_node
              << " time_ms=" << elapsed_ms[k] << std::endl;
  }

  // Verify sub-quadratic scaling: each 4x source increase yields <= 5.0x M2L
  // interactions (vs ~16x quadratic), and M2L interactions per node stay bounded.
  for (size_t k = 1; k < counts.size(); ++k) {
    assert(node_counts[k] > node_counts[k - 1]);
    assert(m2l_counts[k - 1] > 0);
    const double node_ratio =
        static_cast<double>(node_counts[k]) / static_cast<double>(node_counts[k - 1]);
    const double m2l_ratio =
        static_cast<double>(m2l_counts[k]) / static_cast<double>(m2l_counts[k - 1]);
    const double per_node_prev =
        static_cast<double>(m2l_counts[k - 1]) / static_cast<double>(node_counts[k - 1]);
    const double per_node_curr =
        static_cast<double>(m2l_counts[k]) / static_cast<double>(node_counts[k]);
    const double per_node_ratio = per_node_curr / per_node_prev;
    std::cout << "    step " << k << ": node_ratio=" << node_ratio
              << " m2l_ratio=" << m2l_ratio
              << " per_node_ratio=" << per_node_ratio << std::endl;
    assert(m2l_ratio <= 5.0);
    assert(per_node_ratio <= 1.5);
  }
  const double total_node_ratio =
      static_cast<double>(node_counts.back()) / static_cast<double>(node_counts.front());
  const double total_m2l_ratio =
      static_cast<double>(m2l_counts.back()) / static_cast<double>(m2l_counts.front());
  assert(total_m2l_ratio < 0.25 * total_node_ratio * total_node_ratio);

  std::cout << "  test_fmm_subquadratic_m2l_scaling: PASSED" << std::endl;
}

void test_fmm_benchmark_100k_sources() {
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {1.0e-3, 1.0e-3, 100.0e-6};
  dcfg.periodic = {true, true, false};
  dcfg.grid_dx = 4.0e-6;
  dcfg.chemistry_stride = {5, 5, 1};  // 50 x 50 x 25 = 62,500 voxels
  Domain domain;
  domain.init(dcfg);
  assert(domain.nx() == 50 && domain.ny() == 50 && domain.nz() == 25);

  AdvectionConfig acfg;
  acfg.radial_turnover = 1e20;
  acfg.distal_transit_time = 1e20;
  acfg.mucus_thickness = 100.0e-6;
  acfg.distal_length = 1.0e-3;
  AdvectionField adv;
  adv.init(acfg, domain);

  GreensFunction gf;
  gf.init(domain, adv);

  constexpr int kNumSources = 100000;
  RNG rng(159100);
  std::vector<Vec3> positions(kNumSources);
  std::vector<Real> strengths(kNumSources, 1.0e-18);
  GreensFunctionParams base_p;
  base_p.diff_coeff = 4.0e-11;
  base_p.source_rate = 1.0e-18;
  base_p.pI = 7.0;
  base_p.retardation = 5.0;
  base_p.decay_rate = 1.0e-2;  // ~69 s protease half-life (screening length ~28 um)

  for (int i = 0; i < kNumSources; ++i) {
    positions[static_cast<size_t>(i)] = {
        rng.uniform(0.0, 1.0e-3),
        rng.uniform(0.0, 1.0e-3),
        rng.uniform(5.0e-6, 95.0e-6)};
  }

  GreensFunctionParams avg = base_p;
  avg.source_rate = 0.0;
  constexpr Real kCutoff = 60.0e-6;
  constexpr Real kTheta = 0.5;

  const auto t_build_0 = std::chrono::steady_clock::now();
  FMM fmm;
  fmm.build(positions, strengths, domain, 1);
  const auto t_build_1 = std::chrono::steady_clock::now();
  fmm.compute_local_expansions(kTheta, gf, avg, kCutoff);
  const auto t_m2l_1 = std::chrono::steady_clock::now();

  assert(fmm.locals_ready());
  assert(fmm.num_nodes() > 256);

  // Evaluate far-field across the full 50x50x25 = 62,500-voxel grid.
  const Int ncells = domain.ncells();
  std::vector<Real> far_grid(static_cast<size_t>(ncells), 0.0);
  const auto t_eval_0 = std::chrono::steady_clock::now();
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.nx(); ++ix) {
        const Vec3 tgt = domain.cell_center(ix, iy, iz);
        far_grid[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
            fmm.evaluate_far_field(tgt, kTheta, kCutoff, gf, avg);
      }
    }
  }
  const auto t_eval_1 = std::chrono::steady_clock::now();

  // Verify accuracy on 8 sampled target voxels and time brute-force subset.
  constexpr int kSampleCount = 8;
  const std::array<std::array<Int, 3>, kSampleCount> sample_coords = {{
      {5, 5, 12}, {25, 25, 12}, {45, 45, 12}, {2, 48, 8},
      {48, 2, 16}, {15, 35, 6}, {35, 15, 18}, {25, 48, 12},
  }};
  const Real cutoff2 = kCutoff * kCutoff;
  Real err_sq = 0.0;
  Real exact_sq = 0.0;
  const auto t_bf_0 = std::chrono::steady_clock::now();
  for (const auto& coord : sample_coords) {
    const Vec3 tgt = domain.cell_center(coord[0], coord[1], coord[2]);
    Real exact_total = 0.0;
    Real exact_near = 0.0;
    for (int s = 0; s < kNumSources; ++s) {
      const Real c = gf.concentration_bounded(
          positions[static_cast<size_t>(s)], tgt, base_p);
      exact_total += c;
      if (domain.min_image_dist_sq(positions[static_cast<size_t>(s)], tgt) <= cutoff2) {
        exact_near += c;
      }
    }
    const Int cell = domain.cell_index(coord[0], coord[1], coord[2]);
    const Real approx_total = exact_near + far_grid[static_cast<size_t>(cell)];
    assert(std::isfinite(approx_total) && approx_total > 0.0);
    const Real diff = approx_total - exact_total;
    err_sq += diff * diff;
    exact_sq += exact_total * exact_total;
  }
  const auto t_bf_1 = std::chrono::steady_clock::now();

  const Real rel_l2 = std::sqrt(err_sq) / std::sqrt(exact_sq);
  const double build_ms =
      std::chrono::duration<double, std::milli>(t_build_1 - t_build_0).count();
  const double m2l_ms =
      std::chrono::duration<double, std::milli>(t_m2l_1 - t_build_1).count();
  const double eval_ms =
      std::chrono::duration<double, std::milli>(t_eval_1 - t_eval_0).count();
  const double fmm_total_ms = build_ms + m2l_ms + eval_ms;
  const double bf_sample_ms =
      std::chrono::duration<double, std::milli>(t_bf_1 - t_bf_0).count();
  const double bf_extrap_ms =
      bf_sample_ms * (static_cast<double>(ncells) / static_cast<double>(kSampleCount));
  const double speedup = bf_extrap_ms / std::max(fmm_total_ms, 1.0e-6);

  std::cout << "    100k-source FMM benchmark (grid=" << domain.nx() << "x"
            << domain.ny() << "x" << domain.nz() << ", nodes=" << fmm.num_nodes()
            << ", M2L=" << fmm.m2l_interaction_count() << "):\n"
            << "      build_ms=" << build_ms << " m2l_ms=" << m2l_ms
            << " grid_eval_ms=" << eval_ms << " fmm_total_ms=" << fmm_total_ms
            << "\n      bf_extrap_ms=" << bf_extrap_ms
            << " speedup=" << speedup << "x rel_l2=" << rel_l2 << std::endl;

  assert(rel_l2 <= 2.0e-2);
  assert(speedup > 2.0);
  std::cout << "  test_fmm_benchmark_100k_sources: PASSED" << std::endl;
}

namespace {

const std::array<testshards::NamedTest, 12> kTests = {{
    {"coefficient_count", test_fmm_coefficient_count},
    {"build_and_moments", test_fmm_build_and_moments},
    {"degenerate_clusters_match_direct_field",
     test_fmm_degenerate_clusters_match_direct_field},
    {"accuracy_order2_vs_exact", test_fmm_accuracy_order2_vs_exact},
    {"higher_order_vs_monopole",
     test_fmm_higher_order_more_accurate_than_monopole},
    {"local_expansion_nonnegative", test_fmm_local_expansion_nonnegative},
    {"kernel_derivative_accuracy", test_fmm_kernel_derivative_accuracy},
    {"large_tree_locals_ready", test_fmm_large_tree_locals_ready},
    {"config_defaults", test_fmm_config_defaults},
    {"periodic_cross_seam_accuracy", test_fmm_periodic_cross_seam_accuracy},
    {"subquadratic_m2l_scaling", test_fmm_subquadratic_m2l_scaling},
    {"benchmark_100k_sources", test_fmm_benchmark_100k_sources},
}};

}  // namespace

int main(int argc, char** argv) {
  return testshards::shard_main(
      {argv, static_cast<std::size_t>(argc)}, kTests,
      "=== FMM Tests ===\n", "All FMM tests passed.\n");
}
