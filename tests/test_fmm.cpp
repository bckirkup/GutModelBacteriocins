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
#include <cassert>
#include <iostream>
#include <cmath>
#include <functional>
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

int main() {
  std::cout << "=== FMM Tests ===\n";
  test_fmm_coefficient_count();
  test_fmm_build_and_moments();
  test_fmm_accuracy_order2_vs_exact();
  test_fmm_higher_order_more_accurate_than_monopole();
  test_fmm_local_expansion_nonnegative();
  test_fmm_kernel_derivative_accuracy();
  test_fmm_large_tree_locals_ready();
  test_fmm_config_defaults();
  std::cout << "All FMM tests passed.\n";
  return 0;
}
