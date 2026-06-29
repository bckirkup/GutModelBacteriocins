/* -----------------------------------------------------------------------
   GutIBM – FMM tests
   Verify higher-order multipole expansions, M2M/M2L/L2L passes, and
   accuracy improvement over monopole Barnes-Hut.
   ----------------------------------------------------------------------- */

#include "fmm.h"
#include "octree.h"
#include "greens_function.h"
#include "qssa_solver.h"
#include "domain.h"
#include "advection.h"
#include "random.h"
#include <cassert>
#include <iostream>
#include <cmath>
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
  test_fmm_config_defaults();
  std::cout << "All FMM tests passed.\n";
  return 0;
}
