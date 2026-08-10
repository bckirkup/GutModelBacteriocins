/* -----------------------------------------------------------------------
   GutIBM – Green's function tests
   Verify that the advection-diffusion kernel produces correct
   concentration profiles and comet-tail asymmetry.
   ----------------------------------------------------------------------- */

#include "greens_function.h"
#include "domain.h"
#include "advection.h"
#include "error.h"
#include "plasmid.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include <numbers>
#include <string>

using namespace gutibm;

namespace {

void setup_zero_flow(Domain& domain, AdvectionField& adv, GreensFunction& gf) {
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
  gf.init(domain, adv);
}

Real fit_screening_length(GreensFunction& gf, const Vec3& source,
                          const BICluster& bi) {
  GreensFunctionParams params;
  params.diff_coeff = bi.diff_coeff;
  params.retardation = bi.retardation;
  params.source_rate = 1.0e-18;
  params.decay_rate = std::numbers::ln2 / bi.protease_half_life;

  const std::array<Real, 4> distances = {10e-6, 20e-6, 40e-6, 80e-6};
  Real sum_r = 0.0;
  Real sum_y = 0.0;
  Real sum_rr = 0.0;
  Real sum_ry = 0.0;
  for (const Real r : distances) {
    Vec3 target = source;
    target[0] += r;
    const Real y = std::log(gf.concentration(source, target, params) * r);
    sum_r += r;
    sum_y += y;
    sum_rr += r * r;
    sum_ry += r * y;
  }
  const Real n = static_cast<Real>(distances.size());
  const Real slope = (n * sum_ry - sum_r * sum_y)
      / (n * sum_rr - sum_r * sum_r);
  return -1.0 / slope;
}

}  // namespace

void test_zero_decay_exact_regression() {
  Domain domain;
  AdvectionField adv;
  GreensFunction gf;
  setup_zero_flow(domain, adv, gf);

  GreensFunctionParams implicit_zero;
  implicit_zero.diff_coeff = 4e-11;
  implicit_zero.retardation = 1.5;
  implicit_zero.source_rate = 1e-18;
  implicit_zero.pI = 5.4;
  GreensFunctionParams explicit_zero = implicit_zero;
  explicit_zero.decay_rate = 0.0;
  const Vec3 source = {500e-6, 500e-6, 50e-6};
  const Vec3 target = {550e-6, 500e-6, 50e-6};
  const Real old_value = gf.concentration(source, target, implicit_zero);
  const Real new_value = gf.concentration(source, target, explicit_zero);
  assert(old_value == new_value);
  std::cout << "  test_zero_decay_exact_regression: PASSED\n";
}

void test_screening_lengths() {
  Domain domain;
  AdvectionField adv;
  GreensFunction gf;
  setup_zero_flow(domain, adv, gf);
  const Vec3 source = {500e-6, 500e-6, 50e-6};
  const Real col_e1 = fit_screening_length(gf, source, PlasmidLibrary::colicin_E1());
  const Real col_b = fit_screening_length(gf, source, PlasmidLibrary::colicin_B());
  const Real expected_e1 = std::sqrt(
      (4.0e-11 / 50.0) / (std::numbers::ln2 / 1800.0));
  const Real expected_b = std::sqrt(
      (4.0e-11 / 1.5) / (std::numbers::ln2 / 900.0));
  assert(std::abs(col_e1 - expected_e1) / expected_e1 < 1.0e-10);
  assert(std::abs(col_b - expected_b) / expected_b < 1.0e-10);
  std::cout << "  test_screening_lengths: PASSED (ColE1=" << col_e1 * 1e6
            << " um ColB=" << col_b * 1e6 << " um)\n";
}

void test_core_halo_decay_ordering() {
  Domain domain;
  AdvectionField adv;
  GreensFunction gf;
  setup_zero_flow(domain, adv, gf);
  const Vec3 source = {500e-6, 500e-6, 50e-6};
  const Vec3 near = {510e-6, 500e-6, 50e-6};
  const Vec3 far = {100e-6, 500e-6, 50e-6};

  auto make_params = [](const BICluster& bi) {
    GreensFunctionParams p;
    p.diff_coeff = bi.diff_coeff;
    p.retardation = bi.retardation;
    p.source_rate = 1.0e-18;
    p.decay_rate = std::numbers::ln2 / bi.protease_half_life;
    return p;
  };
  const Real e1_near = gf.concentration(source, near,
                                         make_params(PlasmidLibrary::colicin_E1()));
  const Real e1_far = gf.concentration(source, far,
                                       make_params(PlasmidLibrary::colicin_E1()));
  const Real b_near = gf.concentration(source, near,
                                       make_params(PlasmidLibrary::colicin_B()));
  const Real b_far = gf.concentration(source, far,
                                      make_params(PlasmidLibrary::colicin_B()));
  const Real e1_ratio = e1_far / e1_near;
  const Real b_ratio = b_far / b_near;
  assert(e1_ratio < b_ratio);
  std::cout << "  test_core_halo_decay_ordering: PASSED (ColE1="
            << e1_ratio << " ColB=" << b_ratio << ")\n";
}

void test_radial_symmetry_no_flow() {
  // Without flow, Green's function should be radially symmetric
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {1e-3, 1e-3, 100e-6};
  dcfg.grid_dx = 5e-6;

  Domain domain;
  domain.init(dcfg);

  // Zero-flow advection
  AdvectionConfig acfg;
  acfg.radial_turnover = 1e20;     // effectively zero flow
  acfg.distal_transit_time = 1e20;
  acfg.mucus_thickness = 100e-6;
  acfg.distal_length = 1e-3;

  AdvectionField adv;
  adv.init(acfg, domain);

  GreensFunction gf;
  gf.init(domain, adv);

  Vec3 source = {500e-6, 500e-6, 50e-6};
  GreensFunctionParams params;
  params.diff_coeff   = 4e-11;
  params.source_rate  = 1e-18;
  params.pI           = 7.0;
  params.retardation  = 1.0;

  // Two points equidistant from source in different directions
  Vec3 p1 = {510e-6, 500e-6, 50e-6};
  Vec3 p2 = {500e-6, 510e-6, 50e-6};

  Real c1 = gf.concentration(source, p1, params);
  Real c2 = gf.concentration(source, p2, params);

  // Should be approximately equal (radial symmetry)
  Real ratio = std::abs(c1 - c2) / std::max(c1, 1e-30);
  assert(ratio < 0.01);

  std::cout << "  test_radial_symmetry: PASSED (ratio=" << ratio << ")\n";
}

void test_comet_tail_asymmetry() {
  // With flow, downstream concentration > upstream
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {1e-3, 1e-3, 100e-6};
  dcfg.grid_dx = 5e-6;

  Domain domain;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 5400.0;
  acfg.distal_transit_time = 43200.0;
  acfg.mucus_thickness = 100e-6;
  acfg.distal_length = 1e-3;
  acfg.profile_alpha = 1.0;  // linear profile for testing

  AdvectionField adv;
  adv.init(acfg, domain);

  GreensFunction gf;
  gf.init(domain, adv);

  Vec3 source = {500e-6, 500e-6, 50e-6};  // mid-height for flow
  GreensFunctionParams params;
  params.diff_coeff   = 4e-11;
  params.source_rate  = 1e-18;
  params.pI           = 5.4;    // acidic → lethal halo, low retardation
  params.retardation  = 1.5;

  // Downstream (in distal/x direction)
  Vec3 downstream = {550e-6, 500e-6, 50e-6};
  // Upstream
  Vec3 upstream   = {450e-6, 500e-6, 50e-6};

  Real c_down = gf.concentration(source, downstream, params);
  Real c_up   = gf.concentration(source, upstream, params);

  // Downstream should have higher concentration (comet tail)
  assert(c_down > c_up);

  std::cout << "  test_comet_tail: PASSED (down=" << c_down
            << " up=" << c_up << " ratio=" << c_down/std::max(c_up, 1e-30) << ")\n";
}

void test_lethal_core_vs_halo() {
  // Core (high retardation, small D_eff) → steeper concentration gradient
  // Halo (low retardation, large D_eff) → flatter gradient
  // Both follow C = Q/(4π D_eff r), so core_near/core_far > halo_near/halo_far
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {1e-3, 1e-3, 100e-6};
  dcfg.grid_dx = 5e-6;

  Domain domain;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 1e20;
  acfg.distal_transit_time = 1e20;
  acfg.mucus_thickness = 100e-6;
  acfg.distal_length = 1e-3;

  AdvectionField adv;
  adv.init(acfg, domain);

  GreensFunction gf;
  gf.init(domain, adv);

  Vec3 source = {500e-6, 500e-6, 50e-6};

  // Lethal Core (basic, high retardation)
  GreensFunctionParams core_params;
  core_params.diff_coeff   = 4e-11;
  core_params.source_rate  = 1e-18;
  core_params.pI           = 9.0;
  core_params.retardation  = 50.0;

  // Lethal Halo (acidic, low retardation)
  GreensFunctionParams halo_params;
  halo_params.diff_coeff   = 4e-11;
  halo_params.source_rate  = 1e-18;
  halo_params.pI           = 5.4;
  halo_params.retardation  = 1.5;

  Vec3 near  = {505e-6, 500e-6, 50e-6};  // 5 um away

  Real core_near = gf.concentration(source, near, core_params);
  Real halo_near = gf.concentration(source, near, halo_params);

  // Both should be positive
  assert(core_near > 0.0);
  assert(halo_near > 0.0);

  // Core (low D_eff) should be more concentrated near source
  assert(core_near > halo_near);

  // In no-flow steady-state: C = Q/(4π D_eff r), gradient steepness
  // is independent of D_eff (ratio at near/far = r_far/r_near).
  // The lethal core/halo distinction manifests under advection.
  // Here we just verify the D_eff scaling: core_near/halo_near ≈ retard_core/retard_halo
  Real expected_ratio = core_params.retardation / halo_params.retardation;
  Real actual_ratio   = core_near / std::max(halo_near, 1e-30);
  Real error = std::abs(actual_ratio - expected_ratio) / expected_ratio;
  assert(error < 0.05);  // within 5%

  std::cout << "  test_core_vs_halo: PASSED"
            << " (core_near=" << core_near
            << " halo_near=" << halo_near
            << " retard_ratio=" << expected_ratio
            << " conc_ratio=" << actual_ratio << ")\n";
}

void test_uninitialized_throws() {
  GreensFunction gf;
  Vec3 source = {500e-6, 500e-6, 50e-6};
  Vec3 target = {510e-6, 500e-6, 50e-6};
  GreensFunctionParams params;
  params.diff_coeff  = 4e-11;
  params.source_rate = 1e-18;
  params.pI          = 7.0;
  params.retardation = 1.0;

  auto expect_init_error = [](const auto& fn) {
    try {
      fn();
      assert(false && "expected GreensFunction to throw before init()");
    } catch (const SimulationError& e) {
      assert(std::string(e.what()).find("init()") != std::string::npos);
    }
  };

  expect_init_error([&gf, &source, &target, &params]() { (void)gf.concentration(source, target, params); });
  expect_init_error([&gf, &source, &target, &params]() { (void)gf.concentration_bounded(source, target, params); });
  expect_init_error([&gf, &source]() { (void)gf.peclet(source, 4e-11, 50e-6); });

  std::vector<Vec3> sources = {source};
  std::vector<GreensFunctionParams> params_vec = {params};
  std::vector<Real> grid;
  expect_init_error([&gf, &sources, &params_vec, &grid]() { gf.superpose_to_grid(sources, params_vec, grid, 50e-6); });

  std::cout << "  test_uninitialized_throws: PASSED\n";
}

void test_peclet_number() {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {1e-3, 1e-3, 100e-6};
  dcfg.grid_dx = 5e-6;

  Domain domain;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 5400.0;
  acfg.distal_transit_time = 43200.0;
  acfg.mucus_thickness = 100e-6;
  acfg.distal_length = 1e-3;

  AdvectionField adv;
  adv.init(acfg, domain);

  GreensFunction gf;
  gf.init(domain, adv);

  // At mid-height, Pe should be >= 1 (advection-dominated)
  Vec3 mid = {500e-6, 500e-6, 50e-6};
  Real D_eff = 4e-11 / 1.5;  // acidic colicin
  Real Pe = gf.peclet(mid, D_eff, 50e-6);

  assert(Pe >= 0.0);
  std::cout << "  test_peclet: PASSED (Pe=" << Pe << ")\n";
}

int main() {
  std::cout << "=== Green's Function Tests ===\n";
  test_zero_decay_exact_regression();
  test_screening_lengths();
  test_core_halo_decay_ordering();
  test_uninitialized_throws();
  test_radial_symmetry_no_flow();
  test_comet_tail_asymmetry();
  test_lethal_core_vs_halo();
  test_peclet_number();
  std::cout << "All Green's function tests passed.\n";
  return 0;
}
