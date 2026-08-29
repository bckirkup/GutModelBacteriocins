/* -----------------------------------------------------------------------
   GutIBM – Independent Robin lumen-boundary tests
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "domain.h"
#include "greens_function.h"
#include "robin_correction_table.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numbers>

using namespace gutibm;

namespace {

struct TestSystem {
  Domain domain;
  AdvectionField adv;
  GreensFunction gf;
};

TestSystem make_system() {
  TestSystem system;
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {1.0e-3, 1.0e-3, 100.0e-6};
  dcfg.grid_dx = 5.0e-6;
  system.domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 1.0e20;
  acfg.distal_transit_time = 1.0e20;
  acfg.mucus_thickness = 100.0e-6;
  acfg.distal_length = 1.0e-3;
  acfg.taylor_aris_enabled = false;
  system.adv.init(acfg, system.domain);
  system.gf.init(system.domain, system.adv);
  return system;
}

GreensFunctionParams params_for(Real d_eff, Real decay, Real transfer_length) {
  GreensFunctionParams params;
  params.diff_coeff = 4.0e-11;
  params.retardation = params.diff_coeff / d_eff;
  params.source_rate = 1.0e-18;
  params.decay_rate = decay;
  params.lumen_transfer_length = transfer_length;
  params.robin_cutoff = 200.0e-6;
  return params;
}

void require(bool condition, const char* message) {
  if (condition) return;
  std::cerr << message << "\n";
  std::exit(1);
}

void test_flux_residual() {
  auto system = make_system();
  const Real height = system.domain.size()[2];
  const Real d_eff = 2.0e-11;
  const Real d_free = 4.0e-11;
  const Real transfer_length = 100.0e-6;
  const Real decay = 5.0e-5;
  const Real kc = d_free / transfer_length;
  const Real epsilon = 1.0e-9;
  const Real step = 1.0e-9;
  Real maximum = 0.0;

  for (const Real source_fraction : {0.1, 0.35, 0.7, 0.9}) {
    const Vec3 source = {500.0e-6, 500.0e-6, source_fraction * height};
    for (const Real rho : {20.0e-6, 50.0e-6}) {
      for (const int wall : {0, 1}) {
        Vec3 target = {source[0] + rho, source[1],
                       wall == 0 ? epsilon : height - epsilon};
        Vec3 lower = target;
        Vec3 upper = target;
        lower[2] -= step;
        upper[2] += step;
        const Real value = robin::normalized_robin_field(
            source[2], target[2], rho, 0.0, height, d_eff, 4.0e-11,
            decay, transfer_length, 0.0, 0.0, 0.0, robin::kTableModeCount);
        const Real derivative = (
            robin::normalized_robin_field(
                source[2], upper[2], rho, 0.0, height, d_eff, 4.0e-11,
                decay, transfer_length, 0.0, 0.0, 0.0,
                robin::kTableModeCount)
            - robin::normalized_robin_field(
                source[2], lower[2], rho, 0.0, height, d_eff, 4.0e-11,
                decay, transfer_length, 0.0, 0.0, 0.0,
                robin::kTableModeCount))
            / (2.0 * step);
        const Real residual = wall == 0
            ? (-d_eff * derivative) * height / (d_eff * value)
            : (-d_eff * derivative - kc * value) * height
                / (d_eff * value);
        maximum = std::max(maximum, std::abs(residual));
        if (std::abs(residual) > 1.0e-2) {
          std::cerr << "  flux residual source=" << source_fraction
                    << " rho=" << rho << " wall=" << wall
                    << " value=" << value << " residual=" << residual
                    << "\n";
        }
      }
    }
  }
  std::cout << "  test_flux_residual: result (max residual="
            << maximum << ")\n";
  require(maximum < 1.0e-3, "Robin wall flux residual exceeded tolerance");
}

void test_table_against_direct_modes() {
  auto system = make_system();
  const Real z_lo = system.domain.lo()[2];
  const Real z_hi = system.domain.hi()[2];
  const Real d_eff = 2.0e-11;
  const Real decay = 5.0e-5;
  const Real transfer_length = 100.0e-6;
  const auto table = robin::build_table(
      system.adv, z_lo, z_hi, 4.0e-11, d_eff, decay, transfer_length,
      200.0e-6);
  const robin::TableView view{
      table.values.data(), table.z_lo, table.height, table.cutoff};
  Real maximum = 0.0;
  for (const Real source_fraction : {0.13, 0.41, 0.77}) {
    for (const Real target_fraction : {0.21, 0.58, 0.92}) {
      for (const Real rho : {1.0e-6, 7.0e-6, 31.0e-6, 87.0e-6}) {
        const Real source_z = z_lo + source_fraction * (z_hi - z_lo);
        const Real target_z = z_lo + target_fraction * (z_hi - z_lo);
        const Real interpolated = robin::interpolate(
            view, source_z, target_z, rho);
        const Real sealed = robin::normalized_robin_field(
            source_z, target_z, rho, z_lo, z_hi, d_eff, 4.0e-11, decay,
            std::numeric_limits<Real>::infinity(), 0.0, 0.0, 0.0,
            robin::kTableModeCount);
        const Real direct = robin::normalized_robin_field(
            source_z, target_z, rho, z_lo, z_hi, d_eff, 4.0e-11, decay,
            transfer_length, 0.0, 0.0, 0.0, robin::kTableModeCount);
        const Real reconstructed = sealed + interpolated;
        const Real relative = std::abs(reconstructed - direct)
            / std::max(std::abs(direct), 1.0e-30);
        maximum = std::max(maximum, relative);
      }
    }
  }
  require(maximum < 5.0e-3,
          "Robin correction interpolation exceeded specification tolerance");
  std::cout << "  test_table_against_direct_modes: PASSED (max relative error="
            << maximum << ")\n";
}

void test_sealed_limit() {
  auto system = make_system();
  const auto sealed = params_for(
      2.0e-11, 1.0e-4, std::numeric_limits<Real>::infinity());
  const auto zero_transfer = params_for(2.0e-11, 1.0e-4, 0.0);
  const Vec3 source = {500.0e-6, 500.0e-6, 25.0e-6};
  const Vec3 target = {530.0e-6, 500.0e-6, 75.0e-6};
  const Real bounded = system.gf.concentration_bounded(
      source, target, sealed);
  const Real through_robin_path = system.gf.concentration_bounded(
      source, target, zero_transfer);
  require(std::abs(bounded - through_robin_path)
              <= 1.0e-12 * std::max(std::abs(bounded), 1.0e-30),
          "Robin sealed limit differs from sealed Green's function");
  const Real correction = robin::normalized_correction(
      source[2], target[2], 30.0e-6, 0.0, 100.0e-6, 2.0e-11, 4.0e-11,
      1.0e-4, 0.0, 0.0, 0.0, 0.0, robin::kTableModeCount);
  require(correction == 0.0,
          "zero transfer length did not produce zero correction");
  std::cout << "  test_sealed_limit: PASSED\n";
}

void test_sink_limit() {
  auto system = make_system();
  const Real d_eff = 2.0e-11;
  const Real transfer_length = 2.0e-8;
  const auto params = params_for(d_eff, 1.0e-4, transfer_length);
  const Vec3 source = {500.0e-6, 500.0e-6, 50.0e-6};
  const Vec3 wall = {520.0e-6, 500.0e-6, 100.0e-6};
  const Vec3 middle = {520.0e-6, 500.0e-6, 50.0e-6};
  const Real ratio = system.gf.concentration_bounded(source, wall, params)
      / system.gf.concentration_bounded(source, middle, params);
  require(ratio < 1.0e-3, "Robin sink limit did not suppress lumen field");
  std::cout << "  test_sink_limit: PASSED (wall/mid ratio=" << ratio << ")\n";
}

void test_cross_language_anchors() {
  const Real z_lo = 0.0;
  const Real z_hi = 100.0e-6;
  const Real d_eff = 2.0e-11;
  const Real decay = 5.0e-5;
  const Real transfer_length = 100.0e-6;
  struct Anchor {
    Real source_fraction;
    Real target_fraction;
    Real rho;
    Real field;
    Real correction;
  };
  const std::array<Anchor, 5> anchors = {{
      {0.25, 0.25, 2.0e-6, 5.1061945263e5, -3.3149597263e4},
      {0.25, 0.75, 10.0e-6, 1.9964906375e4, -3.7235028089e4},
      {0.50, 0.50, 20.0e-6, 5.0272543392e4, -3.6163301346e4},
      {0.90, 0.99, 5.0e-6, 1.2537873736e5, -7.7533465488e4},
      {0.10, 0.90, 40.0e-6, 1.0721327524e4, -3.7660834321e4},
  }};
  Real maximum = 0.0;
  for (const Anchor& anchor : anchors) {
    const Real source_z = anchor.source_fraction * z_hi;
    const Real target_z = anchor.target_fraction * z_hi;
    const Real field = robin::normalized_robin_field(
        source_z, target_z, anchor.rho, z_lo, z_hi, d_eff, 4.0e-11, decay,
        transfer_length, 0.0, 0.0, 0.0, 3000);
    const Real correction = robin::normalized_correction(
        source_z, target_z, anchor.rho, z_lo, z_hi, d_eff, 4.0e-11, decay,
        transfer_length, 0.0, 0.0, 0.0, 3000);
    maximum = std::max(maximum, std::abs(field - anchor.field)
        / anchor.field);
    maximum = std::max(maximum, std::abs(correction - anchor.correction)
        / std::abs(anchor.correction));
  }
  require(maximum < 1.0e-6,
          "Robin cross-language anchor discrepancy exceeded tolerance");
  std::cout << "  test_cross_language_anchors: PASSED (max relative error="
            << maximum << ")\n";
}

void test_shipped_screening_sealed_series() {
  auto system = make_system();
  const Real d_eff = 2.0e-11;
  const Real source_rate = 1.0e-18;
  const Real decay = 5.0e-5;  // sqrt(lambda/D_eff)*H = 0.158114
  const auto params = params_for(
      d_eff, decay, std::numeric_limits<Real>::infinity());
  const Vec3 source = {500.0e-6, 500.0e-6, 25.0e-6};
  const Vec3 target = {520.0e-6, 500.0e-6, 75.0e-6};
  const Real normalized_image = system.gf.concentration_bounded(
      source, target, params) * (4.0 * std::numbers::pi * d_eff)
      / source_rate;
  const Real normalized_modes = robin::normalized_robin_field(
      source[2], target[2], 20.0e-6, 0.0, 100.0e-6, d_eff, 4.0e-11,
      decay, std::numeric_limits<Real>::infinity(), 0.0, 0.0, 0.0, 2048);
  const Real relative = std::abs(normalized_image - normalized_modes)
      / std::abs(normalized_modes);
  require(relative < 1.0e-6,
          "shipped-screening sealed series is under-resolved");
  std::cout << "  test_shipped_screening_sealed_series: PASSED "
            << "(kH=0.158114, relative error=" << relative << ")\n";
}

}  // namespace

int main() {
  std::cout << "=== Independent Robin Lumen-Boundary Tests ===\n";
  test_flux_residual();
  test_table_against_direct_modes();
  test_sealed_limit();
  test_sink_limit();
  test_cross_language_anchors();
  test_shipped_screening_sealed_series();
  std::cout << "All independent Robin lumen-boundary tests passed.\n";
  return 0;
}
