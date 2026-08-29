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

TestSystem make_system(bool shipped_flow = false) {
  TestSystem system;
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {1.0e-3, 1.0e-3, 100.0e-6};
  dcfg.grid_dx = 5.0e-6;
  system.domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = shipped_flow ? 5400.0 : 1.0e20;
  acfg.distal_transit_time = shipped_flow ? 43200.0 : 1.0e20;
  acfg.distal_length = shipped_flow ? 0.05 : 1.0e-3;
  acfg.mucus_thickness = 100.0e-6;
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
  const Real transfer_length = 100.0e-6;
  const Real decay = 5.0e-5;
  const Real kc = d_eff / transfer_length;
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
  const auto zero_table = robin::build_table(
      system.adv, 0.0, 100.0e-6, 4.0e-11, 2.0e-11, 1.0e-4, 0.0,
      200.0e-6);
  require(std::all_of(zero_table.values.begin(), zero_table.values.end(),
                      [](const Real value) { return value == 0.0; }),
          "zero transfer length produced a nonzero correction table");
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
  };
  const std::array<Anchor, 5> anchors = {{
      {0.25, 0.25, 2.0e-6, 5.048464448396090e5},
      {0.25, 0.75, 10.0e-6, 7.083746627572486e3},
      {0.50, 0.50, 20.0e-6, 5.111459366783190e4},
      {0.90, 0.99, 5.0e-6, 1.227298328345537e5},
      {0.10, 0.90, 40.0e-6, 4.602243835860912e3},
  }};
  const Real flow_x = 0.05 / 43200.0;
  const Real flow_z = 100.0e-6 / 5400.0;
  Real maximum = 0.0;
  for (const Anchor& anchor : anchors) {
    const Real source_z = anchor.source_fraction * z_hi;
    const Real target_z = anchor.target_fraction * z_hi;
    const Real field = robin::normalized_robin_field(
        source_z, target_z, anchor.rho, z_lo, z_hi, d_eff, 4.0e-11, decay,
        transfer_length, flow_x, 0.0, flow_z, 3000,
        robin::TransferBasis::Free);
    maximum = std::max(maximum, std::abs(field - anchor.field)
        / anchor.field);
    std::cout << "  anchor z_s=" << anchor.source_fraction
              << " z_t=" << anchor.target_fraction
              << " rho=" << anchor.rho
              << " field=" << field
              << " field_rel=" << std::abs(field - anchor.field)
                  / anchor.field
              << "\n";
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

void test_shipped_flow_residuals() {
  auto system = make_system(true);
  const Real height = system.domain.size()[2];
  const Real d_eff = 2.0e-11;
  const Real transfer_length = 100.0e-6;
  const Real kc = d_eff / transfer_length;
  const Real step = 1.0e-9;
  Real maximum = 0.0;
  const auto params = params_for(d_eff, 5.0e-5, transfer_length);
  for (const Real source_fraction : {0.25, 0.4}) {
    const Vec3 source = {500.0e-6, 500.0e-6, source_fraction * height};
    const Real rho = 2.0e-6;
    for (const int wall : {0, 1}) {
      const Real wall_z = wall == 0 ? 0.0 : height;
      const Real interior_z = wall == 0 ? step : height - step;
      const Vec3 wall_target = {source[0] + rho, source[1], wall_z};
      const Vec3 interior_target = {
          source[0] + rho, source[1], interior_z};
      const Real wall_value = system.gf.concentration_bounded(
          source, wall_target, params);
      const Real interior_value = system.gf.concentration_bounded(
          source, interior_target, params);
      const Real derivative = wall == 0
          ? (interior_value - wall_value) / step
          : (wall_value - interior_value) / step;
      const Real flow_z = system.adv.velocity(source)[2];
      const Real residual = wall == 0
          ? (-d_eff * derivative + flow_z * wall_value) * height
              / (d_eff * std::max(wall_value, 1.0e-30))
          : (-d_eff * derivative - kc * wall_value) * height
              / (d_eff * std::max(wall_value, 1.0e-30));
      std::cout << "  shipped wall residual source=" << source_fraction
                << " wall=" << wall << " value=" << wall_value
                << " residual=" << residual << "\n";
      maximum = std::max(maximum, std::abs(residual));
    }
  }
  std::cout << "  test_shipped_flow_residuals: max residual="
            << maximum << "\n";
  require(maximum <= 5.0e-3,
          "shipped-flow wall residual exceeded tolerance");
}

void test_shipped_flow_table() {
  auto system = make_system(true);
  const Real height = system.domain.size()[2];
  const Real d_eff = 2.0e-11;
  const Real transfer_length = 100.0e-6;
  const Real source_rate = 1.0e-18;
  const auto table = robin::build_table(
      system.adv, 0.0, height, 4.0e-11, d_eff, 5.0e-5,
      transfer_length, 200.0e-6);
  const robin::TableView view{
      table.values.data(), table.z_lo, table.height, table.cutoff};
  Real maximum = 0.0;
  for (const Real source_fraction : {0.25, 0.4}) {
    const Vec3 source = {500.0e-6, 500.0e-6, source_fraction * height};
    const Vec3 flow = system.adv.velocity(source);
    const Real source_z = source[2];
    for (const Real target_fraction : {0.95, 0.98}) {
      for (const Real rho : {0.5e-6, 2.0e-6, 5.0e-6}) {
        const Real target_z = target_fraction * height;
        const Vec3 target = {source[0] + rho, source[1], target_z};
        const Real base = system.gf.concentration_bounded(
            source, target,
            params_for(d_eff, 5.0e-5,
                       std::numeric_limits<Real>::infinity()));
        const Real correction = robin::interpolate(
            view, source_z, target_z, rho) * std::exp(
                (flow[0] * rho + flow[2] * (target_z - source_z))
                / (2.0 * d_eff));
        const Real reconstructed = base + source_rate
            / (4.0 * std::numbers::pi * d_eff) * correction;
        const Real direct = robin::normalized_robin_field(
            source_z, target_z, rho, 0.0, height, d_eff, 4.0e-11,
            5.0e-5, transfer_length, flow[0], flow[1], flow[2],
            robin::kTableModeCount);
        const Real direct_field = source_rate
            / (4.0 * std::numbers::pi * d_eff) * direct;
        require(reconstructed >= 0.0,
                "shipped-flow reconstructed field became negative");
        const Real relative = std::abs(reconstructed - direct_field)
            / std::max(std::abs(direct_field), 1.0e-30);
        maximum = std::max(maximum, relative);
      }
    }
  }
  std::cout << "  test_shipped_flow_table: max relative error="
            << maximum << "\n";
  require(maximum <= 5.0e-3,
          "shipped-flow table interpolation exceeded tolerance");
}

void test_basis_and_cache() {
  auto system = make_system();
  const Real effective_bi = robin::robin_biot_number(
      4.0e-11, 2.0e-11, 100.0e-6, 100.0e-6,
      robin::TransferBasis::Effective);
  const Real free_bi = robin::robin_biot_number(
      4.0e-11, 2.0e-11, 100.0e-6, 100.0e-6,
      robin::TransferBasis::Free);
  require(std::abs(effective_bi - 1.0) < 1.0e-12,
          "effective transfer basis Bi mismatch");
  require(std::abs(free_bi - 2.0) < 1.0e-12,
          "free transfer basis Bi mismatch");
  auto& cache = robin::global_table_cache();
  const uint64_t built_before = cache.tables_built();
  for (int i = 0; i < 200; ++i) {
    const Real retardation = 2.0 + static_cast<Real>(i) * 1.0e-6;
    cache.get(system.adv, 0.0, 100.0e-6, 4.0e-11,
              4.0e-11 / retardation, 5.0e-5, 100.0e-6,
              200.0e-6, robin::TransferBasis::Effective);
  }
  const uint64_t built_delta = cache.tables_built() - built_before;
  std::cout << "  test_basis_and_cache: cache size=" << cache.size()
            << " tables built=" << cache.tables_built()
            << " evictions=" << cache.table_evictions() << "\n";
  require(cache.size() <= robin::TableCache::kMaximumTables,
          "Robin table cache exceeded capacity");
  require(built_delta < 200,
          "Robin cache built one table per continuous retardation value");
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
  test_shipped_flow_residuals();
  test_shipped_flow_table();
  test_basis_and_cache();
  std::cout << "All independent Robin lumen-boundary tests passed.\n";
  return 0;
}
