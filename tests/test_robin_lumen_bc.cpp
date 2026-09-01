/* -----------------------------------------------------------------------
   GutIBM – Independent Robin lumen-boundary tests
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "domain.h"
#include "greens_function.h"
#include "greens_function_gpu.h"
#include "input_parser.h"
#include "config_json.h"
#include "error.h"
#include "plasmid.h"
#include "robin_correction_table.h"
#include "simulation.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <string>

using namespace gutibm;

namespace {

struct TestSystem {
  Domain domain;
  AdvectionField adv;
  GreensFunction gf;
};

TestSystem make_system(bool shipped_flow = false, bool peristaltic = false) {
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
  acfg.peristaltic_enabled = peristaltic;
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

constexpr uint64_t kFnvOffset = 14695981039346656037ULL;
constexpr uint64_t kFnvPrime = 1099511628211ULL;

void append_hash_bytes(uint64_t& hash, const void* data, size_t size) {
  const auto* bytes = static_cast<const unsigned char*>(data);
  for (size_t index = 0; index < size; ++index) {
    hash ^= bytes[index];
    hash *= kFnvPrime;
  }
}

uint64_t table_identity_hash(const robin::Table& table) {
  uint64_t hash = kFnvOffset;
  for (const int64_t group : table.quantized_key) {
    append_hash_bytes(hash, &group, sizeof(group));
  }
  const int basis = static_cast<int>(table.basis);
  append_hash_bytes(hash, &basis, sizeof(basis));
  append_hash_bytes(hash, &table.z_lo, sizeof(table.z_lo));
  append_hash_bytes(hash, &table.height, sizeof(table.height));
  append_hash_bytes(hash, &table.cutoff, sizeof(table.cutoff));
  for (const double value : table.legacy_values) {
    append_hash_bytes(hash, &value, sizeof(value));
  }
  return hash;
}

SimulationConfig provenance_test_config() {
  SimulationConfig config = InputParser::default_config();
  config.domain.hi = {50.0e-6, 50.0e-6, 25.0e-6};
  config.domain.grid_dx = 5.0e-6;
  config.hdf5.enabled = false;
  config.enabled_fixes.clear();
  config.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 1;
  strain.mu_max = 5.0e-4;
  config.initial_strains.push_back(strain);
  return config;
}

std::shared_ptr<const robin::Table> build_identity_table(
    const TestSystem& system, Real transfer_length) {
  return robin::global_table_cache().get(
      system.adv, 0.0, 100.0e-6, 4.0e-11, 2.0e-11, 5.0e-5,
      transfer_length, 173.0e-6, robin::TransferBasis::Effective);
}

void test_run_scoped_identity_and_eviction() {
  auto first_sim = provenance_test_config();
  Simulation first;
  first.init(first_sim);
  const auto system = make_system();
  const auto& cache = robin::global_table_cache();
  const uint64_t built_before = cache.tables_built();
  const uint64_t evictions_before = cache.table_evictions();
  const uint64_t identity_before = cache.built_identity();

  constexpr Real set_base = 7.123e-6;
  std::array<std::shared_ptr<const robin::Table>, 3> set_tables;
  uint64_t expected_identity = 0;
  for (size_t index = 0; index < set_tables.size(); ++index) {
    const Real transfer_length = set_base
        * std::pow(1.05, static_cast<Real>(index));
    set_tables[index] = build_identity_table(system, transfer_length);
    expected_identity ^= table_identity_hash(*set_tables[index]);
  }
  const uint64_t set_identity = expected_identity;
  const std::string first_hash = first.robin_table_hash();
  require(first.robin_tables_built() == set_tables.size(),
          "first Robin provenance pass did not build exactly three tables");
  require(first_hash == robin::format_identity_hash(set_identity),
          "Robin run identity did not match the first table set");

  constexpr int filler_count = 64;
  for (int index = 0; index < filler_count; ++index) {
    const Real transfer_length = 1.234e-6
        * std::pow(1.05, static_cast<Real>(index));
    const auto filler = build_identity_table(system, transfer_length);
    expected_identity ^= table_identity_hash(*filler);
  }
  require(cache.table_evictions() > evictions_before,
          "Robin cache eviction was not exercised");

  SimulationConfig second_config = provenance_test_config();
  Simulation second;
  second.init(second_config);
  uint64_t reverse_identity = 0;
  for (size_t index = set_tables.size(); index > 0; --index) {
    const Real transfer_length = set_base
        * std::pow(1.05, static_cast<Real>(index - 1));
    const auto table = build_identity_table(system, transfer_length);
    reverse_identity ^= table_identity_hash(*table);
  }

  require(second.robin_tables_built() == set_tables.size(),
          "second Robin provenance pass did not rebuild exactly three tables");
  require(first_hash == second.robin_table_hash(),
          "Robin run identity depended on insertion order");
  expected_identity ^= reverse_identity;
  require((cache.built_identity() ^ identity_before) == expected_identity,
          "Robin identity omitted an evicted table");
  require(cache.tables_built() - built_before == 3 + filler_count + 3,
          "Robin table build count did not match the eviction test workload");
  require(cache.tables_built() > second.robin_tables_built(),
          "Robin run-scoped count was not isolated from process totals");
  std::cout << "  test_run_scoped_identity_and_eviction: 70 builds, "
            << cache.table_evictions() - evictions_before
            << " evictions, order-independent identity passed\n";
}

void test_launch_local_table_mapping() {
  std::vector<std::shared_ptr<const robin::Table>> source_tables;
  source_tables.reserve(128);
  for (int i = 0; i < 128; ++i) {
    auto table = std::make_shared<robin::Table>();
    table->quantized_key[0] = i;
    source_tables.push_back(std::move(table));
  }
  source_tables.push_back(source_tables[7]);
  std::vector<std::shared_ptr<const robin::Table>> launch_tables;
  const auto indices = make_robin_launch_table_indices(
      source_tables, launch_tables);
  require(launch_tables.size() == 128,
          "launch-local Robin table set lost a distinct identity");
  require(indices.size() == source_tables.size(),
          "launch-local Robin index count does not match sources");
  for (size_t source = 0; source < source_tables.size(); ++source) {
    const int index = indices[source];
    require(index >= 0
                && static_cast<size_t>(index) < launch_tables.size(),
            "launch-local Robin index is out of bounds");
    require(launch_tables[static_cast<size_t>(index)].get()
                == source_tables[source].get(),
            "launch-local Robin index resolved to the wrong table");
  }
  source_tables.clear();
  for (int i = 0; i < 257; ++i) {
    auto table = std::make_shared<robin::Table>();
    table->quantized_key[0] = i;
    source_tables.push_back(std::move(table));
  }
  launch_tables.clear();
  bool rejected_overflow = false;
  try {
    (void)make_robin_launch_table_indices(source_tables, launch_tables);
  } catch (const SimulationError&) {
    rejected_overflow = true;
  }
  require(rejected_overflow,
          "launch-local Robin table cap did not reject 257 identities");
  std::cout << "  test_launch_local_table_mapping: 128 identities passed\n";
}

void test_robin_fallback_preflight() {
  const auto system = make_system();
  const std::vector<Vec3> sources = {
      {500.0e-6, 500.0e-6, 4.999e-6},
      {500.0e-6, 500.0e-6, 5.001e-6},
      {500.0e-6, 500.0e-6, 50.0e-6}};
  auto params = params_for(2.0e-11, 5.0e-5, 100.0e-6);
  const std::vector enabled(sources.size(), params);
  auto disabled = enabled;
  for (auto& item : disabled) {
    item.lumen_transfer_length = std::numeric_limits<Real>::infinity();
  }
  const auto fallback = robin_host_fallback_sources(
      system.domain, sources, enabled);
  require(fallback.size() == 1 && fallback.front() == 0,
          "Robin fallback preflight did not honor the cell-radius threshold");
  require(robin_host_fallback_sources(system.domain, sources, disabled).empty(),
          "Robin fallback preflight changed the disabled path");
}

void test_disabled_default_is_inert() {
  const SimulationConfig config = InputParser::default_config();
  require(!robin::transfer_enabled(config.qssa.lumen_transfer_length),
          "Robin transfer must be disabled by default");
  const std::string resolved = ConfigJson::serialize_document(config);
  require(resolved.find("\"toxin.lumen_transfer_length\":\"inf\"")
              != std::string::npos,
          "resolved provenance must serialize disabled Robin transfer");

  auto system = make_system(true);
  const auto params = params_for(2.0e-11, 5.0e-5,
                                 config.qssa.lumen_transfer_length);
  const auto explicit_sealed = params_for(2.0e-11, 5.0e-5, 0.0);
  for (const Real source_fraction : {0.2, 0.5, 0.8}) {
    const Vec3 source = {500.0e-6, 500.0e-6,
                         source_fraction * system.domain.size()[2]};
    for (const Real target_fraction : {0.1, 0.5, 0.9}) {
      const Vec3 target = {530.0e-6, 500.0e-6,
                           target_fraction * system.domain.size()[2]};
      const Real bounded = system.gf.concentration_bounded(
          source, target, params);
      const Real sealed = system.gf.concentration_bounded(
          source, target, explicit_sealed);
      require(std::abs(bounded - sealed) <=
                  1.0e-12 * std::max(std::abs(sealed), 1.0e-30),
              "disabled Robin default changed sealed concentration");
    }
  }
  std::cout << "  test_disabled_default_is_inert: boundary_status=disabled\n";
}

void test_enabled_boundary_impact() {
  auto system = make_system(true);
  const Real height = system.domain.size()[2];
  const Vec3 source = {500.0e-6, 500.0e-6, 0.4 * height};
  const Vec3 near_lumen = {520.0e-6, 500.0e-6, 0.98 * height};
  const Vec3 mid_slab = {520.0e-6, 500.0e-6, 0.5 * height};
  const auto disabled = params_for(
      2.0e-11, 5.0e-5, robin::kZeroTransferLength);
  const auto effective = params_for(2.0e-11, 5.0e-5, 100.0e-6);
  auto free = effective;
  free.lumen_transfer_basis_free = true;
  const Real disabled_near =
      system.gf.concentration_bounded(source, near_lumen, disabled);
  const Real disabled_mid =
      system.gf.concentration_bounded(source, mid_slab, disabled);
  const Real effective_near =
      system.gf.concentration_bounded(source, near_lumen, effective);
  const Real effective_mid =
      system.gf.concentration_bounded(source, mid_slab, effective);
  const Real free_near =
      system.gf.concentration_bounded(source, near_lumen, free);
  const Real free_mid =
      system.gf.concentration_bounded(source, mid_slab, free);
  std::cout << "  enabled impact near-lumen (effective/disabled)="
            << effective_near / disabled_near
            << " mid-slab=" << effective_mid / disabled_mid << "\n";
  std::cout << "  enabled impact near-lumen (free/disabled)="
            << free_near / disabled_near
            << " mid-slab=" << free_mid / disabled_mid << "\n";
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
      table.legacy_values.data(), table.z_lo, table.height, table.cutoff};
  Real maximum = 0.0;
  for (const Real source_fraction : {0.13, 0.41, 0.77}) {
    for (const Real target_fraction : {0.21, 0.58, 0.92}) {
      for (const Real rho : {1.0e-6, 7.0e-6, 31.0e-6, 87.0e-6}) {
        const Real source_z = z_lo + source_fraction * (z_hi - z_lo);
        const Real target_z = z_lo + target_fraction * (z_hi - z_lo);
        const Real interpolated = robin::interpolate_uniform(
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
  const Vec3 source = {500.0e-6, 500.0e-6, 25.0e-6};
  const Vec3 near_lumen = {530.0e-6, 500.0e-6, 99.0e-6};
  const Real sealed_value = system.gf.concentration_bounded(
      source, near_lumen, sealed);

  for (const std::string& invalid : {"0", "-1e-6", "nan", "-inf"}) {
    SimulationConfig config = InputParser::default_config();
    bool rejected = false;
    try {
      InputParser::apply_flat_key(
          config, "toxin.lumen_transfer_length", invalid);
    } catch (const ConfigError&) {
      rejected = true;
    }
    require(rejected, "invalid Robin transfer length was accepted");
  }
  for (const std::string& invalid : {"0", "nan", "-inf"}) {
    SimulationConfig config = InputParser::default_config();
    bool rejected = false;
    try {
      InputParser::apply_flat_key(config, "toxin_cutoff", invalid);
    } catch (const ConfigError&) {
      rejected = true;
    }
    require(rejected, "invalid Robin cutoff was accepted");
  }
  SimulationConfig config = InputParser::default_config();
  require(InputParser::apply_flat_key(
              config, "toxin.lumen_transfer_length", "inf"),
          "positive infinity transfer length was not accepted");
  require(!robin::transfer_enabled(config.qssa.lumen_transfer_length),
          "positive infinity transfer length did not disable Robin");

  const auto broad_layer = params_for(2.0e-11, 1.0e-4, 1.0e-6);
  const auto thin_layer = params_for(2.0e-11, 1.0e-4, 1.0e-8);
  const Real broad_value = system.gf.concentration_bounded(
      source, near_lumen, broad_layer);
  const Real thin_value = system.gf.concentration_bounded(
      source, near_lumen, thin_layer);
  require(broad_value < sealed_value,
          "finite Robin transfer did not reduce near-lumen concentration");
  require(thin_value < broad_value,
          "near-lumen concentration did not decrease as transfer length shrank");
  std::cout << "  test_sealed_limit: PASSED (sealed=" << sealed_value
            << ", delta=1e-6=" << broad_value
            << ", delta=1e-8=" << thin_value << ")\n";
}

void test_robin_mode_residuals() {
  constexpr Real height = 100.0e-6;
  constexpr Real transfer_length = 100.0e-6;
  const auto system = make_system(true);
  const std::array<Real, 3> source_fractions = {0.02, 0.5, 0.98};
  const std::array<robin::TransferBasis, 2> bases = {
      robin::TransferBasis::Effective, robin::TransferBasis::Free};
  for (const auto& entry : PlasmidLibrary::entries()) {
    for (const auto basis : bases) {
      for (const Real source_fraction : source_fractions) {
        const Real z = source_fraction * height;
        const Vec3 flow = system.adv.mean_velocity({0.0, 0.0, z});
        const Real d_eff = entry.cluster.diff_coeff
            / entry.cluster.retardation;
        const Real a = flow[2] / (2.0 * d_eff);
        const Real bi = robin::robin_biot_number(
            entry.cluster.diff_coeff, d_eff, height, transfer_length, basis);
        const Real c_lo = -a;
        const Real c_hi = bi / height + a;
        const auto roots = robin::robin_mode_roots(
            height, c_lo, c_hi, 64);
        require(roots.size() == 64,
                "Robin residual sweep returned the wrong mode count");
        Real previous = -1.0;
        const Real p = c_lo * c_hi * height * height;
        const Real q = (c_hi - c_lo) * height;
        for (const Real beta : roots) {
          const Real x = beta * height;
          require(x > previous,
                  "Robin residual sweep returned duplicate or unordered roots");
          previous = x;
          const Real residual = (x * x + p) * std::sin(x)
              - q * x * std::cos(x);
          const Real scale = (x * x + std::abs(p)) + std::abs(q) * x;
          require(std::abs(residual) / scale <= 1.0e-9,
                  "Robin residual sweep exceeded scaled residual tolerance");
        }
      }
    }
  }
  std::cout << "  test_robin_mode_residuals: PASSED (6 plasmids, 2 bases, "
               "3 source heights, 64 modes)\n";
}

void test_colE1_pole_regression() {
  constexpr Real height = 100.0e-6;
  constexpr Real a_height = 1.1626094919444088;
  const Real c_lo = -a_height / height;
  const Real c_hi = (1.0 + a_height) / height;
  const auto roots = robin::robin_mode_roots(height, c_lo, c_hi, 5);
  const Real first_x = roots.front() * height;
  const Real expected = 1.5763832766320491;
  const Real pole = 1.5856450809382423;
  require(std::abs(first_x - expected) / expected <= 1.0e-9,
          "ColE1 first Robin root does not match the pole-free reference");
  require(std::abs(first_x - pole) / pole >= 1.0e-6,
          "ColE1 first Robin root converged to the rational pole");
  std::cout << "  test_colE1_pole_regression: PASSED (first beta*H="
            << first_x << ")\n";
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

void test_shipped_flow_direct_boundary() {
  auto system = make_system(true);
  const Real height = system.domain.size()[2];
  const Real d_eff = 2.0e-11;
  const Real transfer_length = 100.0e-6;
  const Real kc = d_eff / transfer_length;
  const Real step = 1.0e-9;
  Real maximum_epithelial = 0.0;
  Real maximum_lumen = 0.0;
  for (const Real source_fraction : {0.25, 0.4, 0.6, 0.75}) {
    const Vec3 source = {500.0e-6, 500.0e-6, source_fraction * height};
    const Real rho = 2.0e-6;
    const Vec3 flow = system.adv.velocity(source);
    for (const int wall : {0, 1}) {
      const Real wall_z = wall == 0 ? 0.0 : height;
      const Real interior_z = wall == 0 ? step : height - step;
      const Real wall_value = robin::normalized_robin_field(
          source[2], wall_z, rho, 0.0, height, d_eff, 4.0e-11,
          5.0e-5, transfer_length, flow[0], flow[1], flow[2],
          robin::kTableModeCount);
      const Real interior_value = robin::normalized_robin_field(
          source[2], interior_z, rho, 0.0, height, d_eff, 4.0e-11,
          5.0e-5, transfer_length, flow[0], flow[1], flow[2],
          robin::kTableModeCount);
      const Real derivative = wall == 0
          ? (interior_value - wall_value) / step
          : (wall_value - interior_value) / step;
      const Real flow_z = flow[2];
      const Real residual = wall == 0
          ? (-d_eff * derivative + flow_z * wall_value) * height
              / (d_eff * std::max(wall_value, 1.0e-30))
          : (-d_eff * derivative - kc * wall_value) * height
              / (d_eff * std::max(wall_value, 1.0e-30));
      std::cout << "  direct boundary residual source=" << source_fraction
                << " wall=" << wall << " value=" << wall_value
                << " residual=" << residual << "\n";
      if (wall == 0) {
        maximum_epithelial =
            std::max(maximum_epithelial, std::abs(residual));
      } else {
        maximum_lumen = std::max(maximum_lumen, std::abs(residual));
      }
    }
  }
  std::cout << "  test_shipped_flow_direct_boundary: max epithelial="
            << maximum_epithelial << " max lumen=" << maximum_lumen
            << "\n";
  require(maximum_epithelial <= 5.0e-3,
          "direct epithelial total-flux residual exceeded tolerance");
  require(maximum_lumen <= 1.0e-3,
          "direct lumen flux residual exceeded tolerance");
}

void test_shipped_flow_reconstruction() {
  auto system = make_system(true);
  const Real height = system.domain.size()[2];
  const Real d_eff = 2.0e-11;
  const Real transfer_length = 100.0e-6;
  const Real source_rate = 1.0e-18;
  const auto table = robin::build_table(
      system.adv, 0.0, height, 4.0e-11, d_eff, 5.0e-5,
      transfer_length, 200.0e-6);
  const robin::TableView view{
      table.legacy_values.data(), table.z_lo, table.height, table.cutoff};
  Real maximum = 0.0;
  for (const Real source_fraction : {0.25, 0.4, 0.6, 0.75}) {
    const Vec3 source = {500.0e-6, 500.0e-6, source_fraction * height};
    const Vec3 flow = system.adv.velocity(source);
    const Real source_z = source[2];
    for (const Real target_fraction : {0.95, 0.98, 0.995}) {
      for (const Real rho : {0.5e-6, 2.0e-6, 5.0e-6}) {
        const Real target_z = target_fraction * height;
        const Vec3 target = {source[0] + rho, source[1], target_z};
        const Real base = system.gf.concentration_bounded(
            source, target,
            params_for(d_eff, 5.0e-5,
                       std::numeric_limits<Real>::infinity()));
        const Real correction = robin::interpolate_uniform(
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
        std::cout << "  shipped reconstruction source=" << source_fraction
                  << " target=" << target_fraction << " rho=" << rho
                  << " reconstructed=" << reconstructed
                  << " direct=" << direct_field
                  << " relative=" << relative << "\n";
        maximum = std::max(maximum, relative);
      }
    }
  }
  std::cout << "  test_shipped_flow_reconstruction: max relative error="
            << maximum << "\n";
  require(maximum <= 1.0e-2,
          "shipped-flow table interpolation exceeded tolerance");
}

void test_shipped_flow_interpolated_wall_guard() {
  auto system = make_system(true);
  const Real height = system.domain.size()[2];
  const Real d_eff = 2.0e-11;
  const Real transfer_length = 100.0e-6;
  const Real kc = d_eff / transfer_length;
  const Real step = 1.0e-9;
  Real maximum = 0.0;
  const auto params = params_for(d_eff, 5.0e-5, transfer_length);
  for (const Real source_fraction : {0.25, 0.4, 0.6, 0.75}) {
    const Vec3 source = {500.0e-6, 500.0e-6, source_fraction * height};
    for (const int wall : {0, 1}) {
      const Real wall_z = wall == 0 ? 0.0 : height;
      const Real interior_z = wall == 0 ? step : height - step;
      const Vec3 wall_target = {source[0] + 2.0e-6, source[1], wall_z};
      const Vec3 interior_target = {
          source[0] + 2.0e-6, source[1], interior_z};
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
      std::cout << "  interpolated wall residual source=" << source_fraction
                << " wall=" << wall << " value=" << wall_value
                << " residual=" << residual << "\n";
      maximum = std::max(maximum, std::abs(residual));
    }
  }
  std::cout << "  test_shipped_flow_interpolated_wall_guard: max residual="
            << maximum << "\n";
  // Regression guard only: the direct modal boundary residual is ~1e-4,
  // the inherited image-base inconsistency is ~1e-2, and the remainder is
  // trilinear gradient error inside a 3.1 um cell. Tightening this requires
  // a finer near-wall z grid or a direct-mode wall path, not a looser test.
  require(maximum <= 1.5e-1,
          "interpolated shipped-flow wall residual exceeded guard");
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

  const auto near_wall_params = params_for(
      2.0e-11, 5.0e-5, 100.0e-6);
  const Vec3 near_wall_source = {500.0e-6, 500.0e-6, 2.0e-6};
  const Vec3 near_wall_target = {501.0e-6, 500.0e-6, 3.0e-6};
  (void)system.gf.concentration_bounded(
      near_wall_source, near_wall_target, near_wall_params);
  std::cout << "  test_basis_and_cache: direct evaluations="
            << system.gf.robin_direct_evaluations() << "\n";
  require(system.gf.robin_direct_evaluations() > 0,
          "Robin direct fallback was not exercised");
}

void test_peristaltic_mean_profile() {
  auto system = make_system(true, true);
  const Real height = system.domain.size()[2];
  const Real d_eff = 2.0e-11;
  const Real source_rate = 1.0e-18;
  const Real rho = 2.0e-6;
  const Vec3 source = {500.0e-6, 500.0e-6, 75.0e-6};
  const Vec3 target = {source[0] + rho, source[1], 98.0e-6};
  const auto sealed_params = params_for(
      d_eff, 5.0e-5, std::numeric_limits<Real>::infinity());
  const std::array<Real, 3> factors = {0.5, 1.0, 1.5};
  const std::array<Real, 3> times = {15.0, 0.0, 5.0};
  std::array<Real, 3> corrections{};
  std::array<Real, 3> totals{};
  const robin::Table* first_table = nullptr;
  for (size_t i = 0; i < factors.size(); ++i) {
    system.adv.set_time(times[i]);
    require(std::abs(system.adv.peristaltic_factor(source) - factors[i])
                < 1.0e-12,
            "peristaltic test setup produced the wrong flow factor");
    const auto table = robin::global_table_cache().get(
        system.adv, 0.0, height, 4.0e-11, d_eff, 5.0e-5,
        100.0e-6, 200.0e-6, robin::TransferBasis::Effective);
    if (first_table == nullptr) {
      first_table = table.get();
    } else {
      require(table.get() == first_table,
              "peristaltic modulation changed the Robin cache key");
    }
    const robin::TableView view{
        table->legacy_values.data(), table->z_lo, table->height, table->cutoff};
    const Real correction_base = robin::interpolate_uniform(
        view, source[2], target[2], rho);
    const Vec3 flow = system.adv.velocity(source);
    corrections[i] = source_rate / (4.0 * std::numbers::pi * d_eff)
        * correction_base * std::exp(
            (flow[0] * rho + flow[2] * (target[2] - source[2]))
            / (2.0 * d_eff));
    totals[i] = system.gf.concentration_bounded(
        source, target, sealed_params) + corrections[i];
    require(std::isfinite(corrections[i]) && std::isfinite(totals[i]),
            "peristaltic correction measurement was not finite");
    require(totals[i] >= 0.0,
            "peristaltic reconstructed field became negative");
    std::cout << "  peristaltic factor=" << factors[i]
              << " correction=" << corrections[i]
              << " total=" << totals[i] << "\n";
  }
  Real maximum = 0.0;
  for (size_t i = 0; i < factors.size(); ++i) {
    const Real relative = std::abs(corrections[i] - corrections[1])
        / std::max(std::abs(totals[i]), 1.0e-30);
    maximum = std::max(maximum, relative);
  }
  std::cout << "  test_peristaltic_mean_profile: max relative neglected "
               "correction=" << maximum << "\n";
}

}  // namespace

int main() {
  std::cout << "=== Independent Robin Lumen-Boundary Tests ===\n";
  test_launch_local_table_mapping();
  test_robin_fallback_preflight();
  test_disabled_default_is_inert();
  test_enabled_boundary_impact();
  test_flux_residual();
  test_table_against_direct_modes();
  test_sealed_limit();
  test_robin_mode_residuals();
  test_colE1_pole_regression();
  test_sink_limit();
  test_cross_language_anchors();
  test_shipped_screening_sealed_series();
  test_shipped_flow_direct_boundary();
  test_shipped_flow_reconstruction();
  test_shipped_flow_interpolated_wall_guard();
  test_basis_and_cache();
  test_run_scoped_identity_and_eviction();
  test_peristaltic_mean_profile();
  std::cout << "All independent Robin lumen-boundary tests passed.\n";
  return 0;
}
