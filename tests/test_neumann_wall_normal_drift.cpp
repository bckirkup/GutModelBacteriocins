/* -----------------------------------------------------------------------
   GutIBM – Wall-normal drift validity-envelope tests
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "chemical_field.h"
#include "domain.h"
#include "error.h"
#include "greens_function.h"
#include "neumann_image_series.h"
#include "qssa_solver.h"
#include "robin_correction_table.h"
#include "species_names.h"
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numbers>
#include <string>
#include <string_view>
#include <vector>

using namespace gutibm;

namespace {

constexpr Real kHeight = 100.0e-6;
constexpr Real kDiffusion = 2.0e-11;
constexpr Real kDecay = 5.0e-5;
constexpr Real kSourceRate = 1.0e-18;
constexpr Real kFlowX = 1.0e-3 / 43200.0;

struct Probe {
  Real source_fraction;
  Real target_fraction;
  Real rho;
};

constexpr std::array<Probe, 9> kProbes = {{
    {0.40, 0.98, 5.0e-6},
    {0.40, 0.50, 5.0e-6},
    {0.10, 0.05, 2.0e-6},
    {0.50, 0.50, 20.0e-6},
    {0.90, 0.99, 5.0e-6},
    {0.05, 0.95, 10.0e-6},
    {0.20, 0.20, 50.0e-6},
    {0.02, 0.02, 2.0e-6},
    {0.30, 0.10, 10.0e-6},
}};

struct TestSystem {
  Domain domain;
  AdvectionField adv;
  GreensFunction gf;
};

TestSystem make_system(Real flow_z, Real flow_x = kFlowX,
                       Real grid_dx = 5.0e-6) {
  TestSystem system;
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {1.0e-3, 1.0e-3, kHeight};
  dcfg.grid_dx = grid_dx;
  system.domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.crypts_enabled = false;
  acfg.peristaltic_enabled = false;
  acfg.mucus_thickness = kHeight;
  acfg.profile_alpha = 0.0;
  acfg.distal_length = flow_x == 0.0 ? 0.0 : 1.0e-3;
  acfg.distal_transit_time = flow_x == 0.0 ? 1.0 : acfg.distal_length / flow_x;
  acfg.radial_turnover = flow_z == 0.0
      ? std::numeric_limits<Real>::infinity() : kHeight / flow_z;
  system.adv.init(acfg, system.domain);
  system.gf.init(system.domain, system.adv);
  return system;
}

GreensFunctionParams make_params() {
  GreensFunctionParams params;
  params.diff_coeff = kDiffusion;
  params.retardation = 1.0;
  params.source_rate = kSourceRate;
  params.decay_rate = kDecay;
  params.lumen_transfer_length = robin::kZeroTransferLength;
  params.image_series_relative_tolerance = 1.0e-14;
  params.image_series_max_shells = 4096;
  params.image_series_max_shells_explicit = true;
  params.image_series_legacy_reflections = false;
  return params;
}

Real relative_error(Real actual, Real reference) {
  return std::abs(actual - reference)
      / std::max(std::abs(reference), 1.0e-300);
}

void require(bool condition, const char* message) {
  if (!condition) {
    std::cerr << message << "\n";
    std::exit(1);
  }
}

Real modal_reference(const Probe& probe, Real flow_z, bool physical,
                     Real flow_x = kFlowX) {
  const robin::SealedFieldParams params{
      probe.source_fraction * kHeight, probe.target_fraction * kHeight,
      probe.rho, 0.0, kHeight, kDiffusion, kDecay, flow_x, 0.0, flow_z,
      4000};
  const Real normalized = physical
      ? robin::normalized_sealed_field(params)
      : robin::normalized_image_consistent_sealed_field(params);
  return kSourceRate / (4.0 * std::numbers::pi * kDiffusion) * normalized;
}

Real probe_error(const TestSystem& system,
                 const GreensFunctionParams& params,
                 const Probe& probe, Real flow_x, Real flow_z) {
  const Vec3 source = {500.0e-6, 500.0e-6,
                       probe.source_fraction * kHeight};
  const Vec3 target = {source[0] + probe.rho, source[1],
                       probe.target_fraction * kHeight};
  const Real image = system.gf.concentration_bounded(source, target, params);
  const Real reference = modal_reference(
      probe, flow_z, params.drift_correction, flow_x);
  if (!(std::isfinite(image) && image > 0.0
        && std::isfinite(reference) && reference > 0.0)) {
    std::cerr << "non-positive or non-finite sealed field in probe: flow_x="
              << flow_x << " flow_z=" << flow_z << "\n";
    std::exit(1);
  }
  return relative_error(image, reference);
}

Real worst_case_error(Real flow_x, Real flow_z) {
  auto system = make_system(flow_z, flow_x);
  const auto params = make_params();
  Real worst = 0.0;
  for (const Probe& probe : kProbes) {
    worst = std::max(
        worst, probe_error(system, params, probe, flow_x, flow_z));
  }
  return worst;
}

void test_zero_wall_normal_flow_is_exact() {
  const Real maximum = worst_case_error(kFlowX, 0.0);
  if (!(maximum < 1.0e-8)) {
    std::cerr << "zero-flow sealed image/modal mismatch: " << maximum << "\n";
    std::exit(1);
  }
  std::cout << "  test_zero_wall_normal_flow_is_exact: PASSED (max relative "
               "error=" << maximum << ")\n";
}

void test_wall_normal_flow_sweep_is_sensitive_and_ordered() {
  const std::array<Real, 5> pe_values = {0.0, 0.01, 0.0926, 0.37, 1.5};
  std::array<Real, 5> errors = {};
  for (std::size_t i = 0; i < pe_values.size(); ++i) {
    const Real pe_z = pe_values[i];
    const Real flow_z = pe_z * kDiffusion / kHeight;
    errors[i] = worst_case_error(kFlowX, flow_z);
    if (i > 0 && errors[i] < errors[i - 1]) {
      std::cerr << "wall-normal drift error was not non-decreasing: Pe_z="
                << pe_z << " previous=" << errors[i - 1]
                << " current=" << errors[i] << "\n";
      std::exit(1);
    }
    std::cout << "  Pe_z=" << pe_z << " worst relative error=" << errors[i]
              << "\n";
  }
  if (!(errors[0] <= 1.0e-8
        && errors[4] - errors[2] > 0.1
        && errors[2] >= 0.20 * pe_values[2]
        && errors[2] <= 0.60 * pe_values[2]
        && errors[3] >= 0.20 * pe_values[3]
        && errors[3] <= 0.60 * pe_values[3])) {
    std::cerr << "wall-normal drift sweep outside documented bounds\n";
    std::exit(1);
  }
  // The middle-range lower bounds pin the known documented defect, so a stub
  // or accidental series change cannot pass. An exact treatment must update
  // this test and §§4/6 of the drift document together.
  std::cout << "  test_wall_normal_flow_sweep_is_sensitive_and_ordered: PASSED\n";
}

void test_reflected_flow_reversal_is_detectable() {
  const Probe probe = {0.40, 0.98, 5.0e-6};
  bool valid = true;
  for (const Real pe_z : {0.01, 0.0926, 0.37}) {
    const Real flow_z = pe_z * kDiffusion / kHeight;
    auto system = make_system(flow_z);
    const auto params = make_params();
    const Real error = probe_error(system, params, probe, kFlowX, flow_z);
    if (!(error >= 0.020 * pe_z && error <= 0.070 * pe_z)) {
      std::cerr << "reversal-sensitive asymmetric probe outside bounds: Pe_z="
                << pe_z << " error=" << error << "\n";
      valid = false;
    } else {
      std::cout << "  Pe_z=" << pe_z
                << " asymmetric reversal-sensitive error=" << error << "\n";
    }
  }
  if (!valid) {
    std::exit(1);
  }
  // This asymmetric probe guards the reflected family specifically: the
  // symmetric worst-case probe cannot see image-flow reversal. The two-sided
  // bounds come from docs/NEUMANN_WALL_NORMAL_DRIFT.md §3; flipping
  // image_flow[2] back to same-sign in greens_function.cpp was verified to
  // break the upper bound at all three Pe_z values.
  std::cout << "  test_reflected_flow_reversal_is_detectable: PASSED\n";
}

void test_wall_parallel_flow_is_exact() {
  for (const Real flow_x : {0.0, kFlowX, 10.0 * kFlowX}) {
    const Real error = worst_case_error(flow_x, 0.0);
    if (!(error <= 1.0e-8)) {
      std::cerr << "wall-parallel flow lost sealed exactness: flow_x="
                << flow_x << " worst=" << error << "\n";
      std::exit(1);
    }
    std::cout << "  flow_x=" << flow_x << " wall-parallel worst relative error="
              << error << "\n";
  }
  // The defect is proportional to wall-normal flow and vanishes for
  // wall-parallel flow; this negative control guards that invariant.
  std::cout << "  test_wall_parallel_flow_is_exact: PASSED\n";
}

void test_drift_classifier() {
  if (const Real classified = neumann::wall_normal_peclet(
          0.0926 * kDiffusion / kHeight, kHeight, kDiffusion);
      std::abs(classified - 0.0926) > 1.0e-15) {
    std::cerr << "wall-normal Pe_z classifier was not exact\n";
    std::exit(1);
  }
  if (neumann::drift_envelope_exceeded(
          0.033 * kDiffusion / kHeight, kHeight, kDiffusion)) {
    std::cerr << "in-envelope Pe_z=0.033 was incorrectly rejected\n";
    std::exit(1);
  }
  if (!neumann::drift_envelope_exceeded(
          0.0926 * kDiffusion / kHeight, kHeight, kDiffusion)) {
    std::cerr << "above-envelope Pe_z was not rejected\n";
    std::exit(1);
  }
  std::cout << "  test_drift_classifier: PASSED\n";
}

void test_drift_envelope_policy() {
  auto system = make_system(0.0926 * kDiffusion / kHeight);
  ChemicalSpec bacteriocin;
  bacteriocin.name = species::BACTERIOCIN_BTUB;
  bacteriocin.diff_coeff = kDiffusion;
  bacteriocin.retardation = 1.0;
  bacteriocin.decay_rate = kDecay;
  const std::vector<ChemicalSpec> chemicals = {bacteriocin};

  QSSAConfig cfg;
  cfg.low_screening_policy = "allow";
  cfg.drift_envelope_policy = "error";
  bool rejected = false;
  try {
    QSSASolver qssa;
    qssa.init(cfg, system.domain, system.adv, false, &chemicals);
  } catch (const SimulationError&) {
    rejected = true;
  }
  if (!rejected) {
    std::cerr << "drift-envelope error policy did not reject\n";
    std::exit(1);
  }

  cfg.drift_envelope_policy = "allow";
  {
    QSSASolver qssa;
    qssa.init(cfg, system.domain, system.adv, false, &chemicals);
  }
  cfg.drift_envelope_policy = "warn";
  {
    QSSASolver qssa;
    qssa.init(cfg, system.domain, system.adv, false, &chemicals);
  }
  std::cout << "  test_drift_envelope_policy: PASSED\n";
}

void test_runtime_plasmid_basis_closes_gate_hole() {
  const Real flow_z = 0.033 * kDiffusion / kHeight;
  auto system = make_system(flow_z);
  ChemicalSpec bacteriocin;
  bacteriocin.name = species::BACTERIOCIN_BTUB;
  bacteriocin.diff_coeff = kDiffusion;
  bacteriocin.retardation = 1.0;
  bacteriocin.decay_rate = kDecay;
  const std::vector<ChemicalSpec> chemicals = {bacteriocin};
  const RuntimeDriftEnvelopeBasis high_runtime_basis{
      "ColE1", kDiffusion / 10.0};
  const RuntimeDriftEnvelopeBasis low_runtime_basis{"low", kDiffusion};

  QSSAConfig cfg;
  cfg.low_screening_policy = "allow";
  cfg.drift_envelope_policy = "error";

  bool rejected = false;
  try {
    QSSASolver qssa;
    qssa.init(cfg, system.domain, system.adv, false, &chemicals,
              &high_runtime_basis);
  } catch (const SimulationError& error) {
    const std::string_view message(error.what());
    const std::string_view expected("configured-plasmid basis: ColE1");
    rejected = std::search(
                   message.begin(), message.end(), expected.begin(), expected.end())
               != message.end();
  }
  if (!rejected) {
    std::cerr << "configured-plasmid drift basis did not close gate hole\n";
    std::exit(1);
  }

  {
    QSSASolver qssa;
    qssa.init(cfg, system.domain, system.adv, false, &chemicals,
              &low_runtime_basis);
  }
  {
    QSSASolver qssa;
    qssa.init(cfg, system.domain, system.adv, false, &chemicals);
  }
  // The ChemicalSpec basis has Pe_z=0.033 and passes, while the configured
  // plasmid basis has Pe_z=0.33 and must reject. The absent/low negative
  // controls ensure the new basis, rather than the existing species basis,
  // is what closes this false-negative hole.
  std::cout << "  test_runtime_plasmid_basis_closes_gate_hole: PASSED\n";
}

Real wall_residual(const GreensFunction& gf, const Vec3& source,
                   Real wall, Real h, Real d, Real flow_z,
                   const GreensFunctionParams& params, bool top) {
  std::array<Real, 4> values{};
  for (size_t i = 0; i < values.size(); ++i) {
    const Real offset = h * static_cast<Real>(i);
    const Real z = top ? wall - offset : wall + offset;
    values[i] = gf.concentration_bounded(
        source, {source[0] + 10.0e-6, source[1], z}, params);
  }
  const Real derivative = top
      ? (11.0 * values[0] - 18.0 * values[1] + 9.0 * values[2]
         - 2.0 * values[3]) / (6.0 * h)
      : (-11.0 * values[0] + 18.0 * values[1] - 9.0 * values[2]
         + 2.0 * values[3]) / (6.0 * h);
  const Real diffusive = -d * derivative;
  const Real advective = flow_z * values[0];
  return std::abs(diffusive + advective)
      / (std::abs(diffusive) + std::abs(advective)
         + d * std::abs(values[0]) / kHeight);
}

Real robin_wall_residual(const GreensFunction& gf, const Vec3& source,
                         Real h, Real d, Real transfer_length,
                         const GreensFunctionParams& params) {
  std::array<Real, 4> values{};
  for (size_t i = 0; i < values.size(); ++i) {
    const Real z = kHeight - h * static_cast<Real>(i);
    values[i] = gf.concentration_bounded(
        source, {source[0] + 10.0e-6, source[1], z}, params);
  }
  const Real derivative = (11.0 * values[0] - 18.0 * values[1]
      + 9.0 * values[2] - 2.0 * values[3]) / (6.0 * h);
  const Real diffusive = -d * derivative;
  const Real transfer = d / transfer_length * values[0];
  return std::abs(diffusive - transfer)
      / (std::abs(diffusive) + std::abs(transfer)
         + d * std::abs(values[0]) / kHeight);
}

// The wall law is a statement about the analytic composition, so this probes
// the sub-cell geometries where robin::requires_direct_evaluation() bypasses
// the interpolation table: a 60 um cell radius covers the mid-depth source and
// the near-wall targets below. Interpolated accuracy is bounded separately by
// test_physical_sealed_interpolation().
void test_physical_wall_law_composition() {
  const Real flow_z = 0.822 * kDiffusion / kHeight;
  auto system = make_system(flow_z, kFlowX, 60.0e-6);
  const Vec3 source = {500.0e-6, 500.0e-6, 0.5 * kHeight};
  const Real h = 0.25e-6;
  auto corrected = make_params();
  corrected.drift_correction = true;
  auto uncorrected = corrected;
  uncorrected.drift_correction = false;
  const Real corrected_bottom = wall_residual(
      system.gf, source, 0.0, h, kDiffusion, flow_z, corrected, false);
  const Real corrected_top = wall_residual(
      system.gf, source, kHeight, h, kDiffusion, flow_z, corrected, true);
  const Real uncorrected_bottom = wall_residual(
      system.gf, source, 0.0, h, kDiffusion, flow_z, uncorrected, false);
  const Real uncorrected_top = wall_residual(
      system.gf, source, kHeight, h, kDiffusion, flow_z, uncorrected, true);
  require(corrected_bottom <= 1.0e-5 && corrected_top <= 1.0e-5,
          "physical sealed wall-law residual exceeded tolerance");
  require(uncorrected_bottom >= 0.4 && uncorrected_top >= 0.4,
          "uncorrected sealed wall-law defect was not detected");
  std::cout << "  test_physical_wall_law_composition: PASSED (corrected "
               "bottom/top=" << corrected_bottom << "/" << corrected_top
            << ", uncorrected bottom/top=" << uncorrected_bottom << "/"
            << uncorrected_top << ")\n";
}

void test_zero_drift_correction_invariance() {
  const Real flow_z = 0.033 * kDiffusion / kHeight;
  auto table_system = make_system(flow_z);
  const auto table = robin::build_table(
      table_system.adv, 0.0, kHeight, 4.0e-11, kDiffusion, kDecay,
      robin::kZeroTransferLength, robin::kDefaultCutoff,
      robin::TransferBasis::Effective, 1.0e-14, 1, true, false, true);
  // Verify at a table node that the subtrahend uses exactly the explicit
  // one-shell cap used by the runtime image series.
  const int source_index = 7;
  const int target_index = 19;
  const int rho_index = 11;
  const Real source_z = kHeight * source_index
      / static_cast<Real>(robin::kTableNodes - 1);
  const Real target_z = kHeight * target_index
      / static_cast<Real>(robin::kTableNodes - 1);
  const Real rho = robin::kMinimumTableRho * std::pow(
      robin::kDefaultCutoff / robin::kMinimumTableRho,
      static_cast<Real>(rho_index) / static_cast<Real>(robin::kTableNodes - 1));
  const size_t node = static_cast<size_t>(
      robin::table_index(source_index, target_index, rho_index));
  const Real physical = robin::normalized_sealed_field({
      source_z, target_z, rho, 0.0, kHeight, kDiffusion, kDecay,
      kFlowX, 0.0, flow_z, robin::kTableModeCount});
  const Real gauge = std::exp(
      (kFlowX * rho + flow_z * (target_z - source_z))
      / (2.0 * kDiffusion));
  const Real image_from_table = physical / gauge - table.drift_values[node];
  const Real image_direct = robin::normalized_image_series(
      source_z, target_z, rho, 0.0, kHeight, kDiffusion, kDecay, kFlowX,
      flow_z, 1.0e-14, 1, true, false);
  const Real image_direct_gauge_free = image_direct / gauge;
  require(std::abs(image_from_table - image_direct_gauge_free)
              <= 1.0e-12 * std::max(std::abs(image_direct_gauge_free), 1.0e-30),
          "table image-series subtrahend disagreed with runtime series");
  auto zero_flow_system = make_system(0.0);
  const Vec3 source = {500.0e-6, 500.0e-6, 0.4 * kHeight};
  const Vec3 target = {510.0e-6, 500.0e-6, 0.7 * kHeight};
  auto corrected = make_params();
  corrected.drift_correction = true;
  const auto uncorrected = make_params();
  const Real corrected_value = zero_flow_system.gf.concentration_bounded(
      source, target, corrected);
  const Real uncorrected_value = zero_flow_system.gf.concentration_bounded(
      source, target, uncorrected);
  const Real relative = relative_error(corrected_value, uncorrected_value);
  // With zero drift and transfer disabled, the correction is exactly zero.
  require(relative <= 1.0e-9,
          "zero-drift correction changed the sealed field");
  std::cout << "  test_zero_drift_correction_invariance: PASSED (field="
            << relative << ")\n";
}

void test_physical_sealed_interpolation() {
  const Real flow_z = 0.822 * kDiffusion / kHeight;
  auto system = make_system(flow_z);
  auto params = make_params();
  params.drift_correction = true;
  const std::array<Probe, 9> probes = {{
      {0.137, 0.291, 1.13e-6}, {0.243, 0.814, 3.71e-6},
      {0.367, 0.492, 8.19e-6}, {0.431, 0.713, 14.7e-6},
      {0.559, 0.168, 27.3e-6}, {0.641, 0.887, 44.1e-6},
      {0.726, 0.354, 71.8e-6}, {0.813, 0.629, 113.0e-6},
      {0.914, 0.243, 173.0e-6},
  }};
  std::array<Real, 9> errors{};
  for (size_t i = 0; i < probes.size(); ++i) {
    errors[i] = probe_error(system, params, probes[i], kFlowX, flow_z);
  }
  std::sort(errors.begin(), errors.end());
  const Real median = errors[errors.size() / 2];
  const Real maximum = errors.back();
  require(median <= 1.0e-3 && maximum <= 0.2,
          "physical sealed interpolation exceeded tolerance");
  std::cout << "  test_physical_sealed_interpolation: PASSED (median="
            << median << ", max=" << maximum << ")\n";
}

void test_default_off_identity() {
  auto system = make_system(0.822 * kDiffusion / kHeight);
  const auto params = make_params();
  const std::array<Probe, 5> probes = {{
      {0.137, 0.291, 1.13e-6}, {0.243, 0.814, 8.19e-6},
      {0.500, 0.500, 27.3e-6}, {0.726, 0.354, 71.8e-6},
      {0.914, 0.243, 113.0e-6},
  }};
  const std::array<Real, 10> main_values = {{
      0x1.d0d52ca97fba1p-12, 0x1.d7f536f329f8fp-13,
      0x1.d825310450f75p-13, 0x1.64ef9befbd004p-15,
      0x1.2e5cc26721758p-12, 0x1.2efdc394ed6d3p-13,
      0x1.7a8e9eb585528p-13, 0x1.1331c424b79acp-14,
      0x1.3c5d0c7f8d2bap-13, 0x1.ea8d2b9b7d99ep-15,
  }};
  for (size_t i = 0; i < probes.size(); ++i) {
    const Vec3 source = {500.0e-6, 500.0e-6,
                         probes[i].source_fraction * kHeight};
    const Vec3 target = {500.0e-6 + probes[i].rho, 500.0e-6,
                         probes[i].target_fraction * kHeight};
    const Real sealed = system.gf.concentration_bounded(
        source, target, params);
    auto bounded_params = params;
    bounded_params.lumen_transfer_length = 100.0e-6;
    const Real bounded = system.gf.concentration_bounded(
        source, target, bounded_params);
    require(std::bit_cast<uint64_t>(sealed)
                == std::bit_cast<uint64_t>(main_values[2 * i])
            && std::bit_cast<uint64_t>(bounded)
                == std::bit_cast<uint64_t>(main_values[2 * i + 1]),
            "default-off Green's function changed from origin/main");
  }
  std::cout << "  test_default_off_identity: PASSED (5 fixed geometries)\n";
}

void test_robin_composition() {
  const Real flow_z = 0.822 * kDiffusion / kHeight;
  auto system = make_system(flow_z, kFlowX, 60.0e-6);
  auto params = make_params();
  params.drift_correction = true;
  params.lumen_transfer_length = 100.0e-6;
  const Vec3 source = {500.0e-6, 500.0e-6, 0.5 * kHeight};
  const Real h = 0.25e-6;
  const Real bottom = wall_residual(
      system.gf, source, 0.0, h, kDiffusion, flow_z, params, false);
  const Real top = robin_wall_residual(
      system.gf, source, h, kDiffusion, params.lumen_transfer_length, params);
  require(bottom <= 1.0e-5 && top <= 1.0e-5,
          "corrected Robin composition violated a wall law");
  std::cout << "  test_robin_composition: PASSED (bottom/top=" << bottom
            << "/" << top << ")\n";
}

}  // namespace

int main() {
  std::cout << "=== Wall-Normal Drift Envelope Tests ===\n";
  test_zero_wall_normal_flow_is_exact();
  test_wall_normal_flow_sweep_is_sensitive_and_ordered();
  test_reflected_flow_reversal_is_detectable();
  test_wall_parallel_flow_is_exact();
  test_drift_classifier();
  test_drift_envelope_policy();
  test_runtime_plasmid_basis_closes_gate_hole();
  test_physical_wall_law_composition();
  test_zero_drift_correction_invariance();
  test_physical_sealed_interpolation();
  test_robin_composition();
  test_default_off_identity();
  std::cout << "All wall-normal drift envelope tests passed.\n";
  return 0;
}
