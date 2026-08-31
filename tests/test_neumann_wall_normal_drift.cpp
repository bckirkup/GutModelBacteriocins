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

TestSystem make_system(Real flow_z, Real flow_x = kFlowX) {
  TestSystem system;
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {1.0e-3, 1.0e-3, kHeight};
  dcfg.grid_dx = 5.0e-6;
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

Real modal_reference(const Probe& probe, Real flow_z, Real flow_x = kFlowX) {
  const robin::SealedFieldParams params{
      probe.source_fraction * kHeight, probe.target_fraction * kHeight,
      probe.rho, 0.0, kHeight, kDiffusion, kDecay, flow_x, 0.0, flow_z,
      4000};
  const Real normalized = robin::normalized_sealed_field(params);
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
  const Real reference = modal_reference(probe, flow_z, flow_x);
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
  std::cout << "All wall-normal drift envelope tests passed.\n";
  return 0;
}
