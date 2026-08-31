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
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
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

TestSystem make_system(Real flow_z) {
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
  acfg.distal_length = 1.0e-3;
  acfg.distal_transit_time = 43200.0;
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

Real modal_reference(const Probe& probe, Real flow_z) {
  const Real normalized = robin::normalized_sealed_field(
      probe.source_fraction * kHeight, probe.target_fraction * kHeight,
      probe.rho, 0.0, kHeight, kDiffusion, kDecay, kFlowX, 0.0, flow_z,
      4000);
  return kSourceRate / (4.0 * std::numbers::pi * kDiffusion) * normalized;
}

void test_zero_wall_normal_flow_is_exact() {
  auto system = make_system(0.0);
  const auto params = make_params();
  Real maximum = 0.0;
  for (const Probe& probe : kProbes) {
    const Vec3 source = {500.0e-6, 500.0e-6,
                         probe.source_fraction * kHeight};
    const Vec3 target = {source[0] + probe.rho, source[1],
                         probe.target_fraction * kHeight};
    const Real image = system.gf.concentration_bounded(source, target, params);
    maximum = std::max(maximum, relative_error(
        image, modal_reference(probe, 0.0)));
  }
  if (!(maximum < 1.0e-8)) {
    std::cerr << "zero-flow sealed image/modal mismatch: " << maximum << "\n";
    std::exit(1);
  }
  std::cout << "  test_zero_wall_normal_flow_is_exact: PASSED (max relative "
               "error=" << maximum << ")\n";
}

void test_wall_normal_flow_defect_is_first_order() {
  for (const Real pe_z : {0.0926, 0.37}) {
    const Real flow_z = pe_z * kDiffusion / kHeight;
    auto system = make_system(flow_z);
    const auto params = make_params();
    Real worst = 0.0;
    for (const Probe& probe : kProbes) {
      const Vec3 source = {500.0e-6, 500.0e-6,
                           probe.source_fraction * kHeight};
      const Vec3 target = {source[0] + probe.rho, source[1],
                           probe.target_fraction * kHeight};
      const Real image = system.gf.concentration_bounded(source, target, params);
      worst = std::max(worst, relative_error(
          image, modal_reference(probe, flow_z)));
    }
    if (!(worst <= 0.6 * pe_z && worst >= 0.20 * pe_z)) {
      std::cerr << "wall-normal drift defect outside documented bounds: Pe_z="
                << pe_z << " worst=" << worst << "\n";
      std::exit(1);
    }
    // The lower bound pins the known documented defect, so a stub or an
    // accidental change to the series construction cannot pass. An exact
    // treatment must update this test and the drift document together.
    std::cout << "  Pe_z=" << pe_z << " worst relative error=" << worst
              << "\n";
  }
  std::cout << "  test_wall_normal_flow_defect_is_first_order: PASSED\n";
}

void test_drift_classifier() {
  const Real classified = neumann::wall_normal_peclet(
      0.0926 * kDiffusion / kHeight, kHeight, kDiffusion);
  if (std::abs(classified - 0.0926) > 1.0e-15) {
    std::cerr << "wall-normal Pe_z classifier was not exact\n";
    std::exit(1);
  }
  if (neumann::drift_envelope_exceeded(
          0.033 * kDiffusion / kHeight, kHeight, kDiffusion)) {
    std::cerr << "shipped mid-domain Pe_z was incorrectly rejected\n";
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

}  // namespace

int main() {
  std::cout << "=== Wall-Normal Drift Envelope Tests ===\n";
  test_zero_wall_normal_flow_is_exact();
  test_wall_normal_flow_defect_is_first_order();
  test_drift_classifier();
  test_drift_envelope_policy();
  std::cout << "All wall-normal drift envelope tests passed.\n";
  return 0;
}
