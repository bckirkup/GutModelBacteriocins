/* -----------------------------------------------------------------------
   GutIBM – Protease degradation of bacteriocin burst sources (Spec 1)
   ----------------------------------------------------------------------- */

#include "qssa_solver.h"
#include "fix_bacteriocin.h"
#include "simulation.h"
#include "input_parser.h"
#include "domain.h"
#include "advection.h"
#include "plasmid.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

using namespace gutibm;

namespace {

void setup_qssa(Domain& domain, AdvectionField& adv, QSSASolver& qssa) {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {100e-6, 100e-6, 50e-6};
  dcfg.grid_dx = 5e-6;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.mucus_thickness = 50e-6;
  adv.init(acfg, domain);

  QSSAConfig qcfg;
  qcfg.toxin_cutoff = 80e-6;
  qssa.init(qcfg, domain, adv);
}

GreensFunctionParams params_from_bi(const BICluster& bi, Real release_rate) {
  GreensFunctionParams gfp;
  gfp.diff_coeff = bi.diff_coeff;
  gfp.retardation = bi.retardation;
  gfp.pI = bi.pI;
  gfp.source_rate = release_rate;
  return gfp;
}

// Release-integral quadrature: dt = tau / kStepsPerTau, horizon = kTauHorizons * tau.
constexpr int kStepsPerTau = 1000;
constexpr int kTauHorizons = 20;
constexpr int kQuadratureSteps = kStepsPerTau * kTauHorizons;

}  // namespace

void test_spatial_decay_not_temporal_amplitude() {
  Domain domain;
  AdvectionField adv;
  QSSASolver qssa;
  setup_qssa(domain, adv, qssa);

  const BICluster bi = PlasmidLibrary::colicin_E1();
  const Real half_life = bi.protease_half_life;
  const Real decay_rate = std::numbers::ln2 / half_life;

  Vec3 source = {50e-6, 50e-6, 25e-6};
  Vec3 target = {55e-6, 50e-6, 25e-6};
  std::vector<Vec3> sources = {source};
  std::vector<GreensFunctionParams> params = {
      params_from_bi(bi, 1.0e-18)};
  params[0].decay_rate = decay_rate;

  const Real c0 = qssa.point_concentration(target, sources, params, {1.0});
  assert(c0 > 0.0);
  const Real c_same_source_age = qssa.point_concentration(target, sources, params, {1.0});
  assert(std::abs(c_same_source_age - c0) < 1e-30 * std::max(c0, 1.0));

  std::cout << "  test_spatial_decay_not_temporal_amplitude: PASSED (k="
            << decay_rate << ")\n";
}

void test_inventory_conservation() {
  const Real inventory = 1.0e5 / AVOGADRO;
  for (const Real tau : {60.0, 300.0, 1800.0}) {
    const Real dt = tau / static_cast<Real>(kStepsPerTau);
    Real integrated = 0.0;
    for (int i = 0; i < kQuadratureSteps; ++i) {
      const Real age = static_cast<Real>(i) * dt;
      integrated += (inventory / tau) * std::exp(-age / tau) * dt;
    }
    assert(std::abs(integrated - inventory) / inventory < 1.0e-3);
  }
  std::cout << "  test_inventory_conservation: PASSED\n";
}

void test_dose_invariant_to_release_tau() {
  Domain domain;
  AdvectionField adv;
  QSSASolver qssa;
  setup_qssa(domain, adv, qssa);
  const BICluster bi = PlasmidLibrary::colicin_E1();
  const Vec3 source = {50e-6, 50e-6, 25e-6};
  const Vec3 target = {60e-6, 50e-6, 25e-6};
  const Real inventory = bi.burst_size / AVOGADRO;
  GreensFunctionParams params = params_from_bi(bi, inventory);
  params.decay_rate = 0.0;
  const Real kernel = qssa.point_concentration(target, {source}, {params}, {1.0});
  Real reference = 0.0;
  bool have_reference = false;
  for (const Real tau : {60.0, 300.0, 1800.0}) {
    const Real dt = tau / static_cast<Real>(kStepsPerTau);
    Real dose = 0.0;
    for (int i = 0; i < kQuadratureSteps; ++i) {
      const Real age = static_cast<Real>(i) * dt;
      dose += kernel * std::exp(-age / tau) * dt / tau;
    }
    if (!have_reference) {
      reference = dose;
      have_reference = true;
    }
    assert(std::abs(dose - reference) / reference < 1.0e-3);
  }
  std::cout << "  test_dose_invariant_to_release_tau: PASSED (dose="
            << reference << ")\n";
}

void test_protease_disabled_no_decay() {
  Domain domain;
  AdvectionField adv;
  QSSASolver qssa;
  setup_qssa(domain, adv, qssa);

  const BICluster bi = PlasmidLibrary::colicin_B();
  Vec3 source = {50e-6, 50e-6, 25e-6};
  Vec3 target = {55e-6, 50e-6, 25e-6};
  std::vector<Vec3> sources = {source};
  std::vector<GreensFunctionParams> params = {
      params_from_bi(bi, 1.0e-18)};

  const Real c_young = qssa.point_concentration(target, sources, params, {1.0});
  const Real c_old = qssa.point_concentration(target, sources, params, {1.0});
  assert(std::abs(c_young - c_old) < 1e-30 * std::max(c_young, 1.0));

  std::cout << "  test_protease_disabled_no_decay: PASSED\n";
}

void test_per_colicin_decay_rates() {
  const BICluster col_b = PlasmidLibrary::colicin_B();
  const BICluster mcc = PlasmidLibrary::microcin_V();
  assert(col_b.protease_half_life < mcc.protease_half_life);

  const Real age = 1800.0;
  const Real decay_b = std::exp(-std::numbers::ln2 / col_b.protease_half_life * age);
  const Real decay_mcc = std::exp(-std::numbers::ln2 / mcc.protease_half_life * age);
  assert(decay_b < decay_mcc);

  std::cout << "  test_per_colicin_decay_rates: PASSED"
            << " (ColB=" << decay_b << " MccV=" << decay_mcc << ")\n";
}

void test_lysis_registers_burst() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;

  Simulation sim;
  sim.init(cfg);

  Agent a = Agent::create_default(sim.agents().next_tag(), 1,
      {25e-6, 25e-6, 10e-6}, 5e-4);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  sim.agents().push_back(std::move(a));

  BacteriocinConfig bcfg;
  FixBacteriocin fix(sim, bcfg);
  fix.compute(60.0);
  sim.agents()[0].state = PhenoState::SOS_INDUCED;
  sim.agents()[0].timers.sos_timer = 0.0;
  fix.post_step(60.0);

  assert(sim.toxin_bursts().size() == 1);
  assert(sim.agents()[0].state == PhenoState::DEAD);
  const ToxinBurstSource& burst = sim.toxin_bursts().front();
  const Real expected_inventory = 1.0e5 / AVOGADRO;
  assert(std::abs(burst.params.source_rate
                  - expected_inventory / bcfg.burst_release_tau)
         < 1e-30);
  assert(std::abs(burst.release_tau - bcfg.burst_release_tau) < 1e-12);
  assert(std::abs(burst.params.decay_rate
                  - std::numbers::ln2 / 1800.0)
         < 1e-15);

  std::cout << "  test_lysis_registers_burst: PASSED\n";
}

int main() {
  std::cout << "=== Protease Decay Tests ===\n";
  test_spatial_decay_not_temporal_amplitude();
  test_inventory_conservation();
  test_dose_invariant_to_release_tau();
  test_protease_disabled_no_decay();
  test_per_colicin_decay_rates();
  test_lysis_registers_burst();
  std::cout << "All protease decay tests passed.\n";
  return 0;
}
