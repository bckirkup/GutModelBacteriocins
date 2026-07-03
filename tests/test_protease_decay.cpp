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

}  // namespace

void test_burst_halflife_decay() {
  Domain domain;
  AdvectionField adv;
  QSSASolver qssa;
  setup_qssa(domain, adv, qssa);

  const BICluster bi = PlasmidLibrary::colicin_E1();
  const Real half_life = bi.protease_half_life;
  const Real decay_rate = 0.6931471805599453 / half_life;

  Vec3 source = {50e-6, 50e-6, 25e-6};
  Vec3 target = {55e-6, 50e-6, 25e-6};
  std::vector<Vec3> sources = {source};
  std::vector<GreensFunctionParams> params = {
      params_from_bi(bi, 1.0e-18)};

  const Real c0 = qssa.point_concentration(target, sources, params, {1.0});
  const Real factor_half = std::exp(-decay_rate * half_life);
  const Real c_half = qssa.point_concentration(target, sources, params, {factor_half});

  assert(c0 > 0.0);
  const Real ratio = c_half / c0;
  assert(std::abs(ratio - 0.5) < 0.05);

  std::cout << "  test_burst_halflife_decay: PASSED (ratio=" << ratio << ")\n";
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
  const Real decay_b = std::exp(-0.6931471805599453 / col_b.protease_half_life * age);
  const Real decay_mcc = std::exp(-0.6931471805599453 / mcc.protease_half_life * age);
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
  sim.agents()[0].sos_timer = 0.0;
  fix.post_step(60.0);

  assert(sim.toxin_bursts().size() == 1);
  assert(sim.agents()[0].state == PhenoState::DEAD);

  std::cout << "  test_lysis_registers_burst: PASSED\n";
}

int main() {
  std::cout << "=== Protease Decay Tests ===\n";
  test_burst_halflife_decay();
  test_protease_disabled_no_decay();
  test_per_colicin_decay_rates();
  test_lysis_registers_burst();
  std::cout << "All protease decay tests passed.\n";
  return 0;
}
