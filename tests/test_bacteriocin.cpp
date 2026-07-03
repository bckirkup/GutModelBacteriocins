/* -----------------------------------------------------------------------
   GutIBM – Unit tests for fix_bacteriocin
   ----------------------------------------------------------------------- */

#include "fix_bacteriocin.h"
#include "plasmid.h"
#include "simulation.h"
#include "input_parser.h"
#include "qssa_solver.h"
#include "domain.h"
#include "advection.h"
#include "chemical_field.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

static Simulation make_empty_sim(uint64_t seed = 42) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.seed = seed;

  Simulation sim;
  sim.init(cfg);
  return sim;
}

static Agent make_agent_at_center(Simulation& sim, Int type) {
  Vec3 center = {
    0.5 * (sim.domain().lo()[0] + sim.domain().hi()[0]),
    0.5 * (sim.domain().lo()[1] + sim.domain().hi()[1]),
    0.5 * (sim.domain().lo()[2] + sim.domain().hi()[2]),
  };
  Agent a = Agent::create_default(sim.agents().next_tag(), type, center, 5e-4);
  Int ix;
  Int iy;
  Int iz;
  sim.domain().pos_to_grid(a.x, ix, iy, iz);
  a.grid_cell = sim.domain().cell_index(ix, iy, iz);
  return a;
}

void test_pi_diffusion_classes() {
  auto e1 = PlasmidLibrary::colicin_E1();
  auto col_b = PlasmidLibrary::colicin_B();
  auto mcc = PlasmidLibrary::microcin_V();

  assert(e1.bclass == BacteriocinClass::LETHAL_CORE);
  assert(e1.pI > 8.5);
  assert(col_b.bclass == BacteriocinClass::LETHAL_HALO);
  assert(col_b.pI < 7.0);
  assert(mcc.molecular_weight < 10000.0);

  std::cout << "  test_pi_diffusion_classes: PASSED\n";
}

void test_microcin_mu_penalty() {
  BacteriocinConfig cfg;
  cfg.microcin_mu_penalty = 0.03;

  auto sim = make_empty_sim();
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::microcin_V());
  Real mu_before = a.mu_max;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  Real mu_after = sim.agents()[0].mu_max;
  assert(std::abs(mu_after - mu_before * (1.0 - cfg.microcin_mu_penalty)) < 1e-12);

  std::cout << "  test_microcin_mu_penalty: PASSED\n";
}

void test_sos_induction_requires_bi_loci() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 1.0;

  auto sim = make_empty_sim(1001);
  Agent a = make_agent_at_center(sim, 1);
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::NORMAL);

  std::cout << "  test_sos_induction_requires_bi_loci: PASSED\n";
}

void test_sos_induction_high_basal_rate() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 1.0;

  auto sim = make_empty_sim(2002);
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::SOS_INDUCED);
  assert(sim.agents()[0].sos_timer > 0.0);

  std::cout << "  test_sos_induction_high_basal_rate: PASSED\n";
}

void test_colicin_steady_state_qssa_source() {
  // Regression: ColE1/ColB (MW >= 10 kDa) must contribute continuous QSSA sources.
  const BICluster col_e1 = PlasmidLibrary::colicin_E1();
  assert(col_e1.molecular_weight >= 10000.0);

  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {100e-6, 100e-6, 50e-6};
  dcfg.grid_dx = 5e-6;
  Domain domain;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.mucus_thickness = 50e-6;
  AdvectionField adv;
  adv.init(acfg, domain);

  QSSAConfig qcfg;
  qcfg.toxin_cutoff = 80e-6;
  qcfg.microcin_secretion = 1.0e-20;
  QSSASolver qssa;
  qssa.init(qcfg, domain, adv);

  ChemicalSpec toxin;
  toxin.name = "bacteriocin";
  toxin.diff_coeff = 4.0e-11;
  toxin.retardation = 10.0;
  ChemicalField chem;
  chem.init(domain, {toxin});

  const Vec3 pos = {50e-6, 50e-6, 25e-6};
  AgentPool agents;
  Agent producer = Agent::create_default(1, 1, pos, 5e-4);
  producer.genome.bi_loci.push_back(col_e1);
  agents.push_back(std::move(producer));

  ProteaseConfig protease;
  protease.enabled = false;
  std::vector<ToxinBurstSource> bursts;
  qssa.solve_bacteriocin_field(agents, bursts, 0.0, protease, adv, chem, 0);

  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  domain.pos_to_grid(pos, ix, iy, iz);
  const Int idx = domain.cell_index(ix, iy, iz);
  assert(chem.conc(0, idx) > 0.0);

  std::cout << "  test_colicin_steady_state_qssa_source: PASSED (c="
            << chem.conc(0, idx) << ")\n";
}

void test_sos_lysis_post_step() {
  BacteriocinConfig cfg;

  auto sim = make_empty_sim();
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  a.state = PhenoState::SOS_INDUCED;
  a.sos_timer = 30.0;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.post_step(60.0);

  assert(sim.agents()[0].state == PhenoState::DEAD);

  std::cout << "  test_sos_lysis_post_step: PASSED\n";
}

int main() {
  std::cout << "=== Bacteriocin Fix Tests ===\n";
  test_pi_diffusion_classes();
  test_microcin_mu_penalty();
  test_sos_induction_requires_bi_loci();
  test_sos_induction_high_basal_rate();
  test_colicin_steady_state_qssa_source();
  test_sos_lysis_post_step();
  std::cout << "All bacteriocin fix tests passed.\n";
  return 0;
}
