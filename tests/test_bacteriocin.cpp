/* -----------------------------------------------------------------------
   GutIBM – Unit tests for fix_bacteriocin (Spec 2)
   ----------------------------------------------------------------------- */

#include "fix_bacteriocin.h"
#include "plasmid.h"
#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <cmath>
#include <iostream>

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
  auto e2 = PlasmidLibrary::colicin_E2();
  auto col_b = PlasmidLibrary::colicin_B();
  auto mcc = PlasmidLibrary::microcin_V();

  assert(e1.bclass == BacteriocinClass::LETHAL_CORE);
  assert(e1.pI > 8.5);
  assert(e1.release_mode == ReleaseMode::SOS_LYSIS);
  assert(e1.burst_size == 1.0e5);

  assert(e2.is_nuclease);
  assert(e2.release_mode == ReleaseMode::SOS_LYSIS);

  assert(col_b.bclass == BacteriocinClass::LETHAL_HALO);
  assert(col_b.pI < 7.0);
  assert(col_b.release_mode == ReleaseMode::PHAGE_LYSIS);

  assert(mcc.molecular_weight < 10000.0);
  assert(mcc.release_mode == ReleaseMode::CONTINUOUS);

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

void test_sos_induction_on_division() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 0.0;
  cfg.sos_lysis_prob = 1.0;

  auto sim = make_empty_sim(3003);
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  a.just_divided = true;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::SOS_INDUCED);

  std::cout << "  test_sos_induction_on_division: PASSED\n";
}

void test_no_sos_without_division() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 0.0;
  cfg.sos_lysis_prob = 1.0;

  auto sim = make_empty_sim(4004);
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  a.just_divided = false;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::NORMAL);

  std::cout << "  test_no_sos_without_division: PASSED\n";
}

void test_phage_induction() {
  BacteriocinConfig cfg;

  auto sim = make_empty_sim(5005);
  Agent a = make_agent_at_center(sim, 1);
  BICluster col_b = PlasmidLibrary::colicin_B();
  col_b.phage_induction_rate = 100.0;
  a.genome.bi_loci.push_back(col_b);
  a.mu_realized = 5.0e-4;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::SOS_INDUCED);

  std::cout << "  test_phage_induction: PASSED\n";
}

void test_phage_does_not_trigger_sos_path() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 1.0;

  auto sim = make_empty_sim(6006);
  Agent a = make_agent_at_center(sim, 1);
  BICluster col_b = PlasmidLibrary::colicin_B();
  col_b.phage_induction_rate = 0.0;
  a.genome.bi_loci.push_back(col_b);
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::NORMAL);

  std::cout << "  test_phage_does_not_trigger_sos_path: PASSED\n";
}

void test_microcin_no_lysis() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 1.0;

  auto sim = make_empty_sim(7007);
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::microcin_V());
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::NORMAL);

  std::cout << "  test_microcin_no_lysis: PASSED\n";
}

void test_cross_induction() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 0.0;
  cfg.sos_lysis_prob = 0.0;
  cfg.sos_cross_induction_rate = 1.0e3;

  auto sim = make_empty_sim(8008);
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  sim.agents().push_back(std::move(a));

  Int i_nuc = sim.chemical_field().find("nuclease_toxin");
  assert(i_nuc >= 0);
  sim.chemical_field().conc(i_nuc, sim.agents()[0].grid_cell) = 1.0e-3;

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::SOS_INDUCED);

  std::cout << "  test_cross_induction: PASSED\n";
}

void test_per_colicin_burst_size() {
  BacteriocinConfig cfg;
  cfg.burst_molecules = 1.0e4;

  auto sim = make_empty_sim();
  Agent a = make_agent_at_center(sim, 1);
  BICluster e1 = PlasmidLibrary::colicin_E1();
  assert(e1.burst_size == 1.0e5);
  a.genome.bi_loci.push_back(e1);
  a.state = PhenoState::SOS_INDUCED;
  a.sos_timer = 1.0;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.post_step(2.0);

  assert(sim.toxin_bursts().size() == 1);
  const Real expected = sim.config().qssa.colicin_release_rate * (1.0e5 / 1.0e4);
  assert(std::abs(sim.toxin_bursts()[0].params.source_rate - expected) < 1e-30);
  assert(sim.toxin_bursts()[0].is_nuclease == false);

  std::cout << "  test_per_colicin_burst_size: PASSED\n";
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
  test_sos_induction_on_division();
  test_no_sos_without_division();
  test_phage_induction();
  test_phage_does_not_trigger_sos_path();
  test_microcin_no_lysis();
  test_cross_induction();
  test_per_colicin_burst_size();
  test_sos_lysis_post_step();
  std::cout << "All bacteriocin fix tests passed.\n";
  return 0;
}
