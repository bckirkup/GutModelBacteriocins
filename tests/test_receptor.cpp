/* -----------------------------------------------------------------------
   GutIBM – Unit tests for fix_receptor
   ----------------------------------------------------------------------- */

#include "fix_receptor.h"
#include "plasmid.h"
#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

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

static Agent make_susceptible_agent(Simulation& sim) {
  Vec3 center = {
    0.5 * (sim.domain().lo()[0] + sim.domain().hi()[0]),
    0.5 * (sim.domain().lo()[1] + sim.domain().hi()[1]),
    0.5 * (sim.domain().lo()[2] + sim.domain().hi()[2]),
  };
  Agent a = Agent::create_default(sim.agents().next_tag(), 2, center, 5e-4);
  // Isolate BtuB-mediated killing for colicin E immunity tests.
  a.receptor_expr[to_underlying(ReceptorType::FepA)] = 0.0;
  a.receptor_expr[to_underlying(ReceptorType::CirA)] = 0.0;
  Int ix, iy, iz;
  sim.domain().pos_to_grid(a.x, ix, iy, iz);
  a.grid_cell = sim.domain().cell_index(ix, iy, iz);
  return a;
}

static void set_local_chemistry(Simulation& sim, Int cell,
                                Real tox_conc, Real b12_conc) {
  auto& chem = sim.chemical_field();
  Int i_tox = chem.find("bacteriocin");
  Int i_b12 = chem.find("b12");
  assert(i_tox >= 0 && i_b12 >= 0);
  chem.conc(i_tox, cell) = tox_conc;
  chem.conc(i_b12, cell) = b12_conc;
}

void test_high_toxin_kills_susceptible() {
  ReceptorConfig rcfg;
  rcfg.kill_rate_colicin = 1.0;

  auto sim = make_empty_sim(3003);
  Agent a = make_susceptible_agent(sim);
  Int cell = a.grid_cell;
  set_local_chemistry(sim, cell, 1.0e-4, 0.0);
  sim.agents().push_back(std::move(a));

  FixReceptor fix(sim, rcfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::DEAD);

  std::cout << "  test_high_toxin_kills_susceptible: PASSED\n";
}

void test_immunity_reduces_lethality() {
  ReceptorConfig rcfg;
  rcfg.kill_rate_colicin = 1.0;

  auto sim_immune = make_empty_sim(4004);
  Agent immune = make_susceptible_agent(sim_immune);
  Int cell = immune.grid_cell;
  immune.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  set_local_chemistry(sim_immune, cell, 1.0e-4, 0.0);
  sim_immune.agents().push_back(std::move(immune));

  auto sim_naive = make_empty_sim(4004);
  Agent naive = make_susceptible_agent(sim_naive);
  set_local_chemistry(sim_naive, cell, 1.0e-4, 0.0);
  sim_naive.agents().push_back(std::move(naive));

  FixReceptor fix_immune(sim_immune, rcfg);
  FixReceptor fix_naive(sim_naive, rcfg);
  fix_immune.compute(60.0);
  fix_naive.compute(60.0);

  assert(sim_naive.agents()[0].state == PhenoState::DEAD);
  assert(sim_immune.agents()[0].state != PhenoState::DEAD);

  std::cout << "  test_immunity_reduces_lethality: PASSED\n";
}

void test_ligand_competition_reduces_kill() {
  ReceptorConfig rcfg;
  rcfg.kill_rate_colicin = 1.0;

  auto sim_low_ligand = make_empty_sim(5005);
  Agent a_low = make_susceptible_agent(sim_low_ligand);
  Int cell = a_low.grid_cell;
  set_local_chemistry(sim_low_ligand, cell, 1.0e-4, 0.0);
  sim_low_ligand.agents().push_back(std::move(a_low));

  auto sim_high_ligand = make_empty_sim(5005);
  Agent a_high = make_susceptible_agent(sim_high_ligand);
  set_local_chemistry(sim_high_ligand, cell, 1.0e-4, 1.0);
  sim_high_ligand.agents().push_back(std::move(a_high));

  FixReceptor fix_low(sim_low_ligand, rcfg);
  FixReceptor fix_high(sim_high_ligand, rcfg);
  fix_low.compute(60.0);
  fix_high.compute(60.0);

  assert(sim_low_ligand.agents()[0].state == PhenoState::DEAD);
  assert(sim_high_ligand.agents()[0].state != PhenoState::DEAD);

  std::cout << "  test_ligand_competition_reduces_kill: PASSED\n";
}

void test_partial_resistance_reduces_lethality() {
  ReceptorConfig rcfg;
  rcfg.kill_rate_colicin = 1.0;

  auto sim_resistant = make_empty_sim(6006);
  Agent resistant = make_susceptible_agent(sim_resistant);
  Int cell = resistant.grid_cell;
  resistant.genome.toxin_affinity[to_underlying(ReceptorType::BtuB)] = 1.0e-6;
  set_local_chemistry(sim_resistant, cell, 1.0e-7, 0.0);
  sim_resistant.agents().push_back(std::move(resistant));

  auto sim_wt = make_empty_sim(6006);
  Agent wt = make_susceptible_agent(sim_wt);
  set_local_chemistry(sim_wt, cell, 1.0e-7, 0.0);
  sim_wt.agents().push_back(std::move(wt));

  FixReceptor fix_resistant(sim_resistant, rcfg);
  FixReceptor fix_wt(sim_wt, rcfg);
  fix_resistant.compute(60.0);
  fix_wt.compute(60.0);

  assert(sim_wt.agents()[0].state == PhenoState::DEAD);
  assert(sim_resistant.agents()[0].state != PhenoState::DEAD);

  std::cout << "  test_partial_resistance_reduces_lethality: PASSED\n";
}

int main() {
  std::cout << "=== Receptor Fix Tests ===\n";
  test_high_toxin_kills_susceptible();
  test_immunity_reduces_lethality();
  test_ligand_competition_reduces_kill();
  test_partial_resistance_reduces_lethality();
  std::cout << "All receptor fix tests passed.\n";
  return 0;
}
