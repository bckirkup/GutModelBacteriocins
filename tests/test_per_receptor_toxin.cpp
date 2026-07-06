/* -----------------------------------------------------------------------
   GutIBM – Per-receptor bacteriocin field tests (Spec 4 §1)
   ----------------------------------------------------------------------- */

#include "fix_bacteriocin.h"
#include "fix_receptor.h"
#include "plasmid.h"
#include "simulation.h"
#include "input_parser.h"
#include "qssa_solver.h"
#include "domain.h"
#include "advection.h"
#include "chemical_field.h"
#include "species_names.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

static Simulation make_tiny_sim() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.seed = 4242;

  Simulation sim;
  sim.init(cfg);
  return sim;
}

static Agent agent_at(Simulation& sim, Vec3 pos, Int type) {
  Agent a = Agent::create_default(sim.agents().next_tag(), type, pos, 5e-4);
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  sim.domain().pos_to_grid(a.x, ix, iy, iz);
  a.grid_cell = sim.domain().cell_index(ix, iy, iz);
  return a;
}

void test_cole1_burst_routes_btuB() {
  auto sim = make_tiny_sim();
  Vec3 pos = {50e-6, 50e-6, 25e-6};
  Agent a = agent_at(sim, pos, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  a.state = PhenoState::SOS_INDUCED;
  a.timers.sos_timer = 0.0;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, sim.config().fixes.bacteriocin);
  fix.post_step(1.0);
  sim.step(60.0);

  const auto& chem = sim.chemical_field();
  Int i_btuB = chem.find(species::BACTERIOCIN_BTUB);
  Int i_fepA = chem.find(species::BACTERIOCIN_FEPA);
  assert(i_btuB >= 0 && i_fepA >= 0);

  Int cell = sim.agents().size() == 0 ? 0 : sim.agents()[0].grid_cell;
  (void)cell;
  Real max_btuB = 0.0;
  Real max_fepA = 0.0;
  for (Int c = 0; c < chem.ncells(); ++c) {
    max_btuB = std::max(max_btuB, chem.conc(i_btuB, c));
    max_fepA = std::max(max_fepA, chem.conc(i_fepA, c));
  }
  assert(max_btuB > 0.0);
  assert(max_fepA == 0.0);

  std::cout << "  test_cole1_burst_routes_btuB: PASSED\n";
}

void test_colb_burst_routes_fepA() {
  auto sim = make_tiny_sim();
  Vec3 pos = {50e-6, 50e-6, 25e-6};
  Agent a = agent_at(sim, pos, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_B());
  a.state = PhenoState::SOS_INDUCED;
  a.timers.sos_timer = 0.0;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, sim.config().fixes.bacteriocin);
  fix.post_step(1.0);
  sim.step(60.0);

  const auto& chem = sim.chemical_field();
  Int i_btuB = chem.find(species::BACTERIOCIN_BTUB);
  Int i_fepA = chem.find(species::BACTERIOCIN_FEPA);
  assert(i_btuB >= 0 && i_fepA >= 0);

  Real max_btuB = 0.0;
  Real max_fepA = 0.0;
  for (Int c = 0; c < chem.ncells(); ++c) {
    max_btuB = std::max(max_btuB, chem.conc(i_btuB, c));
    max_fepA = std::max(max_fepA, chem.conc(i_fepA, c));
  }
  assert(max_fepA > 0.0);
  assert(max_btuB == 0.0);

  std::cout << "  test_colb_burst_routes_fepA: PASSED\n";
}

void test_receptor_specific_kill() {
  ReceptorConfig rcfg;
  rcfg.kill_rate_colicin = 1.0;

  auto sim = make_tiny_sim();
  Vec3 pos = {50e-6, 50e-6, 25e-6};
  Agent a = agent_at(sim, pos, 2);
  a.receptor_expr[to_underlying(ReceptorType::BtuB)] = 1.0;
  a.receptor_expr[to_underlying(ReceptorType::FepA)] = 0.0;
  a.receptor_expr[to_underlying(ReceptorType::CirA)] = 0.0;
  a.receptor_expr[to_underlying(ReceptorType::FhuA)] = 0.0;
  Int cell = a.grid_cell;
  sim.agents().push_back(std::move(a));

  auto& chem = sim.chemical_field();
  Int i_btuB = chem.find(species::BACTERIOCIN_BTUB);
  Int i_fepA = chem.find(species::BACTERIOCIN_FEPA);
  chem.conc(i_btuB, cell) = 1.0e-4;
  chem.conc(i_fepA, cell) = 1.0e-4;

  FixReceptor fix_btuB(sim, rcfg);
  fix_btuB.compute(60.0);
  assert(sim.agents()[0].state == PhenoState::DEAD);

  auto sim2 = make_tiny_sim();
  Agent b = agent_at(sim2, pos, 2);
  b.receptor_expr[to_underlying(ReceptorType::BtuB)] = 0.0;
  b.receptor_expr[to_underlying(ReceptorType::FepA)] = 1.0;
  b.receptor_expr[to_underlying(ReceptorType::CirA)] = 0.0;
  b.receptor_expr[to_underlying(ReceptorType::FhuA)] = 0.0;
  Int cell2 = b.grid_cell;
  sim2.agents().push_back(std::move(b));
  auto& chem2 = sim2.chemical_field();
  chem2.conc(chem2.find(species::BACTERIOCIN_BTUB), cell2) = 1.0e-4;
  chem2.conc(chem2.find(species::BACTERIOCIN_FEPA), cell2) = 0.0;

  FixReceptor fix_only_btuB(sim2, rcfg);
  fix_only_btuB.compute(60.0);
  assert(sim2.agents()[0].state != PhenoState::DEAD);

  std::cout << "  test_receptor_specific_kill: PASSED\n";
}

void test_nuclease_cross_induction_uses_btuB_field() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 0.0;
  cfg.sos_lysis_prob = 0.0;
  cfg.sos_cross_induction_rate = 1.0e3;

  auto sim = make_tiny_sim();
  Vec3 pos = {50e-6, 50e-6, 25e-6};
  Agent a = agent_at(sim, pos, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  Int cell = a.grid_cell;
  sim.agents().push_back(std::move(a));

  auto& chem = sim.chemical_field();
  chem.conc(chem.find(species::BACTERIOCIN_BTUB), cell) = 1.0e-3;
  chem.conc(chem.find(species::BACTERIOCIN_FEPA), cell) = 1.0e-3;

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);
  assert(sim.agents()[0].state == PhenoState::SOS_INDUCED);

  auto sim2 = make_tiny_sim();
  Agent a2 = agent_at(sim2, pos, 1);
  a2.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  Int cell2 = a2.grid_cell;
  sim2.agents().push_back(std::move(a2));
  auto& chem2 = sim2.chemical_field();
  chem2.conc(chem2.find(species::BACTERIOCIN_BTUB), cell2) = 0.0;
  chem2.conc(chem2.find(species::BACTERIOCIN_FEPA), cell2) = 1.0e-3;

  FixBacteriocin fix2(sim2, cfg);
  fix2.compute(60.0);
  assert(sim2.agents()[0].state == PhenoState::NORMAL);

  std::cout << "  test_nuclease_cross_induction_uses_btuB_field: PASSED\n";
}

int main() {
  std::cout << "=== Per-Receptor Toxin Tests ===\n";
  test_cole1_burst_routes_btuB();
  test_colb_burst_routes_fepA();
  test_receptor_specific_kill();
  test_nuclease_cross_induction_uses_btuB_field();
  std::cout << "All per-receptor toxin tests passed.\n";
  return 0;
}
