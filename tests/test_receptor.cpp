/* -----------------------------------------------------------------------
   GutIBM – Unit tests for fix_receptor
   ----------------------------------------------------------------------- */

#include "fix_receptor.h"
#include "plasmid.h"
#include "simulation.h"
#include "input_parser.h"
#include "species_names.h"
#include "qssa_solver.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

using namespace gutibm;

static Simulation make_empty_sim(uint64_t seed = 42, bool siderophore = false) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.seed = seed;
  cfg.chem_env.siderophore.enabled = siderophore;

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
  Int ix;
  Int iy;
  Int iz;
  sim.domain().pos_to_grid(a.x, ix, iy, iz);
  a.grid_cell = sim.domain().cell_index(ix, iy, iz);
  return a;
}

static void set_local_chemistry(Simulation& sim, Int cell,
                                Real tox_conc, Real b12_conc) {
  auto& chem = sim.chemical_field();
  Int i_tox = chem.find(species::BACTERIOCIN_BTUB);
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

void test_cira_uses_ferric_enterobactin_ligand() {
  ReceptorConfig rcfg;
  rcfg.kill_rate_microcin = 1.0;

  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.seed = 7010;
  cfg.chem_env.siderophore.enabled = true;
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);

  Agent victim = make_susceptible_agent(sim);
  victim.receptor_expr[to_underlying(ReceptorType::BtuB)] = 0.0;
  victim.receptor_expr[to_underlying(ReceptorType::FepA)] = 0.0;
  victim.receptor_expr[to_underlying(ReceptorType::CirA)] = 1.0;
  victim.receptor_expr[to_underlying(ReceptorType::FhuA)] = 0.0;
  Int cell = victim.grid_cell;

  auto& chem = sim.chemical_field();
  Int i_cira = chem.find(species::BACTERIOCIN_CIRA);
  Int i_ferric_enterobactin = chem.find(species::FERRIC_ENTEROBACTIN);
  assert(i_cira >= 0 && i_ferric_enterobactin >= 0);
  chem.conc(i_cira, cell) = 1.0e-4;
  chem.conc(i_ferric_enterobactin, cell) = 1.0e-3;
  sim.agents().push_back(std::move(victim));

  FixReceptor fix(sim, rcfg);
  fix.compute(60.0);
  assert(sim.agents()[0].state == PhenoState::DEAD);

  std::cout << "  test_cira_uses_ferric_enterobactin_ligand: PASSED\n";
}

void test_fepa_uses_ferric_enterobactin_not_iron() {
  ReceptorConfig rcfg;
  rcfg.kill_rate_colicin = 1.0;

  auto sim_fe = make_empty_sim(7020, true);
  Agent fe_agent = make_susceptible_agent(sim_fe);
  fe_agent.receptor_expr[to_underlying(ReceptorType::BtuB)] = 0.0;
  fe_agent.receptor_expr[to_underlying(ReceptorType::FepA)] = 1.0;
  fe_agent.receptor_expr[to_underlying(ReceptorType::CirA)] = 0.0;
  fe_agent.receptor_expr[to_underlying(ReceptorType::FhuA)] = 0.0;
  const Int fe_cell = fe_agent.grid_cell;
  auto& fe_chem = sim_fe.chemical_field();
  const Int fe_toxin = fe_chem.find(species::BACTERIOCIN_FEPA);
  const Int fe_ligand = fe_chem.find(species::FERRIC_ENTEROBACTIN);
  const Int fe_iron = fe_chem.find(species::IRON);
  assert(fe_toxin >= 0 && fe_ligand >= 0 && fe_iron >= 0);
  fe_chem.conc(fe_toxin, fe_cell) = 1.0;
  fe_chem.conc(fe_ligand, fe_cell) = 1.0e3;
  fe_chem.conc(fe_iron, fe_cell) = 1.0e3;
  sim_fe.agents().push_back(std::move(fe_agent));

  auto sim_iron = make_empty_sim(7020, true);
  Agent iron_agent = make_susceptible_agent(sim_iron);
  iron_agent.receptor_expr[to_underlying(ReceptorType::BtuB)] = 0.0;
  iron_agent.receptor_expr[to_underlying(ReceptorType::FepA)] = 1.0;
  iron_agent.receptor_expr[to_underlying(ReceptorType::CirA)] = 0.0;
  iron_agent.receptor_expr[to_underlying(ReceptorType::FhuA)] = 0.0;
  const Int iron_cell = iron_agent.grid_cell;
  auto& iron_chem = sim_iron.chemical_field();
  const Int iron_toxin = iron_chem.find(species::BACTERIOCIN_FEPA);
  const Int iron_ligand = iron_chem.find(species::FERRIC_ENTEROBACTIN);
  const Int iron_bulk = iron_chem.find(species::IRON);
  assert(iron_toxin >= 0 && iron_ligand >= 0 && iron_bulk >= 0);
  iron_chem.conc(iron_toxin, iron_cell) = 1.0;
  iron_chem.conc(iron_bulk, iron_cell) = 1.0e3;
  sim_iron.agents().push_back(std::move(iron_agent));

  FixReceptor fe_fix(sim_fe, rcfg);
  FixReceptor iron_fix(sim_iron, rcfg);
  fe_fix.compute(60.0);
  iron_fix.compute(60.0);

  assert(sim_fe.agents()[0].state != PhenoState::DEAD);
  assert(sim_iron.agents()[0].state == PhenoState::DEAD);
  std::cout << "  test_fepa_uses_ferric_enterobactin_not_iron: PASSED\n";
}

void test_burst_kills_same_step() {
  // Regression: fix_receptor must see toxin deposited by QSSA in the same step.
  // Before the fix, module_biology ran receptor before module_chemistry, so a
  // newly added burst left the grid at zero during killing.
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.enabled_fixes = {"receptor"};
  cfg.fixes.receptor.kill_rate_colicin = 10.0;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.seed = 6060;

  Simulation sim;
  sim.init(cfg);

  Agent victim = make_susceptible_agent(sim);
  const Vec3 pos = victim.x;
  const Int cell = victim.grid_cell;
  sim.agents().push_back(std::move(victim));

  const BICluster col_e1 = PlasmidLibrary::colicin_E1();
  ToxinBurstSource burst;
  burst.pos = pos;
  burst.params.diff_coeff = col_e1.diff_coeff;
  burst.params.retardation = col_e1.retardation;
  burst.params.pI = col_e1.pI;
  burst.params.source_rate = sim.config().qssa.colicin_release_rate;
  burst.creation_time = sim.time();
  burst.is_nuclease = col_e1.is_nuclease;
  burst.target = col_e1.target;
  sim.add_toxin_burst(burst);

  Int i_tox = sim.chemical_field().find(species::BACTERIOCIN_BTUB);
  assert(i_tox >= 0);
  assert(sim.chemical_field().conc(i_tox, cell) == 0.0);

  sim.step(60.0);

  assert(sim.agents()[0].state == PhenoState::DEAD);

  std::cout << "  test_burst_kills_same_step: PASSED\n";
}

void test_true_unit_receptor_regression() {
  constexpr Real toxin_at_10_um = 4.421248e-6;
  constexpr Real toxin_at_50_um = 3.676568e-7;
  constexpr Real toxin_at_100_um = 6.137581e-8;

  const SimulationConfig defaults = InputParser::default_config();
  const auto b12_spec = std::find_if(
      defaults.chemicals.begin(), defaults.chemicals.end(),
      [](const ChemicalSpec& s) { return s.name == species::B12; });
  assert(b12_spec != defaults.chemicals.end());
  assert(std::abs(b12_spec->initial_conc - 1.0e-3) < 1.0e-15);

  const ReceptorConfig& receptor = defaults.fixes.receptor;
  const Real b12 = 1.0e-3;
  const Real b12_factor = 1.0 + b12 / receptor.kd_b12_btuB;
  const Real apparent_kd = receptor.kd_colicinE_btuB * b12_factor;
  assert(std::abs(b12_factor - 1001.0) < 1.0e-12);
  assert(std::abs(apparent_kd - 5.005e-4) < 1.0e-12);

  const Real occ_10 = toxin_at_10_um / (apparent_kd + toxin_at_10_um);
  const Real occ_50 = toxin_at_50_um / (apparent_kd + toxin_at_50_um);
  const Real occ_100 = toxin_at_100_um / (apparent_kd + toxin_at_100_um);
  assert(std::abs(occ_10 - 0.00875631) < 1.0e-7);
  assert(std::abs(occ_50 - 0.00073404) < 1.0e-7);
  assert(std::abs(occ_100 - 0.00012261) < 1.0e-7);

  const Real fepa_factor = 1.0 + 1.0e-4 / receptor.kd_enterobactin;
  assert(std::abs(fepa_factor - 101.0) < 1.0e-12);
  const Real old_fepa_factor = 1.0 + 1.0e-4 / 1.0e-9;
  assert(std::abs(old_fepa_factor - 100001.0) < 1.0e-9);

  const Real b12_monod = b12 / (receptor.kd_b12_btuB + b12);
  const Real old_b12_monod = 1.0e-6 / (1.0e-9 + 1.0e-6);
  assert(std::abs(b12_monod - old_b12_monod) < 1.0e-12);

  auto corrected = make_empty_sim(8118);
  Agent corrected_agent = make_susceptible_agent(corrected);
  set_local_chemistry(corrected, corrected_agent.grid_cell, toxin_at_50_um, b12);
  corrected.agents().push_back(std::move(corrected_agent));
  ReceptorConfig corrected_cfg;
  corrected_cfg.kill_rate_colicin = 1.0;
  FixReceptor corrected_fix(corrected, corrected_cfg);
  corrected_fix.compute(60.0);
  assert(corrected.agents()[0].state != PhenoState::DEAD);

  auto legacy = make_empty_sim(8118);
  Agent legacy_agent = make_susceptible_agent(legacy);
  set_local_chemistry(legacy, legacy_agent.grid_cell, toxin_at_50_um, 1.0e-6);
  legacy.agents().push_back(std::move(legacy_agent));
  ReceptorConfig legacy_cfg;
  legacy_cfg.kill_rate_colicin = 1.0;
  legacy_cfg.kd_b12_btuB = 1.0e-9;
  legacy_cfg.kd_colicinE_btuB = 5.0e-10;
  FixReceptor legacy_fix(legacy, legacy_cfg);
  legacy_fix.compute(60.0);
  assert(legacy.agents()[0].state == PhenoState::DEAD);

  std::cout << "  test_true_unit_receptor_regression: PASSED\n";
}

int main() {
  std::cout << "=== Receptor Fix Tests ===\n";
  test_high_toxin_kills_susceptible();
  test_immunity_reduces_lethality();
  test_ligand_competition_reduces_kill();
  test_partial_resistance_reduces_lethality();
  test_cira_uses_ferric_enterobactin_ligand();
  test_fepa_uses_ferric_enterobactin_not_iron();
  test_burst_kills_same_step();
  test_true_unit_receptor_regression();
  std::cout << "All receptor fix tests passed.\n";
  return 0;
}
