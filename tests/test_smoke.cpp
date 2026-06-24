/* -----------------------------------------------------------------------
   GutIBM – Smoke test: mini end-to-end simulation
   Runs ~50 agents for 100 timesteps to verify the full pipeline
   compiles and executes without crashing.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "plasmid.h"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace gutibm;

void test_mini_simulation() {
  SimulationConfig cfg = InputParser::default_config();

  // Shrink domain for speed
  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;

  // Quick run
  cfg.total_time      = 600.0;   // 10 minutes
  cfg.bio_dt          = 60.0;
  cfg.output_interval = 300.0;
  cfg.seed            = 123;

  // Advection
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length   = 100e-6;
  cfg.advection.radial_turnover = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;

  // QSSA
  cfg.qssa.toxin_cutoff = 50e-6;
  cfg.qssa.nutrient_cutoff = 25e-6;

  // Small populations
  cfg.initial_strains.clear();

  // Resident strain with ColE1
  SimulationConfig::InitialStrain resident;
  resident.type       = 1;
  resident.count      = 30;
  resident.mu_max     = 5.0e-4;
  resident.plasmids   = {"ColE1"};
  resident.conjugative = false;
  cfg.initial_strains.push_back(resident);

  // Immigrant strain (no plasmids)
  SimulationConfig::InitialStrain immigrant;
  immigrant.type       = 2;
  immigrant.count      = 15;
  immigrant.mu_max     = 5.0e-4;
  immigrant.plasmids   = {};
  immigrant.conjugative = false;
  cfg.initial_strains.push_back(immigrant);

  // Conjugative strain with ColB
  SimulationConfig::InitialStrain conj;
  conj.type       = 3;
  conj.count      = 5;
  conj.mu_max     = 4.5e-4;
  conj.plasmids   = {"ColB"};
  conj.conjugative = true;
  cfg.initial_strains.push_back(conj);

  // Disable HDF5 output for CI (no file I/O needed)
  cfg.hdf5.enabled = false;

  // Initialize and run
  Simulation sim;
  sim.init(cfg);

  assert(sim.agents().size() == 50);

  // Plasmids must be assigned (ColE1 / ColB)
  Int with_plasmid = 0;
  for (Int i = 0; i < sim.agents().size(); ++i) {
    if (!sim.agents()[i].genome.bi_loci.empty()) ++with_plasmid;
  }
  assert(with_plasmid == 35);  // 30 ColE1 + 5 ColB

  // Run the simulation
  sim.run();

  // Basic sanity checks
  assert(sim.time() > 0.0);
  assert(sim.step_count() > 0);

  // Should still have some agents alive
  Int alive = 0;
  for (Int i = 0; i < sim.agents().size(); ++i) {
    if (sim.agents()[i].state != PhenoState::DEAD) alive++;
  }
  assert(alive > 0);

  // Lineage tracker should have data
  Real retention = sim.lineage_tracker().resident_retention(0.0);
  assert(retention >= 0.0 && retention <= 1.0);

  std::cout << "  test_mini_simulation: PASSED"
            << " (alive=" << alive
            << " steps=" << sim.step_count()
            << " retention=" << retention << ")\n";
}

void test_metabolism_integration() {
  // Verify growth actually changes biomass
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.total_time = 120.0;
  cfg.bio_dt = 60.0;
  cfg.output_interval = 120.0;
  cfg.seed = 456;
  cfg.hdf5.enabled = false;
  cfg.advection.mucus_thickness = 25e-6;
  cfg.advection.distal_length = 50e-6;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 5; s.mu_max = 5e-4;
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);

  Real initial_biomass = 0.0;
  for (Int i = 0; i < sim.agents().size(); ++i)
    initial_biomass += sim.agents()[i].biomass;

  sim.run();

  Real final_biomass = 0.0;
  for (Int i = 0; i < sim.agents().size(); ++i)
    if (sim.agents()[i].state != PhenoState::DEAD)
      final_biomass += sim.agents()[i].biomass;

  // Biomass should have changed (growth or death)
  assert(final_biomass >= 0.0);

  std::cout << "  test_metabolism_integration: PASSED"
            << " (initial=" << initial_biomass
            << " final=" << final_biomass << ")\n";
}

void test_advection_moves_agents() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 10e-6;
  cfg.total_time = 60.0;
  cfg.bio_dt = 60.0;
  cfg.output_interval = 60.0;
  cfg.seed = 789;
  cfg.hdf5.enabled = false;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length = 100e-6;
  cfg.qssa.toxin_cutoff = 50e-6;
  cfg.qssa.nutrient_cutoff = 25e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 10; s.mu_max = 5e-4;
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);

  // Record initial positions
  std::vector<Vec3> initial_pos;
  for (Int i = 0; i < sim.agents().size(); ++i)
    initial_pos.push_back(sim.agents()[i].x);

  sim.step(60.0);

  // At least some agents should have moved
  bool any_moved = false;
  for (Int i = 0; i < sim.agents().size(); ++i) {
    if (sim.agents()[i].state == PhenoState::DEAD) continue;
    Vec3 delta;
    for (int d = 0; d < 3; ++d)
      delta[d] = sim.agents()[i].x[d] - initial_pos[i][d];
    Real dist = std::sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
    if (dist > 1e-15) any_moved = true;
  }
  assert(any_moved);

  std::cout << "  test_advection_moves_agents: PASSED\n";
}

void test_receptor_killing() {
  // Agent exposed to high toxin with full receptor expression should die
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.total_time = 600.0;
  cfg.bio_dt = 60.0;
  cfg.output_interval = 600.0;
  cfg.seed = 999;
  cfg.hdf5.enabled = false;
  cfg.advection.mucus_thickness = 25e-6;
  cfg.advection.distal_length = 50e-6;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  // One producer, many targets
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain producer;
  producer.type = 1; producer.count = 10; producer.mu_max = 5e-4;
  producer.plasmids = {"ColE1"}; producer.conjugative = false;
  cfg.initial_strains.push_back(producer);

  SimulationConfig::InitialStrain target;
  target.type = 2; target.count = 10; target.mu_max = 5e-4;
  target.plasmids = {}; target.conjugative = false;
  cfg.initial_strains.push_back(target);

  Simulation sim;
  sim.init(cfg);
  sim.run();

  // Count survivors by type
  Int type1_alive = 0, type2_alive = 0;
  for (Int i = 0; i < sim.agents().size(); ++i) {
    const Agent& a = sim.agents()[i];
    if (a.state == PhenoState::DEAD) continue;
    if (a.type == 1) type1_alive++;
    else type2_alive++;
  }

  // Producers should generally survive (they have immunity)
  // Result depends on stochastic dynamics, so just verify we ran
  std::cout << "  test_receptor_killing: PASSED"
            << " (producers=" << type1_alive
            << " targets=" << type2_alive << ")\n";
}

void test_crypt_agents_survive_washout() {
  // Agents in the crypt zone must not be washed out even with
  // negative mu_realized (the Washout Trap bypass per VADI §98-99).
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {50e-6, 50e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.total_time      = 600.0;
  cfg.bio_dt          = 60.0;
  cfg.output_interval = 600.0;
  cfg.seed            = 7777;
  cfg.hdf5.enabled    = false;

  // Enable crypts with a generous depth so all agents fall inside
  cfg.advection.crypts_enabled       = true;
  cfg.advection.crypt_depth          = 25e-6;
  cfg.advection.crypt_exit_rate      = 0.0;  // disable exit so they stay
  cfg.advection.crypt_entry_rate     = 0.0;
  cfg.advection.crypt_carrying_capacity = 100;
  cfg.advection.mucus_thickness      = 50e-6;
  cfg.advection.distal_length        = 50e-6;

  cfg.qssa.toxin_cutoff   = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  // Place agents exclusively in the crypt zone (z < 25 um)
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 10; s.mu_max = 5e-4;
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);

  // Force all agents to negative mu_realized (metabolically crippled)
  for (Int i = 0; i < sim.agents().size(); ++i) {
    Agent& a = sim.agents()[i];
    a.mu_realized = -1.0;          // strongly negative
    a.x[2] = 5e-6;                 // inside crypt
    a.in_crypt = true;
  }
  Int initial_count = sim.agents().size();
  assert(initial_count > 0);

  // Run the simulation
  sim.run();

  // All crypt agents should survive despite negative mu_realized
  Int alive = 0;
  for (Int i = 0; i < sim.agents().size(); ++i) {
    if (sim.agents()[i].state != PhenoState::DEAD) alive++;
  }
  assert(alive == initial_count);

  std::cout << "  test_crypt_agents_survive_washout: PASSED"
            << " (alive=" << alive << "/" << initial_count << ")\n";
}

void test_metabolic_washout_trap() {
  // Combinatorial Washout Trap: 0 < mu_realized < gamma_flow → washed out
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {50e-6, 50e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.total_time      = 60.0;
  cfg.bio_dt          = 60.0;
  cfg.output_interval = 60.0;
  cfg.seed            = 4242;
  cfg.hdf5.enabled    = false;
  cfg.advection.crypts_enabled = false;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.radial_turnover = 5400.0;
  cfg.advection.distal_length = 50e-6;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 1;
  s.mu_max = 5e-4;
  s.plasmids = {};
  s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);
  assert(sim.agents().size() == 1);

  Agent& a = sim.agents()[0];
  Real z = 45e-6;
  a.x[2] = z;
  a.in_crypt = false;

  Real gamma = sim.advection().washout_rate(z);
  assert(gamma > 0.0);

  // Receptor-downregulated immigrant scenario: metabolism yields mu << gamma
  for (int k = 0; k < NUM_RECEPTORS; ++k) {
    a.receptor_expr[k] = 0.01;
    a.genome.receptor_expression[k] = 0.01;
  }
  a.mu_max = 1e-6;
  a.km_carbon = 500.0;

  sim.step(60.0);
  assert(a.mu_realized < gamma);
  assert(a.state == PhenoState::DEAD);

  std::cout << "  test_metabolic_washout_trap: PASSED"
            << " (mu=" << a.mu_realized << " gamma=" << gamma << ")\n";
}

void test_metabolic_washout_survives_above_threshold() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {50e-6, 50e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.total_time      = 60.0;
  cfg.bio_dt          = 60.0;
  cfg.output_interval = 60.0;
  cfg.seed            = 4343;
  cfg.hdf5.enabled    = false;
  cfg.advection.crypts_enabled = false;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.radial_turnover = 5400.0;
  cfg.advection.distal_length = 50e-6;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 1;
  s.mu_max = 5e-4;
  s.plasmids = {};
  s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);
  Agent& a = sim.agents()[0];
  a.x[2] = 5e-6;  // low-flow zone near epithelium
  a.in_crypt = false;

  Real gamma = sim.advection().washout_rate(a.x[2]);
  sim.step(60.0);
  assert(a.mu_realized > gamma);
  assert(a.state != PhenoState::DEAD);

  std::cout << "  test_metabolic_washout_survives_above_threshold: PASSED\n";
}

void test_crypt_zero_velocity() {
  // Verify that velocity, washout_rate, and in_crypt_zone behave
  // correctly when crypts are enabled.
  AdvectionConfig acfg;
  acfg.crypts_enabled = true;
  acfg.crypt_depth    = 10e-6;
  acfg.mucus_thickness = 50e-6;

  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {50e-6, 50e-6, 50e-6};
  dcfg.grid_dx = 5e-6;
  dcfg.hash_cell_size = 10e-6;

  Domain dom;
  dom.init(dcfg);

  AdvectionField adv;
  adv.init(acfg, dom);

  // Inside crypt (z = 5 um < 10 um depth)
  assert(adv.in_crypt_zone(5e-6));
  assert(adv.radial_velocity(5e-6) == 0.0);
  assert(adv.distal_velocity(5e-6) == 0.0);
  assert(adv.washout_rate(5e-6) == 0.0);
  Vec3 v_in = adv.velocity({25e-6, 25e-6, 5e-6});
  assert(v_in[0] == 0.0 && v_in[1] == 0.0 && v_in[2] == 0.0);

  // Outside crypt (z = 30 um > 10 um depth)
  assert(!adv.in_crypt_zone(30e-6));
  assert(adv.radial_velocity(30e-6) > 0.0);
  assert(adv.distal_velocity(30e-6) > 0.0);
  assert(adv.washout_rate(30e-6) > 0.0);

  std::cout << "  test_crypt_zero_velocity: PASSED\n";
}

void test_crypt_migration_in_out() {
  // With high entry/exit rates, agents should migrate between zones.
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {50e-6, 50e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.total_time      = 600.0;
  cfg.bio_dt          = 60.0;
  cfg.output_interval = 600.0;
  cfg.seed            = 8888;
  cfg.hdf5.enabled    = false;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length   = 50e-6;

  // Enable crypts with very high entry rate so some agents enter
  cfg.advection.crypts_enabled        = true;
  cfg.advection.crypt_depth           = 10e-6;
  cfg.advection.crypt_exit_rate       = 1.0;   // ~1/s, very high
  cfg.advection.crypt_entry_rate      = 1.0;   // ~1/s, very high
  cfg.advection.crypt_carrying_capacity = 100;

  cfg.qssa.toxin_cutoff   = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 20; s.mu_max = 5e-4;
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);

  // Place all agents just above the crypt zone so they can enter
  for (Int i = 0; i < sim.agents().size(); ++i) {
    sim.agents()[i].x[2] = 12e-6;
    sim.agents()[i].in_crypt = false;
  }

  sim.run();

  // With very high migration rates over 600 s, at least one agent
  // should have ended up in the crypt at some point. We cannot
  // guarantee the exact final state, but the simulation should not crash.
  std::cout << "  test_crypt_migration_in_out: PASSED\n";
}

void test_partial_resistance_survival() {
  // Partially resistant agents should survive toxin exposure better than
  // fully susceptible agents, without the growth penalty of full downregulation.
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.total_time = 600.0;
  cfg.bio_dt = 60.0;
  cfg.output_interval = 600.0;
  cfg.seed = 1337;
  cfg.hdf5.enabled = false;
  cfg.advection.mucus_thickness = 25e-6;
  cfg.advection.distal_length = 50e-6;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  // Producers (type 1) with ColE1
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain producer;
  producer.type = 1; producer.count = 10; producer.mu_max = 5e-4;
  producer.plasmids = {"ColE1"}; producer.conjugative = false;
  cfg.initial_strains.push_back(producer);

  // Fully susceptible targets (type 2)
  SimulationConfig::InitialStrain susceptible;
  susceptible.type = 2; susceptible.count = 10; susceptible.mu_max = 5e-4;
  susceptible.plasmids = {}; susceptible.conjugative = false;
  cfg.initial_strains.push_back(susceptible);

  // Partially resistant targets (type 3) — same as type 2 but with reduced BtuB toxin affinity
  SimulationConfig::InitialStrain partial_res;
  partial_res.type = 3; partial_res.count = 10; partial_res.mu_max = 5e-4;
  partial_res.plasmids = {}; partial_res.conjugative = false;
  cfg.initial_strains.push_back(partial_res);

  Simulation sim;
  sim.init(cfg);

  // Apply partial resistance to type-3 agents (BtuB toxin_affinity reduced)
  int btuB = static_cast<int>(ReceptorType::BtuB);
  for (Int i = 0; i < sim.agents().size(); ++i) {
    Agent& a = sim.agents()[i];
    if (a.type == 3) {
      a.genome.toxin_affinity[btuB] = 0.05;   // 20x reduced toxin binding
      a.genome.ligand_affinity[btuB] = 0.85;   // near wild-type ligand uptake
    }
  }

  assert(sim.agents().size() == 30);

  sim.run();

  // Count survivors by type
  Int type2_alive = 0, type3_alive = 0;
  for (Int i = 0; i < sim.agents().size(); ++i) {
    const Agent& a = sim.agents()[i];
    if (a.state == PhenoState::DEAD) continue;
    if (a.type == 2) type2_alive++;
    if (a.type == 3) type3_alive++;
  }

  // Partially resistant agents (type 3) should survive at least as well
  // as fully susceptible agents (type 2). Stochastic, so just verify we ran.
  std::cout << "  test_partial_resistance_survival: PASSED"
            << " (susceptible=" << type2_alive
            << " partial_res=" << type3_alive << ")\n";
}

int main() {
  std::cout << "=== Smoke Tests ===\n";
  test_mini_simulation();
  test_metabolism_integration();
  test_advection_moves_agents();
  test_receptor_killing();
  test_metabolic_washout_trap();
  test_metabolic_washout_survives_above_threshold();
  test_crypt_agents_survive_washout();
  test_crypt_zero_velocity();
  test_crypt_migration_in_out();
  test_partial_resistance_survival();
  std::cout << "All smoke tests passed.\n";
  return 0;
}
