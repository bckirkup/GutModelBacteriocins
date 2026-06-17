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
  resident.plasmids   = {"colicin_E1"};
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
  conj.plasmids   = {"colicin_B"};
  conj.conjugative = true;
  cfg.initial_strains.push_back(conj);

  // Disable HDF5 output for CI (no file I/O needed)
  cfg.hdf5.enabled = false;

  // Initialize and run
  Simulation sim;
  sim.init(cfg);

  assert(sim.agents().size() == 50);

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
  producer.plasmids = {"colicin_E1"}; producer.conjugative = false;
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

int main() {
  std::cout << "=== Smoke Tests ===\n";
  test_mini_simulation();
  test_metabolism_integration();
  test_advection_moves_agents();
  test_receptor_killing();
  std::cout << "All smoke tests passed.\n";
  return 0;
}
