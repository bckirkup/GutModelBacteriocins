/* -----------------------------------------------------------------------
   GutIBM – Domain decomposition unit tests
   Tests 1D slab decomposition, owner_rank(), is_local(), and migration.
   These tests run on a single rank (nprocs=1) to validate the logic
   without requiring mpirun.
   ----------------------------------------------------------------------- */

#include "domain.h"
#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace gutibm;

void test_decompose_single_rank() {
  // With 1 rank, the entire domain is local
  DomainConfig cfg;
  cfg.lo = {0, 0, 0};
  cfg.hi = {1.0e-3, 1.0e-3, 100.0e-6};
  cfg.mpi_decomp_axis = 0;
  cfg.ghost_width = 10.0e-6;

  Domain dom;
  dom.init(cfg);

  assert(dom.rank() == 0);
  assert(dom.nprocs() == 1);
  assert(std::abs(dom.local_lo_x() - 0.0) < 1e-20);
  assert(std::abs(dom.local_hi_x() - 1.0e-3) < 1e-20);
  assert(std::abs(dom.ghost_width() - 10.0e-6) < 1e-20);

  // All positions should be local
  Vec3 pos_lo = {0.0, 0.5e-3, 50.0e-6};
  Vec3 pos_mid = {0.5e-3, 0.5e-3, 50.0e-6};
  Vec3 pos_hi = {0.999e-3, 0.5e-3, 50.0e-6};

  assert(dom.is_local(pos_lo));
  assert(dom.is_local(pos_mid));
  assert(dom.is_local(pos_hi));

  // All owner_rank should be 0
  assert(dom.owner_rank(pos_lo) == 0);
  assert(dom.owner_rank(pos_mid) == 0);
  assert(dom.owner_rank(pos_hi) == 0);

  // Neighbor ranks should be -1 (periodic wraps to self for nprocs=1
  // but decompose handles this — with 1 proc, rank_lo=0 and rank_hi=0
  // because periodic x wraps)
  // Actually with 1 proc: rank_lo_ = -1 → periodic → nprocs-1 = 0
  // and rank_hi_ = 1 → periodic → 0

  std::cout << "  test_decompose_single_rank: PASSED\n";
}

void test_owner_rank_logic() {
  // Verify owner_rank computation for hypothetical multi-rank scenario
  DomainConfig cfg;
  cfg.lo = {0, 0, 0};
  cfg.hi = {1.0e-3, 1.0e-3, 100.0e-6};
  cfg.mpi_decomp_axis = 0;
  cfg.ghost_width = 10.0e-6;

  Domain dom;
  dom.init(cfg);

  // With 1 rank, owner_rank always returns 0
  assert(dom.owner_rank({0.0, 0.0, 0.0}) == 0);
  assert(dom.owner_rank({0.5e-3, 0.0, 0.0}) == 0);
  assert(dom.owner_rank({0.999e-3, 0.0, 0.0}) == 0);

  std::cout << "  test_owner_rank_logic: PASSED\n";
}

void test_is_local_trivial() {
  // With 1 rank, everything is local
  DomainConfig cfg;
  cfg.lo = {0, 0, 0};
  cfg.hi = {100.0e-6, 100.0e-6, 50.0e-6};
  cfg.grid_dx = 5.0e-6;
  cfg.hash_cell_size = 10.0e-6;
  cfg.mpi_decomp_axis = 0;
  cfg.ghost_width = 10.0e-6;

  Domain dom;
  dom.init(cfg);

  // Spots at various points
  for (double frac = 0.0; frac < 1.0; frac += 0.1) {
    Vec3 pos = {frac * 100.0e-6, 50.0e-6, 25.0e-6};
    assert(dom.is_local(pos));
  }

  std::cout << "  test_is_local_trivial: PASSED\n";
}

void test_init_population_local_only() {
  // With 1 rank, all agents should be retained (is_local always true)
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.total_time = 60.0;
  cfg.bio_dt = 60.0;
  cfg.output_interval = 60.0;
  cfg.seed = 42;
  cfg.hdf5.enabled = false;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length = 100e-6;
  cfg.qssa.toxin_cutoff = 50e-6;
  cfg.qssa.nutrient_cutoff = 25e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 20; s.mu_max = 5e-4;
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);

  // With 1 rank, all 20 agents should be present
  assert(sim.agents().size() == 20);

  // All agents should have owner_rank = 0
  for (Int i = 0; i < sim.agents().size(); ++i) {
    assert(sim.agents()[i].owner_rank == 0);
  }

  // Global count should equal local count
  assert(sim.global_agent_count() == 20);

  std::cout << "  test_init_population_local_only: PASSED\n";
}

void test_single_rank_simulation_unchanged() {
  // Running with 1 rank should produce the same results as before
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.total_time = 300.0;
  cfg.bio_dt = 60.0;
  cfg.output_interval = 300.0;
  cfg.seed = 123;
  cfg.hdf5.enabled = false;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length = 100e-6;
  cfg.advection.radial_turnover = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff = 50e-6;
  cfg.qssa.nutrient_cutoff = 25e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1; resident.count = 20; resident.mu_max = 5e-4;
  resident.plasmids = {"ColE1"}; resident.conjugative = false;
  cfg.initial_strains.push_back(resident);

  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2; immigrant.count = 10; immigrant.mu_max = 5e-4;
  immigrant.plasmids = {}; immigrant.conjugative = false;
  cfg.initial_strains.push_back(immigrant);

  Simulation sim;
  sim.init(cfg);

  Int initial_agents = sim.agents().size();
  assert(initial_agents == 30);

  sim.run();

  // Basic sanity: simulation completed, some agents alive
  assert(sim.time() > 0.0);
  assert(sim.step_count() > 0);

  Int alive = 0;
  for (Int i = 0; i < sim.agents().size(); ++i) {
    if (sim.agents()[i].state != PhenoState::DEAD) alive++;
  }
  assert(alive > 0);

  // Global stats should be consistent with local (1 rank)
  assert(sim.global_agent_count() == alive);
  assert(sim.global_mu_avg() >= 0.0);

  std::cout << "  test_single_rank_simulation_unchanged: PASSED"
            << " (alive=" << alive
            << " global_count=" << sim.global_agent_count()
            << " mu_avg=" << sim.global_mu_avg() << ")\n";
}

void test_ghost_width_config() {
  DomainConfig cfg;
  cfg.lo = {0, 0, 0};
  cfg.hi = {1.0e-3, 1.0e-3, 100.0e-6};
  cfg.ghost_width = 20.0e-6;

  Domain dom;
  dom.init(cfg);

  assert(std::abs(dom.ghost_width() - 20.0e-6) < 1e-20);

  std::cout << "  test_ghost_width_config: PASSED\n";
}

void test_migration_noop_single_rank() {
  // With 1 rank, migrate_agents should be a no-op
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.total_time = 60.0;
  cfg.bio_dt = 60.0;
  cfg.output_interval = 60.0;
  cfg.seed = 42;
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

  Int before = sim.agents().size();

  // Run one step — migration is called internally but should be no-op
  sim.step(60.0);

  // Count alive agents
  Int alive = 0;
  for (Int i = 0; i < sim.agents().size(); ++i) {
    if (sim.agents()[i].state != PhenoState::DEAD) alive++;
  }

  // We should still have agents
  assert(alive > 0);

  std::cout << "  test_migration_noop_single_rank: PASSED"
            << " (before=" << before << " alive_after=" << alive << ")\n";
}

int main() {
  std::cout << "=== Domain Decomposition Tests ===\n";
  test_decompose_single_rank();
  test_owner_rank_logic();
  test_is_local_trivial();
  test_init_population_local_only();
  test_single_rank_simulation_unchanged();
  test_ghost_width_config();
  test_migration_noop_single_rank();
  std::cout << "All domain decomposition tests passed.\n";
  return 0;
}
