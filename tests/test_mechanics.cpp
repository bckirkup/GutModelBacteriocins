/* -----------------------------------------------------------------------
   GutIBM – Mechanical repulsion tests (Issue #16)
   Tests Hertzian contact model, force scaling, and EPS adhesion.
   ----------------------------------------------------------------------- */

#include "agent.h"
#include "fix_mechanics.h"
#include "simulation.h"
#include "input_parser.h"

#include <cassert>
#include <iostream>
#include <cmath>

using namespace gutibm;

// Helper: create a minimal simulation with two overlapping agents
static Simulation make_two_agent_sim(Vec3 pos_a, Vec3 pos_b,
                                     MechanicsConfig mcfg = {}) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.domain.hi = {100e-6, 100e-6, 100e-6};
  cfg.domain.grid_dx = 10e-6;
  cfg.domain.hash_cell_size = 20e-6;
  cfg.mechanics = mcfg;
  cfg.hdf5.enabled = false;

  Simulation sim;
  sim.init(cfg);

  Agent a = Agent::create_default(sim.agents().next_tag(), 1, pos_a, 5e-4);
  Agent b = Agent::create_default(sim.agents().next_tag(), 1, pos_b, 5e-4);
  sim.agents().push_back(std::move(a));
  sim.agents().push_back(std::move(b));

  return sim;
}

void test_overlapping_agents_pushed_apart() {
  // Place two cells overlapping by 0.2 um along x-axis
  Real r = CELL_RADIUS_DEFAULT;  // 0.5 um
  Real separation = 2 * r - 0.2e-6;  // 0.8 um apart, overlap = 0.2 um
  Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b = {50e-6 + separation, 50e-6, 50e-6};

  MechanicsConfig mcfg;
  mcfg.hertzian_enabled = true;
  mcfg.hertz_k = 1.0e-6;

  auto sim = make_two_agent_sim(pos_a, pos_b, mcfg);

  // Rebuild spatial hash so neighbors are found
  sim.domain().spatial_hash().clear();
  sim.domain().spatial_hash().insert(0, sim.agents()[0].x);
  sim.domain().spatial_hash().insert(1, sim.agents()[1].x);

  Real x_a_before = sim.agents()[0].x[0];
  Real x_b_before = sim.agents()[1].x[0];
  Real dist_before = x_b_before - x_a_before;

  // Run mechanics fix
  Real dt = 1.0;  // 1 second timestep
  // Directly construct and call the fix
  FixMechanics fix(sim, mcfg);
  fix.compute(dt);

  Real x_a_after = sim.agents()[0].x[0];
  Real x_b_after = sim.agents()[1].x[0];
  Real dist_after = x_b_after - x_a_after;

  // Agents should have moved apart
  assert(dist_after > dist_before);
  // Agent a should have moved left (negative x)
  assert(x_a_after < x_a_before);
  // Agent b should have moved right (positive x)
  assert(x_b_after > x_b_before);

  std::cout << "  test_overlapping_agents_pushed_apart: PASSED\n";
}

void test_hertzian_force_scaling() {
  // Verify F scales as overlap^1.5
  Real r = CELL_RADIUS_DEFAULT;
  Real k = 1.0e-6;

  // Test two different overlaps and verify ratio matches power law
  Real overlap1 = 0.1e-6;
  Real overlap2 = 0.4e-6;

  Real F1 = k * std::pow(overlap1, 1.5);
  Real F2 = k * std::pow(overlap2, 1.5);

  Real ratio_actual = F2 / F1;
  Real ratio_expected = std::pow(overlap2 / overlap1, 1.5);

  assert(std::abs(ratio_actual - ratio_expected) / ratio_expected < 1e-10);

  // Now verify through simulation: larger overlap should give larger displacement
  Real sep1 = 2 * r - overlap1;
  Real sep2 = 2 * r - overlap2;

  Vec3 center = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b1 = {50e-6 + sep1, 50e-6, 50e-6};
  Vec3 pos_b2 = {50e-6 + sep2, 50e-6, 50e-6};

  MechanicsConfig mcfg;
  mcfg.hertzian_enabled = true;
  mcfg.hertz_k = k;

  // Sim 1: small overlap
  auto sim1 = make_two_agent_sim(center, pos_b1, mcfg);
  sim1.domain().spatial_hash().clear();
  sim1.domain().spatial_hash().insert(0, sim1.agents()[0].x);
  sim1.domain().spatial_hash().insert(1, sim1.agents()[1].x);

  Real xa1_before = sim1.agents()[0].x[0];
  FixMechanics fix1(sim1, mcfg);
  fix1.compute(1.0);
  Real displacement1 = xa1_before - sim1.agents()[0].x[0];

  // Sim 2: large overlap
  auto sim2 = make_two_agent_sim(center, pos_b2, mcfg);
  sim2.domain().spatial_hash().clear();
  sim2.domain().spatial_hash().insert(0, sim2.agents()[0].x);
  sim2.domain().spatial_hash().insert(1, sim2.agents()[1].x);

  Real xa2_before = sim2.agents()[0].x[0];
  FixMechanics fix2(sim2, mcfg);
  fix2.compute(1.0);
  Real displacement2 = xa2_before - sim2.agents()[0].x[0];

  // Displacement ratio should approximate the force ratio (overlap^1.5 scaling)
  Real disp_ratio = displacement2 / displacement1;
  // For equal-mass pair: displacement proportional to force
  assert(disp_ratio > 1.0);  // larger overlap → larger displacement
  // Check it's close to the Hertzian ratio (within 5% tolerance for numerical)
  assert(std::abs(disp_ratio - ratio_expected) / ratio_expected < 0.05);

  std::cout << "  test_hertzian_force_scaling: PASSED\n";
}

void test_no_force_without_overlap() {
  // Agents separated by more than sum of radii → no force
  Real r = CELL_RADIUS_DEFAULT;
  Real separation = 2 * r + 1.0e-6;  // 1 um gap
  Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b = {50e-6 + separation, 50e-6, 50e-6};

  MechanicsConfig mcfg;
  mcfg.hertzian_enabled = true;
  mcfg.adhesion_enabled = false;

  auto sim = make_two_agent_sim(pos_a, pos_b, mcfg);
  sim.domain().spatial_hash().clear();
  sim.domain().spatial_hash().insert(0, sim.agents()[0].x);
  sim.domain().spatial_hash().insert(1, sim.agents()[1].x);

  Real x_a_before = sim.agents()[0].x[0];
  Real x_b_before = sim.agents()[1].x[0];

  FixMechanics fix(sim, mcfg);
  fix.compute(1.0);

  // No movement since no overlap and no adhesion
  assert(std::abs(sim.agents()[0].x[0] - x_a_before) < 1e-30);
  assert(std::abs(sim.agents()[1].x[0] - x_b_before) < 1e-30);

  std::cout << "  test_no_force_without_overlap: PASSED\n";
}

void test_adhesion_holds_agents() {
  // Place two agents with a small gap, adhesion should pull them together
  Real r = CELL_RADIUS_DEFAULT;
  Real gap = 0.2e-6;  // 0.2 um gap (within adhesion range)
  Real separation = 2 * r + gap;
  Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b = {50e-6 + separation, 50e-6, 50e-6};

  MechanicsConfig mcfg;
  mcfg.hertzian_enabled = true;
  mcfg.adhesion_enabled = true;
  mcfg.adhesion_strength = 1.0e-12;
  mcfg.adhesion_range = 0.5e-6;  // 0.5 um range

  auto sim = make_two_agent_sim(pos_a, pos_b, mcfg);
  sim.domain().spatial_hash().clear();
  sim.domain().spatial_hash().insert(0, sim.agents()[0].x);
  sim.domain().spatial_hash().insert(1, sim.agents()[1].x);

  Real dist_before = sim.agents()[1].x[0] - sim.agents()[0].x[0];

  FixMechanics fix(sim, mcfg);
  fix.compute(1.0);

  Real dist_after = sim.agents()[1].x[0] - sim.agents()[0].x[0];

  // Adhesion should pull agents closer together
  assert(dist_after < dist_before);

  std::cout << "  test_adhesion_holds_agents: PASSED\n";
}

void test_no_adhesion_beyond_range() {
  // Agents far apart → no adhesion effect
  Real r = CELL_RADIUS_DEFAULT;
  Real gap = 1.0e-6;  // 1 um gap (beyond 0.5 um adhesion range)
  Real separation = 2 * r + gap;
  Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b = {50e-6 + separation, 50e-6, 50e-6};

  MechanicsConfig mcfg;
  mcfg.hertzian_enabled = true;
  mcfg.adhesion_enabled = true;
  mcfg.adhesion_strength = 1.0e-12;
  mcfg.adhesion_range = 0.5e-6;

  auto sim = make_two_agent_sim(pos_a, pos_b, mcfg);
  sim.domain().spatial_hash().clear();
  sim.domain().spatial_hash().insert(0, sim.agents()[0].x);
  sim.domain().spatial_hash().insert(1, sim.agents()[1].x);

  Real x_a_before = sim.agents()[0].x[0];
  Real x_b_before = sim.agents()[1].x[0];

  FixMechanics fix(sim, mcfg);
  fix.compute(1.0);

  // No movement
  assert(std::abs(sim.agents()[0].x[0] - x_a_before) < 1e-30);
  assert(std::abs(sim.agents()[1].x[0] - x_b_before) < 1e-30);

  std::cout << "  test_no_adhesion_beyond_range: PASSED\n";
}

void test_dead_agents_ignored() {
  Real r = CELL_RADIUS_DEFAULT;
  Real separation = 2 * r - 0.2e-6;  // overlapping
  Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b = {50e-6 + separation, 50e-6, 50e-6};

  MechanicsConfig mcfg;
  mcfg.hertzian_enabled = true;

  auto sim = make_two_agent_sim(pos_a, pos_b, mcfg);
  // Mark second agent as dead
  sim.agents()[1].state = PhenoState::DEAD;

  sim.domain().spatial_hash().clear();
  sim.domain().spatial_hash().insert(0, sim.agents()[0].x);
  sim.domain().spatial_hash().insert(1, sim.agents()[1].x);

  Real x_a_before = sim.agents()[0].x[0];

  FixMechanics fix(sim, mcfg);
  fix.compute(1.0);

  // No movement because one agent is dead
  assert(std::abs(sim.agents()[0].x[0] - x_a_before) < 1e-30);

  std::cout << "  test_dead_agents_ignored: PASSED\n";
}

int main() {
  std::cout << "=== Mechanics Tests (Issue #16) ===\n";
  test_overlapping_agents_pushed_apart();
  test_hertzian_force_scaling();
  test_no_force_without_overlap();
  test_adhesion_holds_agents();
  test_no_adhesion_beyond_range();
  test_dead_agents_ignored();
  std::cout << "All mechanics tests passed.\n";
  return 0;
}
