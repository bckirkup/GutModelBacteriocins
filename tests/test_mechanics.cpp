/* -----------------------------------------------------------------------
   GutIBM – Mechanical repulsion tests (Issue #16)
   Tests Hertzian contact model, force scaling, and EPS adhesion.
   ----------------------------------------------------------------------- */

#include "agent.h"
#include "fix_mechanics.h"
#include "simulation.h"
#include "input_parser.h"
#include "mechanics_test_helpers.h"

#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>

using namespace gutibm;
using gutibm::test::make_two_agent_sim;

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

void test_displacement_is_bounded_and_counted() {
  Real r = CELL_RADIUS_DEFAULT;
  Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b = {50e-6 + 2 * r - 0.2e-6, 50e-6, 50e-6};
  MechanicsConfig mcfg;
  mcfg.hertz_k = 1.0e12;
  auto sim = make_two_agent_sim(pos_a, pos_b, mcfg);
  const Real before = sim.agents()[0].x[0];
  FixMechanics fix(sim, mcfg);
  fix.compute(60.0);
  const Real displacement = std::abs(sim.agents()[0].x[0] - before);
  assert(displacement <= 0.1 * r * (1.0 + 1e-12));
  assert(sim.mechanics_step_stats().displacement_clamps > 0);
  std::cout << "  test_displacement_is_bounded_and_counted: PASSED\n";
}

void test_calm_pair_has_no_clamp() {
  Real r = CELL_RADIUS_DEFAULT;
  Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b = {50e-6 + 2 * r + 1e-6, 50e-6, 50e-6};
  MechanicsConfig mcfg;
  auto sim = make_two_agent_sim(pos_a, pos_b, mcfg);
  FixMechanics fix(sim, mcfg);
  fix.compute(60.0);
  assert(sim.mechanics_step_stats().displacement_clamps == 0);
  std::cout << "  test_calm_pair_has_no_clamp: PASSED\n";
}

void test_viscosity_slows_relaxation() {
  Real r = CELL_RADIUS_DEFAULT;
  const Real separation = 2 * r - 0.2e-6;
  Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b = {50e-6 + separation, 50e-6, 50e-6};
  MechanicsConfig mcfg;
  auto low_viscosity = make_two_agent_sim(pos_a, pos_b, mcfg, false, 0.01);
  auto high_viscosity = make_two_agent_sim(pos_a, pos_b, mcfg, false, 0.02);
  FixMechanics low_fix(low_viscosity, mcfg);
  FixMechanics high_fix(high_viscosity, mcfg);
  low_fix.compute(1.0);
  high_fix.compute(1.0);
  const Real low_move = pos_a[0] - low_viscosity.agents()[0].x[0];
  const Real high_move = pos_a[0] - high_viscosity.agents()[0].x[0];
  assert(low_move > high_move);
  assert(low_move > 1.8 * high_move);
  std::cout << "  test_viscosity_slows_relaxation: PASSED\n";
}

void test_floor_is_contained() {
  Real r = CELL_RADIUS_DEFAULT;
  Vec3 pos_a = {50e-6, 50e-6, 0.1e-6};
  Vec3 pos_b = {50e-6, 50e-6, 0.9e-6};
  MechanicsConfig mcfg;
  auto sim = make_two_agent_sim(pos_a, pos_b, mcfg);
  FixMechanics fix(sim, mcfg);
  fix.compute(60.0);
  for (const auto& agent : sim.agents()) {
    assert(agent.x[0] >= sim.domain().lo()[0]);
    assert(agent.x[0] < sim.domain().hi()[0]);
    assert(agent.x[1] >= sim.domain().lo()[1]);
    assert(agent.x[1] < sim.domain().hi()[1]);
    assert(agent.x[2] >= sim.domain().lo()[2]);
  }
  std::cout << "  test_floor_is_contained: PASSED\n";
}

void test_relaxation_is_dissipative() {
  Real r = CELL_RADIUS_DEFAULT;
  const Real separation = 2 * r - 0.1e-6;
  Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b = {50e-6 + separation, 50e-6, 50e-6};
  MechanicsConfig mcfg;
  mcfg.hertz_k = 1.0e-6;
  auto sim = make_two_agent_sim(pos_a, pos_b, mcfg);
  FixMechanics fix(sim, mcfg);
  Real previous = separation;
  for (Int step = 0; step < 100; ++step) {
    sim.domain().spatial_hash().clear();
    sim.domain().spatial_hash().insert(0, sim.agents()[0].x);
    sim.domain().spatial_hash().insert(1, sim.agents()[1].x);
    fix.compute(1.0);
    const Real current = sim.agents()[1].x[0] - sim.agents()[0].x[0];
    assert(current + 1e-18 >= previous);
    assert(current <= 2 * r + 1e-12);
    previous = current;
  }
  std::cout << "  test_relaxation_is_dissipative: PASSED\n";
}

void test_dense_population_stays_contained() {
  constexpr Int kAgentPairs = 200;
  constexpr Real kBioDt = 60.0;
  const Real r = CELL_RADIUS_DEFAULT;
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.domain.hi = {100e-6, 100e-6, 100e-6};
  cfg.domain.grid_dx = 10e-6;
  cfg.domain.hash_cell_size = 20e-6;
  cfg.time.total_time = 0.0;
  cfg.hdf5.enabled = false;
  cfg.vbf.viscosity = 0.01;
  cfg.initial_population.placement = "z_slab";
  cfg.initial_population.z_min = 40e-6;
  cfg.initial_population.z_max = 60e-6;
  cfg.initial_strains.push_back(
      {0, 2 * kAgentPairs, 5.0e-4, {}, false, 0, 0});

  Simulation sim;
  sim.init(cfg);
  assert(sim.agents().size() == 2 * kAgentPairs);
  for (Int i = 0; i < sim.agents().size(); i += 2) {
    const Real x = 10e-6 + (i / 2 % 20) * 4e-6;
    const Real y = 10e-6 + (i / 40) * 4e-6;
    const Real z = 45e-6;
    sim.agents()[i].x = {x, y, z};
    sim.agents()[i + 1].x = {x + 0.1e-6, y, z};
  }

  MechanicsConfig mcfg = cfg.fixes.mechanics;
  FixMechanics fix(sim, mcfg);
  for (Int step = 0; step < 20; ++step) {
    sim.domain().spatial_hash().clear();
    for (Int i = 0; i < sim.agents().size(); ++i) {
      sim.domain().spatial_hash().insert(i, sim.agents()[i].x);
    }
    std::vector<Vec3> before;
    before.reserve(sim.agents().size());
    for (const auto& agent : sim.agents()) before.push_back(agent.x);
    fix.compute(kBioDt);
    for (Int i = 0; i < sim.agents().size(); ++i) {
      const auto& agent = sim.agents()[i];
      const Real dx = agent.x[0] - before[i][0];
      const Real dy = agent.x[1] - before[i][1];
      const Real dz = agent.x[2] - before[i][2];
      const Real displacement = std::sqrt(dx * dx + dy * dy + dz * dz);
      assert(displacement <=
             kMechanicsMaxDisplacementRadiusFraction * r * (1.0 + 1e-12));
      assert(agent.x[0] >= sim.domain().lo()[0]);
      assert(agent.x[0] < sim.domain().hi()[0]);
      assert(agent.x[1] >= sim.domain().lo()[1]);
      assert(agent.x[1] < sim.domain().hi()[1]);
      assert(agent.x[2] >= sim.domain().lo()[2]);
    }
  }
  assert(sim.mechanics_step_stats().displacement_clamps > 0);
  std::cout << "  test_dense_population_stays_contained: PASSED\n";
}

int main() {
  std::cout << "=== Mechanics Tests (Issue #16) ===\n";
  test_overlapping_agents_pushed_apart();
  test_hertzian_force_scaling();
  test_no_force_without_overlap();
  test_adhesion_holds_agents();
  test_no_adhesion_beyond_range();
  test_dead_agents_ignored();
  test_displacement_is_bounded_and_counted();
  test_calm_pair_has_no_clamp();
  test_viscosity_slows_relaxation();
  test_floor_is_contained();
  test_relaxation_is_dissipative();
  test_dense_population_stays_contained();
  std::cout << "All mechanics tests passed.\n";
  return 0;
}
