#include "immigration.h"
#include "input_parser.h"
#include "simulation.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

void test_schedule_and_uniform_placement() {
  ImmigrationConfig cfg;
  cfg.enabled = true;
  cfg.count = 3;
  cfg.schedule = "pulse";
  cfg.step = 2;
  RNG rng(7);
  assert(immigration_event_count(cfg, 1, 60.0, rng) == 0);
  assert(immigration_event_count(cfg, 2, 60.0, rng) == 1);

  cfg.placement = "z_slab";
  cfg.z_min = 2e-6;
  cfg.z_max = 8e-6;
  auto positions = immigration_positions(
      cfg, {0.0, 0.0, 0.0}, {20e-6, 20e-6, 10e-6}, rng, false,
      [](const std::vector<Vec3>&, std::vector<Real>&) {});
  assert(positions.size() == 3);
  for (const Vec3& pos : positions) {
    assert(pos[2] >= cfg.z_min);
    assert(pos[2] < cfg.z_max);
  }
  std::cout << "  test_schedule_and_uniform_placement: PASSED\n";
}

void test_at_distance_selects_candidate() {
  ImmigrationConfig cfg;
  cfg.enabled = true;
  cfg.count = 1;
  cfg.placement = "at_distance";
  cfg.distance = 5e-6;
  cfg.distance_tolerance = 1e-12;
  RNG rng(42);
  int calls = 0;
  auto positions = immigration_positions(
      cfg, {0.0, 0.0, 0.0}, {20e-6, 20e-6, 20e-6}, rng, true,
      [&calls](const std::vector<Vec3>& candidates, std::vector<Real>& out) {
        ++calls;
        for (size_t i = 0; i < candidates.size(); ++i) {
          out[i] = 25e-12;
        }
      });
  assert(calls == 1);
  assert(positions.size() == 1);
  std::cout << "  test_at_distance_selects_candidate: PASSED\n";
}

void test_disabled_is_inert() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {40e-6, 40e-6, 20e-6};
  cfg.time.total_time = 0.0;
  cfg.hdf5.enabled = false;
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 2;
  strain.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(strain);
  Simulation off;
  off.init(cfg);
  cfg.immigration.enabled = true;
  Simulation on;
  on.init(cfg);
  assert(off.agents().size() == on.agents().size());
  std::cout << "  test_disabled_is_inert: PASSED\n";
}

void test_pulse_constructs_full_agents() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {40e-6, 40e-6, 20e-6};
  cfg.time.total_time = 60.0;
  cfg.hdf5.enabled = false;
  cfg.enabled_fixes = {"mechanics"};
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 3;
  resident.plasmids = {};
  cfg.initial_strains.push_back(resident);
  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = 0;
  immigrant.mu_max = 3.2e-4;
  immigrant.plasmids = {"ColE1"};
  immigrant.cdi_type = 7;
  immigrant.cdi_immunity = 9;
  cfg.initial_strains.push_back(immigrant);
  cfg.immigration.enabled = true;
  cfg.immigration.count = 2;
  cfg.immigration.strain_index = 1;
  cfg.immigration.step = 0;

  Simulation sim;
  sim.init(cfg);
  const Int initial_size = sim.agents().size();
  sim.step(60.0);
  assert(sim.agents().size() == initial_size + 2);
  Int formed = 0;
  for (const Agent& agent : sim.agents()) {
    if (agent.identity.type != immigrant.type) continue;
    ++formed;
    assert(agent.mu_max == immigrant.mu_max);
    assert(agent.genome.cdi_type == immigrant.cdi_type);
    assert(agent.genome.cdi_immunity == immigrant.cdi_immunity);
    assert(!agent.genome.bi_loci.empty());
    assert(agent.genome.lineage_id == agent.identity.tag);
  }
  assert(formed == 2);
  std::cout << "  test_pulse_constructs_full_agents: PASSED\n";
}

int main() {
  test_schedule_and_uniform_placement();
  test_at_distance_selects_candidate();
  test_disabled_is_inert();
  test_pulse_constructs_full_agents();
  std::cout << "All immigration tests passed.\n";
  return 0;
}
