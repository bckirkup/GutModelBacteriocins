#include "immigration.h"
#include "input_parser.h"
#include "sim_fingerprint.h"
#include "simulation.h"
#include "species_names.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

using namespace gutibm;

constexpr Real kEncounterDistanceTolerance = 1e-6;

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
      cfg, {0.0, 0.0, 0.0}, {20e-6, 20e-6, 10e-6}, rng,
      {}, false, false,
      [](const std::vector<Vec3>&, std::vector<Real>&) {
        // Intentionally empty: this test isolates placement behavior.
      },
      [](Vec3&) {
        // Intentionally empty: this test isolates placement behavior.
      });
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
      cfg, {0.0, 0.0, 0.0}, {20e-6, 20e-6, 20e-6}, rng,
      {{0.0, 0.0, 0.0}}, true, false,
      [&calls](const std::vector<Vec3>& candidates, std::vector<Real>& out) {
        ++calls;
        for (size_t i = 0; i < candidates.size(); ++i) {
          out[i] = 25e-12;
        }
      },
      [](Vec3&) {
        // Intentionally empty: projection is not under test here.
      });
  assert(calls == 1);
  assert(positions.size() == 1);
  std::cout << "  test_at_distance_selects_candidate: PASSED\n";
}

void test_disabled_is_inert() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {40e-6, 40e-6, 20e-6};
  cfg.time.total_time = 180.0;
  cfg.time.bio_dt = 60.0;
  cfg.hdf5.enabled = false;
  cfg.enabled_fixes = {"mechanics"};
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 2;
  strain.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(strain);
  Simulation off;
  off.init(cfg);
  off.run();
  const uint64_t off_fingerprint = test_util::simulation_fingerprint(off);
  cfg.immigration.enabled = false;
  cfg.immigration.count = 4;
  cfg.immigration.strain_index = 0;
  cfg.immigration.placement = "z_slab";
  cfg.immigration.z_min = 2e-6;
  cfg.immigration.z_max = 8e-6;
  cfg.immigration.schedule = "continuous";
  cfg.immigration.rate = 1.0;
  Simulation on;
  on.init(cfg);
  on.run();
  assert(off_fingerprint == test_util::simulation_fingerprint(on));
  std::cout << "  test_disabled_is_inert: PASSED\n";
}

void test_at_distance_end_to_end_and_empty_fallback() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {80e-6, 80e-6, 30e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.total_time = 60.0;
  cfg.time.bio_dt = 60.0;
  cfg.hdf5.enabled = false;
  cfg.enabled_fixes = {"mechanics"};
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 8;
  cfg.initial_strains.push_back(resident);
  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = 0;
  cfg.initial_strains.push_back(immigrant);
  cfg.immigration.enabled = true;
  cfg.immigration.count = 1;
  cfg.immigration.strain_index = 1;
  cfg.immigration.placement = "at_distance";
  cfg.immigration.distance = 8e-6;
  cfg.immigration.distance_tolerance = 2e-6;
  cfg.immigration.step = 0;

  Simulation sim;
  sim.init(cfg);
  const Vec3 center = {40e-6, 40e-6, 15e-6};
  for (Agent& agent : sim.agents()) {
    agent.x = center;
    agent.flags.in_crypt = true;
  }
  sim.step(60.0);
  const Agent* injected = nullptr;
  for (const Agent& agent : sim.agents()) {
    if (agent.identity.type == immigrant.type) injected = &agent;
  }
  assert(injected != nullptr);
  Real nearest_sq = std::numeric_limits<Real>::max();
  for (const Agent& agent : sim.agents()) {
    if (agent.identity.type == resident.type &&
        agent.state != PhenoState::DEAD) {
      nearest_sq = std::min(nearest_sq,
                            sim.domain().min_image_dist_sq(injected->x, agent.x));
    }
  }
  const Real achieved = std::sqrt(nearest_sq);
  assert(std::abs(achieved - cfg.immigration.distance) <=
         cfg.immigration.distance_tolerance);
  std::cout << "  at_distance achieved error="
            << std::abs(achieved - cfg.immigration.distance) << "\n";

  cfg.initial_strains[0].count = 0;
  cfg.immigration.distance = 5e-6;
  std::ostringstream warning;
  auto* old_buffer = std::cerr.rdbuf(warning.rdbuf());
  Simulation empty;
  empty.init(cfg);
  empty.step(60.0);
  std::cerr.rdbuf(old_buffer);
  assert(empty.agents().size() == 1);
  assert(warning.str().find("no live biomass") != std::string::npos);
  std::cout << "  test_at_distance_empty_fallback: PASSED\n";
}

void test_small_colony_shell_reliability() {
  for (uint64_t seed = 9100; seed < 9108; ++seed) {
    SimulationConfig cfg = InputParser::default_config();
    cfg.seed = seed;
    cfg.domain.hi = {200e-6, 200e-6, 50e-6};
    cfg.domain.grid_dx = 5e-6;
    cfg.time.bio_dt = 60.0;
    cfg.hdf5.enabled = false;
    cfg.enabled_fixes = {"mechanics"};
    cfg.initial_strains.clear();
    SimulationConfig::InitialStrain resident;
    resident.type = 1;
    resident.count = 2;
    cfg.initial_strains.push_back(resident);
    SimulationConfig::InitialStrain immigrant;
    immigrant.type = 2;
    immigrant.count = 0;
    cfg.initial_strains.push_back(immigrant);
    cfg.immigration.enabled = true;
    cfg.immigration.count = 1;
    cfg.immigration.strain_index = 1;
    cfg.immigration.placement = "at_distance";
    cfg.immigration.distance = 5e-6;
    cfg.immigration.distance_tolerance = 1e-6;
    cfg.immigration.step = 0;

    Simulation sim;
    sim.init(cfg);
    for (Agent& agent : sim.agents()) {
      agent.x = {100e-6, 100e-6, 25e-6};
      agent.flags.in_crypt = true;
    }
    sim.step(0.0);
    assert(sim.step_events().immigrations == 1);
    const Agent* injected = nullptr;
    for (const Agent& agent : sim.agents()) {
      if (agent.identity.type == immigrant.type) injected = &agent;
    }
    assert(injected != nullptr);
    Real nearest_sq = std::numeric_limits<Real>::max();
    for (const Agent& agent : sim.agents()) {
      if (agent.identity.type == resident.type) {
        nearest_sq = std::min(
            nearest_sq, sim.domain().min_image_dist_sq(injected->x, agent.x));
      }
    }
    const Real error = std::abs(std::sqrt(nearest_sq) - cfg.immigration.distance);
    assert(error <= cfg.immigration.distance_tolerance);
    std::cout << "  small-colony seed=" << seed
              << " distance error=" << error << "\n";
  }
  std::cout << "  test_small_colony_shell_reliability: PASSED\n";
}

void test_large_colony_shell_regression() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {220e-6, 150e-6, 150e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.bio_dt = 60.0;
  cfg.hdf5.enabled = false;
  cfg.enabled_fixes = {"mechanics"};
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 3000;
  cfg.initial_strains.push_back(resident);
  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = 0;
  cfg.initial_strains.push_back(immigrant);
  cfg.immigration.enabled = true;
  cfg.immigration.count = 1;
  cfg.immigration.strain_index = 1;
  cfg.immigration.placement = "at_distance";
  cfg.immigration.distance = 100e-6;
  cfg.immigration.distance_tolerance = 0.25e-6;
  cfg.immigration.step = 0;

  Simulation sim;
  sim.init(cfg);
  const Vec3 center = {110e-6, 75e-6, 75e-6};
  std::vector<Vec3> colony_positions;
  for (Int x = -10; x <= 10; ++x) {
    for (Int y = -10; y <= 10; ++y) {
      for (Int z = -10; z <= 10; ++z) {
        colony_positions.emplace_back(Vec3{
            {center[0] + static_cast<Real>(x) * 1.2e-6,
             center[1] + static_cast<Real>(y) * 1.2e-6,
             center[2] + static_cast<Real>(z) * 1.2e-6}});
      }
    }
  }
  std::ranges::sort(colony_positions, [center](const Vec3& lhs, const Vec3& rhs) {
    const Real lhs_radius = (lhs[0] - center[0]) * (lhs[0] - center[0])
                          + (lhs[1] - center[1]) * (lhs[1] - center[1])
                          + (lhs[2] - center[2]) * (lhs[2] - center[2]);
    const Real rhs_radius = (rhs[0] - center[0]) * (rhs[0] - center[0])
                          + (rhs[1] - center[1]) * (rhs[1] - center[1])
                          + (rhs[2] - center[2]) * (rhs[2] - center[2]);
    return lhs_radius < rhs_radius;
  });
  for (Int i = 0; i < static_cast<Int>(sim.agents().size()); ++i) {
    sim.agents()[i].x = colony_positions[i];
    sim.agents()[i].flags.in_crypt = true;
  }
  sim.step(0.0);

  const Agent* injected = nullptr;
  for (const Agent& agent : sim.agents()) {
    if (agent.identity.type == immigrant.type) injected = &agent;
  }
  assert(injected != nullptr);
  Real nearest_sq = std::numeric_limits<Real>::max();
  for (const Agent& agent : sim.agents()) {
    if (agent.identity.type == resident.type) {
      nearest_sq = std::min(
          nearest_sq, sim.domain().min_image_dist_sq(injected->x, agent.x));
    }
  }
  const Real error = std::abs(std::sqrt(nearest_sq) -
                              cfg.immigration.distance);
  assert(error <= cfg.immigration.distance_tolerance);
  std::cout << "  test_large_colony_shell_regression: error=" << error
            << " PASSED\n";
}

void test_shell_candidate_rejected_by_global_reducer() {
  ImmigrationConfig cfg;
  cfg.enabled = true;
  cfg.count = 1;
  cfg.placement = "at_distance";
  cfg.distance = 5e-6;
  cfg.distance_tolerance = 1e-12;
  RNG rng(123);
  Int calls = 0;
  const auto positions = immigration_positions(
      cfg, {0.0, 0.0, 0.0}, {20e-6, 20e-6, 20e-6}, rng,
      {{0.0, 0.0, 0.0}}, true, false,
      [&calls](const std::vector<Vec3>&, std::vector<Real>& distances) {
        ++calls;
        // The unrelated cluster is the true nearest biomass for every
        // proposal in this reducer fixture, so the shell geometry alone
        // must not make a candidate acceptable.
        std::ranges::fill(distances, 0.0);
      },
      [](Vec3&) {
        // Intentionally empty: projection is not under test here.
      });
  assert(calls == 2);
  assert(positions.empty());
  std::cout << "  test_shell_candidate_rejected_by_global_reducer: PASSED\n";
}

void test_failed_at_distance_is_fatal() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {20e-6, 20e-6, 20e-6};
  cfg.time.bio_dt = 60.0;
  cfg.hdf5.enabled = false;
  cfg.enabled_fixes = {"mechanics"};
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 1;
  cfg.initial_strains.push_back(resident);
  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = 0;
  cfg.initial_strains.push_back(immigrant);
  cfg.immigration.enabled = true;
  cfg.immigration.count = 1;
  cfg.immigration.strain_index = 1;
  cfg.immigration.placement = "at_distance";
  cfg.immigration.distance = 1e-3;
  cfg.immigration.distance_tolerance = 1e-9;
  cfg.immigration.step = 0;

  Simulation sim;
  sim.init(cfg);
  bool threw = false;
  try {
    sim.step(0.0);
  } catch (const SimulationError&) {
    threw = true;
  }
  assert(threw);
  std::cout << "  test_failed_at_distance_is_fatal: PASSED\n";
}

void test_continuous_schedule_end_to_end() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {40e-6, 40e-6, 20e-6};
  cfg.time.total_time = 600.0;
  cfg.time.bio_dt = 60.0;
  cfg.hdf5.enabled = false;
  cfg.enabled_fixes = {"mechanics"};
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 2;
  cfg.initial_strains.push_back(strain);
  cfg.immigration.enabled = true;
  cfg.immigration.count = 1;
  cfg.immigration.schedule = "continuous";
  cfg.immigration.rate = 1.0 / 60.0;

  Simulation sim;
  sim.init(cfg);
  for (Agent& agent : sim.agents()) agent.flags.in_crypt = true;
  sim.run();
  const Int injected = sim.step_events().immigrations;
  assert(injected > 0);
  assert(injected >= 2 && injected <= 20);
  std::cout << "  test_continuous_schedule_end_to_end: injected="
            << injected << " PASSED\n";
}

struct EncounterResult {
  bool injected = false;
  Int colicin_kills = 0;
  bool alive = false;
  Real distance_error = std::numeric_limits<Real>::max();
};

EncounterResult run_bacteriocin_encounter(uint64_t seed, Real target_distance) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {200e-6, 200e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.total_time = 120.0;
  cfg.time.bio_dt = 60.0;
  cfg.hdf5.enabled = false;
  cfg.seed = seed;
  cfg.enabled_fixes = {"bacteriocin", "receptor", "mechanics"};
  cfg.qssa.toxin_cutoff = 80e-6;
  cfg.advection.radial_turnover = 1.0e30;
  cfg.advection.distal_transit_time = 1.0e30;
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain producer;
  producer.type = 1;
  producer.count = 2;
  producer.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(producer);
  SimulationConfig::InitialStrain sensitive;
  sensitive.type = 2;
  sensitive.count = 0;
  sensitive.plasmids = {};
  cfg.initial_strains.push_back(sensitive);
  cfg.immigration.enabled = true;
  cfg.immigration.count = 1;
  cfg.immigration.strain_index = 1;
  cfg.immigration.placement = "at_distance";
  cfg.immigration.distance = target_distance;
  cfg.immigration.distance_tolerance = kEncounterDistanceTolerance;
  cfg.immigration.step = 1;
  cfg.fixes.bacteriocin.sos_lysis_prob = 1.0;
  cfg.fixes.bacteriocin.burst_release_tau = 300.0;
  cfg.fixes.receptor.kill_rate_colicin = 100.0;

  Simulation sim;
  sim.init(cfg);
  const Vec3 center = {100e-6, 100e-6, 25e-6};
  for (Agent& agent : sim.agents()) {
    agent.x = center;
    agent.flags.in_crypt = true;
    for (BICluster& bi : agent.genome.bi_loci) bi.burst_size = 1.0e12;
  }
  sim.agents()[0].state = PhenoState::SOS_INDUCED;
  sim.agents()[0].timers.sos_timer = 0.0;
  sim.step(60.0);
  const Int immigrations_before = sim.step_events().immigrations;
  sim.step(0.0);
  EncounterResult result;
  result.injected =
      sim.step_events().immigrations == immigrations_before + 1;
  const Agent* sensitive_agent = nullptr;
  for (const Agent& agent : sim.agents()) {
    if (agent.identity.type == sensitive.type) {
      sensitive_agent = &agent;
      break;
    }
  }
  if (!result.injected || sensitive_agent == nullptr) return result;

  Real nearest_sq = std::numeric_limits<Real>::max();
  for (const Agent& agent : sim.agents()) {
    if (agent.identity.type == producer.type &&
        agent.state != PhenoState::DEAD) {
      nearest_sq = std::min(nearest_sq,
                            sim.domain().min_image_dist_sq(
                                sensitive_agent->x, agent.x));
    }
  }
  result.distance_error =
      std::abs(std::sqrt(nearest_sq) - target_distance);
  const Int kills_before = sim.step_events().colicin_kills;
  sim.step(60.0);
  result.colicin_kills = sim.step_events().colicin_kills - kills_before;
  result.alive = std::ranges::any_of(
      sim.agents(), [type = sensitive.type](const Agent& agent) {
        return agent.identity.type == type && agent.state != PhenoState::DEAD;
      });
  return result;
}

void test_near_colony_kill_separation() {
  Int near_kills = 0;
  Int far_kills = 0;
  constexpr Int replicates = 12;
  for (Int i = 0; i < replicates; ++i) {
    const EncounterResult near =
        run_bacteriocin_encounter(7000 + i, 5e-6);
    const EncounterResult far =
        run_bacteriocin_encounter(8000 + i, 120e-6);
    assert(near.injected);
    assert(far.injected);
    assert(near.distance_error <= kEncounterDistanceTolerance);
    assert(far.distance_error <= kEncounterDistanceTolerance);
    assert(far.alive);
    std::cout << std::scientific << std::setprecision(3)
              << "  encounter " << i << " distance errors="
              << near.distance_error << "/" << far.distance_error << "\n";
    near_kills += near.colicin_kills;
    far_kills += far.colicin_kills;
  }
  std::cout << "  near/far colicin kills=" << near_kills << "/"
            << far_kills << "\n";
  assert(near_kills >= 8);
  assert(near_kills >= far_kills + 4);
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
  test_at_distance_end_to_end_and_empty_fallback();
  test_small_colony_shell_reliability();
  test_large_colony_shell_regression();
  test_shell_candidate_rejected_by_global_reducer();
  test_failed_at_distance_is_fatal();
  test_continuous_schedule_end_to_end();
  test_near_colony_kill_separation();
  test_pulse_constructs_full_agents();
  std::cout << "All immigration tests passed.\n";
  return 0;
}
