/* -----------------------------------------------------------------------
   GutIBM – Active cell motility fix implementation
   ----------------------------------------------------------------------- */

#include "fix_motility.h"
#include "species_names.h"
#include "simulation.h"
#include <algorithm>
#include <cmath>

namespace gutibm {

namespace {

constexpr Real kMinRunTimer = 0.01;
constexpr Real kMaxRunTimer = 120.0;
constexpr Real kMinEventDuration = 1.0e-9;

Vec3 random_unit_direction(RNG& rng) {
  const Real theta = rng.uniform(0.0, 2.0 * PI);
  const Real phi = rng.uniform(0.0, PI);
  return {std::sin(phi) * std::cos(theta),
          std::sin(phi) * std::sin(theta),
          std::cos(phi)};
}

Real sample_exponential_duration(RNG& rng, Real mean_duration) {
  if (mean_duration <= 0.0) return kMinEventDuration;
  return std::max(
      rng.exponential(1.0 / mean_duration), kMinEventDuration);
}

}  // namespace

FixMotility::FixMotility(Simulation& sim, const MotilityConfig& cfg)
    : Fix("motility", sim), cfg_(cfg) {}

void FixMotility::init_agent_motility(Agent& agent, const MotilityConfig& cfg,
                                      RNG& rng) {
  agent.motility.swim_direction = random_unit_direction(rng);
  agent.motility.step_displacement = {0.0, 0.0, 0.0};
  agent.motility.swim_speed = cfg.swim_speed;
  agent.motility.run_timer = sample_exponential_duration(
      rng, cfg.run_mean_duration);
  agent.motility.is_stopped = false;
  agent.motility.stop_timer = 0.0;
  agent.motility.prev_carbon = 0.0;
  agent.motility.prev_oxygen = 0.0;
}

void FixMotility::init() {
  if (!cfg_.enabled) return;
  auto& rng = sim_.rng();
  for (Agent& a : sim_.agents()) {
    if (a.state == PhenoState::DEAD) continue;
    init_agent_motility(a, cfg_, rng);
  }
}

void FixMotility::pre_step(Real dt) {
  if (!cfg_.enabled) return;
  for (Agent& a : sim_.agents()) {
    if (a.state == PhenoState::DEAD) continue;
    update_agent(a, dt);
  }
}

void FixMotility::update_agent(Agent& agent, Real dt) {
  auto& mot = agent.motility;
  mot.step_displacement = {0.0, 0.0, 0.0};
  if (dt <= 0.0) return;

  update_chemotaxis(agent, dt);

  Real remaining = dt;
  while (remaining > 0.0) {
    remaining = mot.is_stopped
        ? advance_stopped_interval(mot, remaining)
        : advance_running_interval(agent, remaining);
  }
}

Real FixMotility::advance_stopped_interval(
    Agent::MotilityState& motility, Real remaining) {
  const Real stopped_time = std::min(
      std::max(motility.stop_timer, 0.0), remaining);
  motility.stop_timer -= stopped_time;
  remaining -= stopped_time;
  if (motility.stop_timer > 0.0) return 0.0;
  start_run(motility);
  return remaining;
}

Real FixMotility::advance_running_interval(Agent& agent, Real remaining) {
  auto& mot = agent.motility;
  const Real run_time = std::min(std::max(mot.run_timer, 0.0), remaining);
  for (int d = 0; d < 3; ++d) {
    mot.step_displacement[d] += mot.swim_direction[d] * mot.swim_speed * run_time;
  }
  mot.run_timer -= run_time;
  remaining -= run_time;
  if (mot.run_timer > 0.0) return 0.0;
  complete_run(agent);
  return remaining;
}

void FixMotility::update_chemotaxis(Agent& agent, Real dt) {
  auto& mot = agent.motility;
  if (!cfg_.chemotaxis_enabled || mot.is_stopped || agent.grid_cell < 0) return;

  auto& chem = sim_.chemical_field();
  const Int i_carbon = chem.find(species::CARBON);
  const Int i_oxygen = chem.find(species::OXYGEN);
  const Real carbon = (i_carbon >= 0) ? chem.conc(i_carbon, agent.grid_cell) : 0.0;
  const Real oxygen = (i_oxygen >= 0) ? chem.conc(i_oxygen, agent.grid_cell) : 0.0;
  const Real d_carbon = (carbon - mot.prev_carbon) / dt;
  const Real d_oxygen = (oxygen - mot.prev_oxygen) / dt;
  const Real modifier = std::clamp(
      1.0 + cfg_.chi_carbon * d_carbon + cfg_.chi_oxygen * d_oxygen,
      0.1, 10.0);
  mot.run_timer = std::clamp(
      mot.run_timer * modifier, kMinRunTimer, kMaxRunTimer);
  mot.prev_carbon = carbon;
  mot.prev_oxygen = oxygen;
}

void FixMotility::complete_run(Agent& agent) {
  auto& rng = sim_.rng();
  auto& mot = agent.motility;
  Real reorient_prob = 1.0;
  const auto& hash = sim_.domain().spatial_hash();
  if (auto neighbor_count = static_cast<Int>(hash.query_radius(
          agent.x, cfg_.cluster_suppress_radius).size());
      neighbor_count >= cfg_.cluster_suppress_threshold) {
    reorient_prob *= cfg_.cluster_tumble_factor;
  }

  if (!rng.bernoulli(reorient_prob)) {
    start_run(mot);
    return;
  }

  if (cfg_.chemotaxis_enabled && agent.grid_cell >= 0) {
    auto& chem = sim_.chemical_field();
    Int i_carbon = chem.find(species::CARBON);
    Int i_oxygen = chem.find(species::OXYGEN);
    Real carbon = (i_carbon >= 0) ? chem.conc(i_carbon, agent.grid_cell) : 0.0;
    Real oxygen = (i_oxygen >= 0) ? chem.conc(i_oxygen, agent.grid_cell) : 0.0;
    mot.prev_carbon = carbon;
    mot.prev_oxygen = oxygen;
  }

  if (rng.bernoulli(cfg_.stop_probability)) {
    mot.is_stopped = true;
    mot.stop_timer = sample_duration(cfg_.stop_duration);
    return;
  }

  mot.swim_direction[0] = -mot.swim_direction[0];
  mot.swim_direction[1] = -mot.swim_direction[1];
  mot.swim_direction[2] = -mot.swim_direction[2];
  start_run(mot);
}

void FixMotility::start_run(Agent::MotilityState& motility) {
  motility.is_stopped = false;
  motility.run_timer = sample_duration(cfg_.run_mean_duration);
}

Real FixMotility::sample_duration(Real mean_duration) {
  return sample_exponential_duration(sim_.rng(), mean_duration);
}

}  // namespace gutibm
