/* -----------------------------------------------------------------------
   GutIBM – Active cell motility fix implementation
   ----------------------------------------------------------------------- */

#include "fix_motility.h"
#include "simulation.h"
#include <algorithm>
#include <cmath>

namespace gutibm {

namespace {

constexpr Real kMinRunTimer = 0.01;
constexpr Real kMaxRunTimer = 120.0;

Vec3 random_unit_direction(RNG& rng) {
  const Real theta = rng.uniform(0.0, 2.0 * PI);
  const Real phi = rng.uniform(0.0, PI);
  return {std::sin(phi) * std::cos(theta),
          std::sin(phi) * std::sin(theta),
          std::cos(phi)};
}

}  // namespace

FixMotility::FixMotility(Simulation& sim, const MotilityConfig& cfg)
    : Fix("motility", sim), cfg_(cfg) {}

void FixMotility::init_agent_motility(Agent& agent, const MotilityConfig& cfg,
                                      RNG& rng) {
  agent.motility.swim_direction = random_unit_direction(rng);
  agent.motility.swim_speed = cfg.swim_speed;
  agent.motility.run_timer = rng.exponential(1.0 / cfg.run_mean_duration);
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
  auto& rng = sim_.rng();
  auto& mot = agent.motility;

  if (mot.is_stopped) {
    mot.stop_timer -= dt;
    if (mot.stop_timer <= 0.0) {
      mot.is_stopped = false;
      mot.run_timer = rng.exponential(1.0 / cfg_.run_mean_duration);
    }
    return;
  }

  mot.run_timer -= dt;
  if (mot.run_timer > 0.0) {
    if (cfg_.chemotaxis_enabled && agent.grid_cell >= 0 && dt > 0.0) {
      auto& chem = sim_.chemical_field();
      Int i_carbon = chem.find("carbon");
      Int i_oxygen = chem.find("oxygen");
      Real carbon = (i_carbon >= 0) ? chem.conc(i_carbon, agent.grid_cell) : 0.0;
      Real oxygen = (i_oxygen >= 0) ? chem.conc(i_oxygen, agent.grid_cell) : 0.0;
      const Real d_carbon = (carbon - mot.prev_carbon) / dt;
      const Real d_oxygen = (oxygen - mot.prev_oxygen) / dt;
      Real modifier = 1.0 + cfg_.chi_carbon * d_carbon + cfg_.chi_oxygen * d_oxygen;
      modifier = std::clamp(modifier, 0.1, 10.0);
      mot.run_timer = std::clamp(
          mot.run_timer * modifier, kMinRunTimer, kMaxRunTimer);
      mot.prev_carbon = carbon;
      mot.prev_oxygen = oxygen;
    }
    return;
  }

  Real reorient_prob = 1.0;
  const auto& hash = sim_.domain().spatial_hash();
  if (auto neighbor_count = static_cast<Int>(hash.query_radius(
          agent.x, cfg_.cluster_suppress_radius).size());
      neighbor_count >= cfg_.cluster_suppress_threshold) {
    reorient_prob *= cfg_.cluster_tumble_factor;
  }

  if (!rng.bernoulli(reorient_prob)) {
    mot.run_timer = rng.exponential(1.0 / cfg_.run_mean_duration);
    return;
  }

  if (cfg_.chemotaxis_enabled && agent.grid_cell >= 0) {
    auto& chem = sim_.chemical_field();
    Int i_carbon = chem.find("carbon");
    Int i_oxygen = chem.find("oxygen");
    Real carbon = (i_carbon >= 0) ? chem.conc(i_carbon, agent.grid_cell) : 0.0;
    Real oxygen = (i_oxygen >= 0) ? chem.conc(i_oxygen, agent.grid_cell) : 0.0;
    mot.prev_carbon = carbon;
    mot.prev_oxygen = oxygen;
  }

  if (rng.bernoulli(cfg_.stop_probability)) {
    mot.is_stopped = true;
    mot.stop_timer = rng.exponential(1.0 / cfg_.stop_duration);
    return;
  }

  mot.swim_direction[0] = -mot.swim_direction[0];
  mot.swim_direction[1] = -mot.swim_direction[1];
  mot.swim_direction[2] = -mot.swim_direction[2];
  mot.run_timer = rng.exponential(1.0 / cfg_.run_mean_duration);
}

}  // namespace gutibm
