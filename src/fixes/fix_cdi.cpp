/* -----------------------------------------------------------------------
   GutIBM – Contact-dependent inhibition (CDI) fix implementation
   ----------------------------------------------------------------------- */

#include "fix_cdi.h"
#include "simulation.h"
#include <cmath>

namespace gutibm {

FixCdi::FixCdi(Simulation& sim, const CdiConfig& cfg)
    : Fix("cdi", sim), cfg_(cfg) {}

namespace {

bool is_active_corpse(const Agent& a, Real sim_time, Real corpse_persistence) {
  return a.state == PhenoState::DEAD && a.death_time >= 0.0
      && (sim_time - a.death_time) < corpse_persistence;
}

bool corpse_blocks_cdi(const Agent& attacker, const Agent& victim,
                       const AgentPool& agents, const Domain& domain,
                       Real sim_time, const CdiConfig& cfg) {
  const Real dist_av2 = domain.min_image_dist_sq(attacker.x, victim.x);
  for (const Agent& corpse : agents) {
    if (!is_active_corpse(corpse, sim_time, cfg.corpse_persistence)) continue;
    const Real dist_ac2 = domain.min_image_dist_sq(attacker.x, corpse.x);
    if (dist_ac2 >= dist_av2) continue;
    if (dist_ac2 > cfg.contact_radius * cfg.contact_radius) continue;
    return true;
  }
  return false;
}

bool victim_eligible_for_cdi(const Agent& attacker, const Agent& victim,
                             const AgentPool& agents, const Domain& domain,
                             Real sim_time, const CdiConfig& cfg) {
  if (victim.state == PhenoState::DEAD) return false;
  if (victim.genome.cdi_immunity == attacker.genome.cdi_type) return false;
  if (corpse_blocks_cdi(attacker, victim, agents, domain, sim_time, cfg)) {
    return false;
  }
  if (Real d2 = domain.min_image_dist_sq(attacker.x, victim.x);
      d2 > cfg.contact_radius * cfg.contact_radius) {
    return false;
  }
  return true;
}

void try_cdi_kill(Agent& victim, Real kill_prob, Real sim_time, RNG& rng) {
  if (!rng.bernoulli(kill_prob)) return;
  victim.state = PhenoState::DEAD;
  victim.death_time = sim_time;
}

void process_cdi_neighbors(const Agent& attacker, Int attacker_idx,
                           const std::vector<Int>& neighbors,
                           AgentPool& agents, Simulation& sim,
                           Real kill_prob, Real sim_time, const CdiConfig& cfg) {
  for (Int j : neighbors) {
    if (j == attacker_idx) continue;
    Agent& victim = agents[j];
    if (!victim_eligible_for_cdi(attacker, victim, agents, sim.domain(),
                                 sim_time, cfg)) {
      continue;
    }
    try_cdi_kill(victim, kill_prob, sim_time, sim.rng());
  }
}

}  // namespace

void FixCdi::compute(Real dt) {
  if (!cfg_.enabled) return;

  auto& agents = sim_.agents();
  const auto& hash = sim_.domain().spatial_hash();
  const Real kill_prob = 1.0 - std::exp(-cfg_.kill_rate * dt);
  const Real sim_time = sim_.time();

  for (Int i = 0; i < agents.size(); ++i) {
    const Agent& attacker = agents[i];
    if (attacker.state == PhenoState::DEAD) continue;
    if (attacker.genome.cdi_type == 0) continue;

    auto neighbors = hash.query_radius(attacker.x, cfg_.contact_radius);
    process_cdi_neighbors(attacker, i, neighbors, agents, sim_,
                          kill_prob, sim_time, cfg_);
  }
}

}  // namespace gutibm
