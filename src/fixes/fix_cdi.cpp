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

}  // namespace

void FixCdi::compute(Real dt) {
  if (!cfg_.enabled) return;

  auto& agents = sim_.agents();
  const auto& hash = sim_.domain().spatial_hash();
  auto& rng = sim_.rng();
  const Real kill_prob = 1.0 - std::exp(-cfg_.kill_rate * dt);
  const Real sim_time = sim_.time();

  for (Int i = 0; i < agents.size(); ++i) {
    Agent& attacker = agents[i];
    if (attacker.state == PhenoState::DEAD) continue;
    if (attacker.genome.cdi_type == 0) continue;

    auto neighbors = hash.query_radius(attacker.x, cfg_.contact_radius);
    for (Int j : neighbors) {
      if (j == i) continue;
      Agent& victim = agents[j];
      if (victim.state == PhenoState::DEAD) continue;
      if (victim.genome.cdi_immunity == attacker.genome.cdi_type) continue;
      if (corpse_blocks_cdi(attacker, victim, agents, sim_.domain(),
                            sim_time, cfg_)) {
        continue;
      }

      if (Real d2 = sim_.domain().min_image_dist_sq(attacker.x, victim.x);
          d2 > cfg_.contact_radius * cfg_.contact_radius) {
        continue;
      }

      if (rng.bernoulli(kill_prob)) {
        victim.state = PhenoState::DEAD;
        victim.death_time = sim_.time();
      }
    }
  }
}

}  // namespace gutibm
