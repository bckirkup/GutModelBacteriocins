/* -----------------------------------------------------------------------
   GutIBM – Bacteriocin fix implementation
   ----------------------------------------------------------------------- */

#include "fix_bacteriocin.h"
#include "simulation.h"
#include <cmath>

namespace gutibm {

FixBacteriocin::FixBacteriocin(Simulation& sim, const BacteriocinConfig& cfg)
    : Fix("bacteriocin", sim), cfg_(cfg) {}

void FixBacteriocin::init() {}

void FixBacteriocin::compute(Real dt) {
  auto& agents = sim_.agents();

  for (Int i = 0; i < agents.size(); ++i) {
    Agent& a = agents[i];
    if (a.state == PhenoState::DEAD) continue;

    // Apply microcin secretion penalty to producers
    if (!a.genome.bi_loci.empty() && a.state == PhenoState::NORMAL) {
      // Continuous secretion cost already applied in metabolism via plasmid cost
      // Here we just track the secretion state
    }

    // Check for stochastic SOS induction
    check_sos_induction(a, dt);
  }
}

void FixBacteriocin::post_step(Real dt) {
  // Process SOS-induced cells: they lyse and release toxin burst
  auto& agents = sim_.agents();

  for (Int i = 0; i < agents.size(); ++i) {
    Agent& a = agents[i];
    if (a.state == PhenoState::SOS_INDUCED) {
      a.sos_timer -= dt;
      if (a.sos_timer <= 0.0) {
        lyse_agent(a);
      }
    }
  }
}

void FixBacteriocin::check_sos_induction(Agent& agent, Real dt) {
  if (agent.genome.bi_loci.empty()) return;  // no toxin genes
  if (agent.state == PhenoState::SOS_INDUCED) return;

  auto& rng = sim_.rng();

  // Basal SOS induction rate (spontaneous DNA damage)
  Real p_sos = 1.0 - std::exp(-cfg_.sos_basal_rate * dt);

  if (rng.bernoulli(p_sos)) {
    agent.state     = PhenoState::SOS_INDUCED;
    agent.sos_timer = 300.0;  // 5 min delay before lysis
  }
}

void FixBacteriocin::lyse_agent(Agent& agent) {
  // Cell dies and releases toxin burst
  agent.state = PhenoState::DEAD;

  // The toxin release is handled implicitly:
  // QSSA solver picks up agents with SOS_INDUCED state as burst sources
  // The burst is a large point source at the agent's position

  // Record lysis event for lineage tracking
  sim_.lineage_tracker().record_lysis(agent.tag, agent.x,
                                       agent.genome.lineage_id);
}

Real FixBacteriocin::retardation_for_pI(Real pI) const {
  if (pI > 8.5) return cfg_.retardation_basic;
  if (pI < 6.0) return cfg_.retardation_acidic;
  return cfg_.retardation_neutral;
}

}  // namespace gutibm
