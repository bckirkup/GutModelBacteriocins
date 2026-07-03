/* -----------------------------------------------------------------------
   GutIBM – Bacteriocin fix implementation
   ----------------------------------------------------------------------- */

#include "fix_bacteriocin.h"
#include "simulation.h"
#include "qssa_solver.h"
#include <cmath>

namespace gutibm {

FixBacteriocin::FixBacteriocin(Simulation& sim, const BacteriocinConfig& cfg)
    : Fix("bacteriocin", sim), cfg_(cfg) {}

void FixBacteriocin::init() { /* no-op: configuration is read at construction */ }

void FixBacteriocin::compute(Real dt) {
  auto& agents = sim_.agents();

  for (Agent& a : agents) {
    if (a.state == PhenoState::DEAD) continue;

    // Microcin continuous secretion: small peptides are exported without lysis.
    // Apply mu_max penalty and contribute steady-state toxin source.
    apply_microcin_secretion(a, dt);

    // Check for stochastic SOS induction (colicin-type lysis only)
    check_sos_induction(a, dt);
  }
}

void FixBacteriocin::post_step(Real dt) { // NOLINT(readability-make-member-function-const)
  // Process SOS-induced cells: they lyse and release toxin burst
  auto& agents = sim_.agents();

  for (Agent& a : agents) {
    if (a.state == PhenoState::SOS_INDUCED) {
      a.sos_timer -= dt;
      if (a.sos_timer <= 0.0) {
        lyse_agent(a);
      }
    }
  }
}

void FixBacteriocin::apply_microcin_secretion(Agent& agent, Real /*dt*/) {
  if (agent.genome.bi_loci.empty()) return;
  if (agent.state != PhenoState::NORMAL) return;

  // Identify microcin-class clusters (small peptides, MW < 10kDa)
  bool has_microcin = false;
  for (const auto& bi : agent.genome.bi_loci) {
    if (bi.molecular_weight < 10000.0) {
      has_microcin = true;
      break;
    }
  }
  if (!has_microcin) return;

  // Continuous secretion imposes a growth penalty (2–5% mu_max reduction).
  // This is applied here rather than in fix_metabolism to keep the
  // mechanism clearly associated with its biological source.
  agent.mu_max *= (1.0 - cfg_.microcin_mu_penalty);

  // The secreted microcins contribute to the steady-state toxin field
  // via the QSSA solver (they appear as low-rate continuous point sources
  // rather than burst sources from SOS lysis).
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
  const Real release_rate = sim_.config().qssa.colicin_release_rate;
  const Real creation_time = sim_.time();
  constexpr Real k_ln2 = 0.6931471805599453;

  for (const auto& bi : agent.genome.bi_loci) {
    if (bi.molecular_weight < 10000.0) continue;

    ToxinBurstSource burst;
    burst.pos = agent.x;
    burst.params.diff_coeff = bi.diff_coeff;
    burst.params.retardation = bi.retardation;
    burst.params.pI = bi.pI;
    burst.params.source_rate = release_rate;
    burst.creation_time = creation_time;
    burst.decay_rate = (bi.protease_half_life > 0.0)
        ? k_ln2 / bi.protease_half_life : 0.0;
    if (!sim_.config().protease.enabled) {
      burst.decay_rate = 0.0;
    }
    sim_.add_toxin_burst(burst);
  }

  agent.state = PhenoState::DEAD;

  sim_.lineage_tracker().record_lysis(agent.tag, agent.x,
                                       agent.genome.lineage_id);
}

Real FixBacteriocin::retardation_for_pI(Real pI) const {
  if (pI > 8.5) return cfg_.retardation_basic;
  if (pI < 6.0) return cfg_.retardation_acidic;
  return cfg_.retardation_neutral;
}

}  // namespace gutibm
