/* -----------------------------------------------------------------------
   GutIBM – Bacteriocin fix implementation
   ----------------------------------------------------------------------- */

#include "fix_bacteriocin.h"
#include "simulation.h"
#include "qssa_solver.h"
#include <algorithm>
#include <cmath>
#include <numbers>

namespace gutibm {

namespace {

constexpr Real k_ln2 = std::numbers::ln2;
constexpr Real k_sos_lysis_delay_s = 300.0;
constexpr Real k_phage_lysis_delay_s = 60.0;

}  // namespace

FixBacteriocin::FixBacteriocin(Simulation& sim, const BacteriocinConfig& cfg)
    : Fix("bacteriocin", sim), cfg_(cfg) {}

void FixBacteriocin::init() { /* no-op: configuration is read at construction */ }

bool FixBacteriocin::has_release_mode(const Agent& agent, ReleaseMode mode) const {
  return std::ranges::any_of(agent.genome.bi_loci,
                             [mode](const BICluster& bi) { return bi.release_mode == mode; });
}

void FixBacteriocin::compute(Real dt) {
  auto& agents = sim_.agents();

  for (Agent& a : agents) {
    if (a.state == PhenoState::DEAD) continue;
    if (a.genome.bi_loci.empty()) continue;

    if (has_release_mode(a, ReleaseMode::SOS_LYSIS)) {
      check_sos_induction(a, dt);
    }

    bool has_continuous = false;
    for (const auto& bi : a.genome.bi_loci) {
      if (bi.release_mode == ReleaseMode::PHAGE_LYSIS) {
        check_phage_induction(a, bi, dt);
      } else if (bi.release_mode == ReleaseMode::CONTINUOUS) {
        has_continuous = true;
      }
    }
    if (has_continuous) {
      apply_microcin_secretion(a, dt);
    }
  }
}

void FixBacteriocin::post_step(Real dt) { // NOLINT(readability-make-member-function-const)
  auto& agents = sim_.agents();

  for (Agent& a : agents) {
    if (a.state == PhenoState::SOS_INDUCED) {
      a.timers.sos_timer -= dt;
      if (a.timers.sos_timer <= 0.0) {
        lyse_agent(a);
      }
    }
  }
}

void FixBacteriocin::apply_microcin_secretion(Agent& agent, Real /*dt*/) const {
  if (agent.state != PhenoState::NORMAL) return;
  if (agent.flags.microcin_penalty_applied) return;

  agent.mu_max *= (1.0 - cfg_.microcin_mu_penalty);
  agent.flags.microcin_penalty_applied = true;
}

void FixBacteriocin::check_sos_induction(Agent& agent, Real dt) {
  if (agent.state == PhenoState::SOS_INDUCED) return;

  const Real bio_dt = sim_.config().time.bio_dt;
  Real rate_total = cfg_.sos_basal_rate;

  if (agent.flags.just_divided && bio_dt > 0.0) {
    rate_total += cfg_.sos_lysis_prob / bio_dt;
  }

  const Real nuclease_conc = sim_.local_nuclease_toxin(agent);
  rate_total += cfg_.sos_cross_induction_rate * nuclease_conc;

  rate_total += sim_.ros_induction_rate(agent);

  const Real p_sos = 1.0 - std::exp(-rate_total * dt);

  if (sim_.rng().bernoulli(p_sos)) {
    agent.state     = PhenoState::SOS_INDUCED;
    agent.timers.sos_timer = k_sos_lysis_delay_s;
  }
}

void FixBacteriocin::check_phage_induction(Agent& agent, const BICluster& bi, Real dt) {
  if (agent.state == PhenoState::SOS_INDUCED) return;
  if (bi.phage_induction_rate <= 0.0) return;

  const Real gen_time = (agent.mu_realized > 0.0)
      ? k_ln2 / agent.mu_realized : 1.0e6;
  const Real rate_per_s = bi.phage_induction_rate / gen_time;
  const Real p_induction = 1.0 - std::exp(-rate_per_s * dt);

  if (sim_.rng().bernoulli(p_induction)) {
    agent.state     = PhenoState::SOS_INDUCED;
    agent.timers.sos_timer = k_phage_lysis_delay_s;
  }
}

void FixBacteriocin::lyse_agent(Agent& agent) {
  const Real base_release = sim_.config().qssa.colicin_release_rate;
  const Real creation_time = sim_.time();

  for (const auto& bi : agent.genome.bi_loci) {
    if (bi.release_mode == ReleaseMode::CONTINUOUS) continue;
    if (bi.molecular_weight < 10000.0) continue;

    const Real scale = (bi.burst_size > 0.0)
        ? bi.burst_size / cfg_.burst_molecules : 1.0;

    ToxinBurstSource burst;
    burst.pos = agent.x;
    burst.params.diff_coeff = bi.diff_coeff;
    burst.params.retardation = bi.retardation;
    burst.params.pI = bi.pI;
    burst.params.source_rate = base_release * scale;
    burst.creation_time = creation_time;
    burst.is_nuclease = bi.is_nuclease;
    burst.decay_rate = (bi.protease_half_life > 0.0)
        ? k_ln2 / bi.protease_half_life : 0.0;
    if (!sim_.config().chem_env.protease.enabled) {
      burst.decay_rate = 0.0;
    }
    sim_.add_toxin_burst(burst);
  }

  agent.state = PhenoState::DEAD;

  sim_.lineage_tracker().record_lysis(agent.identity.tag, agent.x,
                                       agent.genome.lineage_id);
}

Real FixBacteriocin::retardation_for_pI(Real pI) const {
  if (pI > 8.5) return cfg_.retardation_basic;
  if (pI < 6.0) return cfg_.retardation_acidic;
  return cfg_.retardation_neutral;
}

}  // namespace gutibm
