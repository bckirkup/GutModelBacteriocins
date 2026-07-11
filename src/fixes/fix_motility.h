/* -----------------------------------------------------------------------
   GutIBM – Active cell motility fix (Spec 3)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_MOTILITY_H
#define GUTIBM_FIX_MOTILITY_H

#include "fix.h"
#include "agent.h"
#include "random.h"
#include "motility_config.h"

namespace gutibm {

class FixMotility : public Fix {
 public:
  FixMotility(Simulation& sim, const MotilityConfig& cfg);

  void init() override;
  void pre_step(Real dt) override;
  void compute(Real /*dt*/) override { /* motility runs in pre_step only */ }

  static void init_agent_motility(Agent& agent, const MotilityConfig& cfg, RNG& rng);

 private:
  void update_agent(Agent& agent, Real dt);
  Real advance_stopped_interval(Agent::MotilityState& motility, Real remaining);
  Real advance_running_interval(Agent& agent, Real remaining);
  void update_chemotaxis(Agent& agent, Real dt);
  void complete_run(Agent& agent);
  void start_run(Agent::MotilityState& motility);
  Real sample_duration(Real mean_duration);

  MotilityConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_MOTILITY_H
