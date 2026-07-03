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
  void compute(Real /*dt*/) override {}

  static void init_agent_motility(Agent& agent, const MotilityConfig& cfg, RNG& rng);

 private:
  void update_agent(Agent& agent, Real dt);

  MotilityConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_MOTILITY_H
