/* -----------------------------------------------------------------------
   GutIBM – AI-2 quorum sensing fix (Spec 11)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_QUORUM_SENSING_H
#define GUTIBM_FIX_QUORUM_SENSING_H

#include "fix.h"
#include "quorum_sensing_config.h"

namespace gutibm {

class FixQuorumSensing : public Fix {
 public:
  FixQuorumSensing(Simulation& sim, const QuorumSensingConfig& cfg);

  void compute(Real dt) override;

 private:
  QuorumSensingConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_QUORUM_SENSING_H
