/* -----------------------------------------------------------------------
   GutIBM – Contact-dependent Horizontal Gene Transfer (HGT)
   
   Conjugation requires physical cell-to-cell contact via F-pili.
   Mating-Pair Stabilization (MPS) probability scales inversely
   with local advective shear:
   
     P(MPS) = P_base * exp(-shear / shear_crit)
   
   Realistic F-pili length: 1–4 um
   HGT is restricted to stable microcolony environments where
   shear is low, preventing "rescue" of immigrants in high-flow.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_CONJUGATION_H
#define GUTIBM_FIX_CONJUGATION_H

#include "fix.h"
#include "agent.h"

namespace gutibm {

struct ConjugationConfig {
  Real pili_length        = 4.0e-6;   // max F-pilus reach (m), used when heterogeneity off
  Real base_transfer_rate = 1.0e-4;   // conjugation events per s per pair
  Real shear_critical     = 10.0;     // critical shear rate (1/s)
  Real plasmid_copy_cost  = 0.02;     // 2% mu penalty per transferred plasmid
  bool pili_heterogeneity = false;    // enable per-event F-pilus length sampling
  Real pili_length_min    = 1.0e-6;   // min F-pilus length (m), uniform lower bound
  Real pili_length_max    = 4.0e-6;   // max F-pilus length (m), uniform upper bound
};

class FixConjugation : public Fix {
 public:
  FixConjugation(Simulation& sim, const ConjugationConfig& cfg);

  void compute(Real dt) override;

 private:
  void attempt_transfer(Agent& donor, Agent& recipient, Real shear, Real dt);

  ConjugationConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_CONJUGATION_H
