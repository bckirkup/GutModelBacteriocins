/* -----------------------------------------------------------------------
   GutIBM – Stochastic SOS-mediated lysis and bacteriocin release
   
   Colicin release:
     - SOS-mediated suicide: 1% probability per division event
     - Upon lysis, cell releases a burst of toxin
     - Modeled as an instantaneous point source for QSSA
   
   Microcin secretion:
     - Continuous secretion with a static mu_max penalty (2–5%)
     - No cell lysis required
   
   pI-dependent diffusion:
     - Basic toxins (pI > 8.5): high retardation → Lethal Core
     - Acidic toxins (pI < 6.0): low retardation → Lethal Halo
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_BACTERIOCIN_H
#define GUTIBM_FIX_BACTERIOCIN_H

#include "fix.h"
#include "agent.h"

namespace gutibm {

struct BacteriocinConfig {
  Real sos_lysis_prob       = 0.01;    // 1% per division
  Real sos_basal_rate       = 1.0e-6;  // spontaneous SOS induction (1/s)

  // Retardation factors based on pI
  Real retardation_basic    = 50.0;    // pI > 8.5, bound to mucin
  Real retardation_acidic   = 1.5;     // pI < 6.0, repelled by mucin
  Real retardation_neutral  = 5.0;     // 6.0 <= pI <= 8.5

  // Free diffusion coefficient for ~50kDa protein (m^2/s)
  Real D_free_colicin       = 4.0e-11;

  // Burst size: molecules released per lysis event
  Real burst_molecules      = 1.0e4;

  // Microcin continuous secretion penalty on mu_max
  Real microcin_mu_penalty  = 0.03;    // 3%
};

class FixBacteriocin : public Fix {
 public:
  FixBacteriocin(Simulation& sim, const BacteriocinConfig& cfg);

  void init() override;
  void compute(Real dt) override;
  void post_step(Real dt) override;

 private:
  void check_sos_induction(Agent& agent, Real dt);
  void lyse_agent(Agent& agent);
  Real retardation_for_pI(Real pI) const;

  BacteriocinConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_BACTERIOCIN_H
