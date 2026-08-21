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

  // Free diffusion coefficient for ~50kDa protein (m^2/s)
  Real D_free_colicin       = 4.0e-11;

  // Exponential release timescale for finite lysis inventory
  Real burst_release_tau    = 300.0;

  // Microcin continuous secretion penalty on mu_max
  Real microcin_mu_penalty  = 0.03;    // 3%

  // Nuclease colicin cross-induction (provoker mechanism)
  Real sos_cross_induction_rate = 1.0e3;  // 1/s per mol/m³ nuclease toxin
};

class FixBacteriocin : public Fix {
 public:
  FixBacteriocin(Simulation& sim, const BacteriocinConfig& cfg);

  void init() override;
  void compute(Real dt) override;
  void post_step(Real dt) override;
  Real nuclease_cross_induction_rate(const Agent& agent,
                                     Int agent_index) const;

 private:
  void check_sos_induction(Agent& agent, Real dt, Int agent_index);
  void check_phage_induction(Agent& agent, const BICluster& bi, Real dt);
  void apply_microcin_secretion(Agent& agent, Real dt) const;
  void lyse_agent(Agent& agent);
  bool has_release_mode(const Agent& agent, ReleaseMode mode) const;

  BacteriocinConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_BACTERIOCIN_H
