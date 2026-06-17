/* -----------------------------------------------------------------------
   GutIBM – Dual-function TBDT receptor competitive binding
   
   Models the "Double-Bind" at TonB-dependent transporters:
     - BtuB: B12 uptake vs. colicin E entry
     - FepA: enterobactin-Fe uptake vs. colicin B/D entry
     - CirA: linearized enterobactin vs. colicin Ia
   
   Competitive occupancy:
     P(kill) = k_tox * [Tox] * receptor_expr / (Kd_tox + [Tox])
                / (1 + [Ligand]/Kd_ligand)
   
   When nutrients are scarce, receptors are unoccupied and
   "unlocked" → maximum toxin sensitivity.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_RECEPTOR_H
#define GUTIBM_FIX_RECEPTOR_H

#include "fix.h"
#include "agent.h"

namespace gutibm {

struct ReceptorConfig {
  // Binding affinities (Kd, mol/m^3)
  Real kd_b12_btuB        = 1.0e-9;   // B12 affinity for BtuB
  Real kd_colicinE_btuB   = 5.0e-10;  // Colicin E affinity for BtuB
  Real kd_enterobactin    = 1.0e-8;   // enterobactin affinity for FepA
  Real kd_colicinB_fepA   = 2.0e-9;   // Colicin B affinity for FepA
  Real kd_lin_enterobactin = 5.0e-8;  // linearized enterobactin for CirA
  Real kd_colicinIa_cirA  = 3.0e-9;   // Colicin Ia affinity for CirA

  // Kill rate constants (1/s per occupied receptor)
  Real kill_rate_colicin  = 1.0e-3;   // single-hit kill rate
  Real kill_rate_microcin = 5.0e-4;   // slower for microcins

  // Immunity protection factor
  Real immunity_factor    = 0.001;    // 1000x protection
};

class FixReceptor : public Fix {
 public:
  FixReceptor(Simulation& sim, const ReceptorConfig& cfg);

  void compute(Real dt) override;

 private:
  // Evaluate kill probability for one agent from local toxin field
  Real compute_kill_prob(const Agent& agent, Real dt) const;

  // Competitive binding fraction: toxin occupancy given ligand competition
  Real toxin_occupancy(Real tox_conc, Real ligand_conc,
                        Real kd_tox, Real kd_ligand,
                        Real receptor_expr) const;

  ReceptorConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_RECEPTOR_H
