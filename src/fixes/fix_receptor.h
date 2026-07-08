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
#include "chemical_field.h"

namespace gutibm {

struct ReceptorConfig {
  // Binding affinities (Kd, mol/m^3)
  // Spec 6 / Receptor Ligand Parameterization: kd_b12_btuB is the BtuB affinity
  // for the dominant corrinoid analog (aka kd_corrinoid_btuB), the key unknown
  // governing colicin-E competition against the ~1 uM corrinoid pool.
  Real kd_b12_btuB        = 1.0e-9;   // corrinoid affinity for BtuB (1 nM)
  Real kd_colicinE_btuB   = 5.0e-10;  // Colicin E affinity for BtuB
  // Spec 6 §4.1 — FepA affinity for ferric enterobactin tightened 10 nM -> 1 nM
  // to match measured high-affinity siderophore uptake.
  Real kd_enterobactin    = 1.0e-9;   // enterobactin affinity for FepA (1 nM)
  Real kd_colicinB_fepA   = 2.0e-9;   // Colicin B affinity for FepA
  Real kd_lin_enterobactin = 5.0e-8;  // linearized enterobactin for CirA
  Real kd_colicinIa_cirA  = 3.0e-9;   // Colicin Ia affinity for CirA
  Real kd_colicinM_fhuA   = 2.5e-9;   // Colicin M affinity for FhuA
  Real kd_ferrichrome     = 1.0e-8;   // ferrichrome affinity for FhuA

  // CirA ligand: fraction of siderophore pool that is linearized enterobactin
  Real cirA_linearized_fraction = 0.3;

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

  Real local_toxin_conc(const ChemicalField& chem, Int cell,
                        const char* species_name) const;

  // Competitive binding fraction: toxin occupancy given ligand competition
  // toxin_aff / ligand_aff scale the respective Kd values for partial resistance
  Real toxin_occupancy(Real tox_conc, Real ligand_conc,
                        Real kd_tox, Real kd_ligand,
                        Real receptor_expr,
                        Real toxin_aff = 1.0,
                        Real ligand_aff = 1.0) const;

  ReceptorConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_RECEPTOR_H
