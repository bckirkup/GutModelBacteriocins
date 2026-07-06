/* -----------------------------------------------------------------------
   GutIBM – Monod kinetics, growth, division, and death
   
   Growth rate:
     mu = mu_max * [S_carbon/(Km_carbon + S_carbon)]
                  * monod_iron
                  * [S_b12/(Km_b12 + S_b12)]
                  - maintenance
   
   Iron uptake uses graded fallback across multiple receptor systems:
     FepA (primary, Km ~10 nM), IroN (salmochelin, 50 nM),
     IutA (aerobactin, 100 nM), Fiu (catecholate, 200 nM).
   When FepA is downregulated, cells switch to secondary receptors
   rather than complete iron starvation.
   
   B12-dependent Km modification:
     Km_b12 = Km_base / receptor_expr[BtuB]
   
   Division: when biomass exceeds 2x initial, cell divides.
   Death: natural decay or if mu_realized < 0 for extended period.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_METABOLISM_H
#define GUTIBM_FIX_METABOLISM_H

#include "fix.h"
#include "agent.h"

namespace gutibm {

struct MetabolismConfig {
  Real division_threshold = 2.0;      // divide at 2x initial biomass
  Real death_threshold    = -0.01;    // net growth below this → death
  Real maintenance_rate   = 1.0e-5;   // maintenance (1/s)
  Real yield_carbon       = 0.5;      // carbon yield coefficient
  Real yield_iron         = 1.0e-6;   // iron yield (mol Fe / kg biomass)
  Real yield_b12          = 1.0e-9;   // B12 yield

  // Iron uptake Km values per receptor system
  Real km_iron_primary    = 10.0e-6;  // FepA: 10 nM in mol/m^3
  Real km_iron_iroN       = 50.0e-6;  // IroN: 50 nM (salmochelin)
  Real km_iron_iutA       = 100.0e-6; // IutA: 100 nM (aerobactin)
  Real km_iron_fiu        = 200.0e-6; // Fiu: 200 nM (catecholate)

  // MetE penalty: 5% proteome cost when BtuB downregulated
  Real metE_penalty       = 0.05;
  // Acetate inhibition of MetE: half-saturation constant (mol/m³)
  Real metE_acetate_km    = 40.0;
  // At saturating acetate, penalty = metE_penalty * metE_acetate_max_factor
  Real metE_acetate_max_factor = 2.5;
  // Ethanolamine utilization: concentration-dependent penalty (Monod)
  Real eut_km             = 0.1e-3;   // half-saturation for eut utilization (mol/m³)
  Real eut_max_penalty    = 0.10;     // max penalty when ethanolamine abundant
};

class FixMetabolism : public Fix {
 public:
  FixMetabolism(Simulation& sim, const MetabolismConfig& cfg);

  void init() override;
  void compute(Real dt) override;
  void post_step(Real dt) override;

 private:
  void compute_growth_rate(Agent& agent);
  void grow_agent(Agent& agent, Real dt);
  void perform_divisions();
  void check_death(Agent& agent);

  MetabolismConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_METABOLISM_H
