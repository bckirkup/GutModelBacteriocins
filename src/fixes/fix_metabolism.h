/* -----------------------------------------------------------------------
   GutIBM – Monod kinetics, growth, division, and death
   
   Growth rate:
     mu = mu_max * [S_carbon/(Km_carbon + S_carbon)]
                  * [S_iron/(Km_iron + S_iron)]
                  * [S_b12/(Km_b12 + S_b12)]
                  - maintenance
   
   Receptor-dependent Km modification:
     Km_iron = Km_base / receptor_expr[FepA]
     Km_b12  = Km_base / receptor_expr[BtuB]
   
   When receptor expression drops, Km increases → growth penalty.
   
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

  // MetE penalty: 5% proteome cost when BtuB downregulated
  Real metE_penalty       = 0.05;
  // Acetate inhibition of MetE: half-saturation constant (mol/m³)
  Real metE_acetate_km    = 40.0;
  // At saturating acetate, penalty = metE_penalty * metE_acetate_max_factor
  Real metE_acetate_max_factor = 2.5;
  // Ethanolamine utilization loss when BtuB downregulated
  Real eut_penalty        = 0.03;
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
  void try_divide(Agent& agent);
  void check_death(Agent& agent);

  MetabolismConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_METABOLISM_H
