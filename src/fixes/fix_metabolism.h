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
#include "delivery_support.h"
#include "uptake_limit.h"
#include <string>
#include <unordered_map>
#include <vector>

namespace gutibm {

struct MetabolismConfig {
  // Derived from chemistry.species_subset; these are not parser keys.
  // They disable optional uptake terms in carbon_only mode.
  bool iron_uptake_enabled = true;
  bool b12_uptake_enabled = true;
  bool eut_enabled = true;
  Real division_threshold = 2.0;      // divide at 2x initial biomass
  Real bacteriostasis_threshold = 1.0e-6; // classify viable non-growing cells
  Real maintenance_rate   = 1.0e-5;   // maintenance (1/s)
  Real carbon_maintenance_rate = 0.0; // non-growth carbon use (mol/(s kg))
  Real yield_carbon       = 0.5;      // carbon yield coefficient
  bool acid_inhibition_enabled = false;
  Real acid_inhibition_max = 0.8;
  Real acid_inhibition_Ki = 50.0;       // mol/m^3 undissociated acetate
  Real acetate_pKa = 4.76;
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

  // Agent-side uptake limitation model.
  std::string uptake_limit = "none";
  // Radius shared by the delivery concentration read and prescribed sink.
  // Zero preserves the historical agent-voxel behavior.
  Real delivery_far_field_radius = 1.0e-5;
  // Resolved from uptake_limit by InputParser::finalize_config.
  UptakeLimitMode uptake_limit_mode = UptakeLimitMode::None;
};

class FixMetabolism : public Fix {
 public:
  FixMetabolism(Simulation& sim, const MetabolismConfig& cfg);

  void init() override;
  void compute(Real dt) override;
  void post_chemistry(Real dt) override;

 private:
  void compute_agent(Agent& agent, Real dt);
  Real delivery_concentration(const Agent& agent, Int species_index) const;
  std::vector<Int> enumerate_delivery_support_cells(const Agent& agent) const;
  void ensure_delivery_support_stencil() const;
  const std::vector<Int>& delivery_support_cells(const Agent& agent) const;
  void prepare_delivery_support_cache();
  void add_delivery_mass(
      Int species_index, const Agent& agent, Real amount) const;
  Real delivery_field_funding(
      Int species_index, const Agent& agent, Real amount,
      const std::vector<Real>& requested_by_cell) const;
  void prepare_delivery_uptake(Agent& agent, Real dt);
  void prepare_delivery_oxygen(Agent& agent, Real dt);
  void commit_delivery_uptake(Real dt);
  void commit_delivery_carbon(Real dt, Int carbon);
  void commit_delivery_oxygen(Real dt, Int oxygen);
  void compute_growth_rate(Agent& agent, Real dt);
  Real uptake_limit_fraction(const Agent& agent, Real d_biomass, Real dt,
                             bool record_diagnostics);
  void charge_carbon_maintenance(Agent& agent, Real dt);
  void prepare_carbon_maintenance();
  void grow_agent(Agent& agent, Real dt);
  void apply_growth_chemistry(const Agent& agent, Real d_biomass, Real dt);
  Real realized_carbon_cost(const Agent& agent) const;
  void apply_siderophore_chemistry(Real dt);
  void apply_siderophore_chelation(Int i_sid, Int i_iron,
                                   Int i_ferric_enterobactin, Int num_cells);
  void apply_siderophore_reimport(Int i_sid, Int i_iron,
                                  Int i_ferric_enterobactin, Int num_cells,
                                  Real cell_volume, Real dt);
  void perform_divisions();

  MetabolismConfig cfg_;
  std::vector<Real> biomass_by_cell_;
  std::vector<Real> carbon_maintenance_available_;
  std::vector<Real> fepA_biomass_by_cell_;
  std::vector<Int> occupancy_by_cell_;
  std::vector<Real> chelation_by_cell_;
  std::vector<Int> touched_cells_;
  mutable std::unordered_map<TagID, std::vector<Int>>
      delivery_support_cache_;
  mutable DeliverySupportStencil delivery_support_stencil_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_METABOLISM_H
