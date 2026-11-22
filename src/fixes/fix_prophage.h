/* -----------------------------------------------------------------------
   GutIBM – Temperate prophage dynamics (Spec 14, GAP-10)

   Implements:
     - Separate genomic prophage carriage (superinfection immunity, including
       cryptic/defective remnants) vs active inducibility
       (docs/SPEC14_PRIOR_REVIEW.md §3).
     - Spontaneous + DNA-damage Hill induction hazard:
         h_ind(D) = k_spont + k_max_damage * (D^hill_n / (K_damage^hill_n + D^hill_n))
       driven by local nuclease colicin exposure (e.g., ColE2 DNase),
       ROS induction signal, and explicit exogenous damage_stimulus.
     - Explicit free virion field (species::VIRION) with implicit diffusion,
       mucus advective washout + decay, receptor-mediated adsorption,
       repressor-mediated superinfection immunity, and lysogenic conversion
       vs lytic infection with exact cell and virion conservation closure.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_PROPHAGE_H
#define GUTIBM_FIX_PROPHAGE_H

#include "agent.h"
#include "fix.h"
#include "types.h"
#include <string>

namespace gutibm {

struct ProphageConfig {
  bool enabled = false;  // default-off for bitwise backward parity

  // Separate genomic prophage carriage vs active inducible subset
  // (docs/SPEC14_PRIOR_REVIEW.md §3):
  //   P(has_prophage_immunity)  = lysogen_prevalence_init
  //   P(has_inducible_prophage) = lysogen_prevalence_init * active_inducible_fraction
  Real lysogen_prevalence_init = 0.8;
  Real active_inducible_fraction = 0.375;

  // Spontaneous + SOS/DNA-damage Hill induction hazard per active lysogen:
  //   h_ind(D) = k_spont + k_max_damage * (D^hill_n / (K_damage^hill_n + D^hill_n))
  Real k_spont = 1.2e-7;       // spontaneous induction hazard (1/s)
  Real k_max_damage = 1.0e-3;  // max DNA-damage induction hazard (1/s)
  Real K_damage = 1.0e-6;      // Hill half-saturation damage signal
  Real hill_n = 2.0;           // Hill cooperativity exponent

  // DNA-damage signal weights:
  //   D = damage_stimulus
  //       + nuclease_damage_weight * local_nuclease_exposure
  //       + ros_damage_weight * ros_induction_rate
  Real damage_stimulus = 0.0;
  Real nuclease_damage_weight = 1.0;
  Real ros_damage_weight = 1.0;

  // Latent lysis delay after induction (s); 0.0 = immediate lysis in step
  Real lysis_delay_s = 0.0;

  // Virion burst, transport, decay/washout, and receptor-mediated adsorption
  Real virion_burst_size = 12.1;      // De Paepe et al. (2016) gut burst size
  Real D_virion = 4.0e-12;            // ~50-100 nm phage Stokes-Einstein in mucus (m^2/s)
  Real virion_retardation = 1.0;      // mucin retardation factor (>= 1.0)
  Real k_virion_decay = 1.0e-4;       // first-order virion decay (1/s)
  Real virion_washout_scale = 1.0;    // coupling to AdvectionField::washout_rate(z)
  Real k_adsorb = 1.0e-15;            // adsorption constant (m^3/s per cell at expr=1)
  Real p_lysogeny = 0.19;             // De Paepe et al. (2016) lysogenization probability
  std::string adsorption_receptor = "none";  // "none" (constitutive) or ReceptorType name
};

// Pure Hill induction hazard evaluation (1/s) for unit tests and CPU/GPU parity.
[[nodiscard]] Real prophage_hill_induction_hazard(Real damage_signal,
                                                  const ProphageConfig& cfg);

// Resolve configured adsorption receptor index (-1 for constitutive "none").
[[nodiscard]] Int prophage_adsorption_receptor_index(
    const ProphageConfig& cfg);

class FixProphage : public Fix {
 public:
  FixProphage(Simulation& sim, const ProphageConfig& cfg);

  void init() override;
  void compute(Real dt) override;

  // Compute local DNA-damage / SOS driver signal D for one agent.
  [[nodiscard]] Real compute_damage_signal(const Agent& agent,
                                           Int agent_index) const;

  // Total free virions currently stored across owned voxels of species::VIRION.
  [[nodiscard]] static Real total_virions_in_field(const Simulation& sim);

 private:
  void resolve_latent_lysis(Int i_virion, Real cell_vol, Real dt);
  void run_hill_induction(Int i_virion, Real cell_vol, Real dt);
  void run_virion_adsorption_and_decay(Int i_virion, Real cell_vol, Real dt);
  void lyse_induced_agent(Agent& agent, Int i_virion, Real cell_vol, Real dt);
  void record_lytic_kill(Agent& agent, Int i_virion, Real cell_vol, Real dt);

  ProphageConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_PROPHAGE_H
