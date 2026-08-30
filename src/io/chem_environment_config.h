/* -----------------------------------------------------------------------
   GutIBM – Chemical environment expansion config (Spec 1)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_CHEM_ENVIRONMENT_CONFIG_H
#define GUTIBM_CHEM_ENVIRONMENT_CONFIG_H

#include "types.h"
#include <string>

namespace gutibm {

enum class RespirationDriver : int {
  Ambient = 0,
  Funded = 1,
};

enum class RosDriver : int {
  Ambient = 0,
  Funded = 1,
};

struct OxygenConfig {
  bool enabled = false;
  Real epithelial_conc = 55.0e-6;   // mol/m^3 (~0.042 mmHg; 42 mmHg ≈ 5.5e-2 mol/m^3)
  bool delivery_uptake_enabled = false;
  std::string respiration_driver = "ambient";
  RespirationDriver respiration_driver_mode = RespirationDriver::Ambient;
  std::string ros_driver = "ambient";
  RosDriver ros_driver_mode = RosDriver::Ambient;
  Real D_free = 2.1e-9;             // m^2/s
  Real Km = 1.0e-6;                 // mol/m^3 Monod half-saturation
  Real boost_max = 2.0;             // max aerobic growth multiplier - 1
  bool metabolic_switch_enabled = false;
  Real mu_crit = 9.7e-5;            // 0.35 h^-1 acetate-overflow onset
  Real aerobic_mu_factor = 1.0;
  Real anaerobic_mu_factor = 0.55;
  // Multipliers on substrate-per-biomass yield (not biomass-per-substrate).
  Real aerobic_carbon_cost_factor = 1.0;
  Real anaerobic_carbon_cost_factor = 4.1;
  Real tau_metabolic_switch = 3600.0; // s
  Real ferm_acid_yield = 1.0;          // acetate per carbon-equivalent consumed
  Real anaerobic_maintenance_factor = 15.0;
  // Pirt-style respiration: OUR_cell = q_consumption * mu_realized (growth-
  // associated) + q_maintenance (basal, density-coupled). The maintenance term
  // is what makes the O2 field track cell *density* rather than only growth:
  // a present-but-non-growing cell (e.g. washing out, mu->0) still respires.
  Real q_consumption = 1.0e-14;     // mol/cell (growth-associated O2 per unit mu)
  Real q_maintenance = 1.0e-18;     // mol/s/cell basal respiration (density-coupled)
  // First-order background O2 uptake RATE CONSTANT (1/s) by the anaerobic
  // majority: reac -= vbf_sink * [O2] (see apply_oxygen_sink). NOT a zero-order
  // mol/m^3/s removal — that form removes O2 that isn't there, hard-zeroing the
  // interior in one bio step and masking per-agent respiration.
  Real vbf_sink = 1.0e-3;           // 1/s first-order background O2 uptake rate
  Real k_ROS = 0.0;                 // opt-in ambient ROS-driven SOS coefficient
  Real k_ROS_respiratory = 0.0;     // kg/mol funded specific-flux coefficient
  Real k_ROS_funded = 0.0;          // mol^-1 funded absolute-flux coefficient
};

struct AcetateConfig {
  bool enabled = false;
  Real D_free = 1.2e-9;
  Real vbf_production = 1.0e-3;     // mol/m^3/s
  Real vbf_consumption = 2.0e-4;    // mol/m^3/s
  Real overflow_threshold = 3.0e-4;   // 1/s
  Real overflow_rate = 1.0e-15;     // mol/s/cell overflow
  Real scavenge_rate = 1.0e-15;     // mol/s/cell max scavenging
  Real scavenge_Km = 5.0;           // mol/m^3
  Real epithelial_uptake = 5.0e-4;  // mol/m^3/s at z=0
};

struct MucinConfig {
  bool enabled = false;
  Real initial_conc = 1.0e-2;       // mol/m^3 bulk mucin polymer
  Real secretion_rate = 1.0e-4;     // mol/m^3/s at epithelium
  Real Km_degradation = 1.0e-3;     // mol/m^3
  Real k_liberation = 5.0e-16;       // mol/(cell·s) per-cell specific liberation rate
  Real D_free = 1.0e-12;            // effectively immobile
  Real retardation = 1000.0;
};

struct ProteaseConfig {
  bool enabled = true;
  Real default_half_life = 1800.0;    // s
  Real dilution_rate = 1.0e-4;        // 1/s fallback when advection washout negligible
};

struct SiderophoreConfig {
  bool enabled = true;
  Real secretion_rate = 1.0e-5;       // mol/(s·kg), constrained estimate
  Real D_free = 1.0e-10;              // m^2/s
  Real chelation_rate = 1.0e3;        // m^3/(mol·s) effective second-order
  Real Km_reimport = 1.0e-6;          // mol/m^3 for FepA-mediated ferric enterobactin reimport
  Real Vmax_reimport = 1.0e-5;        // mol/(s·kg), FepA FeEnt transport capacity
};

struct FerrichromeConfig {
  bool enabled = false;
  Real initial_conc = 0.0;             // ambient ferrichrome concentration (mol/m^3)
  Real boundary_conc = 0.0;            // epithelial/luminal boundary concentration (mol/m^3)
};

}  // namespace gutibm

#endif  // GUTIBM_CHEM_ENVIRONMENT_CONFIG_H
