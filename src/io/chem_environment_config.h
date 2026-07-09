/* -----------------------------------------------------------------------
   GutIBM – Chemical environment expansion config (Spec 1)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_CHEM_ENVIRONMENT_CONFIG_H
#define GUTIBM_CHEM_ENVIRONMENT_CONFIG_H

#include "types.h"

namespace gutibm {

struct OxygenConfig {
  bool enabled = false;
  Real epithelial_conc = 55.0e-6;   // mol/m^3 (~42 mmHg)
  Real D_free = 2.1e-9;             // m^2/s
  Real Km = 1.0e-6;                 // mol/m^3 Monod half-saturation
  Real boost_max = 2.0;             // max aerobic growth multiplier - 1
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
  Real k_ROS = 1.0e2;               // ROS-driven SOS rate coefficient
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
  Real k_liberation = 1.0e-4;       // 1/s rate constant
  Real D_free = 1.0e-12;            // effectively immobile
  Real retardation = 1000.0;
};

struct ProteaseConfig {
  bool enabled = true;
  Real default_half_life = 1800.0;    // s
  Real dilution_rate = 1.0e-4;        // 1/s fallback when advection washout negligible
};

struct SiderophoreConfig {
  bool enabled = false;
  Real secretion_rate = 1.0e-15;      // mol/s/cell
  Real D_free = 1.0e-10;              // m^2/s
  Real chelation_rate = 1.0e3;        // m^3/(mol·s) effective second-order
  Real Km_reimport = 1.0e-9;          // mol/m^3 for FepA-mediated iron recapture
  Real recapture_fraction = 0.5;      // fraction of secreted siderophore iron returned locally
};

}  // namespace gutibm

#endif  // GUTIBM_CHEM_ENVIRONMENT_CONFIG_H
