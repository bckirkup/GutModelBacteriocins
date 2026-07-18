/* -----------------------------------------------------------------------
   GutIBM – Active cell motility config (Spec 3 / Spec 10v2)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_MOTILITY_CONFIG_H
#define GUTIBM_MOTILITY_CONFIG_H

#include "types.h"

namespace gutibm {

struct MotilityConfig {
  bool enabled = true;
  Real swim_speed = 7.76e-6;              // m/s
  Real run_mean_duration = 1.0;           // s
  Real stop_probability = 0.3;
  Real stop_duration = 0.5;               // s

  // Carbon chemotaxis (Tar/Tsr MCP pathway — Weber–Fechner)
  bool chemotaxis_enabled = false;
  Real chi_carbon = 2.0;                  // dimensionless sensitivity
  Real chemotaxis_threshold = 1.0e-6;     // mol/m³ floor for fractional sensing

  // Aerotaxis (Aer receptor — primary directional cue in mucus)
  bool aerotaxis_enabled = true;
  Real aerotaxis_sensitivity = 4.0;       // dimensionless; stronger than carbon

  // Energy taxis: speed reduction under metabolic stress
  bool energy_taxis_enabled = true;
  Real energy_taxis_floor = 0.1;          // minimum speed fraction at mu=0

  // Surface sensing: speed reduction near epithelium (opt-in)
  bool surface_sensing_enabled = false;
  Real surface_sensing_depth = 10.0e-6;   // m
  Real surface_sensing_floor = 0.3;

  // Mucin viscosity drag (opt-in; needs the mucin chemical field enabled)
  bool mucin_drag_enabled = false;
  Real mucin_drag_reference = 1.0e-2;     // mol/m³

  Real cluster_suppress_radius = 10.0e-6; // m
  Int cluster_suppress_threshold = 5;
  Real cluster_tumble_factor = 0.2;
};

}  // namespace gutibm

#endif  // GUTIBM_MOTILITY_CONFIG_H
