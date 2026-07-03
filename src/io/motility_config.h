/* -----------------------------------------------------------------------
   GutIBM – Active cell motility config (Spec 3)
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
  bool chemotaxis_enabled = false;
  Real chi_carbon = 0.1;
  Real chi_oxygen = 0.1;
  Real cluster_suppress_radius = 10.0e-6; // m
  Int cluster_suppress_threshold = 5;
  Real cluster_tumble_factor = 0.2;
};

}  // namespace gutibm

#endif  // GUTIBM_MOTILITY_CONFIG_H
