/* -----------------------------------------------------------------------
   GutIBM – AI-2 quorum sensing config (Spec 11)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_QUORUM_SENSING_CONFIG_H
#define GUTIBM_QUORUM_SENSING_CONFIG_H

#include "types.h"

namespace gutibm {

struct QuorumSensingConfig {
  bool enabled = false;                 // master switch; registers AI-2 species
  Real ai2_basal_rate = 1.0e-20;        // mol/s/cell constitutive LuxS
  Real ai2_growth_coupled = 1.0e-16;    // mol/cell growth-associated
  Real lsr_vmax = 1.0e-18;              // mol/s/cell max Lsr import
  Real lsr_km = 1.0e-7;                 // mol/m³ Lsr half-saturation
  Real ai2_D_free = 5.0e-10;            // m²/s
  Real ai2_decay_rate = 1.0e-4;         // 1/s background degradation

  // Chemotaxis knobs live here; FixMotility reads them for Weber–Fechner AI-2
  bool ai2_chemotaxis_enabled = false;
  Real chi_ai2 = 3.0;
};

}  // namespace gutibm

#endif  // GUTIBM_QUORUM_SENSING_CONFIG_H
