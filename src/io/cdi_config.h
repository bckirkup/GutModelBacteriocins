/* -----------------------------------------------------------------------
   GutIBM – Contact-dependent inhibition config (Spec 3)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_CDI_CONFIG_H
#define GUTIBM_CDI_CONFIG_H

#include "types.h"

namespace gutibm {

struct CdiConfig {
  bool enabled = true;
  Real kill_rate = 5.0e-4;           // 1/s per contact pair
  Real contact_radius = 1.0e-6;      // m
  Real corpse_persistence = 300.0;   // s
};

}  // namespace gutibm

#endif  // GUTIBM_CDI_CONFIG_H
