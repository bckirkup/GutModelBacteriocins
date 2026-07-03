/* -----------------------------------------------------------------------
   GutIBM – Fur-regulated receptor expression config (Spec 3)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FUR_CONFIG_H
#define GUTIBM_FUR_CONFIG_H

#include "types.h"

namespace gutibm {

struct FurConfig {
  bool enabled = true;
  Real Km = 1.0e-5;              // mol/m³ — half-max Fur repression iron
  Real upregulation_max = 4.0;   // max fold-upregulation under starvation
  Real receptor_max = 5.0;       // cap on effective expression
};

}  // namespace gutibm

#endif  // GUTIBM_FUR_CONFIG_H
