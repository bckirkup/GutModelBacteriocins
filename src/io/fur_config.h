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
  // Spec 6 §4.2 — Fur derepression under iron starvation. Raised 4x -> 10x;
  // measured induction of the Fur regulon is 35-56x, so 10x remains a
  // conservative lower bound while restoring meaningful siderophore-receptor
  // upregulation. Effective expression is still capped by receptor_max.
  Real upregulation_max = 10.0;  // max fold-upregulation under starvation
  Real receptor_max = 5.0;       // cap on effective expression
};

}  // namespace gutibm

#endif  // GUTIBM_FUR_CONFIG_H
