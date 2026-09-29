/* -----------------------------------------------------------------------
   GutIBM – Spec 13 Phase 1: tabulated patch growth/capacity

   Tabulated lookup over the carbon-amplitude ladder measured in
   docs/CARBON_LADDER_CAMPAIGN.md (#314, commit a115447): 200x200x100 um
   domain (4.0e-6 mL), 100 founders, metabolism.uptake_limit=delivery,
   oxygen species disabled, oxygen.k_ROS=0.0, dysbiosis_threshold=1e10,
   12 h equilibration. Each arm reports peak cell count and measured
   mu/mu_max; the logistic coefficient uses mu_max = 5e-4 /s.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_LAYER2_PATCH_TABLE_H
#define GUTIBM_LAYER2_PATCH_TABLE_H

#include "types.h"

namespace gutibm {

class PatchTable {
 public:
  PatchTable() = default;

  // Carrying capacity (cells) for a patch of `volume_mL` fed at
  // `supply_mult` times the 1x ladder arm's epithelial carbon amplitude.
  // Capacity scales linearly with volume relative to the 4.0e-6 mL
  // ladder domain; interpolation between arms is log-linear in
  // supply_mult, clamped at the table ends.
  [[nodiscard]] Real capacity_cells(Real supply_mult,
                                    Real volume_mL) const;

  // Logistic growth coefficient (1/s) at `supply_mult`, interpolated
  // log-linear between ladder arms, clamped at the table ends.
  [[nodiscard]] Real growth_rate(Real supply_mult) const;

 private:
  static constexpr Real kLadderVolumeMl = 4.0e-6;
  static constexpr Real kMuMax = 5.0e-4;  // strain mu_max used by the ladder
  static constexpr int kRows = 5;
  // Measured ladder: {supply_mult, peak_agents, mu/mu_max}
  static constexpr Real kSupply[kRows] = {0.25, 0.5, 1.0, 2.0, 4.0};
  static constexpr Real kPeak[kRows] = {100, 100, 122, 1063, 9046};
  static constexpr Real kMuFrac[kRows] = {0.029, 0.083, 0.082, 0.197, 0.189};

  static Real interp_log(Real x, const Real (&ys)[kRows]);
};

}  // namespace gutibm

#endif  // GUTIBM_LAYER2_PATCH_TABLE_H
