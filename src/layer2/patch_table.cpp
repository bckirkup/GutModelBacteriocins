#include "patch_table.h"

#include <algorithm>
#include <cmath>

namespace gutibm {

Real PatchTable::interp_log(Real x, const Real (&ys)[kRows]) {
  const Real lx = std::log(x);
  const Real lo = std::log(kSupply[0]);
  const Real hi = std::log(kSupply[kRows - 1]);
  if (lx <= lo) {
    return ys[0];
  }
  if (lx >= hi) {
    return ys[kRows - 1];
  }
  for (int i = 1; i < kRows; ++i) {
    const Real x1 = std::log(kSupply[i]);
    if (lx <= x1) {
      const Real x0 = std::log(kSupply[i - 1]);
      const Real w = (lx - x0) / (x1 - x0);
      const Real y0 = ys[i - 1];
      const Real y1 = ys[i];
      // Log-linear in supply; linear in the measured quantity (capacity
      // is non-negative everywhere in the ladder, so no log-of-zero).
      return y0 + w * (y1 - y0);
    }
  }
  return ys[kRows - 1];
}

Real PatchTable::capacity_cells(Real supply_mult, Real volume_mL) const {
  const Real density = interp_log(supply_mult, kPeak) / kLadderVolumeMl;
  return std::max(density * volume_mL, 1.0);
}

Real PatchTable::growth_rate(Real supply_mult) const {
  return interp_log(supply_mult, kMuFrac) * kMuMax;
}

}  // namespace gutibm
