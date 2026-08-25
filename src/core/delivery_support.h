/* -----------------------------------------------------------------------
   GutIBM – Physical delivery-support enumeration
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_DELIVERY_SUPPORT_H
#define GUTIBM_DELIVERY_SUPPORT_H

#include "domain.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace gutibm {

inline std::vector<Int> enumerate_physical_delivery_ball(
    const Domain& domain, const Vec3& point, Real radius) {
  if (radius <= 0.0) return {};

  Int center_ix = 0;
  Int center_iy = 0;
  Int center_iz = 0;
  domain.pos_to_grid(point, center_ix, center_iy, center_iz);
  const Int x_span = static_cast<Int>(
      std::ceil(radius / domain.dx_x())) + 1;
  const Int y_span = static_cast<Int>(
      std::ceil(radius / domain.dx_y())) + 1;
  const Int z_span = static_cast<Int>(
      std::ceil(radius / domain.dx_z())) + 1;
  const auto& periodic = domain.config().periodic;
  const bool all_x = 2 * x_span + 1 >= domain.nx();
  const bool all_y = 2 * y_span + 1 >= domain.ny();
  const Real radius_sq = radius * radius;
  std::vector<Int> support;

  const Int x_begin = all_x ? 0 : center_ix - x_span;
  const Int x_end = all_x ? domain.nx() - 1 : center_ix + x_span;
  const Int y_begin = all_y ? 0 : center_iy - y_span;
  const Int y_end = all_y ? domain.ny() - 1 : center_iy + y_span;
  const Int z_begin = std::max(0, center_iz - z_span);
  const Int z_end = std::min(domain.nz() - 1, center_iz + z_span);
  for (Int ix = x_begin; ix <= x_end; ++ix) {
    Int candidate_ix = ix;
    if (periodic[0]) {
      candidate_ix %= domain.nx();
      if (candidate_ix < 0) candidate_ix += domain.nx();
    } else if (candidate_ix < 0 || candidate_ix >= domain.nx()) {
      continue;
    }
    for (Int iy = y_begin; iy <= y_end; ++iy) {
      Int candidate_iy = iy;
      if (periodic[1]) {
        candidate_iy %= domain.ny();
        if (candidate_iy < 0) candidate_iy += domain.ny();
      } else if (candidate_iy < 0 || candidate_iy >= domain.ny()) {
        continue;
      }
      for (Int iz = z_begin; iz <= z_end; ++iz) {
        const Int cell = domain.cell_index(candidate_ix, candidate_iy, iz);
        const Vec3 center = domain.cell_center(
            candidate_ix, candidate_iy, iz);
        if (domain.min_image_dist_sq(point, center) <= radius_sq) {
          support.push_back(cell);
        }
      }
    }
  }
  return support;
}

}  // namespace gutibm

#endif  // GUTIBM_DELIVERY_SUPPORT_H
