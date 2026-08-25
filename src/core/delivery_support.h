/* -----------------------------------------------------------------------
   GutIBM – Physical delivery-support enumeration
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_DELIVERY_SUPPORT_H
#define GUTIBM_DELIVERY_SUPPORT_H

#include "domain.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace gutibm {

struct DeliverySupportOffset {
  Int dx = 0;
  Int dy = 0;
  Int dz = 0;
  bool within_radius = false;
};

struct DeliverySupportStencil {
  Real radius = -1.0;
  Real dx = 0.0;
  Real dy = 0.0;
  Real dz = 0.0;
  Int nx = 0;
  Int ny = 0;
  Int nz = 0;
  bool periodic_x = false;
  bool periodic_y = false;
  std::vector<DeliverySupportOffset> offsets;

  bool matches(const Domain& domain, Real requested_radius) const {
    const auto& periodic = domain.config().periodic;
    return radius == requested_radius
        && dx == domain.dx_x()
        && dy == domain.dx_y()
        && dz == domain.dx_z()
        && nx == domain.nx()
        && ny == domain.ny()
        && nz == domain.nz()
        && periodic_x == periodic[0]
        && periodic_y == periodic[1];
  }
};

inline DeliverySupportStencil make_delivery_support_stencil(
    const Domain& domain, Real radius) {
  DeliverySupportStencil stencil;
  stencil.radius = radius;
  stencil.dx = domain.dx_x();
  stencil.dy = domain.dx_y();
  stencil.dz = domain.dx_z();
  stencil.nx = domain.nx();
  stencil.ny = domain.ny();
  stencil.nz = domain.nz();
  stencil.periodic_x = domain.config().periodic[0];
  stencil.periodic_y = domain.config().periodic[1];
  if (radius <= 0.0) return stencil;

  const Int x_span = static_cast<Int>(
      std::ceil(radius / domain.dx_x())) + 1;
  const Int y_span = static_cast<Int>(
      std::ceil(radius / domain.dx_y())) + 1;
  const Int z_span = static_cast<Int>(
      std::ceil(radius / domain.dx_z())) + 1;
  const bool all_x = stencil.periodic_x
      && 2 * x_span + 1 >= domain.nx();
  const bool all_y = stencil.periodic_y
      && 2 * y_span + 1 >= domain.ny();
  const Int x_begin = all_x ? 0 : -x_span;
  const Int x_end = all_x ? domain.nx() - 1 : x_span;
  const Int y_begin = all_y ? 0 : -y_span;
  const Int y_end = all_y ? domain.ny() - 1 : y_span;
  const Real radius_sq = radius * radius;
  stencil.offsets.reserve(static_cast<size_t>(
      (x_end - x_begin + 1) * (y_end - y_begin + 1)
      * (2 * z_span + 1)));
  const auto periodic_distance = [](Int offset, Int extent) {
    Int wrapped = offset % extent;
    if (wrapped < 0) wrapped += extent;
    return std::min(wrapped, extent - wrapped);
  };
  for (Int offset_x = x_begin; offset_x <= x_end; ++offset_x) {
    const Int image_x = stencil.periodic_x
        ? periodic_distance(offset_x, domain.nx()) : std::abs(offset_x);
    for (Int offset_y = y_begin; offset_y <= y_end; ++offset_y) {
      const Int image_y = stencil.periodic_y
          ? periodic_distance(offset_y, domain.ny())
          : std::abs(offset_y);
      for (Int offset_z = -z_span; offset_z <= z_span; ++offset_z) {
        const Real distance_sq =
            static_cast<Real>(image_x * image_x) * domain.dx_x()
                * domain.dx_x()
            + static_cast<Real>(image_y * image_y) * domain.dx_y()
                * domain.dx_y()
            + static_cast<Real>(offset_z * offset_z) * domain.dx_z()
                * domain.dx_z();
        stencil.offsets.push_back(
            {offset_x, offset_y, offset_z, distance_sq <= radius_sq});
      }
    }
  }
  return stencil;
}

inline bool delivery_support_target(
    const Domain& domain, Int center_ix, Int center_iy, Int center_iz,
    const DeliverySupportOffset& offset, Int& cell) {
  Int ix = center_ix + offset.dx;
  Int iy = center_iy + offset.dy;
  const Int iz = center_iz + offset.dz;
  if (domain.config().periodic[0]) {
    ix %= domain.nx();
    if (ix < 0) ix += domain.nx();
  } else if (ix < 0 || ix >= domain.nx()) {
    return false;
  }
  if (domain.config().periodic[1]) {
    iy %= domain.ny();
    if (iy < 0) iy += domain.ny();
  } else if (iy < 0 || iy >= domain.ny()) {
    return false;
  }
  if (iz < 0 || iz >= domain.nz()) return false;
  cell = domain.cell_index(ix, iy, iz);
  return true;
}

inline void enumerate_physical_delivery_ball(
    const Domain& domain, const Vec3& point,
    const DeliverySupportStencil& stencil, std::vector<Int>& support) {
  support.clear();
  if (stencil.radius <= 0.0) return;
  Int center_ix = 0;
  Int center_iy = 0;
  Int center_iz = 0;
  domain.pos_to_grid(point, center_ix, center_iy, center_iz);
  support.reserve(stencil.offsets.size());
  for (const auto& offset : stencil.offsets) {
    Int cell = -1;
    if (!delivery_support_target(
            domain, center_ix, center_iy, center_iz, offset, cell)) {
      continue;
    }
    const Vec3 center = domain.cell_center(
        cell % domain.nx(), (cell / domain.nx()) % domain.ny(),
        cell / (domain.nx() * domain.ny()));
    if (domain.min_image_dist_sq(point, center)
        <= stencil.radius * stencil.radius) {
      support.push_back(cell);
    }
  }
}

inline std::vector<Int> enumerate_physical_delivery_ball(
    const Domain& domain, const Vec3& point, Real radius) {
  const DeliverySupportStencil stencil =
      make_delivery_support_stencil(domain, radius);
  std::vector<Int> support;
  enumerate_physical_delivery_ball(domain, point, stencil, support);
  return support;
}

}  // namespace gutibm

#endif  // GUTIBM_DELIVERY_SUPPORT_H
