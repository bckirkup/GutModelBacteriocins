/* -----------------------------------------------------------------------
   GutIBM – Advection field implementation
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "domain.h"
#include <cmath>
#include <algorithm>

namespace gutibm {

void AdvectionField::init(const AdvectionConfig& cfg, const Domain& domain) {
  cfg_  = cfg;
  h_    = cfg.mucus_thickness;
  lo_z_ = domain.lo()[2];

  // Max velocities at lumen surface (z = h)
  v_radial_max_ = h_ / cfg.radial_turnover;
  v_distal_max_ = cfg.distal_length / cfg.distal_transit_time;
}

Real AdvectionField::radial_velocity(Real z) const {
  Real zn = std::clamp((z - lo_z_) / h_, 0.0, 1.0);
  return v_radial_max_ * std::pow(zn, cfg_.profile_alpha);
}

Real AdvectionField::distal_velocity(Real z) const {
  Real zn = std::clamp((z - lo_z_) / h_, 0.0, 1.0);
  return v_distal_max_ * std::pow(zn, cfg_.profile_alpha);
}

Vec3 AdvectionField::velocity(const Vec3& pos) const {
  Real vx = distal_velocity(pos[2]);   // distal (x-direction)
  Real vy = 0.0;                       // no lateral flow
  Real vz = radial_velocity(pos[2]);   // radial (z, toward lumen)
  return {vx, vy, vz};
}

Real AdvectionField::shear_rate(const Vec3& pos) const {
  // dv/dz at position
  Real zn = std::clamp((pos[2] - lo_z_) / h_, 0.0, 1.0);
  Real alpha = cfg_.profile_alpha;

  Real dvr_dz = v_radial_max_ * alpha * std::pow(zn, alpha - 1.0) / h_;
  Real dvd_dz = v_distal_max_ * alpha * std::pow(zn, alpha - 1.0) / h_;

  return std::sqrt(dvr_dz * dvr_dz + dvd_dz * dvd_dz);
}

void AdvectionField::advect(Vec3& pos, Real dt) const {
  Vec3 v = velocity(pos);
  pos[0] += v[0] * dt;
  pos[1] += v[1] * dt;
  pos[2] += v[2] * dt;
}

Real AdvectionField::washout_rate(Real z) const {
  return radial_velocity(z) / h_;
}

}  // namespace gutibm
