/* -----------------------------------------------------------------------
   GutIBM – Mucus advection field
   Dual-vector flow: radial (epithelium→lumen, 1–2 h turnover) and
   distal (peristaltic, 8–24 h transit).  Corrected per EARI review.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_ADVECTION_H
#define GUTIBM_ADVECTION_H

#include "types.h"
#include <vector>

namespace gutibm {

class Domain;

struct AdvectionConfig {
  // Radial mucus shedding (z-axis, epithelium→lumen)
  Real radial_turnover     = 5400.0;  // 1.5 h in seconds
  Real mucus_thickness     = 100.0e-6; // 100 um

  // Distal peristaltic flow (x-axis, proximal→distal)
  Real distal_transit_time = 43200.0;  // 12 h in seconds
  Real distal_length       = 1.0e-3;   // domain length (m)

  // Velocity profile: near-zero at epithelium (z=0), max at lumen (z=max)
  // Parabolic profile: v(z) = v_max * (z/h)^alpha
  Real profile_alpha       = 1.5;      // shear profile exponent

  // Taylor-Aris effective dispersion enhancement
  // D_eff = D_mol + (U^2 * h^2) / (210 * D_mol)  for Poiseuille flow
  bool  taylor_aris_enabled = true;

  // Crypt refugia: zero-flow zones near epithelium (VADI §80, §98-99)
  bool  crypts_enabled       = false;
  Real  crypt_depth          = 10.0e-6;  // z < lo_z + crypt_depth is zero-flow (10 um)
  Real  crypt_exit_rate      = 1.0e-4;   // per-second probability of agent exiting crypt
  Real  crypt_entry_rate     = 5.0e-5;   // per-second probability of agent entering crypt
  Int   crypt_carrying_capacity = 50;    // max agents per crypt region
};

class AdvectionField {
 public:
  AdvectionField() = default;

  void init(const AdvectionConfig& cfg, const Domain& domain);

  // Flow velocity at a position (m/s)
  Vec3 velocity(const Vec3& pos) const;

  // Radial velocity magnitude at height z
  Real radial_velocity(Real z) const;

  // Distal velocity magnitude at height z
  Real distal_velocity(Real z) const;

  // Local shear rate magnitude (for conjugation MPS calc)
  Real shear_rate(const Vec3& pos) const;

  // Advect an agent position by dt
  void advect(Vec3& pos, Real dt) const;

  // Washout rate at height z (1/s)
  // gamma_flow = v_radial(z) / mucus_thickness
  Real washout_rate(Real z) const;

  // Taylor-Aris effective dispersion coefficient at height z
  // Accounts for shear-enhanced longitudinal spreading
  Real taylor_aris_D_eff(Real z, Real D_mol) const;

  // Query whether a z-coordinate falls within the crypt zone
  bool in_crypt_zone(Real z) const;

  const AdvectionConfig& config() const { return cfg_; }

 private:
  AdvectionConfig cfg_;
  Real v_radial_max_ = 0.0;
  Real v_distal_max_ = 0.0;
  Real h_ = 100.0e-6;
  Real lo_z_ = 0.0;
};

}  // namespace gutibm

#endif  // GUTIBM_ADVECTION_H
