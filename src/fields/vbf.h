/* -----------------------------------------------------------------------
   GutIBM – Viscoelastic Background Field (VBF)
   Represents the 99% obligate anaerobic microbiota as a continuous
   medium that: (1) exerts physical drag, (2) consumes nutrients as
   a biochemical sink, and (3) liberates monosaccharides from mucin
   to serve as a baseline carbon source.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_VBF_H
#define GUTIBM_VBF_H

#include "types.h"

namespace gutibm {

class Domain;
class ChemicalField;

struct OxygenConfig;
struct AcetateConfig;
struct MucinConfig;

struct VBFConfig {
  Real density          = 1.0e11;   // background cell density (#/m^3)
  Real drag_coeff       = 1.0e-9;   // Stokes-like drag (N·s/m)
  Real nutrient_sink    = 1.0e-4;   // volumetric consumption rate (mol/m^3/s)
  Real mucin_liberation = 5.0e-5;   // monosaccharide release from mucin (mol/m^3/s)
  Real carrying_cap     = 1.0e12;   // local carrying capacity (#/m^3)
  Real viscosity        = 0.01;     // effective viscosity (Pa·s), ~10x water

  // Optional z-dependent mucin liberation scaling with epithelial distance
  bool mucin_z_gradient_enabled = false;
  Real mucin_z_gradient_lambda  = 25.0e-6;  // characteristic decay length (m)

  // When true, carbon source comes from dynamic mucin degradation (Spec 1)
  bool use_dynamic_mucin = false;
};

class VBF {
 public:
  VBF() = default;

  void init(const VBFConfig& cfg, const Domain& domain);

  // Apply VBF nutrient sink/source to chemical field
  void apply_nutrient_coupling(ChemicalField& chem, const Domain& domain,
                                Real dt,
                                const OxygenConfig& oxygen,
                                const AcetateConfig& acetate,
                                const MucinConfig& mucin) const;

  // Compute drag force on an agent at position with velocity
  Vec3 drag_force(const Vec3& agent_vel) const;

  // Local carrying capacity at a grid cell
  Real local_capacity(Int /*cell_idx*/) const { return carrying_cap_; }

  // Effective viscosity at position (for conjugation shear calc)
  Real viscosity() const { return cfg_.viscosity; }

  Real density() const { return cfg_.density; }

  const VBFConfig& config() const { return cfg_; }

  // Mucin liberation rate at a z-position relative to the epithelium
  Real mucin_rate(Real z_rel) const;

 private:
  VBFConfig cfg_;
  Real carrying_cap_ = 0.0;
  Int ncells_ = 0;
};

}  // namespace gutibm

#endif  // GUTIBM_VBF_H
