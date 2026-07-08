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
  // First-order iron uptake RATE CONSTANT (1/s): the sink applied by the VBF is
  // concentration-dependent, reac -= nutrient_sink * [iron] (unsaturated /
  // Monod-at-low-concentration limit), not a constant zero-order mol/m^3/s
  // removal. At the default iron scale (~1e-4 mol/m^3) a zero-order 1e-4
  // mol/m^3/s would deplete iron in ~1 s, which is unphysical.
  Real nutrient_sink    = 1.0e-4;   // first-order iron uptake rate constant (1/s)
  Real mucin_liberation = 5.0e-5;   // monosaccharide release from mucin (mol/m^3/s)
  Real carrying_cap     = 1.0e12;   // local carrying capacity (#/m^3)
  Real viscosity        = 0.01;     // effective viscosity (Pa·s), ~10x water

  // Optional z-dependent mucin liberation scaling with epithelial distance
  bool mucin_z_gradient_enabled = false;
  Real mucin_z_gradient_lambda  = 25.0e-6;  // characteristic decay length (m)

  // When true, carbon source comes from dynamic mucin degradation (Spec 1)
  bool use_dynamic_mucin = false;

  // Spec 5 §1 — VBF carbon sink. The anaerobic majority consumes most of the
  // liberated monosaccharides ("Restaurant Hypothesis"); E. coli scavenges the
  // leftovers. Modeled as a Monod-saturating first-order-at-low-conc sink so
  // carbon reaches a bounded steady state instead of accumulating without
  // bound. Disabled by default (vmax = 0) to preserve legacy behavior; set
  // vbf_carbon_sink_vmax > 0 to activate.
  Real carbon_sink_vmax = 0.0;      // mol/m^3/s max VBF carbon consumption
  Real carbon_sink_km   = 1.0e-3;   // mol/m^3 half-saturation

  // Spec 5 §3 — VBF B12 (cobalamin) production. B12 is consumed by E. coli but
  // has no source, so it drains to zero and forces all cells onto the MetE
  // pathway (a drift artifact, not biology). The anaerobic community produces
  // B12 at a low constant rate that maintains a homeostatic steady state.
  // Disabled by default (0) to preserve legacy behavior.
  Real b12_production = 0.0;         // mol/m^3/s constant B12 source
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
