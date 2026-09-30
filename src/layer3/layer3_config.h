/* -----------------------------------------------------------------------
   GutIBM – Spec 13 Phase 3: Layer-3 colonic gradient configuration

   Three regions (cecum/ascending, transverse, descending/sigmoid), each
   one Layer-2 MucusSegment, coupled by unidirectional luminal transit
   and a per-lineage luminal population compartment. Regional values
   are axial PROFILES (S7): every per-region field is a sourced value
   from the SPEC13_MULTISCALE Layer-3 table — never a free fit —
   collapsed to the transverse row when `uniform_profile` is set (the
   gradient sanity arm). Everything is off by default
   (layer3.enabled = false) and leaves default-mode behaviour untouched.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_LAYER3_CONFIG_H
#define GUTIBM_LAYER3_CONFIG_H

#include "types.h"
#include <string>

namespace gutibm {

// The three evidence-classed regions of the Layer-3 table, in transit
// order. Parser heads: cecum, transverse, descending.
enum class RegionKind : Int {
  CecumAscending = 0,
  Transverse = 1,
  DescendingSigmoid = 2,
};
constexpr Int kRegionCount = 3;
const char* region_name(RegionKind r);

// One row of the Layer-3 axial profile. Each field is a value read
// off the SPEC13_MULTISCALE table (piecewise-constant interpolation —
// the resolution the underlying measurements support); see the table
// for per-row citations.
struct Layer3RegionConfig {
  Real transit_h = 10.0;               // segment transit time (sourced range 5-15 h)
  Real ph = 6.6;                       // luminal pH (6.0-7.0 axial rise)
  Real lumen_carbon_mol_m3 = 5.0e-3;   // available carbon for luminal growth
  Real lumen_volume_L = 0.35;          // luminal content volume of the region
  Real contraction_k0_per_min = 0.4;   // baseline contraction frequency (<=0.6 admissible)
  Real outflow_survival = 0.9;         // fraction of outflow cells reaching downstream
  Real flora_density_cfu_ml = 1.0e7;   // background flora proxy (flatness denominator)
  Real patch_supply_scale = 1.0;       // scalar on per-type supply_mults
};

struct Layer3Config {
  bool enabled = false;
  // Sanity arm: collapse every region to the transverse profile so the
  // run reports the same physics with no axial gradient (S7).
  bool uniform_profile = false;

  // Sourced defaults = transverse row midpoints; per-region values
  // below carry the axial profile.
  Layer3RegionConfig cecum{12.0, 6.25, 8.0e-3, 0.40, 0.50, 0.92, 1.0e7, 1.0};
  Layer3RegionConfig transverse{8.0, 6.60, 5.0e-3, 0.35, 0.40, 0.90, 1.0e7, 1.0};
  Layer3RegionConfig descending{10.0, 6.80, 2.5e-3, 0.30, 0.30, 0.88, 1.0e7, 1.0};

  Int n_patches_per_region = 400;

  // Mucosa->lumen transfer coefficients (spec's three named channels).
  // Channel 1 (mucus turnover) rides the per-type single_cell_loss
  // rates on the Layer-2 config; channel 2 (contractions) is k0 above
  // plus the HAPC boost below; channel 3 is growth-edge shedding.
  Real f_edge = 0.2;               // fraction of mucosal births shed per step
  Real reattach_p0 = 0.05;         // p0 in p_attach = p0*(1 - occupancy_fraction)

  // HAPC schedule: ~6/day, 95% antegrade, 80% in the 12 h day window.
  Real hapc_rate_per_day = 6.0;
  Real hapc_antegrade_fraction = 0.95;
  Real hapc_day_fraction = 0.8;
  Real hapc_day_start_h = 8.0;
  Real alpha_hapc_per_min = 0.2;   // contraction-rate increment during an event

  // Stool observation model (continuous first-order R3 export).
  Real stool_g_per_day = 128.0;    // median daily stool mass for CFU/g

  // Declared flatness bound for the Enterobacteriaceae fraction readout:
  // max/min regional mucosal fraction <= bound. Ahmed 2007 reports
  // 1.3/2.5/0.8 % (P = 0.09, NS) -> observed span ~3.1x; a real axial
  // gradient like crypt occupancy would be >16x.
  Real flatness_bound = 4.0;

  // Luminal seeding + perturbation probe.
  Real lumen_seed_cells = 1.0e3;   // initial cells per region, strain[0] exemplar
  Real purge_at_s = -1.0;          // >=0: one-time luminal flush (recovery probe)
  Real purge_fraction = 0.9;

  std::string timeseries_file;        // CSV path; empty disables
  Int summary_interval_steps = 60;    // observables cadence
  std::string provenance_file;        // JSON run record; empty disables
  std::string checkpoint_file;        // HDF5 checkpoint path; empty disables
  Int checkpoint_interval_steps = 0;  // 0 = no periodic checkpoint
  bool checkpoint_final = false;      // write checkpoint at run end
};

}  // namespace gutibm

#endif  // GUTIBM_LAYER3_CONFIG_H
