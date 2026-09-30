/* -----------------------------------------------------------------------
   GutIBM – Spec 13 Phase 1: Layer-2 mucus-segment configuration

   The Layer-2 segment couples many mucosal patches through a shared
   luminal pool. Parameterization comes from docs/SPEC13_MULTISCALE.md
   (Layer-2 table) and docs/CARBON_LADDER_CAMPAIGN.md (tabulated
   growth/capacity lookups). Everything is off by default
   (layer2.enabled = false) and leaves default-mode behaviour untouched.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_LAYER2_CONFIG_H
#define GUTIBM_LAYER2_CONFIG_H

#include "types.h"
#include <string>

namespace gutibm {

// Patch type is persistent per-patch host identity (review S6). Type
// fractions are anatomical inputs, not fittable parameters.
enum class PatchType : Int { Crypt = 0, ExposedProximal = 1, ExposedDistal = 2 };
constexpr Int kPatchTypeCount = 3;
const char* patch_type_name(PatchType t);

// Per-patch status is distinct from run-level TerminationCause (S5): a
// guarded patch is marked invalid locally without killing the run.
enum class PatchStatus : Int {
  Active = 0,
  BloomHalted = 1,
  SpatialHalted = 2,
};
const char* patch_status_name(PatchStatus s);

// Two loss channels feed the one luminal pool (S1).
enum class PacketKind : Int { Single = 0, Fragment = 1 };
const char* packet_kind_name(PacketKind k);

struct Layer2TypeConfig {
  Real fraction = 0.0;                // anatomical fraction of all patches
  Real disruption_prob = 0.0;         // probability disrupted per contraction
  Real single_cell_loss_per_s = 0.0;  // first-order single-cell export rate
  Real supply_mult = 1.0;             // epithelial carbon amplitude vs 1x ladder arm
};

struct Layer2AuditConfig {
  Int patch_index = -1;        // -1 disables the live-patch audit gate
  Real interval_s = 21600.0;   // discrepancy comparison cadence (6 h)
  Real tolerance = 0.25;       // |N_live - N_shadow| / max(N_shadow, 1) bound
  bool halt_on_exceed = true;  // terminate run vs record-only
};

struct Layer2Config {
  bool enabled = false;

  Int n_patches = 200;
  Real patch_lateral_m = 300.0e-6;   // lateral extent of one patch
  Real patch_depth_m = 150.0e-6;     // mucus depth inside one patch
  Real mucus_thickness_m = 200.0e-6; // segment-level mucus thickness (CFU/cm^2)

  // Defaults at the midpoints of the SPEC13_MULTISCALE Layer-2 table:
  // crypt fraction 5-15%, exposed proximal 40-50%, exposed distal 35-50%;
  // crypt disruption 0.05-0.15, exposed 0.3-0.7; single-cell loss
  // 0.02-0.2 /h crypt-side of spec washout ranges.
  Layer2TypeConfig crypt{0.10, 0.10, 5.6e-6, 2.0};
  Layer2TypeConfig exposed_proximal{0.45, 0.50, 2.8e-5, 1.0};
  Layer2TypeConfig exposed_distal{0.45, 0.50, 5.6e-5, 0.5};

  Real contraction_rate_per_min = 0.5;    // segment contraction frequency
  Real agent_loss_fraction = 0.3;         // fraction of patch population per disruption
  Real transit_half_life_s = 14400.0;     // luminal transit half-life (2-8 h)
  Real reattach_prob_per_transit = 0.05;  // reattachment probability per transit
  Real establish_prob_single = 0.01;      // establishment kernel for single cells
  Real establish_ratio = 10.0;            // fragment:single establishment ratio

  // Spec 13 Phase 3 (Layer-3 coupling): growth-edge shedding sheds a
  // fraction of each step's tabulated births into the pool as Single
  // packets (spec mechanism 3, J = f_edge * mu * N). 0 disables and is
  // the Layer-2 default.
  Real edge_shed_fraction = 0.0;
  // High-bit prefix OR'd into the (patch_index + 1) << 32 cell-tag
  // bases so a chain of segments mints non-colliding tags. The Layer-3
  // driver sets one prefix per region; standalone Layer-2 runs leave
  // this 0.
  TagID tag_prefix = 0;

  Int occupancy_min_agents = 1;      // a patch counts occupied at >= this many cells

  Real initial_patch_fraction = 0.05;  // fraction of patches seeded at t=0
  Int initial_founder_cells = 10;      // founders per seeded patch
  Int initial_pool_cells = 0;          // optional extra single-cell packets at t=0

  // Patch-level guards (S4): relative bloom + spatial capacity at patch
  // level only; the absolute ladder applies to the segment mean below.
  Real bloom_factor = 10.0;            // halted when N > factor x baseline sustained
  Real bloom_sustain_s = 3600.0;       // sustained duration before halt
  Real spatial_packing_fraction = 0.6; // packing limit as fraction of close pack
  Real cell_radius_m = 0.5e-6;         // cell radius for the packing bound

  // Absolute density ladder applied to the segment mean only (S4).
  Real segment_dysbiosis_threshold = 0.0;  // cells/mL; 0 disables
  Real segment_guard_window_s = 1800.0;    // trailing window for the mean check
  Real invalid_patch_fraction_stop = 0.5;  // run stop when invalid fraction exceeds

  Layer2AuditConfig audit;

  std::string timeseries_file;        // CSV path; empty disables
  Int summary_interval_steps = 60;    // observables cadence
  std::string provenance_file;        // JSON run record; empty disables
  std::string checkpoint_file;        // HDF5 checkpoint path; empty disables
  Int checkpoint_interval_steps = 0;  // 0 = no periodic checkpoint
  bool checkpoint_final = false;      // write checkpoint at run end
};

}  // namespace gutibm

#endif  // GUTIBM_LAYER2_CONFIG_H
