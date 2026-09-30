/* -----------------------------------------------------------------------
   GutIBM – Spec 13 Phase 3: per-region luminal population compartment

   The lumen is a POPULATION, not a pool (review S8): per-lineage cell
   counts indexed by BI-locus genotype, carried as a Real count plus a
   single exemplar Agent per lineage that preserves the genotype for
   genotype-correct count->agent instantiation on reattachment. Scalar
   chemistry (carbon, acetate) shares the Layer-1 yield/acid parameters
   through the functions in fixes/metabolic_mode.h — implemented once,
   here, used by every region.

   Two-phase update keeps region/lineage processing order irrelevant
   (G3): compute() derives a step plan from the pre-step snapshot only;
   apply() commits it plus the upstream region's outflow. All draws
   are seeded by stable ids (region, step, lineage).
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_LAYER3_LUMEN_H
#define GUTIBM_LAYER3_LUMEN_H

#include "agent.h"
#include "layer3_config.h"
#include "random.h"
#include "types.h"

#include <map>
#include <vector>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

namespace gutibm {

struct SimulationConfig;

// Luminal physiology resolved once from the shared Layer-1 parameters
// (S8: same yield and acid-inhibition values, read off
// SimulationConfig, never duplicated).
struct LumenPhysio {
  Real yield_carbon = 0.5;            // mol carbon per kg biomass
  Real carbon_cost_factor = 4.1;      // anaerobic multiplier on yield
  Real anaerobic_mu_factor = 0.55;
  Real ferm_acid_yield = 1.0;         // mol acetate per mol carbon
  Real acid_inhibition_max = 0.8;
  Real acid_inhibition_ki = 50.0;     // mol/m^3 undissociated acetate
  Real acetate_pka = 4.76;
  Real diet_carbon_mol_m3 = 8.0e-3;   // upstream boundary for region 0
};

// One luminal lineage: a fractional cell count plus the exemplar that
// carries the BI-locus genotype for reattachment instantiation.
struct LuminalLineage {
  Agent exemplar;
  Real count = 0.0;
};

struct LuminalScalars {
  Real carbon_mol_m3 = 0.0;
  Real acetate_mol_m3 = 0.0;
};

// A computed-but-uncommitted step outcome (see class comment).
struct LumenStepPlan {
  // Per-lineage count after growth and outflow, before reattachment
  // decrements. Populated for every extant lineage.
  std::map<TagID, Real> resident;
  // Cells leaving downstream this step, with the genotype exemplar for
  // downstream lineages / stool bookkeeping.
  std::map<TagID, std::pair<Real, Agent>> outflow;
  // Reattachment attempts drawn this step (sorted by lineage id — the
  // processing order is therefore permutation-invariant by
  // construction).
  std::vector<std::pair<TagID, Int>> reattach_events;
  LuminalScalars scalars;
  Real births = 0.0;          // new cells produced this step
  Real outflow_cells = 0.0;   // total leaving downstream
  Real mean_mu_per_s = 0.0;   // count-weighted realized growth rate
};

class LumenCompartment {
 public:
  void init(const Layer3RegionConfig& rc,
            const std::vector<Agent>& seed_exemplars, Real seed_cells);

  // Mucosal departures routed in by the segment's pool-exit sink.
  void add_cells(const std::vector<Agent>& cells);

  // Phase 1: compute the step from the pre-step snapshot. `occupancy`
  // is the region mucosal occupancy fraction for p_attach.
  [[nodiscard]] LumenStepPlan compute(
      Real dt, const Layer3RegionConfig& rc, const Layer3Config& global,
      const LumenPhysio& physio, const LuminalScalars& upstream_scalars,
      Real occupancy_fraction, Int step_index, int region_index,
      uint64_t seed) const;

  // Phase 2: commit `plan`, subtract `deposited` (cells the chain
  // actually placed into mucosa), and merge the upstream plan's
  // outflow at `upstream_survival`.
  void apply(const LumenStepPlan& plan, const LumenStepPlan* upstream,
             Real upstream_survival,
             const std::map<TagID, Int>& deposited);

  // One-time flush perturbation (recovery-kinematics probe): removes
  // `fraction` of every lineage count; returns cells removed.
  Real purge(Real fraction);

  [[nodiscard]] Real total_cells() const;
  [[nodiscard]] const std::map<TagID, LuminalLineage>& lineages() const {
    return lineages_;
  }
  [[nodiscard]] const LuminalScalars& scalars() const { return scalars_; }

  [[nodiscard]] uint64_t fingerprint() const;

#ifdef GUTIBM_HDF5
  void write_state(hid_t fid, const std::string& prefix) const;
  void read_state(hid_t fid, const std::string& prefix);
#endif

 private:
  std::map<TagID, LuminalLineage> lineages_;  // sorted-key => stable order
  LuminalScalars scalars_;
};

}  // namespace gutibm

#endif  // GUTIBM_LAYER3_LUMEN_H
