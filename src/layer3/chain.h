/* -----------------------------------------------------------------------
   GutIBM – Spec 13 Phase 3: three-region colonic chain driver

   Region r owns one MucusSegment (within-region mechanics are the
   merged Phase-1 Layer 2) and one LumenCompartment (S8 population).
   Coupling is unidirectional proximal->distal:
     mucosa -(3 channels)-> region pool -> pool exit -> region lumen
     lumen -(outflow @ transit)-> downstream lumen -> ... -> stool
     lumen -(reattach p0*(1-occ))-> region mucosa (try_deposit)
     stool = region-2 lumen outflow x outflow_survival

   The chain ledger closes across the mucosa/lumen boundary and across
   region boundaries: a mucosal departure is a transfer, luminal
   division and washout-to-stool are new ledger events.
     stocks == initial + births_seg + births_lumen - stool - deaths - purged
   where births_seg = sum over segments' ledger.births and every
   internal transfer (external_arrivals, distal exits, deposits,
   upstream arrivals) cancels pairwise.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_LAYER3_CHAIN_H
#define GUTIBM_LAYER3_CHAIN_H

#include "agent.h"
#include "layer3_config.h"
#include "lumen.h"
#include "random.h"
#include "segment.h"
#include "termination.h"
#include "types.h"

#include <deque>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace gutibm {

struct SimulationConfig;

// Chain-level population ledger. Mucosal stocks are integer cells;
// luminal stocks are Real counts, so the ledger is Real throughout and
// closure uses a small relative tolerance (fixed-point drift only).
struct ChainLedger {
  Real initial_cells = 0.0;
  Real mucosal_births = 0.0;   // sum of segment ledger.births
  Real luminal_births = 0.0;   // deterministic luminal divisions
  Real stool_exports = 0.0;    // cells leaving region 2 to stool
  Real luminal_deaths = 0.0;   // transfer-survival losses between regions
  Real purged = 0.0;           // one-time luminal flush (recovery probe)

  [[nodiscard]] Real expected_stocks() const {
    return initial_cells + mucosal_births + luminal_births - stool_exports -
           luminal_deaths - purged;
  }
};

// Per-step observables for one region's luminal compartment.
struct LumenObservables {
  Real total_cells = 0.0;
  Int n_lineages = 0;
  Real mean_mu_per_s = 0.0;
  LuminalScalars scalars;
};

// One HAPC event drawn for the step; `regions` are the indices the
// event's propagation path covers this step.
struct HapcEvent {
  bool fired = false;
  bool antegrade = true;
  Int origin = 0;
  std::vector<bool> regions{false, false, false};
};

class ColonicChain {
 public:
  ColonicChain();
  ~ColonicChain();
  ColonicChain(ColonicChain&&) = delete;
  ColonicChain& operator=(ColonicChain&&) = delete;
  ColonicChain(const ColonicChain&) = delete;
  ColonicChain& operator=(const ColonicChain&) = delete;

  void init(const SimulationConfig& cfg);
  void init_from_checkpoint(const SimulationConfig& cfg,
                            const std::string& h5_file);
  int run();

  // One chain timestep of length cfg.time.bio_dt (exposed for tests).
  void step(Real dt);

  // Deterministic full-state fingerprint over segment fingerprints,
  // luminal counts/exemplars, scalars, and ledger counters.
  [[nodiscard]] uint64_t fingerprint() const;

  // Test hook: permute region processing order for a step. Under G3
  // the trajectory must be identical for any permutation.
  void permute_region_order_for_testing(const std::vector<Int>& order);

  struct RegionView {
    const MucusSegment& segment;
    const LumenCompartment& lumen;
    const Layer3RegionConfig& config;
  };
  [[nodiscard]] RegionView region(Int i) const {
    return {*regions_[static_cast<size_t>(i)].segment,
            regions_[static_cast<size_t>(i)].lumen,
            region_cfg(static_cast<RegionKind>(i))};
  }
  [[nodiscard]] const Layer3RegionConfig& region_cfg(RegionKind r) const;
  [[nodiscard]] const ChainLedger& ledger() const { return ledger_; }
  [[nodiscard]] TerminationCause termination_cause() const {
    return termination_cause_;
  }
  [[nodiscard]] const std::string& termination_detail() const {
    return termination_detail_;
  }
  [[nodiscard]] Real time() const { return time_; }
  [[nodiscard]] Int step_count() const { return step_count_; }
  [[nodiscard]] Real stocks() const;
  [[nodiscard]] bool ledger_closed() const;
  [[nodiscard]] LumenObservables lumen_observables(Int i) const;
  [[nodiscard]] Real stool_cfu_per_g() const;

  // Checkpoint I/O (requires GUTIBM_HDF5).
  void write_checkpoint(const std::string& path) const;

 private:
  struct RegionState {
    std::unique_ptr<MucusSegment> segment;
    LumenCompartment lumen;
    TagID next_deposit_tag = 0;  // prefix | (0xFF<<48) | serial
  };

  void validate_config() const;
  // Shared init prelude for init() and init_from_checkpoint(): resolves
  // the uniform-aware region rows, the shared luminal physiology, and the
  // empty regions/order containers.
  void resolve_common_config(const SimulationConfig& cfg);
  // Emplaces and derives one region's SimulationConfig copy.
  SimulationConfig& emplace_region_cfg(const SimulationConfig& cfg, Int r);
  // Attaches the pool-exit sink and the region contraction override.
  void wire_region(Int r);
  [[nodiscard]] std::vector<Agent> seed_lumen_exemplars() const;
  void process_reattach(Int r, const LumenStepPlan& plan,
                        std::map<TagID, Int>& deposited);
  [[nodiscard]] HapcEvent draw_hapc(Real dt, Int step_index) const;
  [[nodiscard]] Real hapc_rate_per_s(Real time_s) const;
  void emit_timeseries_row(std::ofstream& out) const;
  void write_provenance(const std::string& path) const;
  void terminate(TerminationCause cause, std::string detail);

  static uint64_t mix_seed(uint64_t seed, uint64_t a, uint64_t b);
  static uint64_t mix_seed(uint64_t seed, uint64_t a);
  static uint64_t mix_seed(uint64_t seed, const char* tag);

  const SimulationConfig* cfg_ = nullptr;
  LumenPhysio physio_;
  std::vector<Layer3RegionConfig> regions_cfg_;  // resolved (uniform-aware)
  std::vector<RegionState> regions_;             // kRegionCount
  // Per-region SimulationConfig copies backing each segment's cfg_
  // pointer; deque so addresses are stable across the init loop.
  std::deque<SimulationConfig> region_cfgs_;
  std::vector<Int> order_;                     // region processing order
  HapcEvent last_hapc_;

  ChainLedger ledger_;
  Real last_mean_mu_[kRegionCount] = {0.0, 0.0, 0.0};
  bool purged_ = false;

  Real time_ = 0.0;
  Int step_count_ = 0;
  TerminationCause termination_cause_ = TerminationCause::IncompleteUnknown;
  std::string termination_detail_;
};

}  // namespace gutibm

#endif  // GUTIBM_LAYER3_CHAIN_H
