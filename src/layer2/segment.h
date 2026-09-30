/* -----------------------------------------------------------------------
   GutIBM – Spec 13 Phase 1: Layer-2 mucus segment

   One segment couples `n_patches` mucosal patches through a single
   luminal transit pool. Every patch is a lightweight tabulated colony
   whose growth/capacity is interpolated from the CARBON_LADDER_CAMPAIGN
   (#314) ladder; optionally one slot is substituted by a live
   Simulation instance as the audit gate.

   Design contract: docs/SPEC13_IMPLEMENTATION_REVIEW.md
     S1 two loss channels (single-cell washout, contraction fragment)
        feed one luminal pool with distinct reseeding kernels; exactly
        one debit per departing agent.
     S2 a crypt-type patch discounts disruption by type only — no
        agent-level in_crypt double discount exists at Layer 2.
     S5 per-patch status is independent of run TerminationCause.
     G2 contraction = per-step Bernoulli(rate*dt), fired after
        physics/migration and before washout/cleanup.
     G3 all draws use streams seeded by stable ids — results are
        invariant to patch iteration order, rank count, and device.
     G4 pool/colonist serialization reuses agent_transfer format.
     G5 occupancy is type-stratified; segment mean =
        sum_t f_t * occ_t * density_t and is never quoted bare.
     G6 CFU/cm^2 is reported alongside volumetric density.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_LAYER2_SEGMENT_H
#define GUTIBM_LAYER2_SEGMENT_H

#include "agent.h"
#include "layer2_config.h"
#include "random.h"
#include "termination.h"
#include "types.h"

#include <cstdint>
#include <deque>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace gutibm {

struct SimulationConfig;
class Simulation;

// One mucosal patch. `colonists` are real Agent objects so BI-locus
// composition, lineage, and serialization are uniform between the pool
// and patch containers.
struct Patch {
  Int id = 0;
  PatchType type = PatchType::ExposedProximal;
  PatchStatus status = PatchStatus::Active;
  std::vector<Agent> colonists;
  // Cell tags are patch-scoped: (id + 1) << 32 | local serial. This keeps
  // tags (and the pool-packet stream keys derived from them) stable
  // under patch iteration permutations (G3).
  TagID next_cell_tag = 0;
  Int baseline_count = 0;        // post-seed count; bloom guard reference
  Real bloom_excess_since = -1;  // sim time the population first exceeded bound
  bool disrupted = false;        // this step's contraction decision
  bool live = false;             // audit slot: stepped via Simulation
  RNG rng;
};

// A clonal or single-cell packet in luminal transit.
struct LuminalPacket {
  uint64_t id = 0;
  PacketKind kind = PacketKind::Single;
  std::vector<Agent> cells;
  Real entry_time = 0.0;
  RNG rng;
};

// Population ledger for the closure gate: every cell is either in a
// patch, in the pool, exported distally, or never existed. Internal
// patch<->pool transfers are net-zero; closure reads
//   stocks == initial + births - distal_losses.
struct SegmentLedger {
  Int initial_cells = 0;
  Int births = 0;
  Int washout_departures = 0;      // cells patch->pool, single channel
  Int contraction_departures = 0;  // cells patch->pool, fragment channel
  Int reseeds = 0;                 // cells pool->patch
  Int distal_losses = 0;           // cells pool->distal exit

  [[nodiscard]] Int expected_stocks() const {
    return initial_cells + births - distal_losses;
  }
};

struct TypeObservables {
  Int n_patches = 0;
  Int n_occupied = 0;
  Int n_halted = 0;
  Real occupancy = 0.0;
  Real mean_density_occupied = 0.0;  // cells/mL over occupied patches
};

struct SegmentObservables {
  TypeObservables per_type[kPatchTypeCount];
  Real segment_mean_cfu_ml = 0.0;
  Real segment_cfu_cm2 = 0.0;
  Int pool_cells = 0;
  Int pool_packets = 0;
  Int patch_cells = 0;
  Int n_occupied_total = 0;
};

struct AuditRecord {
  Real time = 0.0;
  Int patch_id = -1;
  Int n_live = 0;
  Int n_shadow = 0;
  Real rel_diff = 0.0;
};

// Internal bookkeeping for the live audit slot.
struct LivePatchState {
  std::unique_ptr<Simulation> sim;
  std::vector<Agent> harvested_singles;   // washout exports pending packets
  std::vector<Agent> harvested_fragment;  // contraction exports pending packet
  Real pending_fragment_fraction = 0.0;   // armed before sim.step for the hook
  Patch shadow;                           // tabulated twin for the discrepancy
  Real window_start_time = -1;
  bool shadow_valid = false;
  Int last_live_birth_count = 0;
};

class MucusSegment {
 public:
  // Out-of-line ctor: LivePatchState/table_ hold unique_ptrs to
  // forward-declared types; instantiation must stay in segment.cpp.
  MucusSegment();
  ~MucusSegment();
  // Non-movable: audit hooks registered on live patches capture `this`.
  MucusSegment(MucusSegment&&) = delete;
  MucusSegment& operator=(MucusSegment&&) = delete;
  MucusSegment(const MucusSegment&) = delete;
  MucusSegment& operator=(const MucusSegment&) = delete;

  void init(const SimulationConfig& cfg);
  void init_from_checkpoint(const SimulationConfig& cfg,
                            const std::string& h5_file);
  int run();

  // One segment timestep of length cfg.time.bio_dt (exposed for tests).
  void step(Real dt);

  // Deterministic full-state fingerprint: covers per-patch population,
  // type, status, pool packet kinds/sizes, and ledger counters. Order of
  // `order_` does not enter the hash (keys are stable patch/packet ids).
  [[nodiscard]] uint64_t fingerprint() const;

  // Test hook: permute the per-step patch iteration order. Under G3 the
  // trajectory must be identical for any permutation.
  void permute_patch_order_for_testing(const std::vector<Int>& order);

  [[nodiscard]] const std::vector<Patch>& patches() const { return patches_; }
  [[nodiscard]] const std::deque<LuminalPacket>& pool() const { return pool_; }
  [[nodiscard]] const SegmentLedger& ledger() const { return ledger_; }
  [[nodiscard]] SegmentObservables observables() const;
  [[nodiscard]] TerminationCause termination_cause() const {
    return termination_cause_;
  }
  [[nodiscard]] const std::string& termination_detail() const {
    return termination_detail_;
  }
  [[nodiscard]] Real time() const { return time_; }
  [[nodiscard]] Int step_count() const { return step_count_; }
  [[nodiscard]] const std::vector<AuditRecord>& audit_records() const {
    return audit_records_;
  }
  [[nodiscard]] bool ledger_closed() const;

  // Checkpoint I/O (requires GUTIBM_HDF5; throws ConfigError otherwise).
  void write_checkpoint(const std::string& path) const;

 private:
  // init helpers
  void validate_config() const;
  void assign_patch_types();
  void seed_initial_colonists();
  void make_founder(Patch& patch, Int strain_type, Real mu_max,
                    const std::vector<std::string>& plasmids,
                    bool conjugative);
  void admissibility_check() const;
  void init_audit_slot();
  SimulationConfig live_config_for(const Patch& patch) const;

  // per-step phases
  // `ledger=true` accounts departures into the real pool; the audit
  // shadow passes false so its fictional losses are discarded.
  void grow_tabulated(Patch& patch, Real dt);
  void contract_tabulated(Patch& patch, bool ledger);
  void washout_tabulated(Patch& patch, Real dt, bool ledger);
  void step_live(Patch& patch, Real dt);
  void pool_phase(Real dt);
  void audit_check();
  void shadow_step(Patch& live_patch, Real dt);
  void guard_update(Patch& patch);
  void run_level_guards();
  void emit_timeseries_row(std::ofstream& out) const;
  void write_provenance(const std::string& path) const;
  void write_checkpoint_now() const;

  // helpers
  const Layer2TypeConfig& type_config(PatchType t) const;
  Real patch_volume_mL() const;
  Int patch_count(const Patch& patch) const;
  Int spatial_capacity_cells() const;
  // `stream_key` seeds the packet's RNG; callers pass an order-stable
  // identity (departing cell tag or patch id x step) so results are
  // invariant under patch iteration order (G3).
  void push_packet(PacketKind kind, std::vector<Agent> cells,
                   uint64_t stream_key);
  void terminate(TerminationCause cause, std::string detail);
  static uint64_t mix_seed(uint64_t seed, uint64_t a, uint64_t b);
  static uint64_t mix_seed(uint64_t seed, uint64_t a);
  static uint64_t mix_seed(uint64_t seed, const char* tag);

  const SimulationConfig* cfg_ = nullptr;

  std::vector<Patch> patches_;
  std::vector<Int> order_;          // patch iteration order (G3-invariant)
  std::deque<LuminalPacket> pool_;
  std::unordered_map<Int, LivePatchState> live_;  // by patch id
  std::unique_ptr<class PatchTable> table_;

  SegmentLedger ledger_;
  RNG contraction_stream_;
  RNG init_stream_;

  Real time_ = 0.0;
  Int step_count_ = 0;
  // Pool-origin cell tags live above the patch ranges.
  TagID next_tag_ = 0;
  uint64_t next_packet_id_ = 0;
  std::deque<std::pair<Real, Real>> mean_window_;  // (time, segment mean)
  Real audit_next_time_ = -1;

  TerminationCause termination_cause_ = TerminationCause::IncompleteUnknown;
  std::string termination_detail_;
  std::vector<AuditRecord> audit_records_;
};

}  // namespace gutibm

#endif  // GUTIBM_LAYER2_SEGMENT_H
