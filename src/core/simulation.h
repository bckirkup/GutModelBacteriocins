/* -----------------------------------------------------------------------
   GutIBM – Main simulation engine
   Orchestrates the biological timestep loop inspired by NUFEB's
   nufeb_run.h, with decoupled timescales:
   
   1. Biology module (bio_dt ~60s):
      - Metabolism (Monod growth, division, death)
      - Bacteriocin (SOS lysis scheduling, toxin release)
      - [QSSA bacteriocin deposition before receptor]
      - Receptor (competitive binding, killing)
      - Conjugation (HGT between neighbors)
      - Mutation (stochastic BI locus changes)
   
   2. Chemistry module (instantaneous via QSSA):
      - Green's function superposition for toxin fields
      - Nutrient depletion zones
      - VBF sink/source coupling
   
   3. Physics module (same dt as biology):
      - Advection (mucus flow)
      - VBF drag
      - Cell–cell mechanical repulsion
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_SIMULATION_H
#define GUTIBM_SIMULATION_H

#include "types.h"
#include "step_profiler.h"
#include "dysbiosis_guard.h"
#include "agent.h"
#include "domain.h"
#include "random.h"
#include "chemical_field.h"
#include "advection.h"
#include "vbf.h"
#include "qssa_solver.h"
#include "lineage_tracker.h"
#include "hdf5_writer.h"
#include "step_events.h"
#include "hdf5_reader.h"
#include "input_parser.h"
#include "immigration.h"
#include "fix.h"
#include "chemical_field_gpu.h"
#include "agent_pool_gpu.h"
#include "dispatch.h"
#include "gpu_kernels.h"
#include "termination.h"

#include <memory>
#include <chrono>
#include <vector>

namespace gutibm {

struct ProgressMetrics {
  double pct = 0.0;
  double rate = 0.0;
  double eta_s = 0.0;
};

ProgressMetrics calculate_progress_metrics(Real sim_time,
                                           Real attempt_sim_time,
                                           Real total_time,
                                           double wall_elapsed_s);

class Simulation {
 public:
  struct PopulationStocks {
    Int bacteriostatic_live = 0;
    Int washout_trapped_live = 0;
  };

  Simulation() = default;
  ~Simulation() = default;
  Simulation(Simulation&&) = default;
  Simulation& operator=(Simulation&&) = default;
  Simulation(const Simulation&) = delete;
  Simulation& operator=(const Simulation&) = delete;

  // Initialize from config
  void init(const SimulationConfig& cfg);

  // Initialize domain/modules and restore state from an HDF5 snapshot
  void init_from_checkpoint(const SimulationConfig& cfg,
                            const std::string& h5_file,
                            const std::string& step = "");

  // Run the simulation
  int run();

  // Single timestep
  void step(Real dt);

  // Accessors (const and non-const)
  AgentPool&             agents()          { return agents_; }
  const AgentPool&       agents()    const { return agents_; }

  Domain&                domain()          { return domain_; }
  const Domain&          domain()    const { return domain_; }

  ChemicalField&         chemical_field()       { return chem_; }
  const ChemicalField&   chemical_field() const { return chem_; }

  AdvectionField&        advection()       { return advection_; }
  const AdvectionField&  advection() const { return advection_; }

  VBF&                   vbf()             { return vbf_; }
  const VBF&             vbf()       const { return vbf_; }

  QSSASolver&            qssa()            { return qssa_; }
  const QSSASolver&      qssa()      const { return qssa_; }

  LineageTracker&        lineage_tracker()       { return lineage_; }
  const LineageTracker&  lineage_tracker() const { return lineage_; }

  RNG&                   rng()             { return rng_; }

  Real                   time()      const { return clock_.time; }
  Int                    step_count() const { return clock_.step_count; }

  const SimulationConfig& config() const { return cfg_; }

  const StepEvents& step_events() const { return event_ledger_.step_events; }
  StepEvents&       step_events()       { return event_ledger_.step_events; }
  const StepEvents& cumulative_events() const {
    return event_ledger_.cumulative_events;
  }
  StepEvents& cumulative_events() { return event_ledger_.cumulative_events; }
  const StepEvents& summary_events() const {
    return event_ledger_.summary_events;
  }
  StepEvents& summary_events() { return event_ledger_.summary_events; }
  Int event_window_start_step() const {
    return event_ledger_.window_start_step;
  }
  Real event_window_start_time() const {
    return event_ledger_.window_start_time;
  }
  void set_event_window_start(Int step, Real time) {
    event_ledger_.window_start_step = step;
    event_ledger_.window_start_time = time;
  }
  void prepare_step_events_for_summary();
  void prepare_mechanics_stats_for_summary();
  void prepare_population_stocks_for_summary();
  const PopulationStocks& population_stocks() const {
    return event_ledger_.population_stocks;
  }
  const MechanicsStats& mechanics_stats() const {
    return event_ledger_.mechanics_summary;
  }
  MechanicsStats& mechanics_step_stats() {
    return event_ledger_.mechanics_step;
  }
  const MechanicsStats& mechanics_summary_stats() const {
    return event_ledger_.mechanics_summary;
  }
  const MechanicsStats& mechanics_cumulative_stats() const {
    return event_ledger_.mechanics_cumulative;
  }
  void commit_step_events_after_summary(Int step, Real time) {
    event_ledger_.cumulative_events.add(event_ledger_.summary_events);
    event_ledger_.mechanics_cumulative.add(event_ledger_.mechanics_summary);
    event_ledger_.step_events.reset();
    event_ledger_.summary_events.reset();
    event_ledger_.mechanics_step.reset();
    event_ledger_.mechanics_summary.reset();
    event_ledger_.window_start_step = step + 1;
    event_ledger_.window_start_time = time;
  }
  void reset_step_events_after_summary(Int step, Real time) {
    prepare_step_events_for_summary();
    commit_step_events_after_summary(step, time);
  }
  bool provenance_enabled() const {
    return cfg_.hdf5.enabled && cfg_.hdf5.schedule.provenance > 0;
  }
  const std::vector<KillProvenanceEvent>& kill_provenance() const {
    return event_ledger_.kill_provenance;
  }
  std::vector<KillProvenanceEvent>& kill_provenance() {
    return event_ledger_.kill_provenance;
  }
  void record_kill_provenance(const KillProvenanceEvent& event);
  void clear_kill_provenance() { event_ledger_.kill_provenance.clear(); }

  // Spec 1: local oxygen and ROS induction hook (Spec 2)
  Real local_O2(const Agent& agent) const;
  Real ros_induction_rate(const Agent& agent) const;
  Real local_nuclease_toxin(const Agent& agent, Int agent_index) const;

  // Persistent SOS lysis burst sources (protease decay)
  void add_toxin_burst(const ToxinBurstSource& burst);
  void prune_toxin_bursts(Real current_time);
  void materialize_bacteriocin_fields_for_output();
  const std::vector<ToxinBurstSource>& toxin_bursts() const {
    return bacteriocin_.bursts;
  }

  // Active Fix plugin names in execution order
  std::vector<std::string> fix_names() const;

  // MPI global statistics (valid after allreduce)
  Int  global_agent_count() const {
    return mpi_ghost_.stats.global_agent_count;
  }
  Real global_mu_avg()      const { return mpi_ghost_.stats.global_mu_avg; }

  // Adaptive timestep computation
  Real compute_adaptive_dt() const;

  bool gpu_active() const { return gpu_.active; }
  bool halted_for_dysbiosis() const { return dysbiosis_.halted(); }
  TerminationCause termination_cause() const { return termination_cause_; }
  const std::string& termination_detail() const { return termination_detail_; }
  double termination_wall_seconds() const { return termination_wall_seconds_; }
  Real halt_density_cells_per_mL() const {
    return dysbiosis_.halt_density_cells_per_mL();
  }

  const StepProfile& step_profile() const { return step_profile_; }
  void reset_step_profile() { step_profile_.reset(); }
  void print_step_profile() const;

  ChemicalFieldGpu&       chem_gpu()       { return gpu_.chem; }
  const ChemicalFieldGpu& chem_gpu() const { return gpu_.chem; }
  AgentPoolGpu&           agents_gpu()       { return gpu_.agents; }
  const AgentPoolGpu&     agents_gpu() const { return gpu_.agents; }
  void exchange_ghost_agents();


 private:
  // Initialization helpers
  void init_population(const SimulationConfig& cfg);
  Agent create_strain_agent(const SimulationConfig::InitialStrain& strain,
                            Vec3 pos);
  void apply_checkpoint_snapshot(const HDF5CheckpointSnapshot& snap);
  void update_grid_coupling();
  void rebuild_spatial_hash();
  void remove_dead_agents();
  void check_washout();
  void crypt_migration(Real dt);
  void take_lineage_snapshot();
  void maybe_write_restart();
  void write_restart_now();
  void write_hdf5_step(Real dt);
  void emit_heartbeat(const std::chrono::steady_clock::time_point& wall_start,
                      const std::chrono::steady_clock::time_point& wall_now,
                      std::chrono::steady_clock::time_point& next_heartbeat,
                      bool& heartbeat_emitted) const;
  void emit_progress_if_due(
      Real dt, const std::chrono::steady_clock::time_point& wall_start,
      Real attempt_start_sim_time,
      const std::chrono::steady_clock::time_point& wall_now);
  void update_lineage_snapshot_if_due();
  bool population_stop(int rank) const;
  bool dysbiosis_threshold_exceeded(int rank);
  bool closure_violation(std::string& detail);

  // Module execution (NUFEB-inspired)
  void module_biology(Real dt);
  void module_chemistry(Real dt);
  void module_physics(Real dt);
  void update_bacteriocin_fields();

  // MPI domain decomposition
  void migrate_agents();
  void clear_ghost_agents();
  void allreduce_global_stats();

  // MPI global statistics (valid after allreduce)
  struct MpiStats {
    Int global_agent_count = 0;
    Real global_mu_avg = 0.0;
    Real global_max_abs_mu = 0.0;
    Int global_sos_count = 0;
  };

  struct EventLedger {
    StepEvents step_events;
    StepEvents summary_events;
    StepEvents cumulative_events;
    MechanicsStats mechanics_step;
    MechanicsStats mechanics_summary;
    MechanicsStats mechanics_cumulative;
    PopulationStocks population_stocks;
    Int window_start_step = 1;
    Real window_start_time = 0.0;
    std::vector<KillProvenanceEvent> kill_provenance;
  };

  struct GpuState {
    bool active = false;
    ChemicalFieldGpu chem;
    AgentPoolGpu agents;
  };

  struct MpiGhostState {
    MpiStats stats;
    std::vector<Int> ghost_indices;
  };

  struct ImmigrationState {
    RNG rng;
    Int start_step = 0;
  };

  struct BacteriocinState {
    std::vector<ToxinBurstSource> bursts;
    bool fields_current = false;
    bool grid_materialized = false;
  };

  // Timestep clock and output scheduling
  struct Clock {
    Real time = 0.0;
    Int step_count = 0;
    Real next_output = 0.0;
    Real next_snapshot = 0.0;
  };

  MpiGhostState mpi_ghost_;

  // State
  AgentPool       agents_;
  Domain          domain_;
  ChemicalField   chem_;
  AdvectionField  advection_;
  VBF             vbf_;
  QSSASolver      qssa_;
  LineageTracker  lineage_;
  HDF5Writer      hdf5_;
  RNG             rng_;
  ImmigrationEngine immigration_;

  // Fix modules (mutable: compute() updates simulation state via sim_ reference)
  mutable std::vector<std::unique_ptr<Fix>> fixes_;

  // Config
  SimulationConfig cfg_;

  // Timers
  Clock clock_;

  // Step profiling
  StepProfile step_profile_;

  GpuState gpu_;
  DysbiosisGuard dysbiosis_;

  BacteriocinState bacteriocin_;
  EventLedger event_ledger_;
  TerminationCause termination_cause_ = TerminationCause::IncompleteUnknown;
  std::string termination_detail_ = "run has not completed";
  double termination_wall_seconds_ = 0.0;
  Int zero_realization_steps_ = 0;
};

}  // namespace gutibm

#endif  // GUTIBM_SIMULATION_H
