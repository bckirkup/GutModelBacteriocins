/* -----------------------------------------------------------------------
   GutIBM – Main simulation engine
   Orchestrates the biological timestep loop inspired by NUFEB's
   nufeb_run.h, with decoupled timescales:
   
   1. Biology module (bio_dt ~60s):
      - Metabolism (Monod growth, division, death)
      - Mutation (stochastic BI locus changes)
      - Bacteriocin (SOS lysis, toxin release)
      - Receptor (competitive binding, killing)
      - Conjugation (HGT between neighbors)
   
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
#include "agent.h"
#include "domain.h"
#include "random.h"
#include "chemical_field.h"
#include "advection.h"
#include "vbf.h"
#include "qssa_solver.h"
#include "lineage_tracker.h"
#include "hdf5_writer.h"
#include "hdf5_reader.h"
#include "input_parser.h"
#include "fix.h"
#include "chemical_field_gpu.h"
#include "agent_pool_gpu.h"
#include "spatial_hash_gpu.h"
#include "dispatch.h"
#include "gpu_kernels.h"

#include <memory>
#include <vector>

namespace gutibm {

class Simulation {
 public:
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
  void run();

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

  Real                   time()      const { return time_; }
  Int                    step_count() const { return step_count_; }

  const SimulationConfig& config() const { return cfg_; }

  // Spec 1: local oxygen and ROS induction hook (Spec 2)
  Real local_O2(const Agent& agent) const;
  Real ros_induction_rate(const Agent& agent) const;
  Real local_nuclease_toxin(const Agent& agent) const;

  // Persistent SOS lysis burst sources (protease decay)
  void add_toxin_burst(const ToxinBurstSource& burst);
  void prune_toxin_bursts(Real current_time);
  const std::vector<ToxinBurstSource>& toxin_bursts() const { return toxin_bursts_; }

  // Active Fix plugin names in execution order
  std::vector<std::string> fix_names() const;

  // MPI global statistics (valid after allreduce)
  Int  global_agent_count() const { return global_agent_count_; }
  Real global_mu_avg()      const { return global_mu_avg_; }

  // Adaptive timestep computation
  Real compute_adaptive_dt() const;

  bool gpu_active() const { return gpu_active_; }

  const StepProfile& step_profile() const { return step_profile_; }
  void reset_step_profile() { step_profile_.reset(); }
  void print_step_profile() const;

  ChemicalFieldGpu&       chem_gpu()       { return chem_gpu_; }
  const ChemicalFieldGpu& chem_gpu() const { return chem_gpu_; }
  AgentPoolGpu&           agents_gpu()       { return agents_gpu_; }
  const AgentPoolGpu&     agents_gpu() const { return agents_gpu_; }


 private:
  // Initialization helpers
  void init_population(const SimulationConfig& cfg);
  void apply_checkpoint_snapshot(const HDF5CheckpointSnapshot& snap);
  void update_grid_coupling();
  void rebuild_spatial_hash();
  void remove_dead_agents();
  void check_washout();
  void crypt_migration(Real dt);
  void take_lineage_snapshot();

  // Module execution (NUFEB-inspired)
  void module_biology(Real dt);
  void module_chemistry(Real dt);
  void module_physics(Real dt);

  // MPI domain decomposition
  void migrate_agents();
  void exchange_ghost_agents();
  void clear_ghost_agents();
  void allreduce_global_stats();

  // Ghost agent tracking
  std::vector<Int> ghost_indices_;  // indices of ghost agents in the pool
  Int global_agent_count_ = 0;      // total agents across all ranks
  Real global_mu_avg_ = 0.0;        // global average growth rate

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

  // Fix modules
  std::vector<std::unique_ptr<Fix>> fixes_;

  // Config
  SimulationConfig cfg_;

  // Timers
  Real time_       = 0.0;
  Int  step_count_ = 0;
  Real next_output_ = 0.0;
  Real next_snapshot_ = 0.0;

  // Step profiling
  StepProfile step_profile_;

  // GPU acceleration state
  bool gpu_active_ = false;
  ChemicalFieldGpu chem_gpu_;
  AgentPoolGpu agents_gpu_;
  SpatialHashGpu spatial_hash_gpu_;

  std::vector<ToxinBurstSource> toxin_bursts_;
};

}  // namespace gutibm

#endif  // GUTIBM_SIMULATION_H
