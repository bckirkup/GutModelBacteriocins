/* -----------------------------------------------------------------------
   GutIBM – Simulation engine implementation
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "stop_signal.h"
#include "species_names.h"
#include "agent_transfer.h"
#include "fix_registry.h"
#include "fix_motility.h"
#include "plasmid.h"
#include "dispatch.h"
#include "chemistry_pipeline.h"
#include "gpu_profile.h"
#ifdef GUTIBM_CUDA
#include "device.h"
#include "gpu_kernels.h"
#include "device_memory.h"
#include <cuda_runtime.h>
#endif

#include <format>
#include <iostream>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <numeric>
#include <iomanip>
#include <set>
#include "error.h"
#include <utility>
#ifdef GUTIBM_OPENMP
#include <omp.h>
#endif

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {

namespace {

constexpr uint64_t kImmigrationSeedMix = 0x9e3779b97f4a7c15ULL;

void reject_unsupported_slab_surfaces(
    const SimulationConfig& cfg, const Domain& domain) {
  if (cfg.chemistry_decomposition != "slab") return;
  (void)domain;
}

Real global_density_cells_per_mL(const Domain& domain, Int global_agents) {
  const Vec3 size = domain.size();
  const Real volume_m3 = size[0] * size[1] * size[2];
  constexpr Real kMillilitersPerCubicMeter = 1.0e6;
  const Real volume_mL = volume_m3 * kMillilitersPerCubicMeter;
  return volume_mL > 0.0 ? static_cast<Real>(global_agents) / volume_mL : 0.0;
}

void assign_plasmids(Agent& agent,
                     const std::vector<std::string>& plasmids,
                     const SimulationConfig& cfg,
                     int rank) {
  for (const auto& pname : plasmids) {
    const PlasmidEntry* entry = PlasmidLibrary::find(pname);
    if (!entry) {
      if (rank == 0) {
        std::cerr << "Warning: unknown plasmid '" << pname
                  << "' — agent spawned without BI locus\n";
      }
      continue;
    }
    BICluster cluster = entry->cluster;
    cluster.retardation = retardation_from_pI(
        cluster.pI, cfg.fixes.bacteriocin.mucin_charge);
    if (const auto it = cfg.plasmid_overrides.find(entry->name);
        it != cfg.plasmid_overrides.end()) {
      const auto& values = it->second;
      if (values.retardation.has_value()) {
        cluster.retardation = *values.retardation;
      }
      if (values.diff_coeff.has_value()) {
        cluster.diff_coeff = *values.diff_coeff;
      }
      if (values.burst_size.has_value()) {
        cluster.burst_size = *values.burst_size;
      }
    }
    agent.genome.bi_loci.push_back(cluster);
    if (entry->conjugative) {
      agent.genome.has_conjugative_plasmid = true;
    }
  }
}

void tag_crypt_resident(Agent& agent, const AdvectionField& advection) {
  if (advection.in_crypt_zone(agent.x[2])) {
    agent.flags.in_crypt = true;
  }
}

// Prefer strain mu_max from config when older checkpoints omit /mu_max.
// Never substitute mu_realized — that permanently caps growth after a stressed
// snapshot and triggers a one-step combinatorial washout wipe on resume.
Real fallback_mu_max(const SimulationConfig& cfg, Int type) {
  for (const auto& strain : cfg.initial_strains) {
    if (strain.type == type && strain.mu_max > 0.0) {
      return strain.mu_max;
    }
  }
  return 5.0e-4;
}

std::vector<size_t> build_bi_offsets(const std::vector<Int>& num_bi_loci) {
  const size_t n = num_bi_loci.size();
  std::vector<size_t> offsets(n + 1, 0);
  for (size_t i = 0; i < n; ++i) {
    offsets[i + 1] = offsets[i] + static_cast<size_t>(num_bi_loci[i]);
  }
  return offsets;
}

void validate_checkpoint_genome(const HDF5CheckpointSnapshot& snap,
                              const std::vector<size_t>& bi_offsets) {
  if (!snap.genome.present) return;
  const size_t expected_bi = bi_offsets.back();
  if (snap.genome.bi_toxin_id.size() != expected_bi) {
    throw SimulationError("checkpoint BI locus count mismatch in genome group");
  }
}

void validate_lumped_toxin_species(const SimulationConfig& cfg,
                                   const ChemicalField& chem) {
  if (cfg.qssa.toxin_lumping != "lumped") return;
  if (chem.find(species::BACTERIOCIN_LUMPED) < 0) {
    throw ConfigError(
        "toxin_lumping=lumped requires the bacteriocin_lumped species");
  }
}

bool fix_enabled(const SimulationConfig& cfg, std::string_view name) {
  if (std::find(cfg.disabled_fixes.begin(),
                cfg.disabled_fixes.end(), name)
      != cfg.disabled_fixes.end()) {
    return false;
  }
  return cfg.enabled_fixes.empty()
      || std::find(cfg.enabled_fixes.begin(), cfg.enabled_fixes.end(), name)
          != cfg.enabled_fixes.end();
}

void require_species(const ChemicalField& chem,
                     std::string_view mechanism,
                     const char* name,
                     std::string_view enabling_key) {
  if (chem.find(name) >= 0) return;
  throw ConfigError(
      "mechanism '" + std::string(mechanism) + "' requires species '" + name
      + "'; enable " + std::string(enabling_key)
      + " or repair chemistry.species_subset");
}

void validate_required_species(const SimulationConfig& cfg,
                               const ChemicalField& chem) {
  if (fix_enabled(cfg, "metabolism")) {
    require_species(chem, "fix_metabolism", species::CARBON,
                    "fixes.metabolism");
    if (cfg.fixes.metabolism.iron_uptake_enabled) {
      require_species(chem, "fix_metabolism", species::IRON,
                      "fixes.metabolism.iron_uptake_enabled");
    }
    if (cfg.fixes.metabolism.b12_uptake_enabled) {
      require_species(chem, "fix_metabolism", species::B12,
                      "fixes.metabolism.b12_uptake_enabled");
    }
    if (cfg.chem_env.oxygen.enabled) {
      require_species(chem, "fix_metabolism", species::OXYGEN,
                      "oxygen.enabled");
    }
    if (cfg.chem_env.acetate.enabled) {
      require_species(chem, "fix_metabolism", species::ACETATE,
                      "acetate.enabled");
    }
    if (cfg.fixes.metabolism.eut_enabled) {
      require_species(chem, "fix_metabolism", species::ETHANOLAMINE,
                      "fixes.metabolism.eut_enabled");
    }
    if (cfg.chem_env.siderophore.enabled) {
      require_species(chem, "fix_metabolism", species::SIDEROPHORE,
                      "siderophore.enabled");
      require_species(chem, "fix_metabolism",
                      species::FERRIC_ENTEROBACTIN, "siderophore.enabled");
      require_species(chem, "fix_metabolism", species::IRON,
                      "siderophore.enabled");
    }
  }

  if (fix_enabled(cfg, "receptor")) {
    require_species(chem, "fix_receptor", species::B12, "fixes.receptor");
    require_species(chem, "fix_receptor",
                    species::FERRIC_ENTEROBACTIN, "fixes.receptor");
    if (cfg.chem_env.ferrichrome.enabled) {
      require_species(chem, "fix_receptor", species::FERRICHROME,
                      "ferrichrome.enabled");
    }
    for (ReceptorType target : {ReceptorType::BtuB, ReceptorType::FepA,
                                ReceptorType::CirA, ReceptorType::FhuA}) {
      require_species(chem, "fix_receptor",
                      species::bacteriocin_species_for(
                          target, cfg.qssa.toxin_lumping == "lumped"),
                      "fixes.receptor");
    }
  }

  if (fix_enabled(cfg, "bacteriocin")) {
    for (ReceptorType target : {ReceptorType::BtuB, ReceptorType::FepA,
                                ReceptorType::CirA, ReceptorType::FhuA}) {
      require_species(chem, "fix_bacteriocin",
                      species::bacteriocin_species_for(
                          target, cfg.qssa.toxin_lumping == "lumped"),
                      "fixes.bacteriocin");
    }
  }

  if (fix_enabled(cfg, "motility") && cfg.cell_bio.motility.enabled) {
    require_species(chem, "fix_motility", species::CARBON,
                    "motility.enabled");
    if (cfg.chem_env.oxygen.enabled
        && cfg.cell_bio.motility.aerotaxis_enabled) {
      require_species(chem, "fix_motility", species::OXYGEN,
                      "oxygen.enabled + motility.aerotaxis_enabled");
    }
    if (cfg.cell_bio.motility.mucin_drag_enabled) {
      require_species(chem, "fix_motility", species::MUCIN,
                      "mucin.enabled + motility.mucin_drag_enabled");
    }
    if (cfg.quorum_sensing.enabled
        && cfg.quorum_sensing.ai2_chemotaxis_enabled) {
      require_species(chem, "fix_motility", species::AI2,
                      "quorum_sensing.ai2_chemotaxis_enabled");
    }
  }

  if (fix_enabled(cfg, "quorum_sensing") && cfg.quorum_sensing.enabled) {
    require_species(chem, "fix_quorum_sensing", species::AI2,
                    "quorum_sensing.enabled");
  }

  const bool vbf_carbon_enabled =
      cfg.vbf.mucin_liberation > 0.0
      || cfg.vbf.carbon_sink_vmax > 0.0
      || cfg.vbf.use_dynamic_mucin;
  if (vbf_carbon_enabled) {
    require_species(chem, "VBF", species::CARBON,
                    "vbf.mucin_liberation or vbf.carbon_sink_vmax");
  }
  if (cfg.vbf.nutrient_sink > 0.0) {
    require_species(chem, "VBF", species::IRON, "vbf.nutrient_sink");
  }
  if (cfg.chem_env.oxygen.enabled && cfg.chem_env.oxygen.vbf_sink > 0.0) {
    require_species(chem, "VBF", species::OXYGEN, "oxygen.vbf_sink");
  }
  if (cfg.chem_env.acetate.enabled
      && (cfg.chem_env.acetate.vbf_production > 0.0
          || cfg.chem_env.acetate.vbf_consumption > 0.0
          || cfg.chem_env.acetate.epithelial_uptake > 0.0)) {
    require_species(chem, "VBF", species::ACETATE,
                    "acetate.vbf_production / acetate.vbf_consumption");
  }
  if (cfg.chem_env.mucin.enabled
      && (cfg.chem_env.mucin.secretion_rate > 0.0
          || cfg.vbf.use_dynamic_mucin)) {
    require_species(chem, "VBF", species::MUCIN,
                    "mucin.enabled / vbf.use_dynamic_mucin");
  }
}

void validate_checkpoint_toxin_species(
    const SimulationConfig& cfg,
    const HDF5CheckpointSnapshot& snap) {
  constexpr std::array<const char*, 4> receptor_names = {
      species::BACTERIOCIN_BTUB, species::BACTERIOCIN_FEPA,
      species::BACTERIOCIN_CIRA, species::BACTERIOCIN_FHUA};
  const bool has_lumped =
      snap.grid.species.contains(species::BACTERIOCIN_LUMPED);
  bool has_receptor = true;
  for (const char* name : receptor_names) {
    has_receptor = has_receptor && snap.grid.species.contains(name);
  }
  if (cfg.qssa.toxin_lumping == "lumped" && (!has_lumped || has_receptor)) {
    throw ConfigError(
        "checkpoint bacteriocin species do not match toxin_lumping=lumped");
  }
  if (cfg.qssa.toxin_lumping != "lumped" && (has_lumped || !has_receptor)) {
    throw ConfigError(
        "checkpoint bacteriocin species do not match toxin_lumping=per_receptor");
  }
}

void validate_checkpoint_species_set(
    const SimulationConfig& cfg,
    const HDF5CheckpointSnapshot& snap) {
  std::set<std::string, std::less<>> expected;
  for (const auto& spec : cfg.chemicals) expected.insert(spec.name);
  std::set<std::string, std::less<>> actual;
  for (const auto& [name, values] : snap.grid.species) {
    (void)values;
    actual.insert(name);
  }
  if (expected != actual) {
    throw ConfigError(
        "checkpoint chemical species set does not match current "
        "chemistry.species_subset");
  }
}

void restore_receptor_fields(Agent& agent,
                             const HDF5CheckpointSnapshot& snap,
                             size_t agent_index) {
  const auto& gen = snap.genome;
  Int r = 0;
  for (Real& expr : agent.genome.receptor_expression) {
    const size_t idx = agent_index * NUM_RECEPTORS + static_cast<size_t>(r);
    expr = gen.receptor_expression[idx];
    agent.genome.toxin_affinity[r] = gen.toxin_affinity[idx];
    agent.genome.ligand_affinity[r] = gen.ligand_affinity[idx];
    agent.receptor_expr[r] = gen.receptor_expression[idx];
    agent.receptor_expr_base[r] = gen.receptor_expression[idx];
    ++r;
  }
}

void restore_bi_loci(Agent& agent,
                     const HDF5CheckpointSnapshot& snap,
                     size_t agent_index,
                     const std::vector<size_t>& bi_offsets) {
  const auto& gen = snap.genome;
  const auto& lin = snap.lineage;
  agent.genome.bi_loci.clear();
  agent.genome.bi_loci.reserve(static_cast<size_t>(lin.num_bi_loci[agent_index]));
  for (size_t b = bi_offsets[agent_index]; b < bi_offsets[agent_index + 1]; ++b) {
    BICluster bi;
    bi.toxin_id = static_cast<uint16_t>(gen.bi_toxin_id[b]);
    bi.immunity_id = static_cast<uint16_t>(gen.bi_immunity_id[b]);
    bi.target = static_cast<ReceptorType>(gen.bi_target[b]);
    bi.bclass = static_cast<BacteriocinClass>(gen.bi_bclass[b]);
    bi.pI = gen.bi_pI[b];
    bi.diff_coeff = gen.bi_diff_coeff[b];
    bi.retardation = gen.bi_retardation[b];
    bi.molecular_weight = gen.bi_molecular_weight[b];
    bi.immunity_binding_affinity = gen.bi_immunity_binding_affinity[b];
    agent.genome.bi_loci.push_back(bi);
  }
}

#ifdef GUTIBM_LEGACY_CHECKPOINTS
// Restore genome fields from pre-genome-group checkpoints (no /genome dataset).
// Disabled by default; build with -DGUTIBM_LEGACY_CHECKPOINTS to enable.
void restore_legacy_genome(Agent& agent, const HDF5CheckpointSnapshot& snap, size_t agent_index) {
  Int r = 0;
  for (Real& expr : agent.genome.receptor_expression) {
    expr = agent.receptor_expr[r++];
  }
  agent.genome.bi_loci.clear();
  agent.genome.bi_loci.resize(static_cast<size_t>(snap.lineage.num_bi_loci[agent_index]));
}
#endif  // GUTIBM_LEGACY_CHECKPOINTS

void restore_checkpoint_grid(ChemicalField& chem,
                             const Domain& domain,
                             const HDF5CheckpointSnapshot& snap) {
  for (const auto& [name, values] : snap.grid.species) {
    Int spec = chem.find(name);
    if (spec < 0) {
      if (domain.rank() == 0) {
        std::cerr << "Warning: checkpoint grid species '" << name
                  << "' not in simulation config — skipped\n";
      }
      continue;
    }
    if (static_cast<Int>(values.size()) != chem.global_ncells()) {
      throw SimulationError("checkpoint grid size mismatch for species: " + name);
    }
    if (chem.slab_mode()) {
      for (Int iz = 0; iz < domain.nz(); ++iz) {
        for (Int iy = 0; iy < domain.ny(); ++iy) {
          for (Int ix = chem.owned_global_x_begin();
               ix < chem.owned_global_x_end(); ++ix) {
            const Int global_cell = domain.cell_index(ix, iy, iz);
            chem.conc_global(spec, global_cell) =
                values[static_cast<size_t>(global_cell)];
          }
        }
      }
    } else {
      Int c = 0;
      for (const Real val : values) {
        chem.conc(spec, c++) = val;
      }
    }
  }
  if (chem.slab_mode()) {
    chem.exchange_concentration_halos();
  }
}

void schedule_output_from_time(Real time, Real interval, Real& next_output, Real& next_snapshot) {
  if (time <= 0.0) {
    next_output = 0.0;
    next_snapshot = 0.0;
    return;
  }
  next_output = (std::floor(time / interval) + 1.0) * interval;
  next_snapshot = next_output;
}

bool try_exit_crypt(Agent& agent, Real dt, Real crypt_z, Real epsilon,
                    Real exit_rate, RNG& rng, Int& crypt_pop) {
  if (!agent.flags.in_crypt) return false;
  if (Real p_exit = 1.0 - std::exp(-exit_rate * dt); !rng.bernoulli(p_exit)) return false;
  agent.x[2] = crypt_z + epsilon;
  agent.flags.in_crypt = false;
  --crypt_pop;
  return true;
}

struct CryptEntryParams {
  Real crypt_z;
  Real crypt_depth;
  Real lo_z;
  Real entry_rate;
  Int carrying_capacity;
};

bool try_enter_crypt(Agent& agent, Real dt, const CryptEntryParams& params,
                     RNG& rng, Int& crypt_pop) {
  if (agent.flags.in_crypt) return false;
  if (agent.x[2] >= params.crypt_z + params.crypt_depth) return false;
  if (crypt_pop >= params.carrying_capacity) return false;
  if (Real p_entry = 1.0 - std::exp(-params.entry_rate * dt); !rng.bernoulli(p_entry)) return false;
  agent.x[2] = rng.uniform(params.lo_z, params.crypt_z);
  agent.flags.in_crypt = true;
  ++crypt_pop;
  return true;
}

enum class MigrateSide { None, Lo, Hi };

MigrateSide classify_migration(const Agent& agent, Int my_rank, Int axis,
                               const Domain& domain) {
  using enum MigrateSide;
  if (const Int dest = domain.owner_rank(agent.x); dest == my_rank) {
    return None;
  } else if (dest == domain.rank_lo()) {
    return Lo;
  } else if (dest == domain.rank_hi()) {
    return Hi;
  }
  return (agent.x[axis] < domain.local_lo_x()) ? Lo : Hi;
}

std::string gpu_fallback_reason(const GpuConfig& gpu) {
  if (!gpu.enabled) return {};
#ifndef GUTIBM_CUDA
  return "binary built without CUDA (cmake .. -DGUTIBM_USE_CUDA=ON && make gut_ibm)";
#else
  if (DeviceContext::device_count() <= 0) {
    const std::string err = DeviceContext::last_error();
    return err.empty() ? "no CUDA devices visible (check nvidia-smi)" : err;
  }
  const std::string err = DeviceContext::last_error();
  return err.empty() ? "cudaSetDevice failed" : err;
#endif
}

void print_gpu_status_banner(bool gpu_active, const GpuConfig& gpu) {
  const char* gpu_enabled_str = gpu.enabled ? "true" : "false";
  const std::string gpu_status = gpu_active
      ? std::format(" (device {})", gpu_device().device_id())
      : std::format(" (gpu_enabled={}, device_id={})",
                    gpu_enabled_str,
                    gpu.device_id);
  std::cout << "  GPU: " << (gpu_active ? "ON" : "OFF") << gpu_status << "\n";
  if (!gpu_active && gpu.enabled) {
    std::cerr << "  GPU requested (gpu_enabled) but inactive: "
              << gpu_fallback_reason(gpu) << "\n";
  }
}

}  // namespace

void Simulation::init(const SimulationConfig& cfg) {
  cfg_ = cfg;
  event_ledger_.step_events.reset();
  event_ledger_.summary_events.reset();
  event_ledger_.cumulative_events.reset();
  event_ledger_.mechanics_step.reset();
  event_ledger_.mechanics_summary.reset();
  event_ledger_.mechanics_cumulative.reset();
  event_ledger_.window_start_step = 1;
  event_ledger_.window_start_time = 0.0;
  InputParser::finalize_config(cfg_);
  if (cfg_.gpu.enabled && cfg_.species_subset != "full") {
    throw ConfigError(
        "gpu_enabled=true is not supported with "
        "chemistry.species_subset != full");
  }
  immigration_.validate(cfg_.immigration,
                        static_cast<Int>(cfg_.initial_strains.size()));
  rng_.seed(cfg_.seed);
  immigration_.seed(cfg_.seed ^ kImmigrationSeedMix);
  immigration_.set_start_step(0);

  // Domain
  domain_.init(cfg.domain);
  reject_unsupported_slab_surfaces(cfg_, domain_);

  // Chemical fields
  chem_.init(domain_, cfg_.chemicals, cfg_.chemistry_decomposition);
  validate_lumped_toxin_species(cfg_, chem_);
  validate_required_species(cfg_, chem_);

  // Advection
  advection_.init(cfg.advection, domain_);

  // VBF
  vbf_.init(cfg.vbf, domain_);

  // QSSA solver
  qssa_.init(cfg.qssa, domain_, advection_);

  // Lineage tracker
  lineage_.init(cfg.time.output_interval);

  // HDF5 output
  hdf5_.init(cfg.hdf5, domain_);
  hdf5_.write_run_provenance(*this);

  // GPU acceleration
  gpu_set_config(cfg.gpu);
  gpu_init_for_rank(domain_.rank(), domain_.nprocs());
  gpu_.active = gpu_runtime_enabled();
  if (gpu_.active && qssa_.agent_sampling()) {
    throw ConfigError(
        "chemistry.toxin_evaluation=agents is not supported with GPU execution");
  }
  if (gpu_.active) {
    gpu_.chem.init(chem_);
    gpu_.agents.sync_from_host(agents_);
  }

  // Register biological fixes via plugin registry
  fixes_ = FixRegistry::create_all(*this, cfg_);

  // Initialize fixes
  for (const auto& fix : fixes_) {
    fix->init();
  }

  // Create initial population (only agents local to this rank)
  init_population(cfg);

  // Initial coupling
  rebuild_spatial_hash();
  update_grid_coupling();
  // Slab concentration halos are valid for all biology field reads from this
  // point until chemistry modifies owned concentrations below.
  chem_.exchange_concentration_halos();

  // Initial global stats
  allreduce_global_stats();

  // Timers
  clock_.time          = 0.0;
  clock_.step_count    = 0;
  clock_.next_output   = 0.0;
  clock_.next_snapshot = 0.0;
  dysbiosis_.configure(cfg_.dysbiosis_threshold,
                       cfg_.dysbiosis_sampling_interval,
                       cfg_.dysbiosis_sample_count);
  dysbiosis_.reset(clock_.time);

  int rank = domain_.rank();
  if (rank == 0) {
    if (cfg.dysbiosis_threshold > 0.0 &&
        cfg.dysbiosis_sampling_interval > 0.0 &&
        cfg.dysbiosis_sampling_interval < cfg.time.bio_dt) {
      std::cerr << "Warning: dysbiosis sampling interval ("
                << cfg.dysbiosis_sampling_interval
                << " s) is shorter than bio_dt (" << cfg.time.bio_dt
                << " s); at most one density sample will be taken per step.\n";
    }
    std::cout << "GutIBM initialized:\n"
              << "  Domain: " << domain_.nx() << "x" << domain_.ny()
              << "x" << domain_.nz() << " cells"
              << " (" << domain_.ncells() << " total)\n"
              << "  MPI ranks: " << domain_.nprocs() << "\n"
              << "  Slab [" << domain_.local_lo_x() << ", "
              << domain_.local_hi_x() << ") m\n"
              << "  Local agents: " << agents_.size()
              << "  Global agents: " << mpi_ghost_.stats.global_agent_count << "\n"
              << "  Global density: "
              << global_density_cells_per_mL(domain_, mpi_ghost_.stats.global_agent_count)
              << " cells/mL\n"
              << "  Chemical species: " << chem_.num_species() << "\n"
              << "  Chemistry species subset: " << cfg_.species_subset << "\n";
    if (!cfg_.disabled_mechanisms.empty()) {
      std::cout << "  Disabled mechanisms:";
      for (const auto& mechanism : cfg_.disabled_mechanisms) {
        std::cout << " " << mechanism;
      }
      std::cout << "\n";
    }
    const std::string adaptive_dt_status = cfg.adaptive_dt.enabled
        ? std::format(" [{}s, {}s]", cfg.adaptive_dt.min, cfg.adaptive_dt.max)
        : "";
    std::cout << "  Bio dt: " << cfg.time.bio_dt << " s\n"
              << "  Adaptive dt: " << (cfg.adaptive_dt.enabled ? "ON" : "OFF")
              << adaptive_dt_status << "\n"
              << "  Total time: " << cfg.time.total_time << " s\n";
    print_gpu_status_banner(gpu_.active, cfg.gpu);
    std::cout << std::flush;
  }
}

void Simulation::prepare_step_events_for_summary() {
  event_ledger_.summary_events = event_ledger_.step_events;
#ifdef GUTIBM_MPI
  if (domain_.nprocs() > 1) {
    std::array<Int, 11> local = {
        event_ledger_.step_events.sos_inductions,
        event_ledger_.step_events.phage_inductions,
        event_ledger_.step_events.mortality_colicin,
        event_ledger_.step_events.mortality_cdi,
        event_ledger_.step_events.outflow_washout,
        event_ledger_.step_events.outflow_boundary,
        event_ledger_.step_events.mortality_lysis,
        event_ledger_.step_events.divisions,
        event_ledger_.step_events.conjugation_transfers,
        event_ledger_.step_events.mutations,
        event_ledger_.step_events.immigrations};
    std::array<Int, 11> global{};
    MPI_Allreduce(local.data(), global.data(), static_cast<int>(local.size()),
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    event_ledger_.summary_events = {
        global[0], global[1], global[2], global[3], global[4], global[5],
        global[6], global[7], global[8], global[9], global[10]};
  }
#endif
}

void Simulation::prepare_mechanics_stats_for_summary() {
  event_ledger_.mechanics_summary = event_ledger_.mechanics_step;
#ifdef GUTIBM_MPI
  if (domain_.nprocs() > 1) {
    Int global = 0;
    MPI_Allreduce(&event_ledger_.mechanics_step.displacement_clamps,
                  &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    event_ledger_.mechanics_summary.displacement_clamps = global;
  }
#endif
}

void Simulation::prepare_population_stocks_for_summary() {
  PopulationStocks local;
  for (const Agent& agent : agents_) {
    if (agent.state == PhenoState::DEAD || agent.flags.is_ghost) continue;
    if (agent.mu_realized
        < cfg_.fixes.metabolism.bacteriostasis_threshold) {
      ++local.bacteriostatic_live;
    }
    if (agent.mu_realized < advection_.washout_rate(agent.x[2])) {
      ++local.washout_trapped_live;
    }
  }

#ifdef GUTIBM_MPI
  if (domain_.nprocs() > 1) {
    const std::array<Int, 2> local_values = {
        local.bacteriostatic_live, local.washout_trapped_live};
    std::array<Int, 2> global_values{};
    MPI_Allreduce(local_values.data(), global_values.data(),
                  static_cast<int>(local_values.size()), MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    event_ledger_.population_stocks = {
        global_values[0], global_values[1]};
    return;
  }
#endif
  event_ledger_.population_stocks = local;
}

void Simulation::init_population(const SimulationConfig& cfg) {
  agents_.configure_tags(AgentPool::first_tag_for_rank(domain_.rank(), domain_.nprocs()),
                           AgentPool::tag_stride(domain_.nprocs()));

  const Real z_min = cfg.initial_population.placement == "z_slab"
      ? cfg.initial_population.z_min : domain_.lo()[2];
  const Real z_max = cfg.initial_population.placement == "z_slab"
      ? cfg.initial_population.z_max : domain_.hi()[2] * 0.5;

  for (const auto& strain : cfg.initial_strains) {
    for (Int i = 0; i < strain.count; ++i) {
      Vec3 pos = {
        rng_.uniform(domain_.lo()[0], domain_.hi()[0]),
        rng_.uniform(domain_.lo()[1], domain_.hi()[1]),
        rng_.uniform(z_min, z_max)
      };

      // Only keep agents that belong to this rank's slab
      if (!domain_.is_local(pos)) continue;

      Agent a = create_strain_agent(strain, pos);
      agents_.push_back(std::move(a));
    }
  }
}

void Simulation::init_from_checkpoint(const SimulationConfig& cfg,
                                      const std::string& h5_file,
                                      const std::string& step) {
  cfg_ = cfg;
  event_ledger_.step_events.reset();
  event_ledger_.summary_events.reset();
  event_ledger_.cumulative_events.reset();
  event_ledger_.mechanics_step.reset();
  event_ledger_.mechanics_summary.reset();
  event_ledger_.mechanics_cumulative.reset();
  event_ledger_.window_start_step = 1;
  event_ledger_.window_start_time = 0.0;
  InputParser::finalize_config(cfg_);
  if (cfg_.gpu.enabled && cfg_.species_subset != "full") {
    throw ConfigError(
        "gpu_enabled=true is not supported with "
        "chemistry.species_subset != full");
  }
  immigration_.validate(cfg_.immigration,
                        static_cast<Int>(cfg_.initial_strains.size()));
  rng_.seed(cfg_.seed);
  immigration_.seed(cfg_.seed ^ kImmigrationSeedMix);

  domain_.init(cfg_.domain);
  reject_unsupported_slab_surfaces(cfg_, domain_);
  chem_.init(domain_, cfg_.chemicals, cfg_.chemistry_decomposition);
  validate_lumped_toxin_species(cfg_, chem_);
  validate_required_species(cfg_, chem_);
  advection_.init(cfg.advection, domain_);
  vbf_.init(cfg.vbf, domain_);
  qssa_.init(cfg.qssa, domain_, advection_);
  lineage_.init(cfg.time.output_interval);
  hdf5_.init(cfg.hdf5, domain_);
  hdf5_.write_run_provenance(*this);
  fixes_ = FixRegistry::create_all(*this, cfg_);
  for (const auto& fix : fixes_) {
    fix->init();
  }

  gpu_set_config(cfg.gpu);
  gpu_init_for_rank(domain_.rank(), domain_.nprocs());
  gpu_.active = gpu_runtime_enabled();
  if (gpu_.active && qssa_.agent_sampling()) {
    throw ConfigError(
        "chemistry.toxin_evaluation=agents is not supported with GPU execution");
  }

#ifndef GUTIBM_HDF5
  (void)h5_file;
  (void)step;
  throw SimulationError("checkpoint restart requires HDF5 support");
#else
  HDF5CheckpointSnapshot snap = HDF5Reader::load_snapshot(h5_file, step);
  validate_checkpoint_species_set(cfg_, snap);
  validate_checkpoint_toxin_species(cfg_, snap);
  if (snap.metadata.has_grid_spacing) {
    const Vec3 expected = {
        domain_.dx_x(), domain_.dx_y(), domain_.dx_z()};
    for (Int axis = 0; axis < 3; ++axis) {
      if (std::abs(snap.metadata.grid_spacing[static_cast<size_t>(axis)]
                   - expected[static_cast<size_t>(axis)])
          > 1.0e-12 * std::max(1.0, expected[static_cast<size_t>(axis)])) {
        throw SimulationError(
            "checkpoint chemistry spacing does not match requested stride");
      }
    }
  }
  apply_checkpoint_snapshot(snap);
  dysbiosis_.configure(cfg_.dysbiosis_threshold,
                       cfg_.dysbiosis_sampling_interval,
                       cfg_.dysbiosis_sample_count);
  dysbiosis_.reset(clock_.time);

  // Sync GPU mirrors *after* host restore (agents/grid were empty above).
  if (gpu_.active) {
    gpu_.chem.init(chem_);
    gpu_.agents.sync_from_host(agents_);
  }

  rebuild_spatial_hash();
  update_grid_coupling();
  allreduce_global_stats();

  int rank = domain_.rank();
  if (rank == 0) {
    std::cout << "GutIBM restored from checkpoint:\n"
              << "  File: " << h5_file << "\n"
              << "  Step group: " << snap.step_name << "\n"
              << "  Restored time: " << clock_.time << " s\n"
              << "  Restored step: " << clock_.step_count << "\n"
              << "  Global agents: " << mpi_ghost_.stats.global_agent_count << "\n"
              << "  Local agents: " << agents_.size() << "\n"
              << "  Global density: "
              << global_density_cells_per_mL(domain_, mpi_ghost_.stats.global_agent_count)
              << " cells/mL\n";
    print_gpu_status_banner(gpu_.active, cfg_.gpu);
    std::cout << std::flush;
  }
#endif
}

void Simulation::apply_checkpoint_snapshot(const HDF5CheckpointSnapshot& snap) {
  const auto& atoms = snap.agents;
  const auto& lin   = snap.lineage;
  const auto& gen   = snap.genome;
  const size_t n    = atoms.id.size();

  TagID max_tag = 0;
  agents_.reserve(static_cast<Int>(n));

  const std::vector<size_t> bi_offsets = build_bi_offsets(lin.num_bi_loci);
  validate_checkpoint_genome(snap, bi_offsets);

  size_t i = 0;
  for (int64_t id : atoms.id) {
    (void)id;
    Vec3 pos = {atoms.x[i], atoms.y[i], atoms.z[i]};
    if (!domain_.is_local(pos)) {
      ++i;
      continue;
    }

    Real mu_guess = std::max(atoms.mu[i], 1.0e-8);
    Agent a = Agent::create_default(static_cast<TagID>(atoms.id[i]),
                                    atoms.type[i], pos, mu_guess);
    a.identity.owner_rank   = domain_.rank();
    a.state        = static_cast<PhenoState>(atoms.state[i]);
    a.radius       = atoms.radius[i];
    a.outer_radius = atoms.radius[i] * 1.05;
    a.mass         = sphere_mass(a.radius, CELL_DENSITY_DEFAULT);
    a.biomass      = atoms.biomass[i];
    a.mu_realized  = atoms.mu[i];
    // Intrinsic capacity must survive low mu_realized snapshots. Older files
    // omit /mu_max — fall back to the matching strain (never to mu_realized).
    if (i < atoms.mu_max.size() && atoms.mu_max[i] > 0.0) {
      a.mu_max = atoms.mu_max[i];
    } else {
      a.mu_max = fallback_mu_max(cfg_, atoms.type[i]);
    }
    if (i < atoms.realized_fermentation_fraction.size()) {
      a.realized_fermentation_fraction =
          std::clamp(atoms.realized_fermentation_fraction[i], 0.0, 1.0);
    }

    a.genome.lineage_id = static_cast<TagID>(atoms.lineage[i]);
    a.genome.generation = static_cast<uint32_t>(lin.generation[i]);
    a.receptor_expr[to_underlying(ReceptorType::BtuB)] = lin.btuB_expression[i];
    a.receptor_expr[to_underlying(ReceptorType::FepA)] = lin.fepA_expression[i];

    if (gen.present) {
      a.genome.parent_id = static_cast<TagID>(gen.parent_id[i]);
      a.genome.mutations = static_cast<uint32_t>(gen.mutations[i]);
      a.genome.has_conjugative_plasmid = (gen.has_conjugative_plasmid[i] != 0);
      a.genome.plasmid_cost_amelioration = gen.plasmid_cost_amelioration[i];
      restore_receptor_fields(a, snap, i);
      restore_bi_loci(a, snap, i, bi_offsets);
      if (!gen.cdi_type.empty()) {
        a.genome.cdi_type = static_cast<uint16_t>(gen.cdi_type[i]);
        a.genome.cdi_immunity = static_cast<uint16_t>(gen.cdi_immunity[i]);
      }
    } else {
#ifdef GUTIBM_LEGACY_CHECKPOINTS
      restore_legacy_genome(a, snap, i);
#else
      throw SimulationError(
          "checkpoint has no genome group (legacy format); rebuild with "
          "-DGUTIBM_LEGACY_CHECKPOINTS to load pre-genome snapshots");
#endif
    }

    if (i < atoms.in_crypt.size()) {
      a.flags.in_crypt = (atoms.in_crypt[i] != 0);
    } else {
      tag_crypt_resident(a, advection_);
    }

    max_tag = std::max(max_tag, static_cast<TagID>(atoms.id[i]));
    agents_.push_back(std::move(a));
    ++i;
  }

  agents_.configure_tags(
      AgentPool::next_tag_after_max(max_tag, domain_.rank(), domain_.nprocs()),
      AgentPool::tag_stride(domain_.nprocs()));

  restore_checkpoint_grid(chem_, domain_, snap);

  clock_.time       = snap.metadata.time;
  clock_.step_count = snap.metadata.step;
  event_ledger_.cumulative_events = snap.metadata.cumulative_events;
  if (snap.metadata.flux_accounting.boundary_interval.size()
      == static_cast<size_t>(chem_.num_species())) {
    chem_.flux_accounting() = snap.metadata.flux_accounting;
    auto& flux = chem_.flux_accounting();
    // A closed restart closes the current reporting window at write time.
    // Restore cumulative totals, but never reopen that already-closed window.
    std::ranges::fill(flux.boundary_interval, 0.0);
    std::ranges::fill(flux.boundary_step, 0.0);
    std::ranges::fill(flux.gradient_source_interval, 0.0);
    std::ranges::fill(flux.gradient_source_step, 0.0);
    std::ranges::fill(flux.vbf_source_interval, 0.0);
    std::ranges::fill(flux.vbf_sink_interval, 0.0);
    std::ranges::fill(flux.agent_uptake_interval, 0.0);
    std::ranges::fill(flux.agent_uptake_step, 0.0);
    std::ranges::fill(flux.maintenance_interval, 0.0);
    std::ranges::fill(flux.maintenance_step, 0.0);
    std::ranges::fill(flux.maintenance_shortfall_interval, 0.0);
    std::ranges::fill(flux.maintenance_shortfall_step, 0.0);
    std::ranges::fill(flux.maintenance_limited_agents_interval, 0.0);
    std::ranges::fill(flux.maintenance_limited_agents_step, 0.0);
    std::ranges::fill(flux.uptake_demand_interval, 0.0);
    std::ranges::fill(flux.uptake_demand_step, 0.0);
    std::ranges::fill(flux.uptake_shortfall_interval, 0.0);
    std::ranges::fill(flux.uptake_shortfall_step, 0.0);
    std::ranges::fill(flux.uptake_limited_interval, 0.0);
    std::ranges::fill(flux.uptake_limited_step, 0.0);
    std::ranges::fill(flux.reaction_clip_step, 0.0);
    const auto species_count = static_cast<size_t>(chem_.num_species());
    const auto ensure_sized = [species_count](std::vector<Real>& values) {
      if (values.size() != species_count) values.assign(species_count, 0.0);
    };
    ensure_sized(flux.uptake_demand_interval);
    ensure_sized(flux.gradient_source_interval);
    ensure_sized(flux.gradient_source_step);
    ensure_sized(flux.gradient_source_cumulative);
    ensure_sized(flux.uptake_demand_step);
    ensure_sized(flux.uptake_demand_cumulative);
    ensure_sized(flux.uptake_shortfall_interval);
    ensure_sized(flux.uptake_shortfall_step);
    ensure_sized(flux.uptake_shortfall_cumulative);
    ensure_sized(flux.uptake_limited_interval);
    ensure_sized(flux.uptake_limited_step);
    ensure_sized(flux.uptake_limited_cumulative);
    ensure_sized(flux.agent_uptake_last_step);
    ensure_sized(flux.uptake_demand_last_step);
    ensure_sized(flux.reaction_clip_last_step);
    ensure_sized(flux.maintenance_interval);
    ensure_sized(flux.maintenance_step);
    ensure_sized(flux.maintenance_last_step);
    ensure_sized(flux.maintenance_shortfall_last_step);
    ensure_sized(flux.maintenance_cumulative);
    ensure_sized(flux.maintenance_shortfall_interval);
    ensure_sized(flux.maintenance_shortfall_step);
    ensure_sized(flux.maintenance_shortfall_cumulative);
    ensure_sized(flux.maintenance_limited_agents_interval);
    ensure_sized(flux.maintenance_limited_agents_step);
    ensure_sized(flux.maintenance_limited_agents_cumulative);
    if (flux.reaction_clip_interval.size()
        != static_cast<size_t>(chem_.num_species())) {
      flux.reaction_clip_interval.assign(
          static_cast<size_t>(chem_.num_species()), 0.0);
    }
    if (flux.reaction_clip_cumulative.size()
        != static_cast<size_t>(chem_.num_species())) {
      flux.reaction_clip_cumulative.assign(
          static_cast<size_t>(chem_.num_species()), 0.0);
    }
  }
  event_ledger_.window_start_step = snap.metadata.event_window_end_step > 0
      ? snap.metadata.event_window_end_step + 1
      : clock_.step_count + 1;
  event_ledger_.window_start_time = snap.metadata.event_window_end_step > 0
      ? snap.metadata.event_window_end_time
      : clock_.time;
  immigration_.set_start_step(clock_.step_count);
  schedule_output_from_time(clock_.time, cfg_.time.output_interval, clock_.next_output, clock_.next_snapshot);
}

Agent Simulation::create_strain_agent(
    const SimulationConfig::InitialStrain& strain, Vec3 pos) {
  Agent agent = Agent::create_default(agents_.next_tag(), strain.type, pos,
                                      strain.mu_max);
  agent.identity.owner_rank = domain_.rank();
  assign_plasmids(agent, strain.plasmids, cfg_, domain_.rank());
  for (const auto& [name, expression] : strain.receptor_expression) {
    const auto receptor = receptor_type_from_name(name);
    if (!receptor.has_value()) {
      throw ConfigError("unknown receptor name: " + name);
    }
    const Int index = to_underlying(*receptor);
    agent.receptor_expr_base[index] = expression;
    agent.receptor_expr[index] = expression;
    agent.genome.receptor_expression[index] = expression;
    if (receptor_expression_is_resistant(expression)) {
      agent.state = PhenoState::RESISTANT;
    }
  }
  agent.genome.cdi_type = strain.cdi_type;
  agent.genome.cdi_immunity = strain.cdi_immunity;
  tag_crypt_resident(agent, advection_);
  if (cfg_.cell_bio.motility.enabled) {
    FixMotility::init_agent_motility(agent, cfg_.cell_bio.motility,
                                     rng_);
  }
  return agent;
}

void Simulation::maybe_write_restart() {
  if (!cfg_.restart.enabled || cfg_.restart.interval_steps <= 0) return;
  if (clock_.step_count <= 0) return;
  if (clock_.step_count % cfg_.restart.interval_steps != 0) return;
  write_restart_now();
}

void Simulation::write_restart_now() {
  if (!cfg_.restart.enabled) return;
  if (cfg_.restart.directory.empty()) return;

  const std::string step_name = std::format("step_{:06}", clock_.step_count);
  const std::string path = cfg_.restart.directory + "/" + step_name + ".h5";
  const auto write_t0 = std::chrono::steady_clock::now();
  const bool preserve_event_counters = cfg_.hdf5.enabled
      && cfg_.hdf5.schedule.summary > 0
      && clock_.step_count % cfg_.hdf5.schedule.summary == 0;
  const StepEvents saved_step_events = step_events();
  const StepEvents saved_summary_events = summary_events();
  const StepEvents saved_cumulative_events = cumulative_events();
  const Int saved_window_start_step = event_window_start_step();
  const Real saved_window_start_time = event_window_start_time();
  const bool ok = HDF5Writer::write_closed_restart(
      *this, path, clock_.step_count, clock_.time, cfg_.time.bio_dt,
      preserve_event_counters);
  if (!ok && !preserve_event_counters) {
    step_events() = saved_step_events;
    summary_events() = saved_summary_events;
    cumulative_events() = saved_cumulative_events;
    set_event_window_start(saved_window_start_step, saved_window_start_time);
  }
  const double write_s = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - write_t0).count();
  if (ok && domain_.rank() == 0) {
    const auto old_precision = std::cout.precision();
    std::cout << "Wrote closed restart: " << path
              << "  (t=" << clock_.time << "s step=" << clock_.step_count
              << " write_s=" << std::fixed << std::setprecision(3) << write_s
              << ")\n"
              << std::flush;
    std::cout.unsetf(std::ios_base::floatfield);
    std::cout.precision(old_precision);
  } else if (!ok && domain_.rank() == 0) {
    const auto old_precision = std::cerr.precision();
    std::cerr << "Closed restart failed: " << path
              << "  (write_s=" << std::fixed << std::setprecision(3) << write_s
              << ")\n"
              << std::flush;
    std::cerr.unsetf(std::ios_base::floatfield);
    std::cerr.precision(old_precision);
  }
}

std::vector<std::string> Simulation::fix_names() const {
  std::vector<std::string> names;
  names.reserve(fixes_.size());
  for (const auto& fix : fixes_) {
    names.push_back(fix->name());
  }
  return names;
}

Real Simulation::compute_adaptive_dt() const {
  if (!cfg_.adaptive_dt.enabled) return cfg_.time.bio_dt;

  Real dt = cfg_.adaptive_dt.max;

  // Growth rate constraint: mu_max * dt < growth_limit
  if (mpi_ghost_.stats.global_max_abs_mu > 0) {
    dt = std::min(
        dt, cfg_.adaptive_dt.growth_limit / mpi_ghost_.stats.global_max_abs_mu);
  }

  // SOS cascade constraint: reduce dt during active lysis
  if (mpi_ghost_.stats.global_sos_count > 5)  dt = std::min(dt, 10.0);
  if (mpi_ghost_.stats.global_sos_count > 20) dt = std::min(dt, 2.0);

  // Agent density constraint
  const Real density = global_density_cells_per_mL(
      domain_, mpi_ghost_.stats.global_agent_count);
  // 1e15 cells/m^3 is 1e9 cells/mL.
  if (constexpr Real kDensityLimitCellsPerML = 1.0e9;
      density > kDensityLimitCellsPerML) {
    dt = std::min(dt, 10.0);
  }

  // Apply safety factor and bounds
  dt *= cfg_.adaptive_dt.safety;
  dt = std::clamp(dt, cfg_.adaptive_dt.min, cfg_.adaptive_dt.max);
  return dt;
}

void Simulation::print_step_profile() const {
  if (step_profile_.step_count <= 0) return;

  const int n = step_profile_.step_count;
  const double inv = 1.0 / static_cast<double>(n);
  const double total = step_profile_.total_s() * inv;

  std::cout << "Step profile (mean wall time per step, s):\n"
            << "  ghost_exchange=" << step_profile_.ghost_exchange_s * inv << "\n"
            << "  spatial_hash=" << step_profile_.spatial_hash_s * inv << "\n"
            << "  biology=" << step_profile_.biology_s * inv << "\n"
            << "  chemistry=" << step_profile_.chemistry_s * inv << "\n"
            << "  physics=" << step_profile_.physics_s * inv << "\n"
            << "  mpi_migrate=" << step_profile_.mpi_migrate_s * inv << "\n"
            << "  cleanup=" << step_profile_.cleanup_s * inv << "\n"
            << "  gpu_h2d=" << step_profile_.gpu_h2d_s * inv << "\n"
            << "  gpu_d2h=" << step_profile_.gpu_d2h_s * inv << "\n"
            << "  gpu_slab_x_roundtrip="
            << step_profile_.gpu_slab_x_roundtrip_s * inv << "\n"
            << "  mpi_reaction_reduce=" << step_profile_.mpi_reaction_reduce_s * inv << "\n"
            << "  hdf5=" << step_profile_.hdf5_s * inv << "\n"
            << "  total=" << total << "\n"
            << "PROFILE_CSV steps=" << n
            << " ghost_s=" << step_profile_.ghost_exchange_s * inv
            << " hash_s=" << step_profile_.spatial_hash_s * inv
            << " biology_s=" << step_profile_.biology_s * inv
            << " chemistry_s=" << step_profile_.chemistry_s * inv
            << " physics_s=" << step_profile_.physics_s * inv
            << " mpi_s=" << step_profile_.mpi_migrate_s * inv
            << " cleanup_s=" << step_profile_.cleanup_s * inv
            << " gpu_h2d_s=" << step_profile_.gpu_h2d_s * inv
            << " gpu_d2h_s=" << step_profile_.gpu_d2h_s * inv
            << " gpu_slab_x_roundtrip_s="
            << step_profile_.gpu_slab_x_roundtrip_s * inv
            << " mpi_reaction_reduce_s=" << step_profile_.mpi_reaction_reduce_s * inv
            << " hdf5_s=" << step_profile_.hdf5_s * inv
            << " total_s=" << total
            << " agents=" << mpi_ghost_.stats.global_agent_count
            << " ranks=" << domain_.nprocs()
            << "\n"
            << std::flush;
}

namespace {
constexpr Int kPopulationStopThreshold = 1;

// Progress line fields parsed by deploy/aws/entry.sh into status.json.
}  // namespace

ProgressMetrics calculate_progress_metrics(Real sim_time,
                                           Real attempt_sim_time,
                                           Real total_time,
                                           double wall_elapsed_s) {
  ProgressMetrics metrics;
  metrics.pct = (total_time > 0.0) ? (100.0 * sim_time / total_time) : 0.0;
  metrics.rate = (wall_elapsed_s > 0.0)
                     ? static_cast<double>(attempt_sim_time) / wall_elapsed_s
                     : 0.0;
  if (metrics.rate > 0.0 && total_time > sim_time) {
    metrics.eta_s = static_cast<double>(total_time - sim_time) / metrics.rate;
  }
  return metrics;
}

namespace {
struct ProgressLineData {
  Int step_count;
  Real sim_time;
  Real dt;
  Int global_agents;
  Real global_density;
  Int local_agents;
  Real mu_avg;
  Real total_time;
  Real attempt_sim_time;
  double wall_elapsed_s;
};

void print_progress_line(const ProgressLineData& data) {
  const ProgressMetrics metrics =
      calculate_progress_metrics(data.sim_time, data.attempt_sim_time,
                                 data.total_time, data.wall_elapsed_s);
  std::cout << "Step " << data.step_count
            << "  t=" << data.sim_time << "s"
            << "  dt=" << std::setprecision(3) << data.dt << "s"
            << "  global_agents=" << data.global_agents
            << "  density_cells_per_mL=" << data.global_density
            << "  local_agents=" << data.local_agents
            << "  mu_avg=" << data.mu_avg
            << "  pct=" << std::setprecision(4) << metrics.pct
            << "  rate=" << std::setprecision(4) << metrics.rate
            << "  eta_s=" << std::setprecision(0) << std::fixed
            << metrics.eta_s
            << "\n" << std::flush;
  std::cout.unsetf(std::ios_base::floatfield);
}

void print_heartbeat_line(Int step_count, Real sim_time, double wall_elapsed_s) {
  std::cout << "Heartbeat step=" << step_count
            << "  t=" << sim_time << "s"
            << "  wall_elapsed_s=" << wall_elapsed_s
            << "\n" << std::flush;
}
}  // namespace

void Simulation::write_hdf5_step(Real dt) {
  if (!hdf5_.is_enabled()) return;
  if (qssa_.agent_sampling() && hdf5_.grid_output_due(clock_.step_count)
      && !bacteriocin_.grid_materialized) {
    materialize_bacteriocin_fields_for_output();
  }
  const auto hdf5_t0 = std::chrono::steady_clock::now();
  hdf5_.write_step(*this, clock_.step_count, clock_.time, dt);
  if (cfg_.profile_steps) {
    const auto hdf5_t1 = std::chrono::steady_clock::now();
    step_profile_.hdf5_s +=
        std::chrono::duration<double>(hdf5_t1 - hdf5_t0).count();
  }
}

void Simulation::emit_heartbeat(
    const std::chrono::steady_clock::time_point& wall_start,
    const std::chrono::steady_clock::time_point& wall_now,
    std::chrono::steady_clock::time_point& next_heartbeat,
    bool& heartbeat_emitted) const {
  if (domain_.rank() != 0 ||
      (heartbeat_emitted && wall_now < next_heartbeat)) {
    return;
  }
  const double wall_elapsed_s =
      std::chrono::duration<double>(wall_now - wall_start).count();
  print_heartbeat_line(clock_.step_count, clock_.time, wall_elapsed_s);
  heartbeat_emitted = true;
  next_heartbeat = wall_now + std::chrono::seconds(60);
}

void Simulation::emit_progress_if_due(
    Real dt, const std::chrono::steady_clock::time_point& wall_start,
    Real attempt_start_sim_time,
    const std::chrono::steady_clock::time_point& wall_now) {
  if (clock_.time < clock_.next_output) return;
  if (domain_.rank() == 0) {
    const double wall_elapsed_s =
        std::chrono::duration<double>(wall_now - wall_start).count();
    const ProgressLineData data{
        clock_.step_count,
        clock_.time,
        dt,
        mpi_ghost_.stats.global_agent_count,
        global_density_cells_per_mL(domain_, mpi_ghost_.stats.global_agent_count),
        agents_.size(),
        mpi_ghost_.stats.global_mu_avg,
        cfg_.time.total_time,
        clock_.time - attempt_start_sim_time,
        wall_elapsed_s};
    print_progress_line(data);
  }
  clock_.next_output += cfg_.time.output_interval;
}

void Simulation::update_lineage_snapshot_if_due() {
  if (clock_.time < clock_.next_snapshot) return;
  take_lineage_snapshot();
  clock_.next_snapshot += cfg_.time.output_interval;
}

bool Simulation::population_stop(int rank) const {
  if (mpi_ghost_.stats.global_agent_count > kPopulationStopThreshold) {
    return false;
  }
  if (rank == 0) {
    std::cout << "Population stop: " << mpi_ghost_.stats.global_agent_count
              << " cell(s) — ending simulation.\n"
              << std::flush;
  }
  return true;
}

bool Simulation::dysbiosis_threshold_exceeded(int rank) {
  if (cfg_.dysbiosis_threshold <= 0.0) return false;
  const Real density_cells_per_mL =
      global_density_cells_per_mL(domain_, mpi_ghost_.stats.global_agent_count);
  if (!dysbiosis_.observe(clock_.time, density_cells_per_mL)) return false;
  if (rank == 0) {
    std::cerr << "DYSBIOSIS THRESHOLD EXCEEDED: "
              << density_cells_per_mL
              << " cells/mL > " << cfg_.dysbiosis_threshold
              << " cells/mL, increasing at "
              << dysbiosis_.density_rate_cells_per_mL_per_s()
              << " cells/mL/s with non-decelerating growth"
              << " — halting simulation.\n"
              << std::flush;
  }
  return true;
}

bool Simulation::closure_violation(std::string& detail) {
  const bool delivery_mode =
      cfg_.fixes.metabolism.uptake_limit == "delivery";
  const bool delivery_gate = cfg_.closure.enforce_delivery_realization
      && delivery_mode;
  const bool reaction_gate = cfg_.closure.enforce_reaction_clip;
  if (!delivery_gate && !reaction_gate) {
    return false;
  }
  const Int carbon = chem_.find(species::CARBON);
  bool local_violation = false;
  if (carbon >= 0) {
    const auto& flux = chem_.flux_accounting();
    const Real demand = flux.uptake_demand_for_step(carbon)
        + flux.maintenance_last_step[static_cast<size_t>(carbon)]
        + flux.maintenance_shortfall_last_step[
            static_cast<size_t>(carbon)];
    const Real realized = flux.agent_uptake_for_step(carbon)
        + flux.maintenance_last_step[static_cast<size_t>(carbon)];
    if (delivery_gate) {
      if (demand > 0.0 && realized == 0.0) {
        ++zero_realization_steps_;
      } else {
        zero_realization_steps_ = 0;
      }
      if (zero_realization_steps_ >
          cfg_.closure.zero_realization_grace_steps) {
        local_violation = true;
        detail = std::format(
            "species=carbon step={} demand={} realized_removal={}",
            clock_.step_count, demand, realized);
      }
    } else {
      zero_realization_steps_ = 0;
    }

    if (reaction_gate && !local_violation) {
      const Real clip = flux.reaction_clip_for_step(carbon);
      const Real tolerance =
          cfg_.closure.reaction_clip_tolerance_fraction *
          std::max(demand, std::numeric_limits<Real>::epsilon());
      if (clip > tolerance) {
        local_violation = true;
        detail = std::format(
            "species=carbon step={} reaction_clip={} tolerance={}",
            clock_.step_count, clip, tolerance);
      }
    }
  }
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return local_violation;
  int global_violation = local_violation ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &global_violation, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  if (global_violation != 0 && !local_violation) {
    detail = "closure violation reported by another MPI rank";
  }
  return global_violation != 0;
#else
  return local_violation;
#endif
}

int Simulation::run() {
  int rank = domain_.rank();
  Real last_dt = cfg_.time.bio_dt;
  bool stopped_for_population = false;
  dysbiosis_.reset(clock_.time);
  const auto wall_start = std::chrono::steady_clock::now();
  const Real attempt_start_sim_time = clock_.time;
  auto next_heartbeat = wall_start;
  bool heartbeat_emitted = false;

  dysbiosis_.observe(
      clock_.time,
      global_density_cells_per_mL(domain_, mpi_ghost_.stats.global_agent_count));

  // Initial snapshot only for fresh runs (not checkpoint resume).
  if (hdf5_.is_enabled() && clock_.step_count == 0) {
    hdf5_.write_step(*this, 0, clock_.time, last_dt);
  }
  if (rank == 0) {
    take_lineage_snapshot();
  }

  stopped_for_population = population_stop(rank);
  termination_cause_ = stopped_for_population
      ? TerminationCause::PopulationStop : TerminationCause::IncompleteUnknown;
  termination_detail_ = stopped_for_population
      ? "population stop threshold reached" : "run in progress";

  while (!stopped_for_population && !dysbiosis_.halted() &&
         clock_.time < cfg_.time.total_time) {
    if (gutibm_stop_requested()) {
      termination_cause_ = TerminationCause::StopRequested;
      termination_detail_ = "stop requested before next timestep";
      break;
    }

    Real dt = compute_adaptive_dt();
    last_dt = dt;

    // Clamp so we don't overshoot total_time
    if (clock_.time + dt > cfg_.time.total_time) dt = cfg_.time.total_time - clock_.time;

    step(dt);

    maybe_write_restart();

    // Persist the triggering step before a closure violation breaks the loop.
    write_hdf5_step(last_dt);

    if (closure_violation(termination_detail_)) {
      termination_cause_ = TerminationCause::ClosureViolation;
      break;
    }

    const auto wall_now = std::chrono::steady_clock::now();
    emit_heartbeat(wall_start, wall_now, next_heartbeat, heartbeat_emitted);

    // Console progress and in-memory lineage snapshots use output_interval (seconds).
    emit_progress_if_due(dt, wall_start, attempt_start_sim_time, wall_now);

    // Lineage snapshots
    update_lineage_snapshot_if_due();

    stopped_for_population = population_stop(rank);
    if (!stopped_for_population) {
      const bool halted_for_dysbiosis = dysbiosis_threshold_exceeded(rank);
      if (halted_for_dysbiosis) {
        termination_cause_ = TerminationCause::DysbiosisGuard;
        termination_detail_ = "dysbiosis guard threshold exceeded";
      }
      if (halted_for_dysbiosis && hdf5_.is_enabled()) {
        hdf5_.write_halt_metadata(*this, clock_.step_count);
      }
    } else {
      termination_cause_ = TerminationCause::PopulationStop;
      termination_detail_ = "population stop threshold reached";
    }
  }

  if (termination_cause_ == TerminationCause::IncompleteUnknown) {
    if (clock_.time >= cfg_.time.total_time) {
      termination_cause_ = TerminationCause::HorizonReached;
      termination_detail_ = "requested simulation horizon reached";
    } else if (dysbiosis_.halted()) {
      termination_cause_ = TerminationCause::DysbiosisGuard;
      termination_detail_ = "dysbiosis guard threshold exceeded";
    }
  }

  // Final closed restart so Spot/SIGTERM/early-exit still leaves a usable artifact.
  if (const bool already_checkpointed = cfg_.restart.interval_steps > 0
          && clock_.step_count % cfg_.restart.interval_steps == 0;
      cfg_.restart.enabled && clock_.step_count > 0
      && (dysbiosis_.halted() || !already_checkpointed)) {
    write_restart_now();
  }

  termination_wall_seconds_ = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - wall_start).count();
  if (hdf5_.is_enabled()) {
    hdf5_.write_run_termination(*this, clock_.step_count, clock_.time);
  }
  hdf5_.finalize();

  if (rank == 0) {
    Real retention = lineage_.resident_retention(cfg_.time.total_time * 0.5);
    const auto cause = termination_cause_name(termination_cause_);
    const bool complete = termination_cause_ == TerminationCause::HorizonReached;
    const std::string detail = complete
        ? "" : "  detail=" + termination_detail_;
    std::cout << "\nSimulation " << (complete ? "complete" : "ended early")
              << ": cause=" << cause
              << "  reached_time=" << clock_.time
              << "s  requested_time=" << cfg_.time.total_time
              << "s  step=" << clock_.step_count << detail << "\n"
              << "  Final global agents: " << mpi_ghost_.stats.global_agent_count << "\n"
              << "  Final global density: "
              << global_density_cells_per_mL(domain_, mpi_ghost_.stats.global_agent_count)
              << " cells/mL\n"
              << "  Steps taken: " << clock_.step_count << "\n"
              << "  Resident retention: " << retention * 100.0 << "%\n"
              << "  Dominant lineage: " << lineage_.dominant_lineage() << "\n"
              << std::flush;
    if (cfg_.profile_steps) {
      print_step_profile();
    }
  }
  return termination_cause_ == TerminationCause::ClosureViolation ? 1 : 0;
}

void Simulation::step(Real dt) {
  bacteriocin_.fields_current = false;
  bacteriocin_.grid_materialized = false;
  if (cfg_.profile_steps) {
    gpu_transfer_profile_set_enabled(gpu_.active);
  }

  StepProfiler profiler(cfg_.profile_steps);
  profiler.start();

  for (Agent& a : agents_) {
    a.flags.just_divided = false;
    a.flags.microcin_penalty_applied = false;
  }

  // Update advection time for peristaltic oscillation
  advection_.set_time(clock_.time);

  // Pre-step: clear ghosts from previous step
  clear_ghost_agents();
  const auto& immigration = cfg_.immigration;
  immigration_.inject(
      immigration, clock_.step_count, dt, agents_, domain_,
      [this](const Vec3& position) {
        const auto& strain = cfg_.initial_strains[static_cast<size_t>(
            cfg_.immigration.strain_index)];
        agents_.push_back(create_strain_agent(strain, position));
        ++event_ledger_.step_events.immigrations;
      });
  if (gpu_.active) {
    gpu_.chem.zero_reactions_on_device();
  }
  chem_.zero_reactions();

  // Exchange ghost agents for cross-boundary neighbor queries
  exchange_ghost_agents();
  profiler.lap(step_profile_.ghost_exchange_s);

  // Chemistry changed owned concentrations in the previous step; refresh the
  // halos before any ghost-agent or gradient/stencil read in biology, and
  // before the device mirror is uploaded.
  chem_.exchange_concentration_halos();

  if (gpu_.active) {
    gpu_.agents.sync_from_host(agents_);
    gpu_.chem.sync_to_device(chem_);
#ifdef GUTIBM_CUDA
    gpu::launch_grid_coupling_kernel(
        gpu_.agents.x(), gpu_.agents.y(), gpu_.agents.z(),
        gpu_.agents.grid_cell(), gpu_.agents.state(),
        domain_.lo()[0], domain_.lo()[1], domain_.lo()[2],
        domain_.dx_x(), domain_.dx_y(), domain_.dx_z(),
        domain_.nx(), domain_.ny(), domain_.nz(),
        agents_.size(), gpu_compute_stream());
    gpu_sync_compute();
    gpu_check_error("grid_coupling_kernel");
    gpu_.agents.sync_to_host(agents_);
#endif
  }

  rebuild_spatial_hash();
  update_grid_coupling();
  profiler.lap(step_profile_.spatial_hash_s);

  for (const auto& fix : fixes_) {
    fix->pre_step(dt);
  }

  // 1. Biology module (uses ghost agents for neighbor interactions)
  module_biology(dt);
  profiler.lap(step_profile_.biology_s);

  // Clear ghosts before physics to avoid moving them
  clear_ghost_agents();
  rebuild_spatial_hash();
  update_grid_coupling();

  // 2. Chemistry module (QSSA, instantaneous)
  module_chemistry(dt);
  profiler.lap(step_profile_.chemistry_s);
  for (const auto& fix : fixes_) {
    fix->post_chemistry(dt);
  }
  chem_.debug_report_step(domain_);

  // 3. Physics module (advection + mechanics)
  module_physics(dt);
  profiler.lap(step_profile_.physics_s);

  // Post-step
  for (const auto& fix : fixes_) {
    fix->post_step(dt);
  }

  // Migrate agents that crossed slab boundaries
  migrate_agents();
  profiler.lap(step_profile_.mpi_migrate_s);

  if (gpu_.active) {
    gpu_.agents.sync_from_host(agents_);
  }

  // Cleanup
  check_washout();
  remove_dead_agents();

  // Compute global statistics
  allreduce_global_stats();
  profiler.lap(step_profile_.cleanup_s);

  if (cfg_.profile_steps) {
    step_profile_.step_count++;
    const GpuTransferProfile xfer = gpu_transfer_profile_snapshot();
    step_profile_.gpu_h2d_s += xfer.h2d_s;
    step_profile_.gpu_d2h_s += xfer.d2h_s;
    step_profile_.gpu_slab_x_roundtrip_s += xfer.slab_x_roundtrip_s;
    gpu_transfer_profile_reset();
  }

  clock_.time += dt;
  clock_.step_count++;
}

void Simulation::update_bacteriocin_fields() {
  if (bacteriocin_.fields_current) return;

  prune_toxin_bursts(clock_.time);

  ChemicalFieldGpu* chem_gpu_ptr = gpu_.active ? &gpu_.chem : nullptr;
  const bool materialize_grid =
      hdf5_.grid_output_due(clock_.step_count + 1);
  qssa_.solve_all_bacteriocin_fields(
      agents_, bacteriocin_.bursts, clock_.time, cfg_.chem_env.protease,
      advection_, chem_, chem_gpu_ptr, materialize_grid);
  bacteriocin_.fields_current = true;
  bacteriocin_.grid_materialized = materialize_grid;
}

void Simulation::materialize_bacteriocin_fields_for_output() {
  if (!qssa_.agent_sampling()) return;
  prune_toxin_bursts(clock_.time);
  ChemicalFieldGpu* chem_gpu_ptr = gpu_.active ? &gpu_.chem : nullptr;
  qssa_.solve_all_bacteriocin_fields(
      agents_, bacteriocin_.bursts, clock_.time, cfg_.chem_env.protease,
      advection_, chem_, chem_gpu_ptr, true);
  bacteriocin_.fields_current = true;
  bacteriocin_.grid_materialized = true;
}

void Simulation::module_biology(Real dt) {
  for (const auto& fix : fixes_) {
    // Receptor killing reads the chemical grid; deposit current-step toxin
    // sources (microcin producers + active bursts) before fix_receptor runs.
    if (fix->name() == "receptor") {
      update_bacteriocin_fields();
    }
    fix->compute(dt);
  }
}

void Simulation::module_chemistry(Real dt) {
  update_bacteriocin_fields();

  ChemistryPipelineInput pipeline{
      .gpu_active = gpu_.active,
      .delivery_mode = cfg_.fixes.metabolism.uptake_limit_mode
          == UptakeLimitMode::Delivery,
      .agents_gpu = gpu_.agents,
      .chem_gpu = gpu_.chem,
      .chem = chem_,
      .domain = domain_,
      .vbf = vbf_,
      .qssa = qssa_,
      .agents = agents_,
      .oxygen = cfg_.chem_env.oxygen,
      .acetate = cfg_.chem_env.acetate,
      .mucin = cfg_.chem_env.mucin,
      .num_agents = agents_.size(),
      .flux_accounting = chem_.flux_accounting(),
      .step_profile = cfg_.profile_steps ? &step_profile_ : nullptr,
  };
  (void)run_chemistry_pipeline(pipeline, dt);
}

void Simulation::module_physics(Real dt) {
  // Crypt migration (stochastic entry/exit) before advection
  crypt_migration(dt);

  // Advection pass: each agent independent (read-only on fields)
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Agent& a : agents_) {
    if (a.state == PhenoState::DEAD) continue;

    // Advection: mucus flow carries agent
    advection_.advect(a.x, dt);

    // VBF drag modifies velocity
    Vec3 drag = vbf_.drag_force(a.v);
    Real inv_mass = 1.0 / std::max(a.mass, 1.0e-30);
    a.v[0] += drag[0] * inv_mass * dt;
    a.v[1] += drag[1] * inv_mass * dt;
    a.v[2] += drag[2] * inv_mass * dt;

    // Apply position update from velocity
    a.x[0] += a.v[0] * dt;
    a.x[1] += a.v[1] * dt;
    a.x[2] += a.v[2] * dt;

    if (cfg_.cell_bio.motility.enabled) {
      a.x[0] += a.motility.step_displacement[0];
      a.x[1] += a.motility.step_displacement[1];
      a.x[2] += a.motility.step_displacement[2];
      a.motility.step_displacement = {0.0, 0.0, 0.0};
    }

    // PBC / boundary
    domain_.apply_pbc(a.x);
  }

  // Mechanical repulsion handled by FixMechanics (registered as a fix)
  for (const auto& fix : fixes_) {
    if (fix->name() == "mechanics") {
      fix->compute(dt);
      break;
    }
  }
}

void Simulation::rebuild_spatial_hash() {
  domain_.spatial_hash().clear();
  Int i = 0;
  for (const Agent& a : agents_) {
    const bool live = a.state != PhenoState::DEAD;
    if (const bool corpse = cfg_.cell_bio.cdi.enabled && a.timers.death_time >= 0.0
            && (clock_.time - a.timers.death_time) < cfg_.cell_bio.cdi.corpse_persistence;
        live || corpse) {
      domain_.spatial_hash().insert(i, a.x);
    }
    ++i;
  }
}

void Simulation::update_grid_coupling() {
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Agent& a : agents_) {
    if (a.state == PhenoState::DEAD) continue;

    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    domain_.pos_to_grid(a.x, ix, iy, iz);
    a.grid_cell = domain_.cell_index(ix, iy, iz);
  }
}

void Simulation::check_washout() {
  // Agents that have been advected past the luminal boundary
  // (z > domain_hi) are "washed out"
  Real z_max = domain_.hi()[2];

  for (Agent& a : agents_) {
    if (a.state == PhenoState::DEAD) continue;

    // Agents in crypt refugia bypass washout entirely
    if (a.flags.in_crypt) continue;

    if (a.x[2] >= z_max) {
      a.state = PhenoState::DEAD;
      event_ledger_.step_events.outflow_boundary++;
      if (provenance_enabled()) {
        KillProvenanceEvent event;
        event.victim_id = a.identity.tag;
        event.position = a.x;
        event.strain = a.identity.type;
        event.cause = ProvenanceCause::BOUNDARY;
        record_kill_provenance(event);
      }
      lineage_.record_washout(a.identity.tag, a.genome.lineage_id, a.x);
      continue;
    }

    if (cfg_.advection.washout_trap == WashoutTrapMode::IMPOSED) {
      // The imposed variant removes cells before transport reaches the lumen.
      Real gamma = advection_.washout_rate(a.x[2]);
      if (a.mu_realized < gamma) {
        a.state = PhenoState::DEAD;
        event_ledger_.step_events.outflow_washout++;
        if (provenance_enabled()) {
          KillProvenanceEvent event;
          event.victim_id = a.identity.tag;
          event.position = a.x;
          event.strain = a.identity.type;
          event.cause = ProvenanceCause::WASHOUT;
          record_kill_provenance(event);
        }
        lineage_.record_washout(a.identity.tag, a.genome.lineage_id, a.x);
      }
    }
  }
}

void Simulation::crypt_migration(Real dt) {
  if (!cfg_.advection.crypts_enabled) return;

  Real crypt_z = domain_.lo()[2] + cfg_.advection.crypt_depth;
  Real lo_z    = domain_.lo()[2];
  Real epsilon = cfg_.advection.crypt_depth * 0.01;  // small offset above crypt boundary

  // Count agents currently in the crypt for carrying-capacity enforcement
  Int crypt_pop = 0;
  for (const Agent& a : agents_) {
    if (a.state != PhenoState::DEAD && a.flags.in_crypt)
      ++crypt_pop;
  }

  for (Agent& a : agents_) {
    if (a.state == PhenoState::DEAD) continue;

    if (try_exit_crypt(a, dt, crypt_z, epsilon, cfg_.advection.crypt_exit_rate,
                       rng_, crypt_pop)) {
      continue;
    }
    const CryptEntryParams entry_params{
      crypt_z, cfg_.advection.crypt_depth, lo_z,
      cfg_.advection.crypt_entry_rate, cfg_.advection.crypt_carrying_capacity,
    };
    try_enter_crypt(a, dt, entry_params, rng_, crypt_pop);
  }
}

void Simulation::remove_dead_agents() {
  for (Int i = agents_.size() - 1; i >= 0; --i) {
    if (agents_[i].state != PhenoState::DEAD) continue;
    if (cfg_.cell_bio.cdi.enabled && agents_[i].timers.death_time >= 0.0
        && (clock_.time - agents_[i].timers.death_time) < cfg_.cell_bio.cdi.corpse_persistence) {
      continue;
    }
    agents_.remove(i);
  }
}

void Simulation::take_lineage_snapshot() {
  std::vector<std::pair<TagID, TagID>> lineages;
  for (const Agent& a : agents_) {
    if (a.state != PhenoState::DEAD) {
      lineages.emplace_back(a.identity.tag, a.genome.lineage_id);
    }
  }
  lineage_.take_snapshot(clock_.time, lineages);
}

// ---------------------------------------------------------------------------
//  MPI domain decomposition helpers (serialization in agent_transfer.cpp)
// ---------------------------------------------------------------------------

#ifdef GUTIBM_MPI
namespace {

struct MpiSlabPeers {
  Int rank_lo = -1;
  Int rank_hi = -1;
};

struct MpiDistinctTags {
  int lo_send = 0;
  int lo_recv = 0;
  int hi_send = 0;
  int hi_recv = 0;
};

struct MpiPayloadSizes {
  int send_lo = 0;
  int send_hi = 0;
  int recv_lo = 0;
  int recv_hi = 0;
};

struct MpiBufferXfer {
  const std::vector<char>* send_lo = nullptr;
  const std::vector<char>* send_hi = nullptr;
  std::vector<char>* recv_lo = nullptr;
  std::vector<char>* recv_hi = nullptr;
  int recv_lo_size = 0;
  int recv_hi_size = 0;
};

// Non-blocking neighbor exchange. Sequential MPI_Sendrecv(lo) then
// MPI_Sendrecv(hi) deadlocks on a periodic ring when nprocs > 2: every rank
// waits on its lo neighbor first, forming a wait cycle. Isend/Irecv + Waitall
// posts all outstanding ops before blocking.
void mpi_wait_all(int nreq, MPI_Request* reqs) {
  if (nreq > 0) {
    MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
  }
}

void mpi_exchange_sizes_distinct(const MpiSlabPeers& peers,
                                 const MpiDistinctTags& tags,
                                 MpiPayloadSizes& sizes) {
  std::array<MPI_Request, 4> reqs{};
  int nreq = 0;
  if (peers.rank_lo >= 0) {
    MPI_Isend(&sizes.send_lo, 1, MPI_INT, peers.rank_lo, tags.lo_send,
              MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(&sizes.recv_lo, 1, MPI_INT, peers.rank_lo, tags.lo_recv,
              MPI_COMM_WORLD, &reqs[nreq++]);
  }
  if (peers.rank_hi >= 0) {
    MPI_Isend(&sizes.send_hi, 1, MPI_INT, peers.rank_hi, tags.hi_send,
              MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(&sizes.recv_hi, 1, MPI_INT, peers.rank_hi, tags.hi_recv,
              MPI_COMM_WORLD, &reqs[nreq++]);
  }
  mpi_wait_all(nreq, reqs.data());
}

void mpi_exchange_buffers_distinct(const MpiSlabPeers& peers,
                                   const MpiDistinctTags& tags,
                                   const MpiBufferXfer& xfer) {
  std::array<MPI_Request, 4> reqs{};
  int nreq = 0;
  // Stable address for zero-count ops (MPI does not dereference count==0).
  char sink = 0;
  const void* send_lo_ptr =
      xfer.send_lo->empty() ? static_cast<const void*>(&sink)
                            : static_cast<const void*>(xfer.send_lo->data());
  const void* send_hi_ptr =
      xfer.send_hi->empty() ? static_cast<const void*>(&sink)
                            : static_cast<const void*>(xfer.send_hi->data());
  void* recv_lo_ptr =
      xfer.recv_lo_size > 0 ? static_cast<void*>(xfer.recv_lo->data())
                            : static_cast<void*>(&sink);
  void* recv_hi_ptr =
      xfer.recv_hi_size > 0 ? static_cast<void*>(xfer.recv_hi->data())
                            : static_cast<void*>(&sink);

  if (peers.rank_lo >= 0) {
    MPI_Isend(send_lo_ptr, static_cast<int>(xfer.send_lo->size()), MPI_CHAR,
              peers.rank_lo, tags.lo_send, MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(recv_lo_ptr, xfer.recv_lo_size, MPI_CHAR, peers.rank_lo,
              tags.lo_recv, MPI_COMM_WORLD, &reqs[nreq++]);
  }
  if (peers.rank_hi >= 0) {
    MPI_Isend(send_hi_ptr, static_cast<int>(xfer.send_hi->size()), MPI_CHAR,
              peers.rank_hi, tags.hi_send, MPI_COMM_WORLD, &reqs[nreq++]);
    MPI_Irecv(recv_hi_ptr, xfer.recv_hi_size, MPI_CHAR, peers.rank_hi,
              tags.hi_recv, MPI_COMM_WORLD, &reqs[nreq++]);
  }
  mpi_wait_all(nreq, reqs.data());
}

void mpi_exchange_sizes_collapsed(Int neighbor, int tag, MpiPayloadSizes& sizes) {
  std::array<int, 2> sizes_send = {sizes.send_lo, sizes.send_hi};
  std::array<int, 2> sizes_recv = {0, 0};
  MPI_Sendrecv(sizes_send.data(), 2, MPI_INT, neighbor, tag,
               sizes_recv.data(), 2, MPI_INT, neighbor, tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  sizes.recv_lo = sizes_recv[0];
  sizes.recv_hi = sizes_recv[1];
}

void mpi_exchange_buffers_collapsed(Int neighbor, int tag,
                                    const MpiBufferXfer& xfer) {
  std::vector<char> send_buf;
  send_buf.reserve(xfer.send_lo->size() + xfer.send_hi->size());
  send_buf.insert(send_buf.end(), xfer.send_lo->begin(), xfer.send_lo->end());
  send_buf.insert(send_buf.end(), xfer.send_hi->begin(), xfer.send_hi->end());

  xfer.recv_lo->resize(static_cast<size_t>(xfer.recv_lo_size));
  xfer.recv_hi->resize(static_cast<size_t>(xfer.recv_hi_size));
  std::vector<char> recv_buf(static_cast<size_t>(xfer.recv_lo_size + xfer.recv_hi_size));

  MPI_Sendrecv(send_buf.data(), static_cast<int>(send_buf.size()), MPI_CHAR, neighbor, tag,
               recv_buf.data(), static_cast<int>(recv_buf.size()), MPI_CHAR, neighbor, tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  if (xfer.recv_lo_size > 0) {
    std::memcpy(xfer.recv_lo->data(), recv_buf.data(),
                static_cast<size_t>(xfer.recv_lo_size));
  }
  if (xfer.recv_hi_size > 0) {
    std::memcpy(xfer.recv_hi->data(), recv_buf.data() + xfer.recv_lo_size,
                static_cast<size_t>(xfer.recv_hi_size));
  }
}

}  // namespace
#endif

void Simulation::migrate_agents() {
#ifdef GUTIBM_MPI
  if (domain_.nprocs() <= 1) return;

  Int axis = domain_.config().mpi_decomp_axis;
  Int my_rank = domain_.rank();

  // Collect agents that need to migrate to lo/hi neighbors
  std::vector<Agent> send_lo;
  std::vector<Agent> send_hi;
  std::vector<Int> to_remove;

  Int i = 0;
  for (Agent& a : agents_) {
    if (a.state == PhenoState::DEAD) {
      ++i;
      continue;
    }

    const MigrateSide side = classify_migration(a, my_rank, axis, domain_);
    if (side == MigrateSide::None) {
      ++i;
      continue;
    }

    a.identity.owner_rank = domain_.owner_rank(a.x);
    if (side == MigrateSide::Lo) {
      send_lo.push_back(a);
    } else {
      send_hi.push_back(a);
    }
    to_remove.push_back(i);
    ++i;
  }

  // Remove migrated agents (reverse order)
  std::sort(to_remove.rbegin(), to_remove.rend());
  for (Int idx : to_remove) {
    agents_.remove(idx);
  }

  // Serialize
  std::vector<char> buf_send_lo;
  std::vector<char> buf_send_hi;
  agent_transfer_serialize(send_lo, buf_send_lo);
  agent_transfer_serialize(send_hi, buf_send_hi);

  MpiPayloadSizes sizes;
  sizes.send_lo = static_cast<int>(buf_send_lo.size());
  sizes.send_hi = static_cast<int>(buf_send_hi.size());

  const MpiSlabPeers peers{domain_.rank_lo(), domain_.rank_hi()};
  if (domain_.neighbors_collapsed()) {
    mpi_exchange_sizes_collapsed(domain_.rank_lo(), 0, sizes);
  } else {
    mpi_exchange_sizes_distinct(peers, {0, 1, 1, 0}, sizes);
  }

  std::vector<char> buf_recv_lo(sizes.recv_lo);
  std::vector<char> buf_recv_hi(sizes.recv_hi);
  const MpiBufferXfer xfer{
      &buf_send_lo, &buf_send_hi, &buf_recv_lo, &buf_recv_hi,
      sizes.recv_lo, sizes.recv_hi};

  if (domain_.neighbors_collapsed()) {
    mpi_exchange_buffers_collapsed(domain_.rank_lo(), 2, xfer);
  } else {
    mpi_exchange_buffers_distinct(peers, {2, 3, 3, 2}, xfer);
  }

  // Unpack received agents
  auto recv_lo = agent_transfer_deserialize(buf_recv_lo);
  auto recv_hi = agent_transfer_deserialize(buf_recv_hi);

  for (auto& a : recv_lo) {
    a.identity.owner_rank = my_rank;
    a.flags.is_ghost = false;
    agents_.push_back(std::move(a));
  }
  for (auto& a : recv_hi) {
    a.identity.owner_rank = my_rank;
    a.flags.is_ghost = false;
    agents_.push_back(std::move(a));
  }
#endif
}

void Simulation::exchange_ghost_agents() {
#ifdef GUTIBM_MPI
  if (domain_.nprocs() <= 1) return;

  Int axis = domain_.config().mpi_decomp_axis;
  Real gw = domain_.ghost_width();

  // Collect agents near slab boundaries to send as ghosts
  std::vector<Agent> ghost_lo;
  std::vector<Agent> ghost_hi;

  for (const Agent& a : agents_) {
    if (a.state == PhenoState::DEAD) continue;

    Real pos_ax = a.x[axis];
    if (domain_.rank_lo() >= 0 && pos_ax < domain_.local_lo_x() + gw) {
      ghost_lo.push_back(a);
    }
    if (domain_.rank_hi() >= 0 && pos_ax >= domain_.local_hi_x() - gw) {
      ghost_hi.push_back(a);
    }
  }

  // Serialize
  std::vector<char> buf_send_lo;
  std::vector<char> buf_send_hi;
  agent_transfer_serialize(ghost_lo, buf_send_lo);
  agent_transfer_serialize(ghost_hi, buf_send_hi);

  MpiPayloadSizes sizes;
  sizes.send_lo = static_cast<int>(buf_send_lo.size());
  sizes.send_hi = static_cast<int>(buf_send_hi.size());

  const MpiSlabPeers peers{domain_.rank_lo(), domain_.rank_hi()};
  if (domain_.neighbors_collapsed()) {
    mpi_exchange_sizes_collapsed(domain_.rank_lo(), 10, sizes);
  } else {
    mpi_exchange_sizes_distinct(peers, {10, 11, 11, 10}, sizes);
  }

  std::vector<char> buf_recv_lo(sizes.recv_lo);
  std::vector<char> buf_recv_hi(sizes.recv_hi);
  const MpiBufferXfer xfer{
      &buf_send_lo, &buf_send_hi, &buf_recv_lo, &buf_recv_hi,
      sizes.recv_lo, sizes.recv_hi};

  if (domain_.neighbors_collapsed()) {
    mpi_exchange_buffers_collapsed(domain_.rank_lo(), 12, xfer);
  } else {
    mpi_exchange_buffers_distinct(peers, {12, 13, 13, 12}, xfer);
  }

  // Unpack and add as ghost agents
  auto recv_lo = agent_transfer_deserialize(buf_recv_lo);
  auto recv_hi = agent_transfer_deserialize(buf_recv_hi);

  mpi_ghost_.ghost_indices.clear();
  for (auto& a : recv_lo) {
    Int idx = agents_.size();
    mpi_ghost_.ghost_indices.push_back(idx);
    a.flags.is_ghost = true;
    agents_.push_back(std::move(a));
  }
  for (auto& a : recv_hi) {
    Int idx = agents_.size();
    mpi_ghost_.ghost_indices.push_back(idx);
    a.flags.is_ghost = true;
    agents_.push_back(std::move(a));
  }
#endif
}

void Simulation::clear_ghost_agents() {
#ifdef GUTIBM_MPI
  if (mpi_ghost_.ghost_indices.empty()) return;

  // Remove ghosts in reverse index order
  std::sort(mpi_ghost_.ghost_indices.rbegin(), mpi_ghost_.ghost_indices.rend());
  for (Int idx : mpi_ghost_.ghost_indices) {
    if (idx < agents_.size()) {
      agents_.remove(idx);
    }
  }
  mpi_ghost_.ghost_indices.clear();
#endif
}

void Simulation::allreduce_global_stats() {
  // Compute local stats
  Int local_count = 0;
  Int local_sos_count = 0;
  Real local_mu_sum = 0.0;
  Real local_max_abs_mu = 0.0;
  for (const Agent& a : agents_) {
    if (a.state != PhenoState::DEAD) {
      local_count++;
      local_mu_sum += a.mu_realized;
      // Ghosts duplicate an owned agent, so including one in this max is
      // harmless; the global maximum is unchanged.
      local_max_abs_mu = std::max(local_max_abs_mu, std::abs(a.mu_realized));
      if (a.state == PhenoState::SOS_INDUCED) local_sos_count++;
    }
  }

#ifdef GUTIBM_MPI
  if (domain_.nprocs() > 1) {
    const std::array<Real, 3> local_sums = {
        static_cast<Real>(local_count), static_cast<Real>(local_sos_count),
        local_mu_sum};
    std::array<Real, 3> global_sums = {};
    MPI_Allreduce(local_sums.data(), global_sums.data(), 3, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    Real global_max_abs_mu = 0.0;
    MPI_Allreduce(&local_max_abs_mu, &global_max_abs_mu, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);
    mpi_ghost_.stats.global_agent_count =
        static_cast<Int>(global_sums[0]);
    mpi_ghost_.stats.global_sos_count =
        static_cast<Int>(global_sums[1]);
    mpi_ghost_.stats.global_mu_avg =
        global_sums[0] > 0.0 ? global_sums[2] / global_sums[0] : 0.0;
    mpi_ghost_.stats.global_max_abs_mu = global_max_abs_mu;
    return;
  }
#endif

  mpi_ghost_.stats.global_agent_count = local_count;
  mpi_ghost_.stats.global_sos_count = local_sos_count;
  mpi_ghost_.stats.global_mu_avg = local_count > 0 ? local_mu_sum / local_count : 0.0;
  mpi_ghost_.stats.global_max_abs_mu = local_max_abs_mu;
}

Real Simulation::local_O2(const Agent& agent) const {
  if (!cfg_.chem_env.oxygen.enabled) return 0.0;
  Int i_o2 = chem_.find(species::OXYGEN);
  if (i_o2 < 0 || agent.grid_cell < 0) return 0.0;
  return chem_.conc(i_o2, agent.grid_cell);
}

Real Simulation::ros_induction_rate(const Agent& agent) const {
  if (!cfg_.chem_env.oxygen.enabled) return 0.0;
  return cfg_.chem_env.oxygen.k_ROS * local_O2(agent) * std::max(agent.mu_realized, 0.0);
}

Real Simulation::local_nuclease_toxin(const Agent& /*agent*/,
                                      Int agent_index) const {
  return qssa_.sampled_nuclease_conc(agent_index);
}

void Simulation::add_toxin_burst(const ToxinBurstSource& burst) {
  bacteriocin_.bursts.push_back(burst);
}

void Simulation::record_kill_provenance(const KillProvenanceEvent& event) {
  event_ledger_.kill_provenance.push_back(event);
}

void Simulation::prune_toxin_bursts(Real current_time) {
  if (bacteriocin_.bursts.empty()) return;

  std::vector<ToxinBurstSource> kept;
  kept.reserve(bacteriocin_.bursts.size());

  for (const ToxinBurstSource& burst : bacteriocin_.bursts) {
    // A non-positive release timescale contributes no source (see
    // append_burst_sources), so such a burst is dropped rather than retained.
    if (burst.release_tau <= 0.0) continue;
    const Real age = std::max(0.0, current_time - burst.creation_time);
    const Real max_age = 5.0 * burst.release_tau;
    if (age <= max_age) {
      kept.push_back(burst);
    }
  }

  bacteriocin_.bursts.swap(kept);
}

}  // namespace gutibm
