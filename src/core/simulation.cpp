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
#include <cmath>
#include <cstring>
#include <numeric>
#include <iomanip>
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

void assign_plasmids(Agent& agent, const std::vector<std::string>& plasmids, int rank) {
  for (const auto& pname : plasmids) {
    const PlasmidEntry* entry = PlasmidLibrary::find(pname);
    if (!entry) {
      if (rank == 0) {
        std::cerr << "Warning: unknown plasmid '" << pname
                  << "' — agent spawned without BI locus\n";
      }
      continue;
    }
    agent.genome.bi_loci.push_back(entry->cluster);
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
    if (static_cast<Int>(values.size()) != chem.ncells()) {
      throw SimulationError("checkpoint grid size mismatch for species: " + name);
    }
    Int c = 0;
    for (const Real val : values) {
      chem.conc(spec, c++) = val;
    }
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

}  // namespace

void Simulation::init(const SimulationConfig& cfg) {
  cfg_ = cfg;
  InputParser::finalize_config(cfg_);
  rng_.seed(cfg_.seed);

  // Domain
  domain_.init(cfg.domain);

  // Chemical fields
  chem_.init(domain_, cfg_.chemicals);

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

  // GPU acceleration
  gpu_set_config(cfg.gpu);
  gpu_init_for_rank(domain_.rank(), domain_.nprocs());
  gpu_active_ = gpu_runtime_enabled();
  if (gpu_active_) {
    chem_gpu_.init(chem_);
    agents_gpu_.sync_from_host(agents_);
  }

  // Register biological fixes via plugin registry
  fixes_ = FixRegistry::create_all(*this, cfg);

  // Initialize fixes
  for (const auto& fix : fixes_) {
    fix->init();
  }

  // Create initial population (only agents local to this rank)
  init_population(cfg);

  // Initial coupling
  rebuild_spatial_hash();
  update_grid_coupling();

  // Initial global stats
  allreduce_global_stats();

  // Timers
  clock_.time          = 0.0;
  clock_.step_count    = 0;
  clock_.next_output   = 0.0;
  clock_.next_snapshot = 0.0;

  int rank = domain_.rank();
  if (rank == 0) {
    std::cout << "GutIBM initialized:\n"
              << "  Domain: " << domain_.nx() << "x" << domain_.ny()
              << "x" << domain_.nz() << " cells"
              << " (" << domain_.ncells() << " total)\n"
              << "  MPI ranks: " << domain_.nprocs() << "\n"
              << "  Slab [" << domain_.local_lo_x() << ", "
              << domain_.local_hi_x() << ") m\n"
              << "  Local agents: " << agents_.size()
              << "  Global agents: " << mpi_stats_.global_agent_count << "\n"
              << "  Chemical species: " << chem_.num_species() << "\n";
    const std::string adaptive_dt_status = cfg.adaptive_dt.enabled
        ? std::format(" [{}s, {}s]", cfg.adaptive_dt.min, cfg.adaptive_dt.max)
        : "";
    const std::string gpu_status = gpu_active_
        ? std::format(" (device {})", gpu_device().device_id())
        : std::format(" (gpu_enabled={}, device_id={})",
                      cfg.gpu.enabled ? "true" : "false",
                      cfg.gpu.device_id);
    std::cout << "  Bio dt: " << cfg.time.bio_dt << " s\n"
              << "  Adaptive dt: " << (cfg.adaptive_dt.enabled ? "ON" : "OFF")
              << adaptive_dt_status << "\n"
              << "  Total time: " << cfg.time.total_time << " s\n"
              << "  GPU: " << (gpu_active_ ? "ON" : "OFF") << gpu_status << "\n";
    if (!gpu_active_ && cfg.gpu.enabled) {
      std::cerr << "  GPU requested (gpu_enabled) but inactive: "
                << gpu_fallback_reason(cfg.gpu) << "\n";
    }
    std::cout << std::flush;
  }
}

void Simulation::init_population(const SimulationConfig& cfg) {
  agents_.configure_tags(AgentPool::first_tag_for_rank(domain_.rank(), domain_.nprocs()),
                           AgentPool::tag_stride(domain_.nprocs()));

  for (const auto& strain : cfg.initial_strains) {
    for (Int i = 0; i < strain.count; ++i) {
      Vec3 pos = {
        rng_.uniform(domain_.lo()[0], domain_.hi()[0]),
        rng_.uniform(domain_.lo()[1], domain_.hi()[1]),
        rng_.uniform(domain_.lo()[2], domain_.hi()[2] * 0.5)  // near epithelium
      };

      // Only keep agents that belong to this rank's slab
      if (!domain_.is_local(pos)) continue;

      Agent a = Agent::create_default(agents_.next_tag(), strain.type,
                                       pos, strain.mu_max);
      a.identity.owner_rank = domain_.rank();

      assign_plasmids(a, strain.plasmids, domain_.rank());
      a.genome.cdi_type = strain.cdi_type;
      a.genome.cdi_immunity = strain.cdi_immunity;
      tag_crypt_resident(a, advection_);
      if (cfg_.cell_bio.motility.enabled) {
        FixMotility::init_agent_motility(a, cfg_.cell_bio.motility, rng_);
      }

      agents_.push_back(std::move(a));
    }
  }
}

void Simulation::init_from_checkpoint(const SimulationConfig& cfg,
                                      const std::string& h5_file,
                                      const std::string& step) {
  cfg_ = cfg;
  InputParser::finalize_config(cfg_);
  rng_.seed(cfg_.seed);

  domain_.init(cfg_.domain);
  chem_.init(domain_, cfg_.chemicals);
  advection_.init(cfg.advection, domain_);
  vbf_.init(cfg.vbf, domain_);
  qssa_.init(cfg.qssa, domain_, advection_);
  lineage_.init(cfg.time.output_interval);
  hdf5_.init(cfg.hdf5, domain_);
  fixes_ = FixRegistry::create_all(*this, cfg);
  for (const auto& fix : fixes_) {
    fix->init();
  }

  gpu_set_config(cfg.gpu);
  gpu_init_for_rank(domain_.rank(), domain_.nprocs());
  gpu_active_ = gpu_runtime_enabled();
  if (gpu_active_) {
    chem_gpu_.init(chem_);
    agents_gpu_.sync_from_host(agents_);
  }

#ifndef GUTIBM_HDF5
  (void)h5_file;
  (void)step;
  throw SimulationError("checkpoint restart requires HDF5 support");
#else
  HDF5CheckpointSnapshot snap = HDF5Reader::load_snapshot(h5_file, step);
  apply_checkpoint_snapshot(snap);

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
              << "  Global agents: " << mpi_stats_.global_agent_count << "\n"
              << "  Local agents: " << agents_.size() << "\n"
              << std::flush;
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
    a.mu_max       = std::max(mu_guess, atoms.mu[i]);

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

    tag_crypt_resident(a, advection_);

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
  schedule_output_from_time(clock_.time, cfg_.time.output_interval, clock_.next_output, clock_.next_snapshot);
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
  Real max_mu = 0.0;
  Int sos_count = 0;
  for (const Agent& a : agents_) {
    if (a.state == PhenoState::DEAD) continue;
    max_mu = std::max(max_mu, std::abs(a.mu_realized));
    if (a.state == PhenoState::SOS_INDUCED) sos_count++;
  }
  if (max_mu > 0) dt = std::min(dt, cfg_.adaptive_dt.growth_limit / max_mu);

  // SOS cascade constraint: reduce dt during active lysis
  if (sos_count > 5)  dt = std::min(dt, 10.0);
  if (sos_count > 20) dt = std::min(dt, 2.0);

  // Agent density constraint
  Vec3 sz = domain_.size();
  if (Real volume = sz[0] * sz[1] * sz[2]; volume > 0.0) {
    Real density = static_cast<Real>(agents_.size()) / volume;
    if (density > 1e15) dt = std::min(dt, 10.0);
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
            << " mpi_reaction_reduce_s=" << step_profile_.mpi_reaction_reduce_s * inv
            << " hdf5_s=" << step_profile_.hdf5_s * inv
            << " total_s=" << total
            << " agents=" << mpi_stats_.global_agent_count
            << " ranks=" << domain_.nprocs()
            << "\n"
            << std::flush;
}

namespace {
constexpr Int kPopulationStopThreshold = 1;
}  // namespace

void Simulation::run() {
  int rank = domain_.rank();
  Real last_dt = cfg_.time.bio_dt;
  bool stopped_for_population = false;

  // Initial snapshot (step 0, pre-biology)
  if (hdf5_.is_enabled()) {
    hdf5_.write_step(*this, 0, 0.0, last_dt);
  }
  if (rank == 0) {
    take_lineage_snapshot();
  }

  if (mpi_stats_.global_agent_count <= kPopulationStopThreshold) {
    stopped_for_population = true;
    if (rank == 0) {
      std::cout << "Population stop: " << mpi_stats_.global_agent_count
                << " cell(s) — ending simulation.\n"
                << std::flush;
    }
  }

  while (!stopped_for_population && clock_.time < cfg_.time.total_time) {
    if (gutibm_stop_requested()) break;

    Real dt = compute_adaptive_dt();
    last_dt = dt;

    // Clamp so we don't overshoot total_time
    if (clock_.time + dt > cfg_.time.total_time) dt = cfg_.time.total_time - clock_.time;

    step(dt);

    // HDF5 cadence is controlled solely by hdf5.schedule.* (per-layer intervals).
    if (hdf5_.is_enabled()) {
      const auto hdf5_t0 = std::chrono::steady_clock::now();
      hdf5_.write_step(*this, clock_.step_count, clock_.time, last_dt);
      if (cfg_.profile_steps) {
        const auto hdf5_t1 = std::chrono::steady_clock::now();
        step_profile_.hdf5_s += std::chrono::duration<double>(hdf5_t1 - hdf5_t0).count();
      }
    }

    // Console progress and in-memory lineage snapshots use output_interval (seconds).
    if (clock_.time >= clock_.next_output) {
      if (rank == 0) {
        std::cout << "Step " << clock_.step_count
                  << "  t=" << clock_.time << "s"
                  << "  dt=" << std::setprecision(3) << dt << "s"
                  << "  global_agents=" << mpi_stats_.global_agent_count
                  << "  local_agents=" << agents_.size()
                  << "  mu_avg=" << mpi_stats_.global_mu_avg
                  << "\n" << std::flush;
      }
      clock_.next_output += cfg_.time.output_interval;
    }

    // Lineage snapshots
    if (clock_.time >= clock_.next_snapshot) {
      take_lineage_snapshot();
      clock_.next_snapshot += cfg_.time.output_interval;
    }

    if (mpi_stats_.global_agent_count <= kPopulationStopThreshold) {
      stopped_for_population = true;
      if (rank == 0) {
        std::cout << "Population stop: " << mpi_stats_.global_agent_count
                  << " cell(s) — ending simulation.\n"
                  << std::flush;
      }
      break;
    }

    // Spec 5 §4 — Dysbiosis safety net: halt if density leaves the homeostatic
    // regime the model is calibrated for.
    if (cfg_.dysbiosis_threshold > 0.0) {
      const Vec3 lo = domain_.lo();
      const Vec3 hi = domain_.hi();
      const Real vol_m3 = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
      const Real vol_mL = vol_m3 * 1.0e9;  // 1 m^3 = 1e9 mL
      if (const Real density = vol_mL > 0.0
              ? static_cast<Real>(mpi_stats_.global_agent_count) / vol_mL
              : 0.0;
          density > cfg_.dysbiosis_threshold) {
        stopped_for_population = true;
        if (rank == 0) {
          std::cerr << "DYSBIOSIS THRESHOLD EXCEEDED: " << density
                    << " cells/mL > " << cfg_.dysbiosis_threshold
                    << " cells/mL — halting simulation.\n"
                    << std::flush;
        }
        break;
      }
    }
  }

  hdf5_.finalize();

  if (rank == 0) {
    Real retention = lineage_.resident_retention(cfg_.time.total_time * 0.5);
    std::cout << "\nSimulation complete.\n"
              << "  Final global agents: " << mpi_stats_.global_agent_count << "\n"
              << "  Steps taken: " << clock_.step_count << "\n"
              << "  Resident retention: " << retention * 100.0 << "%\n"
              << "  Dominant lineage: " << lineage_.dominant_lineage() << "\n"
              << std::flush;
    if (cfg_.profile_steps) {
      print_step_profile();
    }
  }
}

void Simulation::step(Real dt) {
  bacteriocin_fields_current_ = false;
  if (cfg_.profile_steps) {
    gpu_transfer_profile_set_enabled(gpu_active_);
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
  if (gpu_active_) {
    chem_gpu_.zero_reactions_on_device();
  }
  chem_.zero_reactions();

  // Exchange ghost agents for cross-boundary neighbor queries
  exchange_ghost_agents();
  profiler.lap(step_profile_.ghost_exchange_s);

  if (gpu_active_) {
    agents_gpu_.sync_from_host(agents_);
    chem_gpu_.sync_to_device(chem_);
#ifdef GUTIBM_CUDA
    gpu::launch_grid_coupling_kernel(
        agents_gpu_.x(), agents_gpu_.y(), agents_gpu_.z(),
        agents_gpu_.grid_cell(), agents_gpu_.state(),
        domain_.lo()[0], domain_.lo()[1], domain_.lo()[2], domain_.dx(),
        domain_.nx(), domain_.ny(), domain_.nz(),
        agents_.size(), gpu_compute_stream());
    gpu_sync_compute();
    gpu_check_error("grid_coupling_kernel");
    agents_gpu_.sync_to_host(agents_);
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

  if (gpu_active_) {
    agents_gpu_.sync_from_host(agents_);
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
    gpu_transfer_profile_reset();
  }

  clock_.time += dt;
  clock_.step_count++;
}

void Simulation::update_bacteriocin_fields() {
  if (bacteriocin_fields_current_) return;

  prune_toxin_bursts(clock_.time);

  ChemicalFieldGpu* chem_gpu_ptr = gpu_active_ ? &chem_gpu_ : nullptr;
  qssa_.solve_all_bacteriocin_fields(agents_, toxin_bursts_, clock_.time,
                                      cfg_.chem_env.protease, advection_, chem_,
                                      chem_gpu_ptr);
  bacteriocin_fields_current_ = true;
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
      .gpu_active = gpu_active_,
      .agents_gpu = agents_gpu_,
      .chem_gpu = chem_gpu_,
      .chem = chem_,
      .domain = domain_,
      .vbf = vbf_,
      .qssa = qssa_,
      .agents = agents_,
      .oxygen = cfg_.chem_env.oxygen,
      .acetate = cfg_.chem_env.acetate,
      .mucin = cfg_.chem_env.mucin,
      .num_agents = static_cast<Int>(agents_.size()),
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

    if (cfg_.cell_bio.motility.enabled && !a.motility.is_stopped) {
      a.x[0] += a.motility.swim_direction[0] * cfg_.cell_bio.motility.swim_speed * dt;
      a.x[1] += a.motility.swim_direction[1] * cfg_.cell_bio.motility.swim_speed * dt;
      a.x[2] += a.motility.swim_direction[2] * cfg_.cell_bio.motility.swim_speed * dt;
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
      step_events_.boundary_deaths++;
      lineage_.record_washout(a.identity.tag, a.genome.lineage_id, a.x);
      continue;
    }

    // Combinatorial Washout Trap: mu_realized < gamma_flow (VADI / EARI)
    Real gamma = advection_.washout_rate(a.x[2]);
    if (a.mu_realized < gamma) {
      a.state = PhenoState::DEAD;
      step_events_.washout_deaths++;
      lineage_.record_washout(a.identity.tag, a.genome.lineage_id, a.x);
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

void mpi_exchange_sizes_distinct(const MpiSlabPeers& peers,
                                 const MpiDistinctTags& tags,
                                 MpiPayloadSizes& sizes) {
  if (peers.rank_lo >= 0) {
    MPI_Sendrecv(&sizes.send_lo, 1, MPI_INT, peers.rank_lo, tags.lo_send,
                 &sizes.recv_lo, 1, MPI_INT, peers.rank_lo, tags.lo_recv,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (peers.rank_hi >= 0) {
    MPI_Sendrecv(&sizes.send_hi, 1, MPI_INT, peers.rank_hi, tags.hi_send,
                 &sizes.recv_hi, 1, MPI_INT, peers.rank_hi, tags.hi_recv,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

void mpi_exchange_buffers_distinct(const MpiSlabPeers& peers,
                                   const MpiDistinctTags& tags,
                                   const MpiBufferXfer& xfer) {
  if (peers.rank_lo >= 0) {
    MPI_Sendrecv(xfer.send_lo->data(), static_cast<int>(xfer.send_lo->size()), MPI_CHAR,
                 peers.rank_lo, tags.lo_send,
                 xfer.recv_lo->data(), xfer.recv_lo_size, MPI_CHAR,
                 peers.rank_lo, tags.lo_recv,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (peers.rank_hi >= 0) {
    MPI_Sendrecv(xfer.send_hi->data(), static_cast<int>(xfer.send_hi->size()), MPI_CHAR,
                 peers.rank_hi, tags.hi_send,
                 xfer.recv_hi->data(), xfer.recv_hi_size, MPI_CHAR,
                 peers.rank_hi, tags.hi_recv,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
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
    agents_.push_back(std::move(a));
  }
  for (auto& a : recv_hi) {
    a.identity.owner_rank = my_rank;
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

  ghost_indices_.clear();
  for (auto& a : recv_lo) {
    Int idx = agents_.size();
    ghost_indices_.push_back(idx);
    agents_.push_back(std::move(a));
  }
  for (auto& a : recv_hi) {
    Int idx = agents_.size();
    ghost_indices_.push_back(idx);
    agents_.push_back(std::move(a));
  }
#endif
}

void Simulation::clear_ghost_agents() {
#ifdef GUTIBM_MPI
  if (ghost_indices_.empty()) return;

  // Remove ghosts in reverse index order
  std::sort(ghost_indices_.rbegin(), ghost_indices_.rend());
  for (Int idx : ghost_indices_) {
    if (idx < agents_.size()) {
      agents_.remove(idx);
    }
  }
  ghost_indices_.clear();
#endif
}

void Simulation::allreduce_global_stats() {
  // Compute local stats
  Int local_count = 0;
  Real local_mu_sum = 0.0;
  for (const Agent& a : agents_) {
    if (a.state != PhenoState::DEAD) {
      local_count++;
      local_mu_sum += a.mu_realized;
    }
  }

#ifdef GUTIBM_MPI
  if (domain_.nprocs() > 1) {
    Int global_count = 0;
    Real global_mu_sum = 0.0;
    MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_mu_sum, &global_mu_sum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    mpi_stats_.global_agent_count = global_count;
    mpi_stats_.global_mu_avg = global_count > 0 ? global_mu_sum / global_count : 0.0;
    return;
  }
#endif

  mpi_stats_.global_agent_count = local_count;
  mpi_stats_.global_mu_avg = local_count > 0 ? local_mu_sum / local_count : 0.0;
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

Real Simulation::local_nuclease_toxin(const Agent& agent) const {
  Int i_btuB = chem_.find(species::BACTERIOCIN_BTUB);
  if (i_btuB < 0 || agent.grid_cell < 0) return 0.0;
  return chem_.conc(i_btuB, agent.grid_cell);
}

void Simulation::add_toxin_burst(const ToxinBurstSource& burst) {
  toxin_bursts_.push_back(burst);
}

void Simulation::prune_toxin_bursts(Real current_time) {
  if (toxin_bursts_.empty()) return;

  std::vector<ToxinBurstSource> kept;
  kept.reserve(toxin_bursts_.size());

  for (const ToxinBurstSource& burst : toxin_bursts_) {
    if (!cfg_.chem_env.protease.enabled || burst.decay_rate <= 0.0) {
      kept.push_back(burst);
      continue;
    }
    const Real age = std::max(0.0, current_time - burst.creation_time);
    const Real max_age = 5.0 / burst.decay_rate;
    if (age <= max_age) {
      kept.push_back(burst);
    }
  }

  toxin_bursts_.swap(kept);
}

}  // namespace gutibm
