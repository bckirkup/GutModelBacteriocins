/* -----------------------------------------------------------------------
   GutIBM – Simulation engine implementation
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "agent_transfer.h"
#include "fix_registry.h"
#include "plasmid.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <iomanip>
#ifdef GUTIBM_OPENMP
#include <omp.h>
#endif

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {

void Simulation::init(const SimulationConfig& cfg) {
  cfg_ = cfg;
  rng_.seed(cfg.seed);

  // Domain
  domain_.init(cfg.domain);

  // Chemical fields
  chem_.init(domain_, cfg.chemicals);

  // Advection
  advection_.init(cfg.advection, domain_);

  // VBF
  vbf_.init(cfg.vbf, domain_);

  // QSSA solver
  qssa_.init(cfg.qssa, domain_, advection_);

  // Lineage tracker
  lineage_.init(cfg.output_interval);

  // HDF5 output
  hdf5_.init(cfg.hdf5);

  // Register biological fixes via plugin registry
  fixes_ = FixRegistry::create_all(*this, cfg);

  // Initialize fixes
  for (auto& fix : fixes_) {
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
  time_          = 0.0;
  step_count_    = 0;
  next_output_   = 0.0;
  next_snapshot_ = 0.0;

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
              << "  Global agents: " << global_agent_count_ << "\n"
              << "  Chemical species: " << chem_.num_species() << "\n"
              << "  Bio dt: " << cfg.bio_dt << " s\n"
              << "  Adaptive dt: " << (cfg.adaptive_dt_enabled ? "ON" : "OFF")
              << (cfg.adaptive_dt_enabled
                    ? (" [" + std::to_string(cfg.dt_min) + "s, "
                       + std::to_string(cfg.dt_max) + "s]")
                    : "")
              << "\n"
              << "  Total time: " << cfg.total_time << " s\n"
              << std::flush;
  }
}

void Simulation::init_population(const SimulationConfig& cfg) {
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
      a.owner_rank = domain_.rank();

      // Assign plasmids (canonical names + legacy aliases via PlasmidLibrary::find)
      for (const auto& pname : strain.plasmids) {
        const PlasmidEntry* entry = PlasmidLibrary::find(pname);
        if (entry) {
          a.genome.bi_loci.push_back(entry->cluster);
          if (entry->conjugative) {
            a.genome.has_conjugative_plasmid = true;
          }
        } else if (domain_.rank() == 0) {
          std::cerr << "Warning: unknown plasmid '" << pname
                    << "' — agent spawned without BI locus\n";
        }
      }

      // Tag agents spawned inside the crypt zone
      if (advection_.in_crypt_zone(a.x[2])) {
        a.in_crypt = true;
      }

      agents_.push_back(std::move(a));
    }
  }
}

Real Simulation::compute_adaptive_dt() const {
  if (!cfg_.adaptive_dt_enabled) return cfg_.bio_dt;

  Real dt = cfg_.dt_max;

  // Growth rate constraint: mu_max * dt < growth_limit
  Real max_mu = 0.0;
  Int sos_count = 0;
  for (Int i = 0; i < agents_.size(); ++i) {
    if (agents_[i].state == PhenoState::DEAD) continue;
    max_mu = std::max(max_mu, std::abs(agents_[i].mu_realized));
    if (agents_[i].state == PhenoState::SOS_INDUCED) sos_count++;
  }
  if (max_mu > 0) dt = std::min(dt, cfg_.dt_growth_limit / max_mu);

  // SOS cascade constraint: reduce dt during active lysis
  if (sos_count > 5)  dt = std::min(dt, 10.0);
  if (sos_count > 20) dt = std::min(dt, 2.0);

  // Agent density constraint
  Vec3 sz = domain_.size();
  Real volume = sz[0] * sz[1] * sz[2];
  if (volume > 0.0) {
    Real density = static_cast<Real>(agents_.size()) / volume;
    if (density > 1e15) dt = std::min(dt, 10.0);
  }

  // Apply safety factor and bounds
  dt *= cfg_.dt_safety;
  dt = std::clamp(dt, cfg_.dt_min, cfg_.dt_max);
  return dt;
}

void Simulation::run() {
  int rank = domain_.rank();

  // Initial output
  hdf5_.write_step(*this, 0, 0.0);
  if (rank == 0) {
    take_lineage_snapshot();
  }

  while (time_ < cfg_.total_time) {
    Real dt = compute_adaptive_dt();

    // Clamp so we don't overshoot total_time
    if (time_ + dt > cfg_.total_time) dt = cfg_.total_time - time_;

    step(dt);

    // Periodic output
    if (time_ >= next_output_) {
      hdf5_.write_step(*this, step_count_, time_);
      if (rank == 0) {
        std::cout << "Step " << step_count_
                  << "  t=" << time_ << "s"
                  << "  dt=" << std::setprecision(3) << dt << "s"
                  << "  global_agents=" << global_agent_count_
                  << "  local_agents=" << agents_.size()
                  << "  mu_avg=" << global_mu_avg_
                  << "\n" << std::flush;
      }
      next_output_ += cfg_.output_interval;
    }

    // Lineage snapshots
    if (time_ >= next_snapshot_) {
      take_lineage_snapshot();
      next_snapshot_ += cfg_.output_interval;
    }
  }

  // Final output
  hdf5_.write_step(*this, step_count_, time_);
  hdf5_.finalize();

  if (rank == 0) {
    Real retention = lineage_.resident_retention(cfg_.total_time * 0.5);
    std::cout << "\nSimulation complete.\n"
              << "  Final global agents: " << global_agent_count_ << "\n"
              << "  Steps taken: " << step_count_ << "\n"
              << "  Resident retention: " << retention * 100.0 << "%\n"
              << "  Dominant lineage: " << lineage_.dominant_lineage() << "\n"
              << std::flush;
  }
}

void Simulation::step(Real dt) {
  // Update advection time for peristaltic oscillation
  advection_.set_time(time_);

  // Pre-step: clear ghosts from previous step
  clear_ghost_agents();
  chem_.zero_reactions();

  // Exchange ghost agents for cross-boundary neighbor queries
  exchange_ghost_agents();

  rebuild_spatial_hash();
  update_grid_coupling();

  for (auto& fix : fixes_) {
    fix->pre_step(dt);
  }

  // 1. Biology module (uses ghost agents for neighbor interactions)
  module_biology(dt);

  // Clear ghosts before physics to avoid moving them
  clear_ghost_agents();

  // 2. Chemistry module (QSSA, instantaneous)
  module_chemistry(dt);

  // 3. Physics module (advection + mechanics)
  module_physics(dt);

  // Post-step
  for (auto& fix : fixes_) {
    fix->post_step(dt);
  }

  // Migrate agents that crossed slab boundaries
  migrate_agents();

  // Cleanup
  check_washout();
  remove_dead_agents();

  // Compute global statistics
  allreduce_global_stats();

  time_ += dt;
  step_count_++;
}

void Simulation::module_biology(Real dt) {
  for (auto& fix : fixes_) {
    fix->compute(dt);
  }
}

void Simulation::module_chemistry(Real dt) {
  // QSSA: compute steady-state toxin field via Green's functions
  Int i_tox = chem_.find("bacteriocin");
  if (i_tox >= 0) {
    qssa_.solve_bacteriocin_field(agents_, chem_, i_tox);
  }

  // Nutrient depletion
  qssa_.solve_nutrient_depletion(agents_, chem_);

  // VBF coupling (nutrient sink/source)
  vbf_.apply_nutrient_coupling(chem_, domain_, dt);

  // Apply reactions to concentrations
  for (Int s = 0; s < chem_.num_species(); ++s) {
    #ifdef GUTIBM_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (Int c = 0; c < chem_.ncells(); ++c) {
      chem_.conc(s, c) += chem_.reac(s, c) * dt;
      chem_.conc(s, c) = std::max(chem_.conc(s, c), 0.0);
    }
  }

  // Boundary conditions
  chem_.apply_boundaries(domain_);
}

void Simulation::module_physics(Real dt) {
  // Crypt migration (stochastic entry/exit) before advection
  crypt_migration(dt);

  // Advection pass: each agent independent (read-only on fields)
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Int i = 0; i < agents_.size(); ++i) {
    Agent& a = agents_[i];
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

    // PBC / boundary
    domain_.apply_pbc(a.x);
  }

  // Mechanical repulsion handled by FixMechanics (registered as a fix)
  for (auto& fix : fixes_) {
    if (fix->name() == "mechanics") {
      fix->compute(dt);
      break;
    }
  }
}

void Simulation::rebuild_spatial_hash() {
  domain_.spatial_hash().clear();
  for (Int i = 0; i < agents_.size(); ++i) {
    if (agents_[i].state != PhenoState::DEAD) {
      domain_.spatial_hash().insert(i, agents_[i].x);
    }
  }
}

void Simulation::update_grid_coupling() {
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Int i = 0; i < agents_.size(); ++i) {
    Agent& a = agents_[i];
    if (a.state == PhenoState::DEAD) continue;

    Int ix, iy, iz;
    domain_.pos_to_grid(a.x, ix, iy, iz);
    a.grid_cell = domain_.cell_index(ix, iy, iz);
  }
}

void Simulation::check_washout() {
  // Agents that have been advected past the luminal boundary
  // (z > domain_hi) are "washed out"
  Real z_max = domain_.hi()[2];

  for (Int i = 0; i < agents_.size(); ++i) {
    Agent& a = agents_[i];
    if (a.state == PhenoState::DEAD) continue;

    // Agents in crypt refugia bypass washout entirely
    if (a.in_crypt) continue;

    if (a.x[2] >= z_max) {
      a.state = PhenoState::DEAD;
      lineage_.record_washout(a.tag, a.genome.lineage_id, a.x);
      continue;
    }

    // Combinatorial Washout Trap: mu_realized < gamma_flow (VADI / EARI)
    Real gamma = advection_.washout_rate(a.x[2]);
    if (a.mu_realized < gamma) {
      a.state = PhenoState::DEAD;
      lineage_.record_washout(a.tag, a.genome.lineage_id, a.x);
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
  for (Int i = 0; i < agents_.size(); ++i) {
    if (agents_[i].state != PhenoState::DEAD && agents_[i].in_crypt)
      ++crypt_pop;
  }

  for (Int i = 0; i < agents_.size(); ++i) {
    Agent& a = agents_[i];
    if (a.state == PhenoState::DEAD) continue;

    if (a.in_crypt) {
      // Stochastic exit from crypt
      Real p_exit = 1.0 - std::exp(-cfg_.advection.crypt_exit_rate * dt);
      if (rng_.bernoulli(p_exit)) {
        a.x[2] = crypt_z + epsilon;
        a.in_crypt = false;
        --crypt_pop;
      }
    } else {
      // Only agents near the crypt boundary can enter
      if (a.x[2] < crypt_z + cfg_.advection.crypt_depth) {
        if (crypt_pop >= cfg_.advection.crypt_carrying_capacity) continue;
        Real p_entry = 1.0 - std::exp(-cfg_.advection.crypt_entry_rate * dt);
        if (rng_.bernoulli(p_entry)) {
          a.x[2] = rng_.uniform(lo_z, crypt_z);
          a.in_crypt = true;
          ++crypt_pop;
        }
      }
    }
  }
}

void Simulation::remove_dead_agents() {
  for (Int i = agents_.size() - 1; i >= 0; --i) {
    if (agents_[i].state == PhenoState::DEAD) {
      agents_.remove(i);
    }
  }
}

void Simulation::take_lineage_snapshot() {
  std::vector<std::pair<TagID, TagID>> lineages;
  for (Int i = 0; i < agents_.size(); ++i) {
    const Agent& a = agents_[i];
    if (a.state != PhenoState::DEAD) {
      lineages.emplace_back(a.tag, a.genome.lineage_id);
    }
  }
  lineage_.take_snapshot(time_, lineages);
}

// ---------------------------------------------------------------------------
//  MPI domain decomposition helpers (serialization in agent_transfer.cpp)
// ---------------------------------------------------------------------------

void Simulation::migrate_agents() {
#ifdef GUTIBM_MPI
  if (domain_.nprocs() <= 1) return;

  Int axis = domain_.config().mpi_decomp_axis;
  Int my_rank = domain_.rank();

  // Collect agents that need to migrate to lo/hi neighbors
  std::vector<Agent> send_lo, send_hi;
  std::vector<Int> to_remove;

  for (Int i = 0; i < agents_.size(); ++i) {
    Agent& a = agents_[i];
    if (a.state == PhenoState::DEAD) continue;

    Int dest = domain_.owner_rank(a.x);
    if (dest != my_rank) {
      a.owner_rank = dest;
      if (dest == domain_.rank_lo()) {
        send_lo.push_back(a);
      } else if (dest == domain_.rank_hi()) {
        send_hi.push_back(a);
      } else {
        // Agent jumped more than one slab (rare, large dt)
        // Send to whichever neighbor is closer in rank
        if (a.x[axis] < domain_.local_lo_x()) {
          send_lo.push_back(a);
        } else {
          send_hi.push_back(a);
        }
      }
      to_remove.push_back(i);
    }
  }

  // Remove migrated agents (reverse order)
  std::sort(to_remove.rbegin(), to_remove.rend());
  for (Int idx : to_remove) {
    agents_.remove(idx);
  }

  // Serialize
  std::vector<char> buf_send_lo, buf_send_hi;
  agent_transfer_serialize(send_lo, buf_send_lo);
  agent_transfer_serialize(send_hi, buf_send_hi);

  // Exchange sizes with neighbors
  int sz_send_lo = static_cast<int>(buf_send_lo.size());
  int sz_send_hi = static_cast<int>(buf_send_hi.size());
  int sz_recv_lo = 0, sz_recv_hi = 0;

  // Sendrecv with lo neighbor
  if (domain_.rank_lo() >= 0) {
    MPI_Sendrecv(&sz_send_lo, 1, MPI_INT, domain_.rank_lo(), 0,
                 &sz_recv_lo, 1, MPI_INT, domain_.rank_lo(), 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  // Sendrecv with hi neighbor
  if (domain_.rank_hi() >= 0) {
    MPI_Sendrecv(&sz_send_hi, 1, MPI_INT, domain_.rank_hi(), 1,
                 &sz_recv_hi, 1, MPI_INT, domain_.rank_hi(), 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  // Exchange agent data
  std::vector<char> buf_recv_lo(sz_recv_lo), buf_recv_hi(sz_recv_hi);

  if (domain_.rank_lo() >= 0) {
    MPI_Sendrecv(buf_send_lo.data(), sz_send_lo, MPI_CHAR,
                 domain_.rank_lo(), 2,
                 buf_recv_lo.data(), sz_recv_lo, MPI_CHAR,
                 domain_.rank_lo(), 3,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (domain_.rank_hi() >= 0) {
    MPI_Sendrecv(buf_send_hi.data(), sz_send_hi, MPI_CHAR,
                 domain_.rank_hi(), 3,
                 buf_recv_hi.data(), sz_recv_hi, MPI_CHAR,
                 domain_.rank_hi(), 2,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  // Unpack received agents
  auto recv_lo = agent_transfer_deserialize(buf_recv_lo);
  auto recv_hi = agent_transfer_deserialize(buf_recv_hi);

  for (auto& a : recv_lo) {
    a.owner_rank = my_rank;
    agents_.push_back(std::move(a));
  }
  for (auto& a : recv_hi) {
    a.owner_rank = my_rank;
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
  std::vector<Agent> ghost_lo, ghost_hi;

  for (Int i = 0; i < agents_.size(); ++i) {
    const Agent& a = agents_[i];
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
  std::vector<char> buf_send_lo, buf_send_hi;
  agent_transfer_serialize(ghost_lo, buf_send_lo);
  agent_transfer_serialize(ghost_hi, buf_send_hi);

  // Exchange sizes
  int sz_send_lo = static_cast<int>(buf_send_lo.size());
  int sz_send_hi = static_cast<int>(buf_send_hi.size());
  int sz_recv_lo = 0, sz_recv_hi = 0;

  if (domain_.rank_lo() >= 0) {
    MPI_Sendrecv(&sz_send_lo, 1, MPI_INT, domain_.rank_lo(), 10,
                 &sz_recv_lo, 1, MPI_INT, domain_.rank_lo(), 11,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (domain_.rank_hi() >= 0) {
    MPI_Sendrecv(&sz_send_hi, 1, MPI_INT, domain_.rank_hi(), 11,
                 &sz_recv_hi, 1, MPI_INT, domain_.rank_hi(), 10,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  // Exchange data
  std::vector<char> buf_recv_lo(sz_recv_lo), buf_recv_hi(sz_recv_hi);

  if (domain_.rank_lo() >= 0) {
    MPI_Sendrecv(buf_send_lo.data(), sz_send_lo, MPI_CHAR,
                 domain_.rank_lo(), 12,
                 buf_recv_lo.data(), sz_recv_lo, MPI_CHAR,
                 domain_.rank_lo(), 13,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (domain_.rank_hi() >= 0) {
    MPI_Sendrecv(buf_send_hi.data(), sz_send_hi, MPI_CHAR,
                 domain_.rank_hi(), 13,
                 buf_recv_hi.data(), sz_recv_hi, MPI_CHAR,
                 domain_.rank_hi(), 12,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
  for (Int i = 0; i < agents_.size(); ++i) {
    if (agents_[i].state != PhenoState::DEAD) {
      local_count++;
      local_mu_sum += agents_[i].mu_realized;
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
    global_agent_count_ = global_count;
    global_mu_avg_ = global_count > 0 ? global_mu_sum / global_count : 0.0;
    return;
  }
#endif

  global_agent_count_ = local_count;
  global_mu_avg_ = local_count > 0 ? local_mu_sum / local_count : 0.0;
}

}  // namespace gutibm
