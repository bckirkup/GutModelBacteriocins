/* -----------------------------------------------------------------------
   GutIBM – Simulation engine implementation
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "fix_metabolism.h"
#include "fix_receptor.h"
#include "fix_bacteriocin.h"
#include "fix_conjugation.h"
#include "fix_mutation.h"
#include "fix_mechanics.h"
#include "plasmid.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
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

  // Register biological fixes
  fixes_.push_back(std::make_unique<FixMetabolism>(*this, cfg.metabolism));
  fixes_.push_back(std::make_unique<FixBacteriocin>(*this, cfg.bacteriocin));
  fixes_.push_back(std::make_unique<FixReceptor>(*this, cfg.receptor));
  fixes_.push_back(std::make_unique<FixConjugation>(*this, cfg.conjugation));
  fixes_.push_back(std::make_unique<FixMutation>(*this, cfg.mutation));
  fixes_.push_back(std::make_unique<FixMechanics>(*this, cfg.mechanics));

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
              << "  Total time: " << cfg.total_time << " s\n"
              << std::flush;
  }
}

void Simulation::init_population(const SimulationConfig& cfg) {
  const auto& lib = PlasmidLibrary::entries();

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

      // Assign plasmids
      for (const auto& pname : strain.plasmids) {
        for (const auto& entry : lib) {
          if (entry.name == pname) {
            a.genome.bi_loci.push_back(entry.cluster);
            if (entry.conjugative) {
              a.genome.has_conjugative_plasmid = true;
            }
            break;
          }
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

void Simulation::run() {
  Real dt = cfg_.bio_dt;
  Int total_steps = static_cast<Int>(std::ceil(cfg_.total_time / dt));

  int rank = domain_.rank();

  // Initial output
  hdf5_.write_step(*this, 0, 0.0);
  if (rank == 0) {
    take_lineage_snapshot();
  }

  for (Int s = 1; s <= total_steps; ++s) {
    step(dt);

    // Periodic output
    if (time_ >= next_output_) {
      hdf5_.write_step(*this, s, time_);
      if (rank == 0) {
        std::cout << "Step " << s << "/" << total_steps
                  << "  t=" << time_ << "s"
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
  hdf5_.write_step(*this, total_steps, time_);
  hdf5_.finalize();

  if (rank == 0) {
    Real retention = lineage_.resident_retention(cfg_.total_time * 0.5);
    std::cout << "\nSimulation complete.\n"
              << "  Final global agents: " << global_agent_count_ << "\n"
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
  module_chemistry();

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

void Simulation::module_chemistry() {
  // QSSA: compute steady-state toxin field via Green's functions
  Int i_tox = chem_.find("bacteriocin");
  if (i_tox >= 0) {
    qssa_.solve_bacteriocin_field(agents_, chem_, i_tox);
  }

  // Nutrient depletion
  qssa_.solve_nutrient_depletion(agents_, chem_);

  // VBF coupling (nutrient sink/source)
  vbf_.apply_nutrient_coupling(chem_, domain_, cfg_.bio_dt);

  // Apply reactions to concentrations
  for (Int s = 0; s < chem_.num_species(); ++s) {
    #ifdef GUTIBM_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (Int c = 0; c < chem_.ncells(); ++c) {
      chem_.conc(s, c) += chem_.reac(s, c) * cfg_.bio_dt;
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
    }

    // Also check washout criterion: mu_realized < gamma_flow
    Real gamma = advection_.washout_rate(a.x[2]);
    if (a.mu_realized < gamma && a.mu_realized < 0.0) {
      // Metabolically crippled in high-flow zone → washout
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
//  MPI domain decomposition helpers
// ---------------------------------------------------------------------------

// Flat struct for MPI transfer of agent data
namespace {
struct AgentTransferData {
  int64_t  tag;
  int32_t  type;
  int32_t  owner_rank;
  double   x[3];
  double   v[3];
  double   radius;
  double   mass;
  double   outer_radius;
  double   mu_max;
  double   mu_realized;
  double   biomass;
  double   maintenance;
  double   receptor_expr[NUM_RECEPTORS];
  double   km_iron;
  double   km_b12;
  double   km_carbon;
  int32_t  state;
  double   age;
  double   sos_timer;
  int32_t  grid_cell;
  // Genome (simplified for transfer)
  int64_t  lineage_id;
  int64_t  parent_id;
  uint32_t generation;
  uint32_t mutations;
  int32_t  has_conjugative_plasmid;
  double   plasmid_cost_amelioration;
  double   genome_receptor_expr[NUM_RECEPTORS];
  int32_t  num_bi_loci;  // how many BIClusters follow
};

struct BIClusterTransferData {
  uint16_t toxin_id;
  uint16_t immunity_id;
  int32_t  target;   // ReceptorType as int
  int32_t  bclass;   // BacteriocinClass as int
  double   pI;
  double   diff_coeff;
  double   retardation;
  double   molecular_weight;
};

void pack_agent(const Agent& a, AgentTransferData& d) {
  d.tag = a.tag;
  d.type = a.type;
  d.owner_rank = a.owner_rank;
  for (int k = 0; k < 3; ++k) { d.x[k] = a.x[k]; d.v[k] = a.v[k]; }
  d.radius = a.radius;
  d.mass = a.mass;
  d.outer_radius = a.outer_radius;
  d.mu_max = a.mu_max;
  d.mu_realized = a.mu_realized;
  d.biomass = a.biomass;
  d.maintenance = a.maintenance;
  for (int k = 0; k < NUM_RECEPTORS; ++k)
    d.receptor_expr[k] = a.receptor_expr[k];
  d.km_iron = a.km_iron;
  d.km_b12 = a.km_b12;
  d.km_carbon = a.km_carbon;
  d.state = static_cast<int32_t>(a.state);
  d.age = a.age;
  d.sos_timer = a.sos_timer;
  d.grid_cell = a.grid_cell;
  d.lineage_id = a.genome.lineage_id;
  d.parent_id = a.genome.parent_id;
  d.generation = a.genome.generation;
  d.mutations = a.genome.mutations;
  d.has_conjugative_plasmid = a.genome.has_conjugative_plasmid ? 1 : 0;
  d.plasmid_cost_amelioration = a.genome.plasmid_cost_amelioration;
  for (int k = 0; k < NUM_RECEPTORS; ++k)
    d.genome_receptor_expr[k] = a.genome.receptor_expression[k];
  d.num_bi_loci = static_cast<int32_t>(a.genome.bi_loci.size());
}

void pack_bi_cluster(const BICluster& c, BIClusterTransferData& d) {
  d.toxin_id = c.toxin_id;
  d.immunity_id = c.immunity_id;
  d.target = static_cast<int32_t>(c.target);
  d.bclass = static_cast<int32_t>(c.bclass);
  d.pI = c.pI;
  d.diff_coeff = c.diff_coeff;
  d.retardation = c.retardation;
  d.molecular_weight = c.molecular_weight;
}

Agent unpack_agent(const AgentTransferData& d,
                   const BIClusterTransferData* bi_data) {
  Agent a{};
  a.tag = d.tag;
  a.type = d.type;
  a.owner_rank = d.owner_rank;
  for (int k = 0; k < 3; ++k) { a.x[k] = d.x[k]; a.v[k] = d.v[k]; }
  a.radius = d.radius;
  a.mass = d.mass;
  a.outer_radius = d.outer_radius;
  a.mu_max = d.mu_max;
  a.mu_realized = d.mu_realized;
  a.biomass = d.biomass;
  a.maintenance = d.maintenance;
  for (int k = 0; k < NUM_RECEPTORS; ++k)
    a.receptor_expr[k] = d.receptor_expr[k];
  a.km_iron = d.km_iron;
  a.km_b12 = d.km_b12;
  a.km_carbon = d.km_carbon;
  a.state = static_cast<PhenoState>(d.state);
  a.age = d.age;
  a.sos_timer = d.sos_timer;
  a.grid_cell = d.grid_cell;
  a.genome.lineage_id = d.lineage_id;
  a.genome.parent_id = d.parent_id;
  a.genome.generation = d.generation;
  a.genome.mutations = d.mutations;
  a.genome.has_conjugative_plasmid = (d.has_conjugative_plasmid != 0);
  a.genome.plasmid_cost_amelioration = d.plasmid_cost_amelioration;
  for (int k = 0; k < NUM_RECEPTORS; ++k)
    a.genome.receptor_expression[k] = d.genome_receptor_expr[k];

  a.genome.bi_loci.resize(d.num_bi_loci);
  for (int k = 0; k < d.num_bi_loci; ++k) {
    const auto& bd = bi_data[k];
    a.genome.bi_loci[k].toxin_id = bd.toxin_id;
    a.genome.bi_loci[k].immunity_id = bd.immunity_id;
    a.genome.bi_loci[k].target = static_cast<ReceptorType>(bd.target);
    a.genome.bi_loci[k].bclass = static_cast<BacteriocinClass>(bd.bclass);
    a.genome.bi_loci[k].pI = bd.pI;
    a.genome.bi_loci[k].diff_coeff = bd.diff_coeff;
    a.genome.bi_loci[k].retardation = bd.retardation;
    a.genome.bi_loci[k].molecular_weight = bd.molecular_weight;
  }
  return a;
}

// Serialize a list of agents into byte buffers
void serialize_agents(const std::vector<Agent>& agents,
                      std::vector<char>& buf) {
  buf.clear();
  int32_t count = static_cast<int32_t>(agents.size());
  buf.insert(buf.end(), reinterpret_cast<const char*>(&count),
             reinterpret_cast<const char*>(&count) + sizeof(count));
  for (const auto& a : agents) {
    AgentTransferData d;
    pack_agent(a, d);
    buf.insert(buf.end(), reinterpret_cast<const char*>(&d),
               reinterpret_cast<const char*>(&d) + sizeof(d));
    for (const auto& c : a.genome.bi_loci) {
      BIClusterTransferData bd;
      pack_bi_cluster(c, bd);
      buf.insert(buf.end(), reinterpret_cast<const char*>(&bd),
                 reinterpret_cast<const char*>(&bd) + sizeof(bd));
    }
  }
}

// Deserialize agents from a byte buffer
std::vector<Agent> deserialize_agents(const std::vector<char>& buf) {
  std::vector<Agent> result;
  if (buf.size() < sizeof(int32_t)) return result;
  size_t offset = 0;
  int32_t count;
  std::memcpy(&count, buf.data() + offset, sizeof(count));
  offset += sizeof(count);

  result.reserve(count);
  for (int32_t i = 0; i < count; ++i) {
    AgentTransferData d;
    std::memcpy(&d, buf.data() + offset, sizeof(d));
    offset += sizeof(d);

    std::vector<BIClusterTransferData> bi_data(d.num_bi_loci);
    for (int32_t k = 0; k < d.num_bi_loci; ++k) {
      std::memcpy(&bi_data[k], buf.data() + offset, sizeof(bi_data[k]));
      offset += sizeof(bi_data[k]);
    }
    result.push_back(unpack_agent(d, bi_data.data()));
  }
  return result;
}
}  // anonymous namespace

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
  serialize_agents(send_lo, buf_send_lo);
  serialize_agents(send_hi, buf_send_hi);

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
  auto recv_lo = deserialize_agents(buf_recv_lo);
  auto recv_hi = deserialize_agents(buf_recv_hi);

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
  serialize_agents(ghost_lo, buf_send_lo);
  serialize_agents(ghost_hi, buf_send_hi);

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
  auto recv_lo = deserialize_agents(buf_recv_lo);
  auto recv_hi = deserialize_agents(buf_recv_hi);

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
