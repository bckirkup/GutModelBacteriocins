/* -----------------------------------------------------------------------
   GutIBM – Simulation engine implementation
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "fix_metabolism.h"
#include "fix_receptor.h"
#include "fix_bacteriocin.h"
#include "fix_conjugation.h"
#include "fix_mutation.h"
#include "plasmid.h"

#include <iostream>
#include <algorithm>
#include <cmath>

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

  // Initialize fixes
  for (auto& fix : fixes_) {
    fix->init();
  }

  // Create initial population
  init_population(cfg);

  // Initial coupling
  rebuild_spatial_hash();
  update_grid_coupling();

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
              << "  Agents: " << agents_.size() << "\n"
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

      Agent a = Agent::create_default(agents_.next_tag(), strain.type,
                                       pos, strain.mu_max);

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

      agents_.push_back(std::move(a));
    }
  }
}

void Simulation::run() {
  Real dt = cfg_.bio_dt;
  Int total_steps = static_cast<Int>(std::ceil(cfg_.total_time / dt));

  int rank = domain_.rank();

  // Initial output
  if (rank == 0) {
    hdf5_.write_step(*this, 0, 0.0);
    take_lineage_snapshot();
  }

  for (Int s = 1; s <= total_steps; ++s) {
    step(dt);

    // Periodic output
    if (time_ >= next_output_) {
      if (rank == 0) {
        hdf5_.write_step(*this, s, time_);
        std::cout << "Step " << s << "/" << total_steps
                  << "  t=" << time_ << "s"
                  << "  agents=" << agents_.size()
                  << "  mu_avg=" << [&]() {
                       Real sum = 0; Int cnt = 0;
                       for (Int i = 0; i < agents_.size(); ++i) {
                         if (agents_[i].state != PhenoState::DEAD) {
                           sum += agents_[i].mu_realized; cnt++;
                         }
                       }
                       return cnt > 0 ? sum / cnt : 0.0;
                     }()
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
  if (rank == 0) {
    hdf5_.write_step(*this, total_steps, time_);
    hdf5_.finalize();

    Real retention = lineage_.resident_retention(cfg_.total_time * 0.5);
    std::cout << "\nSimulation complete.\n"
              << "  Final agents: " << agents_.size() << "\n"
              << "  Resident retention: " << retention * 100.0 << "%\n"
              << "  Dominant lineage: " << lineage_.dominant_lineage() << "\n"
              << std::flush;
  }
}

void Simulation::step(Real dt) {
  // Pre-step
  chem_.zero_reactions();
  rebuild_spatial_hash();
  update_grid_coupling();

  for (auto& fix : fixes_) {
    fix->pre_step(dt);
  }

  // 1. Biology module
  module_biology(dt);

  // 2. Chemistry module (QSSA, instantaneous)
  module_chemistry();

  // 3. Physics module (advection + mechanics)
  module_physics(dt);

  // Post-step
  for (auto& fix : fixes_) {
    fix->post_step(dt);
  }

  // Cleanup
  check_washout();
  remove_dead_agents();

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
  vbf_.apply_nutrient_coupling(chem_, cfg_.bio_dt);

  // Apply reactions to concentrations
  for (Int s = 0; s < chem_.num_species(); ++s) {
    for (Int c = 0; c < chem_.ncells(); ++c) {
      chem_.conc(s, c) += chem_.reac(s, c) * cfg_.bio_dt;
      chem_.conc(s, c) = std::max(chem_.conc(s, c), 0.0);
    }
  }

  // Boundary conditions
  chem_.apply_boundaries(domain_);
}

void Simulation::module_physics(Real dt) {
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

  // Simple mechanical repulsion between overlapping cells
  auto& hash = domain_.spatial_hash();
  for (Int i = 0; i < agents_.size(); ++i) {
    Agent& ai = agents_[i];
    if (ai.state == PhenoState::DEAD) continue;

    auto neighbors = hash.query_neighbors(ai.x);
    for (Int j : neighbors) {
      if (j <= i) continue;
      Agent& aj = agents_[j];
      if (aj.state == PhenoState::DEAD) continue;

      Vec3 delta = domain_.min_image_delta(ai.x, aj.x);
      Real d2 = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
      Real sum_r = ai.radius + aj.radius;

      if (d2 < sum_r * sum_r && d2 > 0.0) {
        Real d = std::sqrt(d2);
        Real overlap = sum_r - d;
        Real force_mag = 1.0e-10 * overlap;  // Hertzian-like repulsion

        Vec3 n = {delta[0]/d, delta[1]/d, delta[2]/d};
        Real push = force_mag * dt / std::max(ai.mass, 1.0e-30);

        ai.x[0] -= n[0] * push * 0.5;
        ai.x[1] -= n[1] * push * 0.5;
        ai.x[2] -= n[2] * push * 0.5;
        aj.x[0] += n[0] * push * 0.5;
        aj.x[1] += n[1] * push * 0.5;
        aj.x[2] += n[2] * push * 0.5;
      }
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

}  // namespace gutibm
