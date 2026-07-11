/* -----------------------------------------------------------------------
   GutIBM – Shared MPI integration test helpers
   ----------------------------------------------------------------------- */

#pragma once

#include "input_parser.h"
#include "simulation.h"

#include <algorithm>
#include <cassert>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {
namespace test {

inline SimulationConfig make_mpi_config(unsigned seed, int agent_count) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.domain.ghost_width = 10e-6;
  cfg.domain.periodic = {false, true, false};
  cfg.time.total_time = 300.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 300.0;
  cfg.seed = seed;
  cfg.hdf5.enabled = false;
  cfg.cell_bio.fur.enabled = false;
  cfg.cell_bio.cdi.enabled = false;
  cfg.cell_bio.motility.enabled = false;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length = 100e-6;
  cfg.advection.radial_turnover = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff = 50e-6;
  cfg.qssa.nutrient_cutoff = 25e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = agent_count;
  s.mu_max = 5e-4;
  s.plasmids = {};
  s.conjugative = false;
  cfg.initial_strains.push_back(s);
  return cfg;
}

#ifdef GUTIBM_MPI

inline std::vector<TagID> gather_live_tags_flat(const Simulation& sim) {
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  std::vector<TagID> local;
  for (const Agent& a : sim.agents()) {
    if (a.state != PhenoState::DEAD) {
      local.push_back(a.identity.tag);
    }
  }

  auto local_n = static_cast<int>(local.size());
  std::vector<int> counts(nprocs, 0);
  MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  int total = 0;
  std::vector<int> displ(nprocs, 0);
  size_t r = 0;
  for (int count : counts) {
    displ[r++] = total;
    total += count;
  }

  std::vector<TagID> all(static_cast<size_t>(total));
  MPI_Allgatherv(local.data(), local_n, MPI_INT64_T,
                 all.data(), counts.data(), displ.data(), MPI_INT64_T,
                 MPI_COMM_WORLD);
  return all;
}

inline void assert_unique_tags(const std::vector<TagID>& tags) {
  std::vector<TagID> sorted = tags;
  std::ranges::sort(sorted);
  assert(std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end());
}

inline void require_mpi_ranks(int expected) {
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  assert(nprocs == expected);
}

#endif  // GUTIBM_MPI

}  // namespace test
}  // namespace gutibm
