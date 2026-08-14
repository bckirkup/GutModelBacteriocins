/* -----------------------------------------------------------------------
   GutIBM – Shared MPI integration test helpers
   ----------------------------------------------------------------------- */

#pragma once

#include "input_parser.h"
#include "simulation.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm::test {

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
  std::vector counts(nprocs, 0);
  MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  int total = 0;
  std::vector displ(nprocs, 0);
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

inline void assert_cell_aligned_slab_contract(const Domain& dom) {
  const Int rank = dom.rank();
  const Int nprocs = dom.nprocs();
  const Int nx = dom.nx();
  const auto [grid_begin, grid_end] =
      Domain::grid_x_range_for_rank(nx, nprocs, rank);
  const Real expected_lo =
      dom.lo()[0] + static_cast<Real>(grid_begin) * dom.dx();
  const Real expected_hi =
      rank == nprocs - 1
          ? dom.hi()[0]
          : dom.lo()[0] + static_cast<Real>(grid_end) * dom.dx();
  const Real boundary_tolerance = 1e-15;
  assert(std::abs(dom.local_lo_x() - expected_lo) < boundary_tolerance);
  assert(std::abs(dom.local_hi_x() - expected_hi) < boundary_tolerance);

  std::vector<Real> all_local_lo(static_cast<std::size_t>(nprocs));
  std::vector<Real> all_local_hi(static_cast<std::size_t>(nprocs));
  std::vector<Int> all_cell_counts(static_cast<std::size_t>(nprocs));
  const Real local_lo_x = dom.local_lo_x();
  const Real local_hi_x = dom.local_hi_x();
  const Int local_cell_count = grid_end - grid_begin;
  MPI_Allgather(&local_lo_x, 1, MPI_DOUBLE, all_local_lo.data(), 1,
                MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(&local_hi_x, 1, MPI_DOUBLE, all_local_hi.data(), 1,
                MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(&local_cell_count, 1, MPI_INT, all_cell_counts.data(), 1,
                MPI_INT, MPI_COMM_WORLD);

  assert(std::abs(all_local_lo.front() - dom.lo()[0]) < boundary_tolerance);
  assert(std::abs(all_local_hi.back() - dom.hi()[0]) < boundary_tolerance);
  for (Int owner = 1; owner < nprocs; ++owner) {
    assert(std::abs(all_local_lo[static_cast<std::size_t>(owner)]
                    - all_local_hi[static_cast<std::size_t>(owner - 1)])
           < boundary_tolerance);
  }
  for (const Real boundary : all_local_lo) {
    const Real scaled_boundary = (boundary - dom.lo()[0]) / dom.dx();
    assert(std::abs(scaled_boundary - std::round(scaled_boundary))
           < 1e-12);
  }

  Int minimum_cells = all_cell_counts.front();
  Int maximum_cells = all_cell_counts.front();
  Int total_cells = 0;
  for (const Int count : all_cell_counts) {
    minimum_cells = std::min(minimum_cells, count);
    maximum_cells = std::max(maximum_cells, count);
    total_cells += count;
  }
  assert(maximum_cells - minimum_cells <= 1);
  assert(total_cells == nx);

  for (Int ix = 0; ix < nx; ++ix) {
    const Vec3 center = dom.cell_center(ix, 0, 0);
    const Int expected_owner =
        Domain::grid_owner_rank_for_cell(nx, nprocs, ix);
    assert(dom.owner_rank(center) == expected_owner);
    assert(dom.is_local(center) == (expected_owner == rank));
  }
}

#endif  // GUTIBM_MPI

}  // namespace gutibm::test
