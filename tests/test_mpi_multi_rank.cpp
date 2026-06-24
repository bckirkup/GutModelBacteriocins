/* -----------------------------------------------------------------------
   GutIBM – Multi-rank MPI integration tests (issue #43)
   Run with: mpirun -np 2 test_mpi_multi_rank
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

using namespace gutibm;

namespace {

SimulationConfig make_mpi_config() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.domain.ghost_width = 10e-6;
  // Non-periodic x: with 2 ranks, each slab has a unique neighbor (avoids
  // duplicate Sendrecv to the same rank when rank_lo == rank_hi).
  cfg.domain.periodic = {false, true, false};
  cfg.total_time = 300.0;
  cfg.bio_dt = 60.0;
  cfg.output_interval = 300.0;
  cfg.seed = 4242;
  cfg.hdf5.enabled = false;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length = 100e-6;
  cfg.advection.radial_turnover = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff = 50e-6;
  cfg.qssa.nutrient_cutoff = 25e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 40;
  s.mu_max = 5e-4;
  s.plasmids = {};
  s.conjugative = false;
  cfg.initial_strains.push_back(s);
  return cfg;
}

#ifdef GUTIBM_MPI

std::vector<TagID> gather_live_tags_flat(const Simulation& sim) {
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  std::vector<TagID> local;
  for (Int i = 0; i < sim.agents().size(); ++i) {
    const Agent& a = sim.agents()[i];
    if (a.state != PhenoState::DEAD) {
      local.push_back(a.tag);
    }
  }

  int local_n = static_cast<int>(local.size());
  std::vector<int> counts(nprocs, 0);
  MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  int total = 0;
  std::vector<int> displ(nprocs, 0);
  for (int r = 0; r < nprocs; ++r) {
    displ[r] = total;
    total += counts[r];
  }

  std::vector<TagID> all(static_cast<size_t>(total));
  MPI_Allgatherv(local.data(), local_n, MPI_INT64_T,
                  all.data(), counts.data(), displ.data(), MPI_INT64_T,
                  MPI_COMM_WORLD);
  return all;
}

void require_two_ranks() {
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  assert(nprocs == 2);
}

void test_slab_decomposition() {
  require_two_ranks();

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  DomainConfig cfg;
  cfg.lo = {0, 0, 0};
  cfg.hi = {100e-6, 100e-6, 50e-6};
  cfg.mpi_decomp_axis = 0;
  cfg.ghost_width = 10e-6;
  cfg.periodic = {false, true, false};

  Domain dom;
  dom.init(cfg);

  if (rank == 0) {
    assert(dom.rank_lo() == -1);
    assert(dom.rank_hi() == 1);
    assert(std::abs(dom.local_lo_x() - 0.0) < 1e-15);
    assert(std::abs(dom.local_hi_x() - 50e-6) < 1e-15);
    assert(dom.is_local({25e-6, 25e-6, 25e-6}));
    assert(!dom.is_local({75e-6, 25e-6, 25e-6}));
    assert(dom.owner_rank({25e-6, 25e-6, 25e-6}) == 0);
    assert(dom.owner_rank({75e-6, 25e-6, 25e-6}) == 1);
  } else {
    assert(dom.rank_lo() == 0);
    assert(dom.rank_hi() == -1);
    assert(std::abs(dom.local_lo_x() - 50e-6) < 1e-15);
    assert(std::abs(dom.local_hi_x() - 100e-6) < 1e-15);
    assert(!dom.is_local({25e-6, 25e-6, 25e-6}));
    assert(dom.is_local({75e-6, 25e-6, 25e-6}));
  }

  if (rank == 0) {
    std::cout << "  test_slab_decomposition: PASSED\n";
  }
}

void test_init_population_partitioned() {
  require_two_ranks();

  SimulationConfig cfg = make_mpi_config();
  Simulation sim;
  sim.init(cfg);

  assert(sim.global_agent_count() == 40);

  for (Int i = 0; i < sim.agents().size(); ++i) {
    const Agent& a = sim.agents()[i];
    assert(sim.domain().is_local(a.x));
    assert(a.owner_rank == sim.domain().rank());
  }

  auto tags = gather_live_tags_flat(sim);
  assert(static_cast<Int>(tags.size()) == sim.global_agent_count());

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "  test_init_population_partitioned: PASSED\n";
  }
}

void test_migration_preserves_global_count() {
  require_two_ranks();

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_config();
  cfg.advection.distal_transit_time = 1e12;
  Simulation sim;
  sim.init(cfg);

  const Int initial_global = sim.global_agent_count();
  assert(initial_global == 40);

  Real moved_x = 0.0;
  if (rank == 0) {
    for (Int i = 0; i < sim.agents().size(); ++i) {
      Agent& a = sim.agents()[i];
      if (a.state != PhenoState::DEAD) {
        moved_x = sim.domain().local_hi_x() + 10e-6;
        a.x[0] = moved_x;
        a.x[1] = 50e-6;
        a.x[2] = 25e-6;
        a.mu_realized = a.mu_max;
        break;
      }
    }
  }
  MPI_Bcast(&moved_x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  assert(moved_x > 0.0);

  sim.step(cfg.bio_dt);

  int found_on_rank1 = 0;
  if (rank == 1) {
    for (Int i = 0; i < sim.agents().size(); ++i) {
      const Agent& a = sim.agents()[i];
      if (a.state != PhenoState::DEAD &&
          std::abs(a.x[0] - moved_x) < 2e-6) {
        found_on_rank1 = 1;
        break;
      }
    }
  }
  int global_found = 0;
  MPI_Allreduce(&found_on_rank1, &global_found, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  assert(global_found == 1);

  if (rank == 0) {
    std::cout << "  test_migration_preserves_global_count: PASSED\n";
  }
}

void test_boundary_ghost_exchange_runs() {
  require_two_ranks();

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_config();
  Simulation sim;
  sim.init(cfg);

  const Real gw = sim.domain().ghost_width();
  for (Int i = 0; i < sim.agents().size(); ++i) {
    Agent& a = sim.agents()[i];
    if (a.state == PhenoState::DEAD) continue;
    if (rank == 0) {
      a.x[0] = sim.domain().local_hi_x() - gw * 0.5;
    } else {
      a.x[0] = sim.domain().local_lo_x() + gw * 0.5;
    }
    break;
  }

  const Int before = sim.global_agent_count();
  sim.step(cfg.bio_dt);
  assert(sim.global_agent_count() > 0);
  assert(sim.global_agent_count() <= before);

  if (rank == 0) {
    std::cout << "  test_boundary_ghost_exchange_runs: PASSED\n";
  }
}

void test_multirank_simulation_steps() {
  require_two_ranks();

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_config();
  cfg.total_time = 120.0;
  Simulation sim;
  sim.init(cfg);

  const Int initial_global = sim.global_agent_count();
  sim.run();

  assert(sim.global_agent_count() > 0);
  assert(sim.global_agent_count() <= initial_global);

  if (rank == 0) {
    std::cout << "  test_multirank_simulation_steps: PASSED"
              << " (global_agents=" << sim.global_agent_count() << ")\n";
  }
}

#endif  // GUTIBM_MPI

}  // namespace

int main(int argc, char** argv) {
#ifdef GUTIBM_MPI
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "=== MPI Multi-Rank Tests (np=2) ===\n";
  }

  test_slab_decomposition();
  test_init_population_partitioned();
  test_migration_preserves_global_count();
  test_boundary_ghost_exchange_runs();
  test_multirank_simulation_steps();

  if (rank == 0) {
    std::cout << "All MPI multi-rank tests passed.\n";
  }

  MPI_Finalize();
#else
  (void)argc;
  (void)argv;
  std::cout << "MPI disabled at build time — skipping multi-rank tests.\n";
#endif
  return 0;
}
