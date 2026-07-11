/* -----------------------------------------------------------------------
   GutIBM – Multi-rank MPI integration tests (issue #43)
   Run with: mpirun -np 2 test_mpi_multi_rank
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "mpi_test_helpers.h"

#include <cassert>
#include <cmath>
#include <iostream>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

using namespace gutibm;
using gutibm::test::assert_unique_tags;
using gutibm::test::gather_live_tags_flat;
using gutibm::test::make_mpi_config;
using gutibm::test::require_mpi_ranks;

namespace {

#ifdef GUTIBM_MPI

SimulationConfig make_mpi_periodic_config() {
  SimulationConfig cfg = make_mpi_config(4242, 40);
  cfg.domain.periodic = {true, true, false};
  return cfg;
}

void test_reaction_sum_and_diffusion_are_rank_identical() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  DomainConfig domain_cfg;
  domain_cfg.lo = {0.0, 0.0, 0.0};
  domain_cfg.hi = {20e-6, 15e-6, 15e-6};
  domain_cfg.grid_dx = 5e-6;
  domain_cfg.periodic = {true, true, false};
  Domain domain;
  domain.init(domain_cfg);

  ChemicalSpec oxygen;
  oxygen.name = "oxygen";
  oxygen.diff_coeff = 2.1e-9;
  oxygen.retardation = 1.0;
  oxygen.initial_conc = 0.0;
  oxygen.boundary_conc = 1.0;
  oxygen.diffusion_enabled = true;
  ChemicalField chem;
  chem.init(domain, {oxygen});

  const Int reaction_cell = domain.cell_index(1, 1, 1);
  chem.reac(0, reaction_cell) = static_cast<Real>(rank + 1);
  chem.sum_reactions_across_ranks();
  assert(std::abs(chem.reac(0, reaction_cell) - 3.0) < 1e-15);

  chem.conc(0, reaction_cell) += chem.reac(0, reaction_cell);
  chem.apply_diffusion(domain, 60.0);

  Real checksum = 0.0;
  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    checksum += chem.conc(0, cell) * static_cast<Real>(cell + 1);
  }
  Real minimum = 0.0;
  Real maximum = 0.0;
  MPI_Allreduce(&checksum, &minimum, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&checksum, &maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  assert(std::abs(maximum - minimum) < 1e-12);

  if (rank == 0) {
    std::cout << "  test_reaction_sum_and_diffusion_are_rank_identical: PASSED\n";
  }
}

void test_slab_decomposition() {
  require_mpi_ranks(2);

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

void test_slab_decomposition_periodic_x() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  DomainConfig cfg;
  cfg.lo = {0, 0, 0};
  cfg.hi = {100e-6, 100e-6, 50e-6};
  cfg.mpi_decomp_axis = 0;
  cfg.ghost_width = 10e-6;
  cfg.periodic = {true, true, false};

  Domain dom;
  dom.init(cfg);

  assert(dom.neighbors_collapsed());
  if (rank == 0) {
    assert(dom.rank_lo() == 1);
    assert(dom.rank_hi() == 1);
  } else {
    assert(dom.rank_lo() == 0);
    assert(dom.rank_hi() == 0);
  }

  if (rank == 0) {
    std::cout << "  test_slab_decomposition_periodic_x: PASSED\n";
  }
}

void test_init_population_partitioned() {
  require_mpi_ranks(2);

  SimulationConfig cfg = make_mpi_config(4242, 40);
  Simulation sim;
  sim.init(cfg);

  assert(sim.global_agent_count() == 40);

  for (const Agent& a : sim.agents()) {
    assert(sim.domain().is_local(a.x));
    assert(a.identity.owner_rank == sim.domain().rank());
  }

  auto tags = gather_live_tags_flat(sim);
  assert(static_cast<Int>(tags.size()) == sim.global_agent_count());
  assert_unique_tags(tags);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "  test_init_population_partitioned: PASSED\n";
  }
}

void test_migration_preserves_global_count() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_config(4242, 40);
  cfg.advection.distal_transit_time = 1e12;
  Simulation sim;
  sim.init(cfg);

  const Int initial_global = sim.global_agent_count();
  assert(initial_global == 40);

  Real moved_x = 0.0;
  if (rank == 0) {
    for (Agent& a : sim.agents()) {
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

  sim.step(cfg.time.bio_dt);

  auto tags = gather_live_tags_flat(sim);
  assert_unique_tags(tags);

  int found_on_rank1 = 0;
  if (rank == 1) {
    for (const Agent& a : sim.agents()) {
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
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_config(4242, 40);
  Simulation sim;
  sim.init(cfg);

  const Real gw = sim.domain().ghost_width();
  for (Agent& a : sim.agents()) {
    if (a.state == PhenoState::DEAD) continue;
    if (rank == 0) {
      a.x[0] = sim.domain().local_hi_x() - gw * 0.5;
    } else {
      a.x[0] = sim.domain().local_lo_x() + gw * 0.5;
    }
    break;
  }

  const Int before = sim.global_agent_count();
  sim.step(cfg.time.bio_dt);
  assert(sim.global_agent_count() > 0);
  assert(sim.global_agent_count() <= before);

  if (rank == 0) {
    std::cout << "  test_boundary_ghost_exchange_runs: PASSED\n";
  }
}

void test_periodic_x_ghost_and_migration() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_periodic_config();
  cfg.advection.distal_transit_time = 1e12;
  cfg.advection.crypts_enabled = true;
  cfg.advection.crypt_exit_rate = 0.0;
  cfg.advection.crypt_entry_rate = 0.0;
  Simulation sim;
  sim.init(cfg);

  assert(sim.domain().neighbors_collapsed());
  assert_unique_tags(gather_live_tags_flat(sim));

  const Real gw = sim.domain().ghost_width();
  Real moved_x = 0.0;
  if (rank == 0) {
    for (Agent& a : sim.agents()) {
      if (a.state == PhenoState::DEAD) continue;
      a.x[0] = sim.domain().local_hi_x() - gw * 0.5;
      moved_x = sim.domain().local_hi_x() + 5e-6;
      a.mu_realized = a.mu_max;
      break;
    }
  } else {
    for (Agent& a : sim.agents()) {
      if (a.state == PhenoState::DEAD) continue;
      a.x[0] = sim.domain().local_lo_x() + gw * 0.5;
      break;
    }
  }

  // Crypt refugia bypass washout so this test isolates periodic MPI exchange.
  for (Agent& a : sim.agents()) {
    if (a.state != PhenoState::DEAD) {
      a.flags.in_crypt = true;
    }
  }

  sim.step(cfg.time.bio_dt);
  assert(sim.global_agent_count() == 40);
  assert_unique_tags(gather_live_tags_flat(sim));

  if (rank == 0) {
    moved_x = sim.domain().local_hi_x() + 5e-6;
    for (Agent& a : sim.agents()) {
      if (a.state != PhenoState::DEAD) {
        a.x[0] = moved_x;
        a.x[1] = 50e-6;
        a.x[2] = 25e-6;
        a.mu_realized = a.mu_max;
        break;
      }
    }
  }
  MPI_Bcast(&moved_x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  sim.step(cfg.time.bio_dt);
  assert(sim.global_agent_count() == 40);
  assert_unique_tags(gather_live_tags_flat(sim));

  if (rank == 0) {
    std::cout << "  test_periodic_x_ghost_and_migration: PASSED\n";
  }
}

void test_multirank_simulation_steps() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_config(4242, 40);
  cfg.time.total_time = 120.0;
  Simulation sim;
  sim.init(cfg);

  const Int initial_global = sim.global_agent_count();
  sim.run();

  assert(sim.global_agent_count() > 0);
  assert(sim.global_agent_count() <= initial_global);
  assert_unique_tags(gather_live_tags_flat(sim));

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

  test_reaction_sum_and_diffusion_are_rank_identical();
  test_slab_decomposition();
  test_slab_decomposition_periodic_x();
  test_init_population_partitioned();
  test_migration_preserves_global_count();
  test_boundary_ghost_exchange_runs();
  test_periodic_x_ghost_and_migration();
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
