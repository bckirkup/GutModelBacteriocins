/* -----------------------------------------------------------------------
   GutIBM – Four-rank MPI integration tests (issue #154)
   Run with: mpirun -np 4 test_mpi_four_rank
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

void test_four_rank_slab_decomposition() {
  require_mpi_ranks(4);

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

  const Real slab = 25e-6;
  assert(std::abs(dom.local_lo_x() - static_cast<Real>(rank) * slab) < 1e-15);
  assert(std::abs(dom.local_hi_x() - static_cast<Real>(rank + 1) * slab) < 1e-15);

  const Real mid_x = (static_cast<Real>(rank) + 0.5) * slab;
  assert(dom.is_local({mid_x, 25e-6, 25e-6}));
  assert(dom.owner_rank({mid_x, 25e-6, 25e-6}) == rank);

  if (rank == 0) {
    assert(dom.rank_lo() == -1);
    assert(dom.rank_hi() == 1);
  } else if (rank == 3) {
    assert(dom.rank_lo() == 2);
    assert(dom.rank_hi() == -1);
  } else {
    assert(dom.rank_lo() == rank - 1);
    assert(dom.rank_hi() == rank + 1);
  }

  if (rank == 0) {
    std::cout << "  test_four_rank_slab_decomposition: PASSED\n";
  }
}

void test_reaction_sum_four_ranks() {
  require_mpi_ranks(4);

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
  assert(std::abs(chem.reac(0, reaction_cell) - 10.0) < 1e-15);

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
    std::cout << "  test_reaction_sum_four_ranks: PASSED\n";
  }
}

void test_init_population_partitioned_four_ranks() {
  require_mpi_ranks(4);

  SimulationConfig cfg = make_mpi_config(4404, 80);
  Simulation sim;
  sim.init(cfg);

  assert(sim.global_agent_count() == 80);

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
    std::cout << "  test_init_population_partitioned_four_ranks: PASSED\n";
  }
}

void test_migration_across_four_ranks() {
  require_mpi_ranks(4);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_config(4404, 80);
  cfg.advection.distal_transit_time = 1e12;
  Simulation sim;
  sim.init(cfg);

  assert(sim.global_agent_count() == 80);

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
  assert_unique_tags(gather_live_tags_flat(sim));

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
    std::cout << "  test_migration_across_four_ranks: PASSED\n";
  }
}

void test_multirank_simulation_steps_four_ranks() {
  require_mpi_ranks(4);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_config(4404, 80);
  cfg.time.total_time = 120.0;
  Simulation sim;
  sim.init(cfg);

  const Int initial_global = sim.global_agent_count();
  sim.run();

  assert(sim.global_agent_count() > 0);
  assert(sim.global_agent_count() <= initial_global);
  assert_unique_tags(gather_live_tags_flat(sim));

  if (rank == 0) {
    std::cout << "  test_multirank_simulation_steps_four_ranks: PASSED"
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
    std::cout << "=== MPI Four-Rank Tests (np=4) ===\n";
  }

  test_four_rank_slab_decomposition();
  test_reaction_sum_four_ranks();
  test_init_population_partitioned_four_ranks();
  test_migration_across_four_ranks();
  test_multirank_simulation_steps_four_ranks();

  if (rank == 0) {
    std::cout << "All MPI four-rank tests passed.\n";
  }

  MPI_Finalize();
#else
  (void)argc;
  (void)argv;
  std::cout << "MPI disabled at build time — skipping four-rank tests.\n";
#endif
  return 0;
}
