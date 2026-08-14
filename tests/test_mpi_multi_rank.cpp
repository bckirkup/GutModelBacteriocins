/* -----------------------------------------------------------------------
   GutIBM – Multi-rank MPI integration tests (issue #43)
   Run with: mpirun -np 2 test_mpi_multi_rank
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "mpi_test_helpers.h"
#include "species_names.h"
#include "plasmid.h"
#include "greens_function.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>

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

void test_chemical_field_layout_mapping() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  DomainConfig cfg;
  cfg.lo = {0.0, 0.0, 0.0};
  cfg.hi = {100e-6, 20e-6, 20e-6};
  cfg.grid_dx = 5e-6;
  cfg.ghost_width = 10e-6;
  cfg.grid_halo_width = 2;
  Domain domain;
  domain.init(cfg);

  ChemicalSpec spec;
  spec.name = "carbon";
  ChemicalField chem;
  chem.init(domain, {spec}, "slab");

  Int owned = 0;
  Int halo = 0;
  for (Int cell = 0; cell < domain.ncells(); ++cell) {
    const Int local_x = domain.global_to_local_grid_x(cell % domain.nx());
    if (local_x < 0) {
      continue;
    }
    assert(chem.storage_to_global_cell(
               chem.global_to_storage_cell(cell)) == cell);
    if (chem.owns_global_cell(cell)) {
      ++owned;
    } else {
      assert(chem.global_cell_in_halo(cell));
      ++halo;
    }
  }
  assert(owned > 0);
  assert(halo > 0);
  if (domain.local_grid_x_begin() == 0) {
    assert(chem.owned_global_x_begin() == 0);
  }
  if (rank == 0) {
    std::cout << "  test_chemical_field_layout_mapping: PASSED\n";
  }
}

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

void test_cross_rank_bacteriocin_source_exchange() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_config(4243, 40);
  cfg.qssa.toxin_cutoff = 80e-6;
  cfg.enabled_fixes = {"receptor"};
  Simulation sim;
  sim.init(cfg);

  for (Agent& agent : sim.agents()) {
    agent.genome.bi_loci.clear();
  }
  if (rank == 0) {
    Agent& producer = sim.agents()[0];
    producer.x = {49.5e-6, 50e-6, 25e-6};
    producer.grid_cell = sim.domain().cell_index(9, 10, 5);
    producer.genome.bi_loci.push_back(PlasmidLibrary::microcin_V());
  }

  const Int target_ix = 11;
  const Int target_iy = 10;
  const Int target_iz = 5;
  const Vec3 target_position =
      sim.domain().cell_center(target_ix, target_iy, target_iz);
  assert(sim.domain().owner_rank(target_position) == 1);

  Int local_sources = rank == 0 ? 1 : 0;
  Int global_sources = 0;
  MPI_Allreduce(&local_sources, &global_sources, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  assert(global_sources == 1);

  const Int species_idx =
      sim.chemical_field().find(species::BACTERIOCIN_CIRA);
  assert(species_idx >= 0);
  sim.step(cfg.time.bio_dt);

  const Int target_cell =
      sim.domain().cell_index(target_ix, target_iy, target_iz);
  const Real local_concentration =
      sim.chemical_field().conc(species_idx, target_cell);
  Real minimum = 0.0;
  Real maximum = 0.0;
  MPI_Allreduce(&local_concentration, &minimum, 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(&local_concentration, &maximum, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  assert(minimum > 0.0);
  assert(std::abs(maximum - minimum)
         <= 1.0e-12 * std::max(1.0, maximum));

  GreensFunction greens;
  greens.init(sim.domain(), sim.advection());
  GreensFunctionParams params;
  const BICluster source_cluster = PlasmidLibrary::microcin_V();
  params.diff_coeff = source_cluster.diff_coeff;
  params.retardation = source_cluster.retardation;
  params.pI = source_cluster.pI;
  params.source_rate = cfg.qssa.microcin_secretion;
  const Real protease_decay =
      cfg.chem_env.protease.enabled && source_cluster.protease_half_life > 0.0
          ? std::numbers::ln2 / source_cluster.protease_half_life
          : 0.0;
  const Real dilution_decay = std::max(
      sim.advection().washout_rate(25e-6),
      cfg.chem_env.protease.dilution_rate);
  params.decay_rate = protease_decay + dilution_decay;
  const Real expected = greens.concentration_bounded(
      {49.5e-6, 50e-6, 25e-6}, target_position, params);
  assert(std::abs(minimum - expected)
         <= 1.0e-12 * std::max(std::abs(expected), 1.0e-30));

  if (rank == 0) {
    std::cout << "  test_cross_rank_bacteriocin_source_exchange: PASSED\n";
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

  gutibm::test::assert_cell_aligned_slab_contract(dom);

  if (rank == 0) {
    assert(dom.rank_lo() == -1);
    assert(dom.rank_hi() == 1);
    const Real mid_x = 0.5 * (dom.local_lo_x() + dom.local_hi_x());
    assert(dom.is_local({mid_x, 25e-6, 25e-6}));
    assert(dom.owner_rank({mid_x, 25e-6, 25e-6}) == rank);
  } else {
    assert(dom.rank_lo() == 0);
    assert(dom.rank_hi() == -1);
    const Real mid_x = 0.5 * (dom.local_lo_x() + dom.local_hi_x());
    assert(dom.is_local({mid_x, 25e-6, 25e-6}));
    assert(dom.owner_rank({mid_x, 25e-6, 25e-6}) == rank);
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

void test_adaptive_dt_is_rank_identical() {
  require_mpi_ranks(2);

  SimulationConfig cfg = make_mpi_config(8080, 40);
  cfg.adaptive_dt.enabled = true;
  cfg.adaptive_dt.min = 1.0;
  cfg.adaptive_dt.max = 300.0;
  cfg.adaptive_dt.safety = 0.8;
  cfg.adaptive_dt.growth_limit = 0.1;
  cfg.initial_strains.clear();

  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 40;
  resident.mu_max = 1.0e-8;
  cfg.initial_strains.push_back(resident);

  SimulationConfig::InitialStrain fast;
  fast.type = 2;
  fast.count = 1;
  fast.mu_max = 5.0e-2;
  cfg.initial_strains.push_back(fast);

  Simulation sim;
  sim.init(cfg);

  const Real dt = sim.compute_adaptive_dt();
  Real local_max_abs_mu = 0.0;
  for (const Agent& agent : sim.agents()) {
    local_max_abs_mu = std::max(local_max_abs_mu,
                                std::abs(agent.mu_realized));
  }
  Real minimum_local_max_abs_mu = 0.0;
  Real maximum_local_max_abs_mu = 0.0;
  MPI_Allreduce(&local_max_abs_mu, &minimum_local_max_abs_mu, 1, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_max_abs_mu, &maximum_local_max_abs_mu, 1, MPI_DOUBLE,
                MPI_MAX, MPI_COMM_WORLD);
  assert(minimum_local_max_abs_mu < 1.0e-6);
  assert(maximum_local_max_abs_mu > 1.0e-2);
  Real minimum = 0.0;
  Real maximum = 0.0;
  MPI_Allreduce(&dt, &minimum, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&dt, &maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  assert(std::abs(maximum - minimum) < 1e-12);
  assert(dt < cfg.adaptive_dt.max);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "  test_adaptive_dt_is_rank_identical: PASSED"
              << " (dt=" << dt << ")\n";
  }
}

void test_multirank_immigration_ownership() {
  require_mpi_ranks(2);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_config(9876, 20);
  cfg.enabled_fixes = {"mechanics"};
  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = 0;
  immigrant.mu_max = 3.5e-4;
  immigrant.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(immigrant);
  cfg.immigration.enabled = true;
  cfg.immigration.count = 3;
  cfg.immigration.strain_index = 1;
  cfg.immigration.step = 0;

  Simulation sim;
  sim.init(cfg);
  sim.step(cfg.time.bio_dt);

  const std::vector<TagID> tags = gather_live_tags_flat(sim);
  assert_unique_tags(tags);
  Int local_immigrants = 0;
  for (const Agent& agent : sim.agents()) {
    if (agent.identity.type == immigrant.type) ++local_immigrants;
  }
  Int global_immigrants = 0;
  MPI_Allreduce(&local_immigrants, &global_immigrants, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  assert(global_immigrants == 3);
  if (rank == 0) {
    std::cout << "  test_multirank_immigration_ownership: PASSED\n";
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
  test_chemical_field_layout_mapping();
  test_cross_rank_bacteriocin_source_exchange();
  test_slab_decomposition();
  test_slab_decomposition_periodic_x();
  test_init_population_partitioned();
  test_migration_preserves_global_count();
  test_boundary_ghost_exchange_runs();
  test_periodic_x_ghost_and_migration();
  test_multirank_simulation_steps();
  test_adaptive_dt_is_rank_identical();
  test_multirank_immigration_ownership();

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
