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
#include "error.h"
#include "hdf5_test_helpers.h"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <format>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

using namespace gutibm;
using gutibm::test::assert_unique_tags;
using gutibm::test::gather_live_tags_flat;
using gutibm::test::make_mpi_config;
using gutibm::test::require_mpi_ranks;
using gutibm::test::resolve_shared_test_h5_path;

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
  DomainConfig narrow_cfg = cfg;
  narrow_cfg.hi[0] = 10e-6;
  Domain narrow_domain;
  narrow_domain.init(narrow_cfg);
  bool rejected = false;
  try {
    ChemicalField narrow_field;
    narrow_field.init(narrow_domain, {spec}, "slab");
  } catch (const ConfigError& error) {
    rejected = std::string(error.what()).find("owned x-slab") != std::string::npos;
  }
  assert(rejected);
  if (rank == 0) {
    std::cout << "  test_chemical_field_layout_mapping: PASSED\n";
  }
}

void test_slab_delivery_boundary_accounting() {
  require_mpi_ranks(2);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  DomainConfig cfg;
  cfg.lo = {0.0, 0.0, 0.0};
  cfg.hi = {40e-6, 15e-6, 20e-6};
  cfg.grid_dx = 5e-6;
  cfg.grid_halo_width = 2;
  Domain domain;
  domain.init(cfg);

  ChemicalSpec spec;
  spec.name = "carbon";
  spec.diff_coeff = 2.1e-9;
  spec.initial_conc = 0.0;
  spec.boundary_conc = 1.0;
  spec.diffusion_enabled = true;
  spec.epithelial_boundary_mode = EpithelialBoundaryMode::Flux;
  spec.epithelial_flux = 1.0e-10;
  ChemicalField chem;
  chem.init(domain, {spec}, "slab");
  chem.apply_diffusion(domain, 2.0);
  chem.flux_accounting().commit_boundary_and_reaction_step();
  chem.flux_accounting().close_interval();

  const Real local_exchange =
      chem.flux_accounting().boundary_cumulative.front();
  Real global_exchange = 0.0;
  MPI_Allreduce(&local_exchange, &global_exchange, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  const Real expected = spec.epithelial_flux * domain.nx() * domain.ny()
      * domain.dx_x() * domain.dx_y() * 2.0;
  assert(std::abs(global_exchange - expected) < 1.0e-24);
  if (rank == 0) {
    std::cout << "  test_slab_delivery_boundary_accounting: PASSED\n";
  }
}

#ifdef GUTIBM_HDF5

int32_t read_event(hid_t file, const std::string& path) {
  hid_t dataset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dataset >= 0);
  int32_t value = 0;
  assert(H5Dread(dataset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, &value) >= 0);
  H5Dclose(dataset);
  return value;
}

void set_all_events(Simulation& sim, Int value) {
  StepEvents events;
  events.sos_inductions = value;
  events.phage_inductions = value;
  events.mortality_colicin = value;
  events.mortality_cdi = value;
  events.outflow_washout = value;
  events.outflow_boundary = value;
  events.mortality_lysis = value;
  events.divisions = value;
  events.conjugation_transfers = value;
  events.mutations = value;
  events.immigrations = value;
  sim.step_events() = events;
}

void test_global_event_counter_reduction() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  SimulationConfig cfg = make_mpi_config(51001, 2);
  cfg.domain.hi = {20e-6, 20e-6, 20e-6};
  cfg.hdf5.enabled = false;
  Simulation sim;
  sim.init(cfg);
  const Int local_increment = rank + 1;
  set_all_events(sim, local_increment);
  Int expected = 0;
  MPI_Allreduce(&local_increment, &expected, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);

  const std::string path = resolve_shared_test_h5_path(
      "GUTIBM_MPI_EVENT_COUNTERS_H5", "mpi_event_counters");
  HDF5Config hdf5_cfg;
  hdf5_cfg.enabled = true;
  hdf5_cfg.filename = path;
  hdf5_cfg.schedule.summary = 1;
  HDF5Writer writer;
  writer.init(hdf5_cfg, sim.domain());
  assert(writer.is_enabled());
  writer.write_step(sim, 1, 60.0, 60.0);
  set_all_events(sim, local_increment);
  writer.write_step(sim, 2, 120.0, 60.0);
  writer.finalize();
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    assert(file >= 0);
    const std::string prefix = "summary/step_000001/events/";
    assert(read_event(file, prefix + "sos_inductions") == expected);
    assert(read_event(file, prefix + "phage_inductions") == expected);
    assert(read_event(file, prefix + "mortality_colicin") == expected);
    assert(read_event(file, prefix + "mortality_cdi") == expected);
    assert(read_event(file, prefix + "outflow_washout") == expected);
    assert(read_event(file, prefix + "outflow_boundary") == expected);
    assert(read_event(file, prefix + "mortality_lysis") == expected);
    assert(read_event(file, prefix + "divisions") == expected);
    assert(read_event(file, prefix + "conjugation_transfers") == expected);
    assert(read_event(file, prefix + "mutations") == expected);
    assert(read_event(file, prefix + "immigrations") == expected);
    assert(read_event(file, prefix + "cumulative_outflow_washout") == expected);
    assert(read_event(file, prefix + "cumulative_mortality_lysis") == expected);
    const std::string second_prefix = "summary/step_000002/events/";
    assert(read_event(file, second_prefix + "outflow_washout") == expected);
    assert(read_event(file, second_prefix + "cumulative_outflow_washout")
           == 2 * expected);
    assert(read_event(file, second_prefix + "cumulative_mortality_lysis")
           == 2 * expected);
    H5Fclose(file);
    std::cout << "  test_global_event_counter_reduction: PASSED\n";
  }
}

void test_population_ledger_closure() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  SimulationConfig cfg = make_mpi_config(51002, 160);
  cfg.time.total_time = 6.0 * cfg.time.bio_dt;
  cfg.time.output_interval = cfg.time.total_time;
  cfg.enabled_fixes = {"bacteriocin"};
  cfg.fixes.bacteriocin.sos_basal_rate = 1.0;
  cfg.initial_strains[0].plasmids = {"ColE1"};
  cfg.advection.radial_turnover = 1.0e12;
  cfg.advection.distal_transit_time = 1.0e12;
  SimulationConfig::InitialStrain survivor = cfg.initial_strains[0];
  survivor.plasmids.clear();
  survivor.count = 80;
  cfg.initial_strains[0].count = 80;
  cfg.initial_strains.push_back(survivor);
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = resolve_shared_test_h5_path(
      "GUTIBM_MPI_LEDGER_H5", "mpi_population_ledger");
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 0;
  cfg.hdf5.schedule.grid = 0;
  cfg.hdf5.schedule.lineage = 0;
  cfg.hdf5.schedule.genome = 0;
  cfg.hdf5.schedule.provenance = 0;

  Simulation sim;
  sim.init(cfg);
  const Int initial = sim.global_agent_count();
  sim.run();
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    hid_t file = H5Fopen(cfg.hdf5.filename.c_str(), H5F_ACC_RDONLY,
                         H5P_DEFAULT);
    assert(file >= 0);
    const std::string step_name =
        std::format("step_{:06}", sim.step_count());
    const auto ledger = gutibm::test::read_population_ledger(file, step_name);
    gutibm::test::assert_population_ledger_closure(ledger, initial);
    H5Fclose(file);
    std::cout << "  test_population_ledger_closure: PASSED"
              << " (mortality_lysis=" << ledger.lysis << ")\n";
  }
}

void test_bacteriocin_ghost_accounting() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  SimulationConfig cfg = make_mpi_config(51003, 2);
  cfg.time.total_time = 6.0 * cfg.time.bio_dt;
  cfg.time.output_interval = cfg.time.total_time;
  cfg.enabled_fixes = {"bacteriocin"};
  cfg.fixes.bacteriocin.sos_basal_rate = 1.0;
  cfg.advection.radial_turnover = 1.0e12;
  cfg.advection.distal_transit_time = 1.0e12;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = resolve_shared_test_h5_path(
      "GUTIBM_MPI_BACTERIOCIN_GHOST_H5", "mpi_bacteriocin_ghost");

  Simulation sim;
  sim.init(cfg);
  while (sim.agents().size() > 0) {
    sim.agents().remove(sim.agents().size() - 1);
  }
  if (rank == 0) {
    Agent producer = Agent::create_default(
        sim.agents().next_tag(), 1, {49.5e-6, 50.0e-6, 25.0e-6}, 5.0e-4);
    const PlasmidEntry* plasmid = PlasmidLibrary::find("ColE1");
    assert(plasmid != nullptr);
    producer.genome.bi_loci.push_back(plasmid->cluster);
    producer.identity.owner_rank = rank;
    sim.agents().push_back(std::move(producer));
  } else {
    Agent survivor = Agent::create_default(
        sim.agents().next_tag(), 2, {75.0e-6, 50.0e-6, 25.0e-6}, 5.0e-4);
    survivor.identity.owner_rank = rank;
    sim.agents().push_back(std::move(survivor));
  }
  sim.run();

  if (rank == 0) {
    hid_t file = H5Fopen(cfg.hdf5.filename.c_str(), H5F_ACC_RDONLY,
                         H5P_DEFAULT);
    assert(file >= 0);
    const std::string step_name =
        std::format("step_{:06}", sim.step_count());
    const std::string prefix = "summary/" + step_name + "/events/";
    const Int sos = read_event(file, prefix + "cumulative_sos_inductions");
    const Int lysis = read_event(file, prefix + "cumulative_mortality_lysis");
    assert(sos == 1);
    assert(lysis == 1);
    H5Fclose(file);
    std::cout << "  test_bacteriocin_ghost_accounting: PASSED"
              << " (sos_inductions=" << sos
              << ", mortality_lysis=" << lysis << ")\n";
  }
}

void test_population_stocks_reduce_owned_agents() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  SimulationConfig cfg = make_mpi_config(51004, 4);
  cfg.time.total_time = cfg.time.bio_dt;
  cfg.time.output_interval = cfg.time.total_time;
  cfg.fixes.metabolism.bacteriostasis_threshold = 1.0e-3;
  cfg.advection.radial_turnover = 1.0e12;
  cfg.advection.distal_transit_time = 1.0e12;
  cfg.hdf5.enabled = false;

  Simulation sim;
  sim.init(cfg);
  sim.step(cfg.time.bio_dt);
  sim.prepare_population_stocks_for_summary();

  const Int expected = sim.global_agent_count();
  assert(sim.population_stocks().bacteriostatic_live == expected);
  assert(sim.population_stocks().washout_trapped_live == 0);
  if (rank == 0) {
    std::cout << "  test_population_stocks_reduce_owned_agents: PASSED"
              << " (starving=" << expected << ")\n";
  }
}

#endif

uint64_t hash_chemical_owned_cells(const ChemicalField& field,
                                   const Domain& domain);
uint64_t hash_simulation_chemistry(const Simulation& simulation);
void assert_equal_ledgers(const Simulation& slab,
                          const Simulation& replicated);

void test_slab_chemistry_transpose_halos_and_ledger() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  DomainConfig cfg;
  cfg.lo = {0.0, 0.0, 0.0};
  cfg.hi = {40e-6, 10e-6, 10e-6};
  cfg.grid_dx = 5e-6;
  cfg.ghost_width = 10e-6;
  cfg.grid_halo_width = 2;
  cfg.periodic = {true, true, false};
  Domain domain;
  domain.init(cfg);

  ChemicalSpec spec;
  spec.name = "carbon";
  spec.diff_coeff = 2.0e-9;
  spec.diffusion_enabled = true;
  spec.initial_conc = 0.0;
  ChemicalField slab;
  slab.init(domain, {spec}, "slab");
  ChemicalField replicated;
  replicated.init(domain, {spec}, "replicated");

  for (Int iz = 0; iz < domain.nz(); ++iz) {
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.nx(); ++ix) {
        const Int global = domain.cell_index(ix, iy, iz);
        const Real value = static_cast<Real>(1 + ix + 10 * iy + 100 * iz);
        replicated.conc(0, global) = value;
        if (slab.owns_global_cell(global)) {
          slab.conc_global(0, global) = value;
        }
      }
    }
  }
  slab.exchange_concentration_halos();

  const Int seam_x = rank == 0 ? domain.nx() - 1 : 0;
  const Int seam_cell = domain.cell_index(seam_x, 0, 0);
  const Real expected_seam = static_cast<Real>(1 + seam_x);
  assert(slab.global_cell_in_halo(seam_cell));
  assert(slab.conc_global(0, seam_cell) == expected_seam);

  slab.apply_diffusion(domain, 60.0);
  replicated.apply_diffusion(domain, 60.0);
  assert(hash_chemical_owned_cells(slab, domain)
         == hash_chemical_owned_cells(replicated, domain));

  SimulationConfig simulation_cfg = make_mpi_config(24680, 8);
  simulation_cfg.time.total_time = 3.0 * simulation_cfg.time.bio_dt;
  simulation_cfg.time.output_interval = simulation_cfg.time.total_time;
  simulation_cfg.domain.periodic = {true, true, false};
  simulation_cfg.domain.grid_halo_width = 2;
  simulation_cfg.chemistry_decomposition = "slab";
  Simulation slab_simulation;
  slab_simulation.init(simulation_cfg);
  slab_simulation.run();

  simulation_cfg.chemistry_decomposition = "replicated";
  Simulation replicated_simulation;
  replicated_simulation.init(simulation_cfg);
  replicated_simulation.run();

  const uint64_t slab_hash = hash_simulation_chemistry(slab_simulation);
  const uint64_t replicated_hash =
      hash_simulation_chemistry(replicated_simulation);
  if (slab_hash != replicated_hash) {
    if (rank == 0) {
      std::cerr << "slab chemistry hash " << slab_hash
                << ", replicated chemistry hash " << replicated_hash << '\n';
    }
    const auto& slab_chem = slab_simulation.chemical_field();
    const auto& replicated_chem = replicated_simulation.chemical_field();
    for (Int species = 0; species < slab_chem.num_species(); ++species) {
      Real slab_sum = 0.0;
      Real replicated_sum = 0.0;
      for (Int iz = 0; iz < slab_simulation.domain().nz(); ++iz) {
        for (Int iy = 0; iy < slab_simulation.domain().ny(); ++iy) {
          for (Int ix = slab_simulation.domain().local_grid_x_begin();
               ix < slab_simulation.domain().local_grid_x_end(); ++ix) {
            const Int cell = slab_simulation.domain().cell_index(ix, iy, iz);
            slab_sum += slab_chem.conc_global(species, cell);
            replicated_sum += replicated_chem.conc_global(species, cell);
          }
        }
      }
      Real global_slab_sum = 0.0;
      Real global_replicated_sum = 0.0;
      MPI_Allreduce(&slab_sum, &global_slab_sum, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(&replicated_sum, &global_replicated_sum, 1, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
      if (rank == 0) {
        std::cerr << "species " << species << " sums "
                  << global_slab_sum << " vs " << global_replicated_sum
                  << '\n';
      }
    }
  }
  assert(slab_hash == replicated_hash);
  assert_equal_ledgers(slab_simulation, replicated_simulation);

  if (rank == 0) {
    std::cout << "  test_slab_chemistry_transpose_halos_and_ledger: PASSED\n";
  }
}

SimulationConfig make_mpi_periodic_config() {
  SimulationConfig cfg = make_mpi_config(4242, 40);
  cfg.domain.periodic = {true, true, false};
  return cfg;
}

uint64_t hash_chemical_owned_cells(const ChemicalField& field,
                                   const Domain& domain) {
  uint64_t hash = 0;
  for (Int species = 0; species < field.num_species(); ++species) {
    for (Int iz = 0; iz < domain.nz(); ++iz) {
      for (Int iy = 0; iy < domain.ny(); ++iy) {
        for (Int ix = domain.local_grid_x_begin();
             ix < domain.local_grid_x_end(); ++ix) {
          const Int global = domain.cell_index(ix, iy, iz);
          const auto quantized = static_cast<uint64_t>(
              std::llround(field.conc_global(species, global) * 1.0e12));
          hash ^= quantized + 0x9e3779b97f4a7c15ULL
              + (hash << 6) + (hash >> 2);
        }
      }
    }
  }
#ifdef GUTIBM_MPI
  uint64_t reduced_hash = 0;
  MPI_Allreduce(&hash, &reduced_hash, 1, MPI_UINT64_T, MPI_BXOR,
                MPI_COMM_WORLD);
  hash = reduced_hash;
#endif
  return hash;
}

uint64_t hash_simulation_chemistry(const Simulation& simulation) {
  return hash_chemical_owned_cells(
      simulation.chemical_field(), simulation.domain());
}

void assert_equal_ledgers(const Simulation& slab,
                          const Simulation& replicated) {
  const auto& slab_flux = slab.chemical_field().flux_accounting();
  const auto& replicated_flux =
      replicated.chemical_field().flux_accounting();
  const auto compare = [](const char* name,
                          const std::vector<Real>& slab_values,
                          const std::vector<Real>& replicated_values,
                          Real absolute_tolerance = 0.0) {
    assert(slab_values.size() == replicated_values.size());
    for (size_t species = 0; species < slab_values.size(); ++species) {
      Real slab_min = 0.0;
      Real slab_max = 0.0;
      MPI_Allreduce(&slab_values[species], &slab_min, 1, MPI_DOUBLE,
                    MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&slab_values[species], &slab_max, 1, MPI_DOUBLE,
                    MPI_MAX, MPI_COMM_WORLD);
      Real replicated_min = 0.0;
      Real replicated_max = 0.0;
      MPI_Allreduce(&replicated_values[species], &replicated_min, 1,
                    MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&replicated_values[species], &replicated_max, 1,
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      assert(slab_min == slab_max);
      assert(replicated_min == replicated_max);
      // Ledger totals sum identical cell contributions in different orders:
      // slab reduces rank-local partials, while replicated sweeps the grid.
      const Real scale = std::max(
          std::max(std::abs(slab_min), std::abs(replicated_min)), 1.0e-300);
      if (std::abs(slab_min - replicated_min)
          > std::max(1.0e-12 * scale, absolute_tolerance)) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
          std::cerr << std::setprecision(17);
          std::cerr << "ledger mismatch " << name << '[' << species << "] "
                    << slab_min << " vs " << replicated_min << '\n';
        }
      }
      assert(std::abs(slab_min - replicated_min)
             <= std::max(1.0e-12 * scale, absolute_tolerance));
    }
  };
  compare("boundary_interval", slab_flux.boundary_interval,
          replicated_flux.boundary_interval);
  compare("boundary_cumulative", slab_flux.boundary_cumulative,
          replicated_flux.boundary_cumulative);
  compare("gradient_source_interval", slab_flux.gradient_source_interval,
          replicated_flux.gradient_source_interval, 1.0e-28);
  compare("gradient_source_cumulative", slab_flux.gradient_source_cumulative,
          replicated_flux.gradient_source_cumulative, 1.0e-28);
  compare("vbf_source_interval", slab_flux.vbf_source_interval,
          replicated_flux.vbf_source_interval);
  compare("vbf_source_cumulative", slab_flux.vbf_source_cumulative,
          replicated_flux.vbf_source_cumulative);
  compare("vbf_sink_interval", slab_flux.vbf_sink_interval,
          replicated_flux.vbf_sink_interval);
  compare("vbf_sink_cumulative", slab_flux.vbf_sink_cumulative,
          replicated_flux.vbf_sink_cumulative);
  compare("agent_uptake_interval", slab_flux.agent_uptake_interval,
          replicated_flux.agent_uptake_interval);
  compare("agent_uptake_step", slab_flux.agent_uptake_step,
          replicated_flux.agent_uptake_step);
  compare("agent_uptake_cumulative", slab_flux.agent_uptake_cumulative,
          replicated_flux.agent_uptake_cumulative);
  compare("maintenance_interval", slab_flux.maintenance_interval,
          replicated_flux.maintenance_interval);
  compare("maintenance_cumulative", slab_flux.maintenance_cumulative,
          replicated_flux.maintenance_cumulative);
  compare("uptake_demand_interval", slab_flux.uptake_demand_interval,
          replicated_flux.uptake_demand_interval);
  compare("uptake_demand_cumulative", slab_flux.uptake_demand_cumulative,
          replicated_flux.uptake_demand_cumulative);
  compare("uptake_shortfall_interval", slab_flux.uptake_shortfall_interval,
          replicated_flux.uptake_shortfall_interval);
  compare("uptake_shortfall_cumulative", slab_flux.uptake_shortfall_cumulative,
          replicated_flux.uptake_shortfall_cumulative);
  compare("reaction_clip_interval", slab_flux.reaction_clip_interval,
          replicated_flux.reaction_clip_interval);
  compare("reaction_clip_cumulative", slab_flux.reaction_clip_cumulative,
          replicated_flux.reaction_clip_cumulative);
}

void test_delivery_sink_slab_matches_replicated() {
  require_mpi_ranks(2);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  DomainConfig cfg;
  cfg.lo = {0.0, 0.0, 0.0};
  cfg.hi = {40e-6, 15e-6, 15e-6};
  cfg.grid_dx = 5e-6;
  cfg.periodic = {true, true, false};
  cfg.grid_halo_width = 2;
  Domain domain;
  domain.init(cfg);
  ChemicalSpec spec;
  spec.name = species::CARBON;
  spec.diff_coeff = 2.1e-9;
  spec.initial_conc = 1.0e-6;
  spec.boundary_conc = 1.0e-6;
  spec.diffusion_enabled = true;
  ChemicalField slab;
  ChemicalField replicated;
  slab.init(domain, {spec}, "slab");
  replicated.init(domain, {spec}, "replicated");
  const Int target = domain.cell_index(1, 1, 1);
  const Int carbon = 0;
  for (Int cell = 0; cell < domain.ncells(); ++cell) {
    replicated.conc_global(carbon, cell) = spec.initial_conc;
    if (slab.owns_global_cell(cell)) {
      slab.conc_global(carbon, cell) = spec.initial_conc;
    }
  }
  replicated.add_sink_rate_global(target, 0.2);
  if (domain.owner_rank(domain.cell_center(1, 1, 1)) == rank) {
    slab.add_sink_rate_global(target, 0.2);
  }
  slab.apply_diffusion(domain, 60.0);
  replicated.apply_diffusion(domain, 60.0);
  for (Int cell = 0; cell < domain.ncells(); ++cell) {
    if (!slab.owns_global_cell(cell)) continue;
    const Real slab_value = slab.conc_global(carbon, cell);
    const Real replicated_value = replicated.conc_global(carbon, cell);
    assert(std::abs(slab_value - replicated_value)
           <= 1.0e-12 * std::max(1.0, std::abs(replicated_value)));
  }
  Real slab_removed = 0.0;
  for (Int cell = 0; cell < domain.ncells(); ++cell) {
    if (slab.owns_global_cell(cell)) {
      slab_removed += slab.sink_realized_global(cell);
    }
  }
  Real global_slab_removed = 0.0;
  MPI_Allreduce(&slab_removed, &global_slab_removed, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  assert(std::abs(global_slab_removed
                  - replicated.sink_realized_global(target))
         <= 1.0e-12 * std::max(
             1.0, std::abs(replicated.sink_realized_global(target))));

  ChemicalSpec gradient_spec = spec;
  gradient_spec.z_gradient_enabled = true;
  gradient_spec.z_gradient_lambda = 10.0e-6;
  gradient_spec.delivery_enabled = true;
  ChemicalField gradient_slab;
  ChemicalField gradient_replicated;
  gradient_slab.init(domain, {gradient_spec}, "slab");
  gradient_replicated.init(domain, {gradient_spec}, "replicated");
  gradient_replicated.add_sink_rate_global(target, 0.2);
  if (domain.owner_rank(domain.cell_center(1, 1, 1)) == rank) {
    gradient_slab.add_sink_rate_global(target, 0.2);
  }
  gradient_slab.apply_diffusion(domain, 60.0);
  gradient_replicated.apply_diffusion(domain, 60.0);
  for (Int cell = 0; cell < domain.ncells(); ++cell) {
    if (!gradient_slab.owns_global_cell(cell)) continue;
    const Real slab_value = gradient_slab.conc_global(carbon, cell);
    const Real replicated_value =
        gradient_replicated.conc_global(carbon, cell);
    assert(std::abs(slab_value - replicated_value)
           <= 1.0e-12 * std::max(1.0, std::abs(replicated_value)));
  }
  Real gradient_slab_removed = 0.0;
  for (Int cell = 0; cell < domain.ncells(); ++cell) {
    if (gradient_slab.owns_global_cell(cell)) {
      gradient_slab_removed += gradient_slab.sink_realized_global(cell);
    }
  }
  Real global_gradient_slab_removed = 0.0;
  MPI_Allreduce(&gradient_slab_removed, &global_gradient_slab_removed, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const Real local_gradient_source =
      gradient_slab.flux_accounting().gradient_source_step.front();
  Real global_gradient_source = 0.0;
  MPI_Allreduce(&local_gradient_source, &global_gradient_source, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  assert(std::abs(global_gradient_source
                  - gradient_replicated.flux_accounting()
                      .gradient_source_step.front())
         <= 1.0e-12 * std::max(
             1.0, std::abs(
                 gradient_replicated.flux_accounting()
                     .gradient_source_step.front())));
  assert(global_gradient_slab_removed > 0.0);
  assert(std::abs(global_gradient_slab_removed
                  - gradient_replicated.sink_realized_global(target))
         <= 1.0e-12 * std::max(
             1.0, std::abs(
                 gradient_replicated.sink_realized_global(target))));

  ChemicalSpec oxygen_spec = spec;
  oxygen_spec.name = species::OXYGEN;
  oxygen_spec.delivery_enabled = true;
  ChemicalField oxygen_slab;
  ChemicalField oxygen_replicated;
  oxygen_slab.init(domain, {oxygen_spec}, "slab");
  oxygen_replicated.init(domain, {oxygen_spec}, "replicated");
  oxygen_replicated.add_sink_rate_global(0, target, 0.2);
  if (domain.owner_rank(domain.cell_center(1, 1, 1)) == rank) {
    oxygen_slab.add_sink_rate_global(0, target, 0.2);
  }
  oxygen_slab.apply_diffusion(domain, 60.0);
  oxygen_replicated.apply_diffusion(domain, 60.0);
  Real oxygen_slab_removed = 0.0;
  for (Int cell = 0; cell < domain.ncells(); ++cell) {
    if (oxygen_slab.owns_global_cell(cell)) {
      oxygen_slab_removed += oxygen_slab.sink_realized_global(0, cell);
    }
  }
  Real global_oxygen_slab_removed = 0.0;
  MPI_Allreduce(&oxygen_slab_removed, &global_oxygen_slab_removed, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  assert(global_oxygen_slab_removed > 0.0);
  assert(std::abs(global_oxygen_slab_removed
                  - oxygen_replicated.sink_realized_global(0, target))
         <= 1.0e-12 * std::max(
             1.0, std::abs(
                 oxygen_replicated.sink_realized_global(0, target))));
  if (rank == 0) {
    std::cout << "  test_delivery_sink_slab_matches_replicated: PASSED"
              << " (gradient_removed=" << global_gradient_slab_removed
              << ")\n";
  }
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

void test_cross_rank_agent_toxin_sampling() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  SimulationConfig cfg = make_mpi_config(4244, 40);
  cfg.qssa.toxin_evaluation = "agents";
  cfg.qssa.toxin_cutoff = 80e-6;
  Simulation sim;
  sim.init(cfg);

  for (Agent& agent : sim.agents()) {
    agent.genome.bi_loci.clear();
  }
  const Vec3 source_position = {49.5e-6, 50e-6, 25e-6};
  const Vec3 target_position = {57.5e-6, 50e-6, 25e-6};
  if (rank == 0) {
    sim.agents()[0].x = source_position;
    sim.agents()[0].grid_cell = sim.domain().cell_index(9, 10, 5);
    sim.agents()[0].genome.bi_loci.push_back(PlasmidLibrary::microcin_V());
  } else {
    for (Agent& agent : sim.agents()) {
      agent.x = target_position;
      agent.grid_cell = sim.domain().cell_index(11, 10, 5);
    }
  }

  sim.qssa().solve_all_bacteriocin_fields(
      sim.agents(), {}, 0.0, cfg.chem_env.protease, sim.advection(),
      sim.chemical_field(), nullptr, false);
  const Int species_idx =
      sim.chemical_field().find(species::BACTERIOCIN_CIRA);
  assert(species_idx >= 0);
  const Real local_sample =
      sim.qssa().sampled_toxin_conc(0, species_idx);
  if (rank == 1) {
    assert(local_sample > 0.0);
  }
  Real minimum = 0.0;
  MPI_Allreduce(&local_sample, &minimum, 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  assert(minimum > 0.0);

  ToxinBurstSource nuclease_burst;
  nuclease_burst.pos = source_position;
  const BICluster col_e2 = PlasmidLibrary::colicin_E2();
  nuclease_burst.params.diff_coeff = col_e2.diff_coeff;
  nuclease_burst.params.retardation = col_e2.retardation;
  nuclease_burst.params.pI = col_e2.pI;
  nuclease_burst.params.source_rate = 1.0e-18;
  nuclease_burst.is_nuclease = true;
  nuclease_burst.target = ReceptorType::BtuB;
  const std::vector<ToxinBurstSource> bursts =
      rank == 0 ? std::vector<ToxinBurstSource>{nuclease_burst}
                : std::vector<ToxinBurstSource>{};
  sim.qssa().solve_all_bacteriocin_fields(
      sim.agents(), bursts, 0.0, cfg.chem_env.protease, sim.advection(),
      sim.chemical_field(), nullptr, false);
  const Real nuclease_sample = sim.qssa().sampled_nuclease_conc(0);
  if (rank == 1) {
    assert(nuclease_sample > 0.0);
  }
  Real minimum_nuclease = 0.0;
  MPI_Allreduce(&nuclease_sample, &minimum_nuclease, 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  assert(minimum_nuclease > 0.0);

  if (rank == 1) {
    sim.agents()[0].flags.is_ghost = true;
    sim.agents()[0].genome.bi_loci.push_back(col_e2);
    BacteriocinConfig bacteriocin_cfg;
    bacteriocin_cfg.sos_basal_rate = 0.0;
    bacteriocin_cfg.sos_lysis_prob = 0.0;
    bacteriocin_cfg.sos_cross_induction_rate = 1.0e9;
    FixBacteriocin bacteriocin_fix(sim, bacteriocin_cfg);
    bacteriocin_fix.compute(60.0);
    assert(sim.agents()[0].state == PhenoState::NORMAL);
    assert(sim.step_events().sos_inductions == 0);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    std::cout << "  test_cross_rank_agent_toxin_sampling: PASSED\n";
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
  cfg.initial_population.placement = "z_slab";
  cfg.initial_population.z_min = 0.0;
  cfg.initial_population.z_max = 25e-6;
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

void test_ghost_agents_do_not_double_count_biology() {
  require_mpi_ranks(2);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  SimulationConfig cfg = make_mpi_config(4244, 1);
  cfg.enabled_fixes = {"metabolism"};
  cfg.time.bio_dt = 1.0;
  cfg.chem_env.oxygen.enabled = false;
  cfg.chem_env.acetate.enabled = false;
  cfg.chem_env.siderophore.enabled = false;
  cfg.quorum_sensing.enabled = false;
  cfg.advection.radial_turnover = 1.0e12;
  cfg.advection.distal_transit_time = 1.0e12;

  Simulation sim;
  sim.init(cfg);

  while (sim.agents().size() > 0) {
    sim.agents().remove(sim.agents().size() - 1);
  }
  sim.agents().configure_tags(
      AgentPool::first_tag_for_rank(rank, 2), AgentPool::tag_stride(2));
  if (rank == 0) {
    const Vec3 position = {
        sim.domain().local_hi_x() - 0.5 * sim.domain().ghost_width(),
        50e-6,
        25e-6};
    Agent agent = Agent::create_default(
        sim.agents().next_tag(), 1, position, cfg.initial_strains[0].mu_max);
    agent.identity.owner_rank = rank;
    sim.agents().push_back(std::move(agent));
  }

  const Real initial_mass = sphere_mass(
      CELL_RADIUS_DEFAULT, CELL_DENSITY_DEFAULT);
  TagID parent_tag = 0;
  Real initial_biomass = 0.0;
  for (Agent& agent : sim.agents()) {
    if (agent.state == PhenoState::DEAD) continue;
    parent_tag = agent.identity.tag;
    initial_biomass = 2.5 * initial_mass;
    agent.biomass = initial_biomass;
    agent.mass = agent.biomass;
    agent.radius = std::cbrt(
        3.0 * agent.biomass / (4.0 * PI * CELL_DENSITY_DEFAULT));
  }

  Int local_initial = parent_tag != 0 ? 1 : 0;
  Int global_initial = 0;
  MPI_Allreduce(&local_initial, &global_initial, 1, MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);
  assert(global_initial == 1);

  sim.step(cfg.time.bio_dt);

  assert(sim.global_agent_count() == 2);

  Real expected_uptake = 0.0;
  for (const Agent& agent : sim.agents()) {
    if (agent.identity.tag == parent_tag) {
      expected_uptake = agent.mu_realized * initial_biomass
          * cfg.time.bio_dt * cfg.fixes.metabolism.yield_carbon;
    }
  }
  Real global_expected_uptake = 0.0;
  MPI_Allreduce(&expected_uptake, &global_expected_uptake, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);

  const Int carbon = sim.chemical_field().find(species::CARBON);
  assert(carbon >= 0);
  const Real actual_uptake =
      sim.chemical_field().flux_accounting().agent_uptake_interval[
          static_cast<size_t>(carbon)];
  assert(std::abs(actual_uptake - global_expected_uptake)
         <= 1.0e-12 * std::max(1.0, std::abs(global_expected_uptake)));

  if (rank == 0) {
    std::cout << "  test_ghost_agents_do_not_double_count_biology: PASSED\n";
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

void test_multirank_closure_violation_is_synchronized() {
  require_mpi_ranks(2);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  SimulationConfig cfg = make_mpi_config(9191, 2);
  cfg.time.total_time = 600.0;
  cfg.time.output_interval = 600.0;
  cfg.hdf5.enabled = false;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.carbon_maintenance_rate = 1.0e-3;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.closure.zero_realization_grace_steps = 1;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.z_gradient_enabled = false;
      chemical.initial_conc = 0.0;
      chemical.boundary_conc = 0.0;
    }
  }
  Simulation sim;
  sim.init(cfg);
  const int status = sim.run();
  int min_status = status;
  int max_status = status;
  MPI_Allreduce(&status, &min_status, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&status, &max_status, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  const Int local_step = sim.step_count();
  Int min_step = local_step;
  Int max_step = local_step;
  MPI_Allreduce(&local_step, &min_step, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_step, &max_step, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  assert(min_status == 1 && max_status == 1);
  assert(min_step == 2 && max_step == 2);
  assert(sim.termination_cause() == TerminationCause::ClosureViolation);
  if (rank == 0) {
    std::cout << "  test_multirank_closure_violation_is_synchronized: PASSED\n";
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
  test_slab_delivery_boundary_accounting();
  test_slab_chemistry_transpose_halos_and_ledger();
  test_delivery_sink_slab_matches_replicated();
  test_cross_rank_bacteriocin_source_exchange();
  test_cross_rank_agent_toxin_sampling();
  test_slab_decomposition();
  test_slab_decomposition_periodic_x();
  test_init_population_partitioned();
  test_migration_preserves_global_count();
  test_boundary_ghost_exchange_runs();
  test_ghost_agents_do_not_double_count_biology();
  test_periodic_x_ghost_and_migration();
  test_multirank_simulation_steps();
  test_adaptive_dt_is_rank_identical();
  test_multirank_immigration_ownership();
  test_multirank_closure_violation_is_synchronized();
#ifdef GUTIBM_HDF5
  test_global_event_counter_reduction();
  test_population_ledger_closure();
  test_bacteriocin_ghost_accounting();
  test_population_stocks_reduce_owned_agents();
#endif

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
