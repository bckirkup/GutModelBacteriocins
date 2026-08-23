/* -----------------------------------------------------------------------
   GutIBM – HDF5 write/read round-trip tests (issue #52)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "plasmid.h"
#include "hdf5_reader.h"
#include "config_json.h"
#include "path_utils.h"
#include "hdf5_test_helpers.h"
#include "gutibm_git_sha.h"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#ifdef GUTIBM_HDF5

std::string read_string_dataset(hid_t file, const char* path) {
  hid_t dataset = H5Dopen2(file, path, H5P_DEFAULT);
  assert(dataset >= 0);
  hid_t type = H5Dget_type(dataset);
  char* value = nullptr;
  assert(H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value) >= 0);
  const std::string result = value == nullptr ? "" : value;
  H5free_memory(value);
  H5Tclose(type);
  H5Dclose(dataset);
  return result;
}

void assert_run_provenance(hid_t file,
                           const gutibm::SimulationConfig& expected) {
  assert(gutibm::test::hdf5_dataset_exists(
      file, "run_provenance/resolved_config"));
  assert(gutibm::test::hdf5_dataset_exists(file, "run_provenance/git_sha"));
  assert(std::string(GUTIBM_GIT_SHA).size() > 0);
  assert(gutibm::test::hdf5_dataset_exists(file, "run_provenance/version"));
  assert(gutibm::test::hdf5_dataset_exists(
      file, "run_provenance/mpi_rank_count"));
  assert(gutibm::test::hdf5_dataset_exists(
      file, "run_provenance/termination_cause_code"));
  assert(gutibm::test::hdf5_dataset_exists(
      file, "run_provenance/termination_cause"));
  assert(gutibm::test::hdf5_dataset_exists(
      file, "run_provenance/termination_detail"));
  assert(gutibm::test::hdf5_dataset_exists(
      file, "run_provenance/termination_wall_seconds"));
  assert(gutibm::test::hdf5_read_scalar<int32_t>(
             file, "run_provenance/termination_cause_code",
             H5T_NATIVE_INT32) == 0);
  assert(read_string_dataset(file, "run_provenance/termination_cause") ==
         "horizon_reached");
  assert(read_string_dataset(file, "run_provenance/git_sha") ==
         std::string(GUTIBM_GIT_SHA));
  assert(gutibm::test::hdf5_read_scalar<double>(
             file, "run_provenance/termination_wall_seconds",
             H5T_NATIVE_DOUBLE) >= 0.0);
  const std::string content =
      read_string_dataset(file, "run_provenance/resolved_config");
  gutibm::SimulationConfig restored = gutibm::InputParser::default_config();
  assert(gutibm::ConfigJson::parse_document(restored, content));
  gutibm::InputParser::finalize_config(restored);
  assert(restored.domain.hi == expected.domain.hi);
  assert(restored.domain.grid_dx == expected.domain.grid_dx);
  assert(restored.time.bio_dt == expected.time.bio_dt);
  assert(restored.seed == expected.seed);
  assert(restored.chemistry_decomposition == expected.chemistry_decomposition);
  assert(restored.chem_env.oxygen.delivery_uptake_enabled
         == expected.chem_env.oxygen.delivery_uptake_enabled);
  assert(restored.initial_strains.size() == expected.initial_strains.size());
}
extern "C" {
#include <hdf5.h>
}
#endif

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

using namespace gutibm;
using gutibm::test::collect_agent_snapshots;
using gutibm::test::compare_agent_snapshots;
using gutibm::test::hdf5_dataset_exists;
using gutibm::test::hdf5_read_scalar;
using gutibm::test::hdf5_read_dataset_1d;
using gutibm::test::kAgentSnapshotTol;
using gutibm::test::read_agent_snapshots;
using gutibm::test::resolve_shared_test_h5_path;

namespace {

constexpr Real kTol = kAgentSnapshotTol;

#ifdef GUTIBM_HDF5

SimulationConfig make_roundtrip_config(std::string_view filename, bool parallel) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.total_time = 120.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 60.0;
  cfg.seed = 24680;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = filename;
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 1;
  cfg.hdf5.schedule.grid = 1;
  cfg.hdf5.schedule.lineage = 1;
  cfg.hdf5.schedule.genome = 1;
  cfg.hdf5.schedule.grid_species = {"all"};
  cfg.hdf5.parallel = parallel;
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.oxygen.delivery_uptake_enabled = true;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.advection.mucus_thickness = 25e-6;
  cfg.advection.distal_length = 50e-6;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 8;
  resident.mu_max = 5e-4;
  resident.plasmids = {"ColE1"};
  resident.conjugative = false;
  cfg.initial_strains.push_back(resident);

  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = 4;
  immigrant.mu_max = 5e-4;
  immigrant.plasmids = {};
  immigrant.conjugative = false;
  cfg.initial_strains.push_back(immigrant);

  return cfg;
}

void assert_schema(hid_t file, const std::string& step) {
  assert(hdf5_dataset_exists(file, "agents/" + step + "/id"));
  assert(hdf5_dataset_exists(file, "agents/" + step + "/type"));
  assert(hdf5_dataset_exists(file, "agents/" + step + "/state"));
  assert(hdf5_dataset_exists(file, "agents/" + step + "/x"));
  assert(hdf5_dataset_exists(file, "agents/" + step + "/y"));
  assert(hdf5_dataset_exists(file, "agents/" + step + "/z"));
  assert(hdf5_dataset_exists(file, "agents/" + step + "/radius"));
  assert(hdf5_dataset_exists(file, "agents/" + step + "/biomass"));
  assert(hdf5_dataset_exists(file, "agents/" + step + "/mu_realized"));
  assert(hdf5_dataset_exists(file, "agents/" + step + "/lineage_id"));

  assert(hdf5_dataset_exists(file, "grid/" + step + "/bacteriocin_BtuB"));
  assert(hdf5_dataset_exists(file, "grid/" + step + "/carbon"));
  assert(hdf5_dataset_exists(file, "summary/" + step + "/time"));
  assert(hdf5_dataset_exists(file, "summary/" + step + "/step"));
  assert(hdf5_dataset_exists(file, "summary/" + step + "/num_agents"));
  assert(hdf5_dataset_exists(file, "summary/" + step + "/num_lineages"));
  assert(hdf5_dataset_exists(file,
                            "summary/" + step + "/mechanics/displacement_clamps"));
  assert(hdf5_dataset_exists(
      file, "summary/" + step + "/mechanics/cumulative_displacement_clamps"));
  assert(hdf5_dataset_exists(file, "lineage/" + step + "/btuB_expression"));
  assert(hdf5_dataset_exists(file, "lineage/" + step + "/num_bi_loci"));
  assert(hdf5_dataset_exists(file, "genome/" + step + "/parent_id"));
  assert(hdf5_dataset_exists(file, "genome/" + step + "/id"));
  assert(hdf5_dataset_exists(file, "genome/" + step + "/bi_offset"));
  assert(hdf5_dataset_exists(file, "genome/" + step + "/bi_count"));
  assert(hdf5_dataset_exists(file, "genome/" + step + "/bi_toxin_id"));
  assert(hdf5_dataset_exists(file, "genome/" + step + "/bi_pI"));
}

void validate_step_schema(hid_t file, const std::string& step) {
  assert_schema(file, step);
}

void validate_step_metadata(hid_t file,
                          const std::string& step,
                          Int expected_step,
                          Real expected_time,
                          Int expected_agents) {
  auto meta_agents = hdf5_read_scalar<int32_t>(file, "summary/" + step + "/num_agents",
                                               H5T_NATIVE_INT32);
  auto meta_step = hdf5_read_scalar<int32_t>(file, "summary/" + step + "/step",
                                             H5T_NATIVE_INT32);
  auto meta_time = hdf5_read_scalar<double>(file, "summary/" + step + "/time",
                                            H5T_NATIVE_DOUBLE);

  assert(meta_step == expected_step);
  assert(std::abs(meta_time - expected_time) < kTol);
  assert(meta_agents == expected_agents);

  auto file_agents = read_agent_snapshots(file, step);
  assert(static_cast<int32_t>(file_agents.size()) == meta_agents);
}

void validate_step_agents_match_sim(hid_t file,
                                    const std::string& step,
                                    const Simulation& sim) {
  assert_schema(file, step);

  auto meta_agents = hdf5_read_scalar<int32_t>(file, "summary/" + step + "/num_agents",
                                               H5T_NATIVE_INT32);
  assert(meta_agents == static_cast<int32_t>(sim.global_agent_count()));

  auto file_agents = read_agent_snapshots(file, step);
  assert(static_cast<int32_t>(file_agents.size()) == meta_agents);

  auto local_agents = collect_agent_snapshots(sim);
#ifdef GUTIBM_MPI
  if (sim.domain().nprocs() == 1) {
    compare_agent_snapshots(local_agents, file_agents);
  } else {
    assert(!file_agents.empty());
    assert(file_agents.size() == static_cast<size_t>(sim.global_agent_count()));
    for (const auto& a : file_agents) {
      assert(a.radius > 0.0);
      assert(std::isfinite(a.x));
      assert(std::isfinite(a.y));
      assert(std::isfinite(a.z));
    }
    (void)local_agents;
  }
#else
  compare_agent_snapshots(local_agents, file_agents);
#endif

  auto grid_btuB_dset = H5Dopen2(file, ("grid/" + step + "/bacteriocin_BtuB").c_str(), H5P_DEFAULT);
  assert(grid_btuB_dset >= 0);
  hid_t grid_space = H5Dget_space(grid_btuB_dset);
  std::array<hsize_t, 3> dims{0, 0, 0};
  H5Sget_simple_extent_dims(grid_space, dims.data(), nullptr);
  assert(dims[0] == static_cast<hsize_t>(sim.domain().nz()));
  assert(dims[1] == static_cast<hsize_t>(sim.domain().ny()));
  assert(dims[2] == static_cast<hsize_t>(sim.domain().nx()));
  const size_t grid_elems = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1])
      * static_cast<size_t>(dims[2]);
  H5Sclose(grid_space);
  H5Dclose(grid_btuB_dset);
  assert(grid_elems == static_cast<size_t>(sim.chemical_field().global_ncells()));
}

void run_slab_grid_pattern() {
  const std::string filename =
      resolve_test_h5_path("GUTIBM_SLAB_GRID_PATTERN_H5", "slab_grid_pattern");
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = 97531;
  cfg.chemistry_decomposition = "slab";
  cfg.domain.hi = {20e-6, 15e-6, 10e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.grid_halo_width = static_cast<Int>(
      std::ceil(cfg.domain.ghost_width / cfg.domain.grid_dx));
  cfg.time.total_time = 0.0;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = filename;
  cfg.hdf5.schedule.summary = 0;
  cfg.hdf5.schedule.agents = 0;
  cfg.hdf5.schedule.lineage = 0;
  cfg.hdf5.schedule.genome = 0;
  cfg.hdf5.schedule.grid = 1;
  cfg.hdf5.schedule.grid_species = {"carbon"};
  cfg.initial_strains.clear();

  Simulation sim;
  sim.init(cfg);
  const Int carbon = sim.chemical_field().find("carbon");
  assert(carbon >= 0);
  for (Int iz = 0; iz < sim.domain().nz(); ++iz) {
    for (Int iy = 0; iy < sim.domain().ny(); ++iy) {
      for (Int ix = sim.chemical_field().owned_global_x_begin();
           ix < sim.chemical_field().owned_global_x_end(); ++ix) {
        const Int cell = sim.domain().cell_index(ix, iy, iz);
        sim.chemical_field().conc_global(carbon, cell) =
            static_cast<Real>(1000 * iz + 100 * iy + ix);
      }
    }
  }
  sim.run();

  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  hid_t dataset =
      H5Dopen2(file, "grid/step_000000/carbon", H5P_DEFAULT);
  assert(dataset >= 0);
  hid_t space = H5Dget_space(dataset);
  std::array<hsize_t, 3> dims{0, 0, 0};
  H5Sget_simple_extent_dims(space, dims.data(), nullptr);
  assert(dims[0] == static_cast<hsize_t>(sim.domain().nz()));
  assert(dims[1] == static_cast<hsize_t>(sim.domain().ny()));
  assert(dims[2] == static_cast<hsize_t>(sim.domain().nx()));
  std::vector<double> values(static_cast<size_t>(
      sim.domain().nx() * sim.domain().ny() * sim.domain().nz()));
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          values.data());
  for (Int iz = 0; iz < sim.domain().nz(); ++iz) {
    for (Int iy = 0; iy < sim.domain().ny(); ++iy) {
      for (Int ix = 0; ix < sim.domain().nx(); ++ix) {
        const size_t index = static_cast<size_t>(
            iz * sim.domain().ny() * sim.domain().nx()
            + iy * sim.domain().nx() + ix);
        assert(values[index] == static_cast<double>(1000 * iz + 100 * iy + ix));
      }
    }
  }
  H5Sclose(space);
  H5Dclose(dataset);
  H5Fclose(file);
}

void validate_genome_slices(hid_t file, const std::string& step,
                            const Simulation& sim) {
  const std::string prefix = "genome/" + step + "/";
  const auto ids = hdf5_read_dataset_1d<int64_t>(
      file, prefix + "id", H5T_NATIVE_INT64);
  const auto offsets = hdf5_read_dataset_1d<int64_t>(
      file, prefix + "bi_offset", H5T_NATIVE_INT64);
  const auto counts = hdf5_read_dataset_1d<int32_t>(
      file, prefix + "bi_count", H5T_NATIVE_INT32);
  const auto toxin_ids = hdf5_read_dataset_1d<int32_t>(
      file, prefix + "bi_toxin_id", H5T_NATIVE_INT32);
  const auto immunity_ids = hdf5_read_dataset_1d<int32_t>(
      file, prefix + "bi_immunity_id", H5T_NATIVE_INT32);
  assert(ids.size() == offsets.size());
  assert(ids.size() == counts.size());
  assert(toxin_ids.size() == immunity_ids.size());

  int64_t previous_end = 0;
  for (size_t i = 0; i < ids.size(); ++i) {
    assert(offsets[i] >= previous_end);
    assert(counts[i] >= 0);
    const auto start = static_cast<size_t>(offsets[i]);
    const auto count = static_cast<size_t>(counts[i]);
    assert(start + count <= toxin_ids.size());
    previous_end = offsets[i] + counts[i];
  }
  assert(previous_end == static_cast<int64_t>(toxin_ids.size()));

  for (const Agent& agent : sim.agents()) {
    size_t index = 0;
    for (; index < ids.size(); ++index) {
      if (ids[index] == agent.identity.tag) break;
    }
    assert(index < ids.size());
    assert(counts[index] == static_cast<int32_t>(agent.genome.bi_loci.size()));
    const auto start = static_cast<size_t>(offsets[index]);
    for (size_t locus = 0; locus < agent.genome.bi_loci.size(); ++locus) {
      assert(toxin_ids[start + locus] ==
             static_cast<int32_t>(agent.genome.bi_loci[locus].toxin_id));
      assert(immunity_ids[start + locus] ==
             static_cast<int32_t>(agent.genome.bi_loci[locus].immunity_id));
    }
  }
}

void validate_step_genome(hid_t file, const std::string& step,
                          const Simulation& sim) {
  const BICluster ref = PlasmidLibrary::colicin_E1();
  int with_bi = 0;
  for (const Agent& a : sim.agents()) {
    if (a.genome.bi_loci.empty()) continue;
    ++with_bi;
    const BICluster& bi = a.genome.bi_loci[0];
    assert(bi.toxin_id == ref.toxin_id);
    assert(bi.immunity_id == ref.immunity_id);
    assert(bi.target == ref.target);
    assert(std::abs(bi.pI - ref.pI) < kTol);
    assert(std::abs(bi.diff_coeff - ref.diff_coeff) < kTol);
    assert(std::abs(bi.retardation - ref.retardation) < kTol);
  }
#ifdef GUTIBM_MPI
  int global_with_bi = 0;
  MPI_Allreduce(&with_bi, &global_with_bi, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  assert(global_with_bi > 0);
#else
  assert(with_bi > 0);
#endif
  if (file >= 0) {
    validate_genome_slices(file, step, sim);
  }
}

void validate_step(hid_t file,
                   const std::string& step,
                   const Simulation& sim,
                   Int expected_step,
                   Real expected_time) {
  validate_step_schema(file, step);
  validate_step_metadata(file, step, expected_step, expected_time,
                         sim.global_agent_count());
  validate_step_agents_match_sim(file, step, sim);
  validate_step_genome(file, step, sim);
}

void validate_parallel_roundtrip(const Simulation& sim, const std::string& filename) {
  HDF5CheckpointSnapshot snap = HDF5Reader::load_snapshot(filename, "step_000002");
  assert(snap.genome.present);
  assert(snap.metadata.step == 2);
  assert(std::abs(snap.metadata.time - 120.0) < kTol);
  assert(snap.metadata.num_agents == sim.global_agent_count());
  assert(static_cast<Int>(snap.agents.id.size()) == sim.global_agent_count());

  auto grid_it = snap.grid.species.find("bacteriocin_BtuB");
  assert(grid_it != snap.grid.species.end());
  assert(grid_it->second.size() == static_cast<size_t>(sim.chemical_field().ncells()));

  for (const Agent& a : sim.agents()) {
    size_t j = 0;
    for (; j < snap.agents.id.size(); ++j) {
      if (snap.agents.id[j] == a.identity.tag) break;
    }
    assert(j < snap.agents.id.size());
    assert(std::abs(snap.agents.x[j] - a.x[0]) < kTol);
  }

#ifdef GUTIBM_MPI
  int with_bi = 0;
  for (const Agent& a : sim.agents()) {
    if (!a.genome.bi_loci.empty()) ++with_bi;
  }
  int global_with_bi = 0;
  MPI_Allreduce(&with_bi, &global_with_bi, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  assert(global_with_bi > 0);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    assert(file >= 0);
    validate_genome_slices(file, "step_000002", sim);
    H5Fclose(file);
  }
#else
  validate_step_genome(static_cast<hid_t>(-1), "step_000002", sim);
#endif
}

void run_roundtrip(bool parallel_io) {
  std::string filename = resolve_shared_test_h5_path(
      "GUTIBM_ROUNDTRIP_H5",
      parallel_io ? "roundtrip_parallel" : "roundtrip_serial");

  SimulationConfig cfg = make_roundtrip_config(filename, parallel_io);
  Simulation sim;
  sim.init(cfg);
  sim.run();

#ifdef GUTIBM_MPI
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs > 1) {
    MPI_Barrier(MPI_COMM_WORLD);
    validate_parallel_roundtrip(sim, filename);
    return;
  }
#endif

  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);

  validate_step_schema(file, "step_000000");
  validate_step_metadata(file, "step_000000", 0, 0.0, 12);
  validate_step(file, "step_000002", sim, 2, 120.0);
  assert_run_provenance(file, cfg);
  const std::string original_config =
      read_string_dataset(file, "run_provenance/resolved_config");

  H5Fclose(file);
  SimulationConfig distinct_cfg = make_roundtrip_config(filename + ".distinct",
                                                         false);
  distinct_cfg.seed += 1;
  Simulation distinct;
  distinct.init(distinct_cfg);
  distinct.run();
  hid_t distinct_file =
      H5Fopen(distinct_cfg.hdf5.filename.c_str(), H5F_ACC_RDONLY,
              H5P_DEFAULT);
  assert(distinct_file >= 0);
  assert_run_provenance(distinct_file, distinct_cfg);
  assert(original_config
         != read_string_dataset(distinct_file, "run_provenance/resolved_config"));
  H5Fclose(distinct_file);
}

#ifdef GUTIBM_MPI
void run_parallel_slab_grid_equivalence() {
  const std::string slab_file = resolve_shared_test_h5_path(
      "GUTIBM_SLAB_ROUNDTRIP_H5", "roundtrip_slab");
  const std::string replicated_file = resolve_shared_test_h5_path(
      "GUTIBM_REPLICATED_ROUNDTRIP_H5", "roundtrip_replicated");

  SimulationConfig slab_cfg = make_roundtrip_config(slab_file, false);
  slab_cfg.chemistry_decomposition = "slab";
  slab_cfg.domain.grid_halo_width = 2;
  {
    Simulation slab;
    slab.init(slab_cfg);
    slab.run();
  }

  SimulationConfig replicated_cfg = make_roundtrip_config(replicated_file, false);
  replicated_cfg.chemistry_decomposition = "replicated";
  {
    Simulation replicated;
    replicated.init(replicated_cfg);
    replicated.run();
  }

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) return;
  hid_t slab_file_id = H5Fopen(slab_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t replicated_file_id =
      H5Fopen(replicated_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(slab_file_id >= 0);
  assert(replicated_file_id >= 0);
  assert_run_provenance(slab_file_id, slab_cfg);
  assert_run_provenance(replicated_file_id, replicated_cfg);
  assert(read_string_dataset(slab_file_id, "run_provenance/resolved_config")
         != read_string_dataset(replicated_file_id,
                                "run_provenance/resolved_config"));
  const std::string path = "grid/step_000000/carbon";
  hid_t slab_dataset = H5Dopen2(slab_file_id, path.c_str(), H5P_DEFAULT);
  hid_t replicated_dataset =
      H5Dopen2(replicated_file_id, path.c_str(), H5P_DEFAULT);
  assert(slab_dataset >= 0);
  assert(replicated_dataset >= 0);
  hid_t slab_space = H5Dget_space(slab_dataset);
  hid_t replicated_space = H5Dget_space(replicated_dataset);
  std::array<hsize_t, 3> slab_dims{0, 0, 0};
  std::array<hsize_t, 3> replicated_dims{0, 0, 0};
  H5Sget_simple_extent_dims(slab_space, slab_dims.data(), nullptr);
  H5Sget_simple_extent_dims(replicated_space, replicated_dims.data(), nullptr);
  assert(slab_dims == replicated_dims);
  const size_t count = static_cast<size_t>(slab_dims[0] * slab_dims[1]
                                            * slab_dims[2]);
  std::vector<double> slab_values(count);
  std::vector<double> replicated_values(count);
  H5Dread(slab_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          slab_values.data());
  H5Dread(replicated_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          replicated_values.data());
  for (size_t i = 0; i < count; ++i) {
    assert(slab_values[i] == replicated_values[i]);
  }
  H5Sclose(slab_space);
  H5Sclose(replicated_space);
  H5Dclose(slab_dataset);
  H5Dclose(replicated_dataset);
  H5Fclose(slab_file_id);
  H5Fclose(replicated_file_id);
}
#endif

#endif  // GUTIBM_HDF5

}  // namespace

int main(int argc, char** argv) {
#ifdef GUTIBM_MPI
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int rank = 0;
#endif

  (void)argc;
  (void)argv;

#ifndef GUTIBM_HDF5
  if (rank == 0) {
    std::cout << "HDF5 disabled at build time — skipping round-trip tests.\n";
  }
#else
#ifdef GUTIBM_MPI
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if (nprocs == 1) {
    if (rank == 0) std::cout << "=== HDF5 Serial Round-Trip Tests ===\n";
    run_roundtrip(false);
    run_slab_grid_pattern();
    if (rank == 0) {
      std::cout << "  test_serial_roundtrip: PASSED\n";
      std::cout << "  test_serial_slab_grid_pattern: PASSED\n";
      std::cout << "All HDF5 round-trip tests passed.\n";
    }
  } else {
    if (rank == 0) std::cout << "=== HDF5 Parallel Round-Trip Tests ===\n";
    run_roundtrip(true);
    run_parallel_slab_grid_equivalence();
    if (rank == 0) {
      std::cout << "  test_parallel_roundtrip: PASSED\n";
      std::cout << "All HDF5 round-trip tests passed.\n";
    }
  }
#else
  if (rank == 0) std::cout << "=== HDF5 Serial Round-Trip Tests ===\n";
  run_roundtrip(false);
  run_slab_grid_pattern();
  if (rank == 0) {
    std::cout << "  test_serial_roundtrip: PASSED\n";
    std::cout << "  test_serial_slab_grid_pattern: PASSED\n";
    std::cout << "All HDF5 round-trip tests passed.\n";
  }
#endif
#endif

#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif
  return 0;
}
