/* -----------------------------------------------------------------------
   GutIBM – HDF5 write/read round-trip tests (issue #52)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "plasmid.h"
#include "hdf5_reader.h"
#include "path_utils.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

using namespace gutibm;

namespace {

constexpr Real kTol = 1e-12;

#ifdef GUTIBM_HDF5

bool dataset_exists(hid_t file, const std::string& path) {
  return H5Lexists(file, path.c_str(), H5P_DEFAULT) > 0;
}

template <typename T>
std::vector<T> read_dataset_1d(hid_t file, const std::string& path, hid_t h5_type) {
  hid_t dset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dset >= 0);

  hid_t space = H5Dget_space(dset);
  int ndims = H5Sget_simple_extent_ndims(space);
  assert(ndims == 1);
  hsize_t len = 0;
  H5Sget_simple_extent_dims(space, &len, nullptr);

  std::vector<T> data(static_cast<size_t>(len));
  H5Dread(dset, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

  H5Sclose(space);
  H5Dclose(dset);
  return data;
}

template <typename T>
T read_scalar(hid_t file, const std::string& path, hid_t h5_type) {
  auto data = read_dataset_1d<T>(file, path, h5_type);
  assert(data.size() == 1);
  return data[0];
}

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

struct AgentSnapshot {
  int64_t id;
  int32_t type;
  int32_t state;
  double x;
  double y;
  double z;
  double radius;
  double biomass;
  double mu;
  int64_t lineage;
};

std::vector<AgentSnapshot> collect_all_agents(const Simulation& sim) {
  std::vector<AgentSnapshot> out;
  out.reserve(sim.agents().size());
  for (const Agent& a : sim.agents()) {
    out.emplace_back(
      a.identity.tag,
      a.identity.type,
      static_cast<int32_t>(to_underlying(a.state)),
      a.x[0], a.x[1], a.x[2],
      a.radius,
      a.biomass,
      a.mu_realized,
      a.genome.lineage_id);
  }
  std::ranges::sort(out, [](const AgentSnapshot& lhs, const AgentSnapshot& rhs) {
              return lhs.id < rhs.id;
            });
  return out;
}

std::vector<AgentSnapshot> read_agent_snapshots(hid_t file, const std::string& step) {
  const std::string prefix = "agents/" + step + "/";
  auto ids      = read_dataset_1d<int64_t>(file, prefix + "id",      H5T_NATIVE_INT64);
  auto types    = read_dataset_1d<int32_t>(file, prefix + "type",    H5T_NATIVE_INT32);
  auto states   = read_dataset_1d<int32_t>(file, prefix + "state",   H5T_NATIVE_INT32);
  auto xs       = read_dataset_1d<double>(file, prefix + "x",        H5T_NATIVE_DOUBLE);
  auto ys       = read_dataset_1d<double>(file, prefix + "y",        H5T_NATIVE_DOUBLE);
  auto zs       = read_dataset_1d<double>(file, prefix + "z",        H5T_NATIVE_DOUBLE);
  auto radii    = read_dataset_1d<double>(file, prefix + "radius",   H5T_NATIVE_DOUBLE);
  auto biomass  = read_dataset_1d<double>(file, prefix + "biomass",  H5T_NATIVE_DOUBLE);
  auto mus      = read_dataset_1d<double>(file, prefix + "mu_realized", H5T_NATIVE_DOUBLE);
  auto lineages = read_dataset_1d<int64_t>(file, prefix + "lineage_id", H5T_NATIVE_INT64);

  const size_t n = ids.size();
  assert(types.size() == n);
  assert(states.size() == n);
  assert(xs.size() == n);
  assert(ys.size() == n);
  assert(zs.size() == n);
  assert(radii.size() == n);
  assert(biomass.size() == n);
  assert(mus.size() == n);
  assert(lineages.size() == n);

  std::vector<AgentSnapshot> out(n);
  size_t i = 0;
  for (int64_t id : ids) {
    (void)id;
    out[i] = AgentSnapshot{
      ids[i], types[i], states[i],
      xs[i], ys[i], zs[i],
      radii[i], biomass[i], mus[i], lineages[i],
    };
    ++i;
  }
  std::ranges::sort(out, [](const AgentSnapshot& lhs, const AgentSnapshot& rhs) {
              return lhs.id < rhs.id;
            });
  return out;
}

void assert_schema(hid_t file, const std::string& step) {
  assert(dataset_exists(file, "agents/" + step + "/id"));
  assert(dataset_exists(file, "agents/" + step + "/type"));
  assert(dataset_exists(file, "agents/" + step + "/state"));
  assert(dataset_exists(file, "agents/" + step + "/x"));
  assert(dataset_exists(file, "agents/" + step + "/y"));
  assert(dataset_exists(file, "agents/" + step + "/z"));
  assert(dataset_exists(file, "agents/" + step + "/radius"));
  assert(dataset_exists(file, "agents/" + step + "/biomass"));
  assert(dataset_exists(file, "agents/" + step + "/mu_realized"));
  assert(dataset_exists(file, "agents/" + step + "/lineage_id"));

  assert(dataset_exists(file, "grid/" + step + "/bacteriocin_BtuB"));
  assert(dataset_exists(file, "grid/" + step + "/carbon"));
  assert(dataset_exists(file, "summary/" + step + "/time"));
  assert(dataset_exists(file, "summary/" + step + "/step"));
  assert(dataset_exists(file, "summary/" + step + "/num_agents"));
  assert(dataset_exists(file, "summary/" + step + "/num_lineages"));
  assert(dataset_exists(file, "lineage/" + step + "/btuB_expression"));
  assert(dataset_exists(file, "lineage/" + step + "/num_bi_loci"));
  assert(dataset_exists(file, "genome/" + step + "/parent_id"));
  assert(dataset_exists(file, "genome/" + step + "/bi_toxin_id"));
  assert(dataset_exists(file, "genome/" + step + "/bi_pI"));
}

void compare_agent_snapshots(const std::vector<AgentSnapshot>& expected,
                             const std::vector<AgentSnapshot>& actual) {
  assert(expected.size() == actual.size());
  auto it_actual = actual.begin();
  for (const AgentSnapshot& exp : expected) {
    const AgentSnapshot& act = *it_actual++;
    assert(exp.id == act.id);
    assert(exp.type == act.type);
    assert(exp.state == act.state);
    assert(std::abs(exp.x - act.x) < kTol);
    assert(std::abs(exp.y - act.y) < kTol);
    assert(std::abs(exp.z - act.z) < kTol);
    assert(std::abs(exp.radius - act.radius) < kTol);
    assert(std::abs(exp.biomass - act.biomass) < kTol);
    assert(std::abs(exp.mu - act.mu) < kTol);
    assert(exp.lineage == act.lineage);
  }
}

void validate_step_schema(hid_t file, const std::string& step) {
  assert_schema(file, step);
}

void validate_step_metadata(hid_t file,
                          const std::string& step,
                          Int expected_step,
                          Real expected_time,
                          Int expected_agents) {
  auto meta_agents = read_scalar<int32_t>(file, "summary/" + step + "/num_agents",
                                          H5T_NATIVE_INT32);
  auto meta_step   = read_scalar<int32_t>(file, "summary/" + step + "/step",
                                          H5T_NATIVE_INT32);
  auto meta_time   = read_scalar<double>(file, "summary/" + step + "/time",
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

  auto meta_agents = read_scalar<int32_t>(file, "summary/" + step + "/num_agents",
                                          H5T_NATIVE_INT32);
  assert(meta_agents == static_cast<int32_t>(sim.global_agent_count()));

  auto file_agents = read_agent_snapshots(file, step);
  assert(static_cast<int32_t>(file_agents.size()) == meta_agents);

  auto local_agents = collect_all_agents(sim);
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
  hsize_t dims[3] = {0, 0, 0};
  H5Sget_simple_extent_dims(grid_space, dims, nullptr);
  const size_t grid_elems = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1])
      * static_cast<size_t>(dims[2]);
  H5Sclose(grid_space);
  H5Dclose(grid_btuB_dset);
  assert(grid_elems == static_cast<size_t>(sim.chemical_field().ncells()));
}

void validate_step_genome(hid_t /*file*/, const std::string& /*step*/,
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
#else
  validate_step_genome(static_cast<hid_t>(-1), "step_000002", sim);
#endif
}

std::string resolve_shared_test_h5_path(const char* env_var, const std::string& tag) {
#ifdef GUTIBM_MPI
  int rank = 0;
  int nprocs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  std::string filename;
  if (rank == 0) {
    filename = resolve_test_h5_path(env_var, tag);
  }

  if (nprocs > 1) {
    int len = 0;
    if (rank == 0) len = static_cast<int>(filename.size());
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) filename.resize(static_cast<size_t>(len));
    if (len > 0) {
      MPI_Bcast(filename.data(), len, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
  }

  return filename;
#else
  return resolve_test_h5_path(env_var, tag);
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

  H5Fclose(file);
}

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
    if (rank == 0) {
      std::cout << "  test_serial_roundtrip: PASSED\n";
      std::cout << "All HDF5 round-trip tests passed.\n";
    }
  } else {
    if (rank == 0) std::cout << "=== HDF5 Parallel Round-Trip Tests ===\n";
    run_roundtrip(true);
    if (rank == 0) {
      std::cout << "  test_parallel_roundtrip: PASSED\n";
      std::cout << "All HDF5 round-trip tests passed.\n";
    }
  }
#else
  if (rank == 0) std::cout << "=== HDF5 Serial Round-Trip Tests ===\n";
  run_roundtrip(false);
  if (rank == 0) {
    std::cout << "  test_serial_roundtrip: PASSED\n";
    std::cout << "All HDF5 round-trip tests passed.\n";
  }
#endif
#endif

#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif
  return 0;
}
