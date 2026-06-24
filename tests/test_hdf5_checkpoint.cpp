/* -----------------------------------------------------------------------
   GutIBM – HDF5 checkpoint restart tests (issues #44, #59)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "hdf5_reader.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

using namespace gutibm;

namespace {

constexpr Real kTol = 1e-12;

SimulationConfig make_checkpoint_config(const std::string& filename) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.total_time = 180.0;
  cfg.bio_dt = 60.0;
  cfg.output_interval = 60.0;
  cfg.seed = 13579;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = filename;
  cfg.hdf5.dump_every = 1;
  cfg.hdf5.parallel = false;
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

std::vector<AgentSnapshot> collect_agents(const Simulation& sim) {
  std::vector<AgentSnapshot> out;
  out.reserve(sim.agents().size());
  for (Int i = 0; i < sim.agents().size(); ++i) {
    const Agent& a = sim.agents()[i];
    out.push_back(AgentSnapshot{
      a.tag,
      a.type,
      static_cast<int32_t>(a.state),
      a.x[0], a.x[1], a.x[2],
      a.radius,
      a.biomass,
      a.mu_realized,
      a.genome.lineage_id,
    });
  }
  std::sort(out.begin(), out.end(),
            [](const AgentSnapshot& lhs, const AgentSnapshot& rhs) {
              return lhs.id < rhs.id;
            });
  return out;
}

std::vector<AgentSnapshot> snapshot_from_hdf5(const HDF5CheckpointAgents& atoms) {
  const size_t n = atoms.id.size();
  std::vector<AgentSnapshot> out(n);
  for (size_t i = 0; i < n; ++i) {
    out[i] = AgentSnapshot{
      atoms.id[i], atoms.type[i], atoms.state[i],
      atoms.x[i], atoms.y[i], atoms.z[i],
      atoms.radius[i], atoms.biomass[i], atoms.mu[i], atoms.lineage[i],
    };
  }
  std::sort(out.begin(), out.end(),
            [](const AgentSnapshot& lhs, const AgentSnapshot& rhs) {
              return lhs.id < rhs.id;
            });
  return out;
}

void compare_snapshots(const std::vector<AgentSnapshot>& expected,
                       const std::vector<AgentSnapshot>& actual) {
  assert(expected.size() == actual.size());
  for (size_t i = 0; i < expected.size(); ++i) {
    assert(expected[i].id == actual[i].id);
    assert(expected[i].type == actual[i].type);
    assert(expected[i].state == actual[i].state);
    assert(std::abs(expected[i].x - actual[i].x) < kTol);
    assert(std::abs(expected[i].y - actual[i].y) < kTol);
    assert(std::abs(expected[i].z - actual[i].z) < kTol);
    assert(std::abs(expected[i].radius - actual[i].radius) < kTol);
    assert(std::abs(expected[i].biomass - actual[i].biomass) < kTol);
    assert(std::abs(expected[i].mu - actual[i].mu) < kTol);
    assert(expected[i].lineage == actual[i].lineage);
  }
}

#ifdef GUTIBM_HDF5

void test_hdf5_reader_api(const std::string& filename) {
  HDF5Reader reader;
  assert(reader.open(filename));

  auto steps = reader.list_steps();
  assert(steps.size() >= 3);
  assert(steps.front() == "step_000000");
  assert(reader.latest_step() == steps.back());

  auto snap = reader.load_step("step_000002");
  assert(snap.step_name == "step_000002");
  assert(snap.metadata.step == 2);
  assert(std::abs(snap.metadata.time - 120.0) < kTol);
  assert(snap.metadata.num_agents == static_cast<Int>(snap.agents.id.size()));
  assert(snap.metadata.num_agents == static_cast<Int>(snap.lineage.generation.size()));
  assert(snap.grid.species.count("carbon") > 0);
  assert(snap.grid.species.at("carbon").size() > 0);

  reader.close();
  assert(!reader.is_open());
}

void test_checkpoint_restart(const std::string& filename) {
  SimulationConfig run_cfg = make_checkpoint_config(filename);
  Simulation baseline;
  baseline.init(run_cfg);
  baseline.run();

  test_hdf5_reader_api(filename);

  HDF5CheckpointSnapshot ckpt = HDF5Reader::load_snapshot(filename, "step_000000");
  assert(ckpt.metadata.num_agents > 0);
  auto expected_agents = snapshot_from_hdf5(ckpt.agents);

  SimulationConfig resume_cfg = run_cfg;
  resume_cfg.hdf5.enabled = false;
  resume_cfg.total_time = ckpt.metadata.time + 120.0;
  resume_cfg.initial_strains.clear();

  Simulation resumed;
  resumed.init_from_checkpoint(resume_cfg, filename, "step_000000");

  assert(resumed.time() == ckpt.metadata.time);
  assert(resumed.step_count() == ckpt.metadata.step);
  assert(resumed.global_agent_count() == ckpt.metadata.num_agents);
  assert(resumed.global_agent_count() > 0);

  auto restored_agents = collect_agents(resumed);
  assert(static_cast<Int>(restored_agents.size()) == ckpt.metadata.num_agents);
  compare_snapshots(expected_agents, restored_agents);

  resumed.run();
  assert(resumed.time() > ckpt.metadata.time);
  assert(resumed.step_count() > ckpt.metadata.step);
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
    std::cout << "HDF5 disabled at build time — skipping checkpoint tests.\n";
  }
#else
  const char* tmpdir = std::getenv("TMPDIR");
  const char* env_path = std::getenv("GUTIBM_CHECKPOINT_H5");
  std::string base = tmpdir ? tmpdir : "/tmp";
  std::string filename = env_path ? env_path : base + "/gutibm_checkpoint_test.h5";

  if (rank == 0) std::cout << "=== HDF5 Checkpoint Restart Tests ===\n";
  test_checkpoint_restart(filename);
  if (rank == 0) {
    std::cout << "  test_hdf5_reader_api: PASSED\n";
    std::cout << "  test_checkpoint_restart: PASSED\n";
    std::cout << "All HDF5 checkpoint tests passed.\n";
  }
#endif

#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif
  return 0;
}
