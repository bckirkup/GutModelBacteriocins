/* -----------------------------------------------------------------------
   GutIBM – HDF5 checkpoint restart tests (issues #44, #59)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "hdf5_reader.h"
#include "path_utils.h"
#include "plasmid.h"
#include "error.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

using namespace gutibm;

namespace {

constexpr Real kTol = 1e-12;

SimulationConfig make_checkpoint_config(std::string_view filename) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.total_time = 180.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 60.0;
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

std::vector<AgentSnapshot> snapshot_from_hdf5(const HDF5CheckpointAgents& atoms) {
  const size_t n = atoms.id.size();
  std::vector<AgentSnapshot> out(n);
  size_t i = 0;
  for (int64_t id : atoms.id) {
    (void)id;
    out[i] = AgentSnapshot{
      atoms.id[i], atoms.type[i], atoms.state[i],
      atoms.x[i], atoms.y[i], atoms.z[i],
      atoms.radius[i], atoms.biomass[i], atoms.mu[i], atoms.lineage[i],
    };
    ++i;
  }
  std::ranges::sort(out, [](const AgentSnapshot& lhs, const AgentSnapshot& rhs) {
              return lhs.id < rhs.id;
            });
  return out;
}

void compare_snapshots(const std::vector<AgentSnapshot>& expected,
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

void assert_genome_bi_identity(const Simulation& sim) {
  const BICluster ref = PlasmidLibrary::colicin_E1();
  int with_bi = 0;
  for (const Agent& a : sim.agents()) {
    if (a.genome.bi_loci.empty()) continue;
    ++with_bi;
    assert(a.genome.bi_loci.size() > 0);
    const BICluster& bi = a.genome.bi_loci[0];
    assert(bi.toxin_id == ref.toxin_id);
    assert(bi.immunity_id == ref.immunity_id);
    assert(bi.target == ref.target);
    assert(std::abs(bi.pI - ref.pI) < kTol);
    assert(std::abs(bi.diff_coeff - ref.diff_coeff) < kTol);
    assert(std::abs(bi.retardation - ref.retardation) < kTol);
    assert(std::abs(bi.molecular_weight - ref.molecular_weight) < kTol);
  }
  assert(with_bi > 0);
}

void assert_genome_matches_snapshot(const Simulation& sim,
                                    const HDF5CheckpointSnapshot& snap) {
  assert(snap.genome.present);
  const size_t n = snap.agents.id.size();
  std::vector<size_t> bi_offsets(n + 1, 0);
  bi_offsets[0] = 0;
  size_t i = 0;
  for (Int num_bi : snap.lineage.num_bi_loci) {
    bi_offsets[i + 1] = bi_offsets[i] + static_cast<size_t>(num_bi);
    ++i;
  }

  int matched = 0;
  for (const Agent& a : sim.agents()) {
    size_t gi = 0;
    for (; gi < n; ++gi) {
      if (static_cast<TagID>(snap.agents.id[gi]) == a.identity.tag) break;
    }
    assert(gi < n);
    ++matched;

    assert(a.genome.parent_id == static_cast<TagID>(snap.genome.parent_id[gi]));
    assert(a.genome.mutations == static_cast<uint32_t>(snap.genome.mutations[gi]));
    assert(a.genome.has_conjugative_plasmid ==
           (snap.genome.has_conjugative_plasmid[gi] != 0));
    assert(std::abs(a.genome.plasmid_cost_amelioration -
                    snap.genome.plasmid_cost_amelioration[gi]) < kTol);
    assert(static_cast<Int>(a.genome.bi_loci.size()) == snap.lineage.num_bi_loci[gi]);

    for (const BICluster& bi : a.genome.bi_loci) {
      const size_t flat = bi_offsets[gi] + (&bi - a.genome.bi_loci.data());
      assert(bi.toxin_id == static_cast<uint16_t>(snap.genome.bi_toxin_id[flat]));
      assert(bi.immunity_id == static_cast<uint16_t>(snap.genome.bi_immunity_id[flat]));
      assert(static_cast<int32_t>(bi.target) == snap.genome.bi_target[flat]);
      assert(std::abs(bi.pI - snap.genome.bi_pI[flat]) < kTol);
      assert(std::abs(bi.immunity_binding_affinity -
                      snap.genome.bi_immunity_binding_affinity[flat]) < kTol);
    }
  }
  assert(matched > 0);
}

#ifdef GUTIBM_HDF5

void test_hdf5_reader_api(const std::string& filename) {
  HDF5Reader reader;
  if (!reader.open(filename)) {
    throw HDF5Error("HDF5Reader::open failed for " + filename);
  }

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
  assert(snap.genome.present);
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
  resume_cfg.time.total_time = ckpt.metadata.time + 120.0;
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
  assert_genome_bi_identity(resumed);
  assert_genome_matches_snapshot(resumed, ckpt);

  resumed.run();
  assert(resumed.time() > ckpt.metadata.time);
  assert(resumed.step_count() > ckpt.metadata.step);
}

void test_split_run_matches_uninterrupted(const std::string& filename) {
  SimulationConfig split_cfg = make_checkpoint_config(filename + ".split.h5");
  split_cfg.time.total_time = 60.0;
  Simulation first;
  first.init(split_cfg);
  first.run();

  HDF5CheckpointSnapshot mid =
      HDF5Reader::load_snapshot(split_cfg.hdf5.filename, "step_000001");
  assert(mid.genome.present);
  assert(mid.metadata.step == 1);
  assert(std::abs(mid.metadata.time - 60.0) < kTol);

  SimulationConfig resume_cfg = split_cfg;
  resume_cfg.hdf5.enabled = false;
  resume_cfg.time.total_time = 120.0;
  resume_cfg.initial_strains.clear();

  Simulation resumed;
  resumed.init_from_checkpoint(resume_cfg, split_cfg.hdf5.filename, "step_000001");
  assert(resumed.time() == mid.metadata.time);
  assert(resumed.step_count() == mid.metadata.step);
  assert(resumed.global_agent_count() == mid.metadata.num_agents);
  assert_genome_matches_snapshot(resumed, mid);

  resumed.run();
  assert(std::abs(resumed.time() - 120.0) < kTol);
  assert(resumed.step_count() == 2);
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
  std::string filename = resolve_test_h5_path("GUTIBM_CHECKPOINT_H5", "checkpoint");

  if (rank == 0) std::cout << "=== HDF5 Checkpoint Restart Tests ===\n";
  test_checkpoint_restart(filename);
  test_split_run_matches_uninterrupted(filename);
  if (rank == 0) {
    std::cout << "  test_hdf5_reader_api: PASSED\n";
    std::cout << "  test_checkpoint_restart: PASSED\n";
    std::cout << "  test_split_run_matches_uninterrupted: PASSED\n";
    std::cout << "All HDF5 checkpoint tests passed.\n";
  }
#endif

#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif
  return 0;
}
