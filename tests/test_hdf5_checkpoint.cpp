/* -----------------------------------------------------------------------
   GutIBM – HDF5 checkpoint restart tests (issues #44, #59)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "hdf5_reader.h"
#include "path_utils.h"
#include "plasmid.h"
#include "error.h"
#include "hdf5_test_helpers.h"

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
using gutibm::test::agent_snapshots_from_checkpoint;
using gutibm::test::collect_agent_snapshots;
using gutibm::test::compare_agent_snapshots;
using gutibm::test::kAgentSnapshotTol;
using gutibm::test::resolve_shared_test_h5_path;

namespace {

constexpr Real kTol = kAgentSnapshotTol;

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
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 1;
  cfg.hdf5.schedule.grid = 1;
  cfg.hdf5.schedule.lineage = 1;
  cfg.hdf5.schedule.genome = 1;
  cfg.hdf5.schedule.grid_species = {"all"};
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

#ifdef GUTIBM_MPI
void assert_checkpoint_agent_partition(
    const Simulation& resumed,
    const std::vector<gutibm::test::AgentSnapshot>& expected_agents,
    const std::vector<gutibm::test::AgentSnapshot>& local_agents,
    Int expected_count) {
  const int local_count = static_cast<int>(local_agents.size());
  int global_count = 0;
  MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  assert(global_count == expected_count);

  std::vector<int> counts(static_cast<size_t>(resumed.domain().nprocs()));
  MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                MPI_COMM_WORLD);
  std::vector<int> displacements(counts.size(), 0);
  int total_count = 0;
  for (size_t rank = 0; rank < counts.size(); ++rank) {
    displacements[rank] = total_count;
    total_count += counts[rank];
  }

  std::vector<int> byte_counts(counts.size());
  std::vector<int> byte_displacements(counts.size());
  for (size_t rank = 0; rank < counts.size(); ++rank) {
    byte_counts[rank] =
        counts[rank] * static_cast<int>(sizeof(gutibm::test::AgentSnapshot));
    byte_displacements[rank] =
        displacements[rank] * static_cast<int>(sizeof(gutibm::test::AgentSnapshot));
  }
  std::vector<gutibm::test::AgentSnapshot> all_agents(
      static_cast<size_t>(total_count));
  MPI_Allgatherv(
      local_agents.data(),
      local_count * static_cast<int>(sizeof(gutibm::test::AgentSnapshot)),
      MPI_BYTE, all_agents.data(), byte_counts.data(), byte_displacements.data(),
      MPI_BYTE, MPI_COMM_WORLD);

  for (const Agent& agent : resumed.agents()) {
    assert(resumed.domain().owner_rank(agent.x) == resumed.domain().rank());
  }
  std::ranges::sort(all_agents, [](const auto& lhs, const auto& rhs) {
    return lhs.id < rhs.id;
  });
  for (size_t i = 1; i < all_agents.size(); ++i) {
    assert(all_agents[i - 1].id != all_agents[i].id);
  }
  compare_agent_snapshots(expected_agents, all_agents);
}
#endif

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
  auto expected_agents = agent_snapshots_from_checkpoint(ckpt.agents);

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

  auto restored_agents = collect_agent_snapshots(resumed);
#ifdef GUTIBM_MPI
  assert_checkpoint_agent_partition(resumed, expected_agents, restored_agents,
                                    ckpt.metadata.num_agents);
#else
  assert(static_cast<Int>(restored_agents.size()) == ckpt.metadata.num_agents);
  compare_agent_snapshots(expected_agents, restored_agents);
#endif
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

void test_toxin_lumping_checkpoint_mismatch(const std::string& per_filename,
                                            const std::string& lumped_filename) {
  auto write_checkpoint = [](const std::string& filename,
                             std::string_view lumping) {
    SimulationConfig cfg = make_checkpoint_config(filename);
    cfg.time.total_time = 60.0;
    cfg.qssa.toxin_lumping = std::string(lumping);
    Simulation sim;
    sim.init(cfg);
    sim.run();
  };
  auto expect_mismatch = [](const std::string& filename,
                            std::string_view lumping) {
    SimulationConfig cfg = make_checkpoint_config(filename);
    cfg.qssa.toxin_lumping = std::string(lumping);
    cfg.hdf5.enabled = false;
    cfg.initial_strains.clear();
    Simulation resumed;
    bool rejected = false;
    try {
      resumed.init_from_checkpoint(cfg, filename, "step_000001");
    } catch (const ConfigError&) {
      rejected = true;
    }
    assert(rejected);
  };

  write_checkpoint(per_filename, "per_receptor");
  expect_mismatch(per_filename, "lumped");
  write_checkpoint(lumped_filename, "lumped");
  expect_mismatch(lumped_filename, "per_receptor");
}

void test_slab_checkpoint_matches_uninterrupted(const std::string& filename) {
  SimulationConfig split_cfg = make_checkpoint_config(filename + ".split.h5");
  split_cfg.chemistry_decomposition = "slab";
  split_cfg.domain.grid_halo_width = 2;
  split_cfg.time.total_time = 60.0;

  Simulation first;
  first.init(split_cfg);
  first.run();
  const HDF5CheckpointSnapshot mid =
      HDF5Reader::load_snapshot(split_cfg.hdf5.filename, "step_000001");

  SimulationConfig baseline_cfg = split_cfg;
  baseline_cfg.hdf5.enabled = false;
  baseline_cfg.time.total_time = 120.0;
  Simulation baseline;
  baseline.init(baseline_cfg);
  baseline.run();

  SimulationConfig resume_cfg = split_cfg;
  resume_cfg.hdf5.enabled = false;
  resume_cfg.time.total_time = 120.0;
  resume_cfg.initial_strains.clear();
  Simulation resumed;
  resumed.init_from_checkpoint(resume_cfg, split_cfg.hdf5.filename,
                               "step_000001");
  resumed.run();

  const auto& baseline_chem = baseline.chemical_field();
  const auto& resumed_chem = resumed.chemical_field();
  assert(baseline_chem.slab_mode());
  assert(resumed_chem.slab_mode());
  int local_chem_mismatch = 0;
  for (Int species = 0; species < baseline_chem.num_species(); ++species) {
    for (Int iz = 0; iz < baseline.domain().nz(); ++iz) {
      for (Int iy = 0; iy < baseline.domain().ny(); ++iy) {
        for (Int ix = baseline_chem.owned_global_x_begin();
             ix < baseline_chem.owned_global_x_end(); ++ix) {
          const Int cell = baseline.domain().cell_index(ix, iy, iz);
          if (baseline_chem.conc_global(species, cell) !=
              resumed_chem.conc_global(species, cell)) {
            ++local_chem_mismatch;
          }
        }
      }
    }
  }
#ifdef GUTIBM_MPI
  int global_chem_mismatch = 0;
  MPI_Allreduce(&local_chem_mismatch, &global_chem_mismatch, 1, MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);
  assert(global_chem_mismatch == 0);
#else
  assert(local_chem_mismatch == 0);
#endif
  assert(mid.metadata.step == 1);
  assert(resumed.step_count() == baseline.step_count());
  assert(resumed.time() == baseline.time());
}

void test_checkpoint_fork_immigration(const std::string& filename) {
  SimulationConfig run_cfg = make_checkpoint_config(filename);
  run_cfg.time.total_time = 60.0;
  run_cfg.initial_strains[1].count = 0;
  run_cfg.hdf5.schedule.summary = 1;
  run_cfg.hdf5.schedule.agents = 1;
  run_cfg.hdf5.schedule.genome = 1;
  run_cfg.hdf5.schedule.lineage = 1;
  Simulation baseline;
  baseline.init(run_cfg);
  baseline.run();

  const HDF5CheckpointSnapshot checkpoint =
      HDF5Reader::load_snapshot(filename, "step_000000");
  std::vector<TagID> checkpoint_tags;
  checkpoint_tags.reserve(checkpoint.agents.id.size());
  for (const int64_t id : checkpoint.agents.id) {
    checkpoint_tags.push_back(static_cast<TagID>(id));
  }

  SimulationConfig fork_cfg = run_cfg;
  fork_cfg.hdf5.enabled = false;
  fork_cfg.restart.enabled = false;
  fork_cfg.enabled_fixes = {"mechanics"};
  fork_cfg.initial_strains[1].count = 0;
  fork_cfg.initial_strains[1].plasmids = {"ColE1"};
  fork_cfg.immigration.enabled = true;
  fork_cfg.immigration.count = 2;
  fork_cfg.immigration.strain_index = 1;
  fork_cfg.immigration.step = 0;

  SimulationConfig control_cfg = fork_cfg;
  control_cfg.immigration.enabled = false;
  Simulation control;
  control.init_from_checkpoint(control_cfg, filename, "step_000000");
  control.step(60.0);

  Simulation fork;
  fork.init_from_checkpoint(fork_cfg, filename, "step_000000");
  fork.step(60.0);

  Int immigrants = 0;
  for (const Agent& agent : fork.agents()) {
    if (agent.identity.type != fork_cfg.initial_strains[1].type) continue;
    ++immigrants;
    assert(std::ranges::find(checkpoint_tags, agent.identity.tag) ==
           checkpoint_tags.end());
    assert(!agent.genome.bi_loci.empty());
    assert(std::abs(agent.mu_max - fork_cfg.initial_strains[1].mu_max) < kTol);
    assert(agent.receptor_expr[to_underlying(ReceptorType::BtuB)] > 0.0);
  }
  Int control_immigrants = 0;
  for (const Agent& agent : control.agents()) {
    if (agent.identity.type == control_cfg.initial_strains[1].type) {
      ++control_immigrants;
    }
  }
#ifdef GUTIBM_MPI
  Int global_immigrants = 0;
  Int global_control_immigrants = 0;
  MPI_Allreduce(&immigrants, &global_immigrants, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&control_immigrants, &global_control_immigrants, 1, MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);
  assert(global_immigrants ==
         global_control_immigrants + fork_cfg.immigration.count);
#else
  assert(immigrants == control_immigrants + fork_cfg.immigration.count);
#endif
  std::cout << "  test_checkpoint_fork_immigration: PASSED\n";
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
  std::string filename =
      resolve_shared_test_h5_path("GUTIBM_CHECKPOINT_H5", "checkpoint");

  if (rank == 0) std::cout << "=== HDF5 Checkpoint Restart Tests ===\n";
  test_checkpoint_restart(filename);
  test_split_run_matches_uninterrupted(filename);
  const std::string toxin_mismatch_per =
      resolve_shared_test_h5_path("GUTIBM_TOXIN_MISMATCH_PER_H5",
                                  "toxin_mismatch_per");
  const std::string toxin_mismatch_lumped =
      resolve_shared_test_h5_path("GUTIBM_TOXIN_MISMATCH_LUMPED_H5",
                                  "toxin_mismatch_lumped");
  test_toxin_lumping_checkpoint_mismatch(toxin_mismatch_per,
                                         toxin_mismatch_lumped);
  const std::string slab_filename =
      resolve_shared_test_h5_path("GUTIBM_SLAB_CHECKPOINT_H5",
                                  "slab_checkpoint");
  test_slab_checkpoint_matches_uninterrupted(slab_filename);
  const std::string immigration_filename =
      resolve_shared_test_h5_path("GUTIBM_CHECKPOINT_IMMIGRATION_H5",
                                  "checkpoint_immigration");
  test_checkpoint_fork_immigration(immigration_filename);
  if (rank == 0) {
    std::cout << "  test_hdf5_reader_api: PASSED\n";
    std::cout << "  test_checkpoint_restart: PASSED\n";
    std::cout << "  test_split_run_matches_uninterrupted: PASSED\n";
    std::cout << "  test_toxin_lumping_checkpoint_mismatch: PASSED\n";
    std::cout << "  test_slab_checkpoint_matches_uninterrupted: PASSED\n";
    std::cout << "All HDF5 checkpoint tests passed.\n";
  }
#endif

#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif
  return 0;
}
