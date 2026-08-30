/* -----------------------------------------------------------------------
   GutIBM – Closed midstream restart artifacts (Tier 2: agents + grid)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "hdf5_reader.h"
#include "hdf5_writer.h"
#include "path_utils.h"
#include "error.h"

#include <cassert>
#include <cmath>
#include <filesystem>
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

namespace fs = std::filesystem;

SimulationConfig make_restart_config(std::string_view out_h5,
                                     std::string_view restart_dir) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.total_time = 180.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 60.0;
  cfg.seed = 424242;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = std::string(out_h5);
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 1;
  cfg.hdf5.schedule.grid = 0;  // analysis trail may omit grid
  cfg.hdf5.schedule.lineage = 1;
  cfg.hdf5.schedule.genome = 1;
  cfg.advection.mucus_thickness = 25e-6;
  cfg.advection.distal_length = 50e-6;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  cfg.restart.enabled = true;
  cfg.restart.directory = std::string(restart_dir);
  cfg.restart.interval_steps = 2;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 6;
  resident.mu_max = 5e-4;
  resident.plasmids = {"ColE1"};
  resident.conjugative = false;
  cfg.initial_strains.push_back(resident);

  return cfg;
}

double grid_l2(const HDF5CheckpointGrid& a, const HDF5CheckpointGrid& b) {
  assert(!a.species.empty());
  assert(a.species.size() == b.species.size());
  double sum = 0.0;
  for (const auto& [name, vals] : a.species) {
    const auto it = b.species.find(name);
    assert(it != b.species.end());
    assert(vals.size() == it->second.size());
    for (size_t i = 0; i < vals.size(); ++i) {
      const double d = vals[i] - it->second[i];
      sum += d * d;
    }
  }
  return std::sqrt(sum);
}

void test_closed_restart_roundtrip() {
  const std::string out = resolve_test_h5_path("GUTIBM_TEST_RESTART_H5", "restart_out");
  const fs::path scratch = fs::path(out).parent_path() / "gutibm_restart_artifacts";
  const fs::path restart_dir = scratch / "a";
  fs::create_directories(restart_dir);

  SimulationConfig cfg = make_restart_config(out, restart_dir.string());
  Simulation sim;
  sim.init(cfg);
  sim.run();

  const fs::path step2 = restart_dir / "step_000002.h5";
  const fs::path step3 = restart_dir / "step_000003.h5";  // final
  assert(fs::exists(step2));
  assert(fs::exists(step3));
  assert(fs::file_size(step2) > 4096);

  HDF5CheckpointSnapshot snap = HDF5Reader::load_snapshot(step2.string(), "");
  assert(snap.step_name == "step_000002");
  assert(std::abs(snap.metadata.time - 120.0) < 1e-9);
  assert(snap.metadata.step == 2);
  assert(snap.genome.present);
  assert(!snap.agents.id.empty());
  assert(!snap.grid.species.empty());

  const std::string resume_out =
      resolve_test_h5_path("GUTIBM_TEST_RESTART_RESUME_H5", "restart_resume");
  SimulationConfig resume_cfg =
      make_restart_config(resume_out, (scratch / "resume").string());
  resume_cfg.hdf5.enabled = false;
  resume_cfg.restart.enabled = false;
  resume_cfg.time.total_time = 120.0;
  resume_cfg.checkpoint.file = step2.string();

  Simulation resumed;
  resumed.init_from_checkpoint(resume_cfg, step2.string(), "");
  assert(std::abs(resumed.time() - 120.0) < 1e-9);
  assert(resumed.step_count() == 2);
  assert(resumed.global_agent_count() == snap.metadata.num_agents);

  HDF5CheckpointSnapshot snap_reload =
      HDF5Reader::load_snapshot(step2.string(), "step_000002");
  assert(grid_l2(snap.grid, snap_reload.grid) < 1e-12);

  // Sensitivity: changing restart.interval_steps must change which files exist.
  const fs::path restart_dir_b = scratch / "b";
  fs::create_directories(restart_dir_b);
  SimulationConfig cfg_b = make_restart_config(
      resolve_test_h5_path("GUTIBM_TEST_RESTART_H5_B", "restart_out_b"),
      restart_dir_b.string());
  cfg_b.restart.interval_steps = 3;
  Simulation sim_b;
  sim_b.init(cfg_b);
  sim_b.run();
  assert(fs::exists(restart_dir_b / "step_000003.h5"));
  assert(!fs::exists(restart_dir_b / "step_000002.h5"));

  std::cout << "PASS: closed restart round-trip (agents+grid+time)\n";
}

void test_resume_preserves_mu_max_and_in_crypt() {
  // Reproduces the AWS baseline wipe: stressed mu_realized with intact mu_max
  // must survive checkpoint restore; in_crypt must come from the file, not z.
  const std::string out =
      resolve_test_h5_path("GUTIBM_TEST_RESTART_STRESS_H5", "restart_stress_out");
  const fs::path scratch = fs::path(out).parent_path() / "gutibm_restart_stress";
  const fs::path restart_dir = scratch / "ckpt";
  fs::create_directories(restart_dir);

  SimulationConfig cfg = make_restart_config(out, restart_dir.string());
  cfg.hdf5.enabled = false;
  cfg.restart.enabled = false;
  cfg.time.total_time = 60.0;
  cfg.initial_strains[0].mu_max = 5.5e-4;
  cfg.advection.crypts_enabled = true;
  cfg.advection.crypt_depth = 5e-6;

  Simulation sim;
  sim.init(cfg);
  assert(sim.agents().size() > 0);

  constexpr Real kStressedMu = 5.0e-6;
  constexpr Real kTrueMuMax = 5.5e-4;
  constexpr Real kMicrocinPenalty = 0.03;
  // Outside crypt zone: z-retag would NOT set in_crypt, and washout_rate at
  // zn≈0.6 exceeds kStressedMu — crushed mu_max would wipe the population.
  const Real outside_crypt_z = sim.domain().lo()[2] + 15e-6;
  for (Agent& a : sim.agents()) {
    a.genome.bi_loci.push_back(PlasmidLibrary::microcin_V());
    a.mu_max = kTrueMuMax * (1.0 - kMicrocinPenalty);
    a.mu_realized = kStressedMu;
    a.x[2] = outside_crypt_z;
    a.flags.in_crypt = true;  // written to HDF5; must round-trip
    a.flags.microcin_penalty_applied = true;
  }

  const fs::path ckpt_crypt = restart_dir / "step_000010.h5";
  assert(HDF5Writer::write_closed_restart(
      sim, ckpt_crypt.string(), 10, 600.0, 60.0));

  HDF5CheckpointSnapshot snap = HDF5Reader::load_snapshot(ckpt_crypt.string(), "");
  assert(!snap.agents.mu_max.empty());
  assert(!snap.agents.in_crypt.empty());
  assert(snap.agents.mu_max.size() == snap.agents.id.size());
  assert(snap.genome.bi_release_mode.size() == snap.genome.bi_toxin_id.size());
  for (size_t i = 0; i < snap.agents.id.size(); ++i) {
    assert(std::abs(snap.agents.mu_max[i]
                    - kTrueMuMax * (1.0 - kMicrocinPenalty)) < 1e-12);
    assert(std::abs(snap.agents.mu[i] - kStressedMu) < 1e-12);
    assert(snap.agents.in_crypt[i] != 0);
    assert(snap.agents.microcin_penalty_applied[i] != 0);
  }
  const auto continuous_count = std::count(
      snap.genome.bi_release_mode.begin(), snap.genome.bi_release_mode.end(),
      static_cast<int32_t>(to_underlying(ReleaseMode::CONTINUOUS)));
  assert(continuous_count == static_cast<int>(snap.agents.id.size()));

  const std::string resume_out =
      resolve_test_h5_path("GUTIBM_TEST_RESTART_STRESS_RESUME_H5",
                           "restart_stress_resume");
  SimulationConfig resume_cfg =
      make_restart_config(resume_out, (scratch / "resume").string());
  resume_cfg.hdf5.enabled = false;
  resume_cfg.restart.enabled = false;
  resume_cfg.time.total_time = 720.0;
  resume_cfg.initial_strains[0].mu_max = kTrueMuMax;
  resume_cfg.advection.crypts_enabled = true;
  resume_cfg.advection.crypt_depth = 5e-6;
  resume_cfg.checkpoint.file = ckpt_crypt.string();

  Simulation resumed_crypt;
  resumed_crypt.init_from_checkpoint(resume_cfg, ckpt_crypt.string(), "");
  assert(resumed_crypt.global_agent_count() == snap.metadata.num_agents);
  for (const Agent& a : resumed_crypt.agents()) {
    assert(std::abs(a.mu_max
                    - kTrueMuMax * (1.0 - kMicrocinPenalty)) < 1e-12);
    assert(std::abs(a.mu_realized - kStressedMu) < 1e-12);
    assert(a.flags.in_crypt);
    assert(a.flags.microcin_penalty_applied);
  }
  resumed_crypt.step(60.0);
  for (const Agent& a : resumed_crypt.agents()) {
    assert(std::abs(a.mu_max
                    - kTrueMuMax * (1.0 - kMicrocinPenalty)) < 1e-12);
  }

  // Legacy checkpoint: without the guard dataset, a continuous BI locus
  // implies that the already-penalized mu_max must not be charged again.
  const fs::path legacy_ckpt = restart_dir / "step_000010_legacy.h5";
  fs::copy_file(ckpt_crypt, legacy_ckpt,
                fs::copy_options::overwrite_existing);
  hid_t legacy_file = H5Fopen(legacy_ckpt.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(legacy_file >= 0);
  assert(H5Ldelete(legacy_file, "agents/step_000010/microcin_penalty_applied",
                   H5P_DEFAULT) >= 0);
  assert(H5Fclose(legacy_file) >= 0);

  HDF5CheckpointSnapshot legacy_snap =
      HDF5Reader::load_snapshot(legacy_ckpt.string(), "");
  assert(legacy_snap.agents.microcin_penalty_applied.empty());
  Simulation legacy_resumed;
  legacy_resumed.init_from_checkpoint(resume_cfg, legacy_ckpt.string(), "");
  for (const Agent& a : legacy_resumed.agents()) {
    assert(std::abs(a.mu_max
                    - kTrueMuMax * (1.0 - kMicrocinPenalty)) < 1e-12);
    assert(a.flags.microcin_penalty_applied);
  }
  legacy_resumed.step(60.0);
  for (const Agent& a : legacy_resumed.agents()) {
    assert(std::abs(a.mu_max
                    - kTrueMuMax * (1.0 - kMicrocinPenalty)) < 1e-12);
  }

  // Second artifact: same stress, but NOT crypt-protected — proves mu_max
  // (not in_crypt) is what prevents the one-step washout wipe.
  for (Agent& a : sim.agents()) {
    a.mu_max = kTrueMuMax;
    a.mu_realized = kStressedMu;
    a.x[2] = outside_crypt_z;
    a.flags.in_crypt = false;
  }
  const fs::path ckpt_open = restart_dir / "step_000011.h5";
  assert(HDF5Writer::write_closed_restart(
      sim, ckpt_open.string(), 11, 660.0, 60.0));

  SimulationConfig open_cfg = resume_cfg;
  open_cfg.checkpoint.file = ckpt_open.string();
  Simulation resumed_open;
  resumed_open.init_from_checkpoint(open_cfg, ckpt_open.string(), "");
  const Int n_before = resumed_open.global_agent_count();
  assert(n_before > 1);
  for (const Agent& a : resumed_open.agents()) {
    assert(std::abs(a.mu_max - kTrueMuMax) < 1e-12);
    assert(!a.flags.in_crypt);
  }
  resumed_open.step(60.0);
  assert(resumed_open.global_agent_count() > 1);
  assert(resumed_open.global_agent_count() >= n_before / 2);

  std::cout << "PASS: resume preserves mu_max + in_crypt (no one-step wipe)\n";
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GUTIBM_MPI
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  try {
    test_closed_restart_roundtrip();
    test_resume_preserves_mu_max_and_in_crypt();
  } catch (const std::exception& ex) {
    std::cerr << "FAIL: " << ex.what() << "\n";
#ifdef GUTIBM_MPI
    MPI_Finalize();
#endif
    return 1;
  }
#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif
  return 0;
}
