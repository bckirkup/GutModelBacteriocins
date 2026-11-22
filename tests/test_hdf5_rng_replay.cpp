/* -----------------------------------------------------------------------
   GutIBM – Bitwise Layer-1 HDF5 RNG checkpoint/restart replay (Issue #472)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "h5_util.h"
#include "hdf5_reader.h"
#include "hdf5_writer.h"
#include "path_utils.h"

#include <cassert>
#include <unordered_set>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>
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

bool bitwise_equal(Real a, Real b) {
  return std::memcmp(&a, &b, sizeof(Real)) == 0;
}

bool bitwise_equal_vec3(const Vec3& a, const Vec3& b) {
  return bitwise_equal(a[0], b[0])
      && bitwise_equal(a[1], b[1])
      && bitwise_equal(a[2], b[2]);
}

void assert_step_events_bitwise_equal(const StepEvents& lhs,
                                      const StepEvents& rhs) {
  assert(lhs.sos_inductions == rhs.sos_inductions);
  assert(lhs.phage_inductions == rhs.phage_inductions);
  assert(lhs.mortality_colicin == rhs.mortality_colicin);
  assert(lhs.mortality_cdi == rhs.mortality_cdi);
  assert(lhs.outflow_washout == rhs.outflow_washout);
  assert(lhs.outflow_boundary == rhs.outflow_boundary);
  assert(lhs.mortality_lysis == rhs.mortality_lysis);
  assert(lhs.divisions == rhs.divisions);
  for (size_t i = 0; i < lhs.divisions_by_type.size(); ++i) {
    assert(lhs.divisions_by_type[i] == rhs.divisions_by_type[i]);
  }
  assert(lhs.conjugation_transfers == rhs.conjugation_transfers);
  assert(lhs.mutations == rhs.mutations);
  assert(lhs.immigrations == rhs.immigrations);
  assert(bitwise_equal(lhs.sos_basal_rate, rhs.sos_basal_rate));
  assert(bitwise_equal(lhs.sos_post_division_rate, rhs.sos_post_division_rate));
  assert(bitwise_equal(lhs.sos_nuclease_cross_induction_rate,
                       rhs.sos_nuclease_cross_induction_rate));
  assert(bitwise_equal(lhs.sos_ros_rate, rhs.sos_ros_rate));
}

void assert_agents_bitwise_equal(const AgentPool& lhs, const AgentPool& rhs) {
  assert(lhs.size() == rhs.size());
  for (Int i = 0; i < lhs.size(); ++i) {
    const Agent& a = lhs[i];
    const Agent& b = rhs[i];
    assert(a.identity.tag == b.identity.tag);
    assert(a.identity.type == b.identity.type);
    assert(a.identity.owner_rank == b.identity.owner_rank);
    assert(a.state == b.state);
    assert(a.flags.is_ghost == b.flags.is_ghost);
    assert(a.flags.in_crypt == b.flags.in_crypt);
    assert(a.flags.just_divided == b.flags.just_divided);
    assert(a.flags.microcin_penalty_applied == b.flags.microcin_penalty_applied);
    assert(bitwise_equal_vec3(a.x, b.x));
    assert(bitwise_equal(a.radius, b.radius));
    assert(bitwise_equal(a.outer_radius, b.outer_radius));
    assert(bitwise_equal(a.mass, b.mass));
    assert(bitwise_equal(a.biomass, b.biomass));
    assert(bitwise_equal(a.mu_realized, b.mu_realized));
    assert(bitwise_equal(a.realized_fermentation_fraction,
                         b.realized_fermentation_fraction));
    assert(bitwise_equal(a.respired_oxygen_rate, b.respired_oxygen_rate));
    assert(bitwise_equal(a.timers.age, b.timers.age));
    assert(bitwise_equal(a.timers.sos_timer, b.timers.sos_timer));
    assert(bitwise_equal(a.timers.death_time, b.timers.death_time));
    assert(bitwise_equal(a.pending_growth_chemistry_biomass,
                         b.pending_growth_chemistry_biomass));
    for (Int r = 0; r < NUM_RECEPTORS; ++r) {
      assert(bitwise_equal(a.receptor_expr[r], b.receptor_expr[r]));
      assert(bitwise_equal(a.receptor_expr_base[r], b.receptor_expr_base[r]));
    }
    assert(bitwise_equal_vec3(a.motility.swim_direction,
                              b.motility.swim_direction));
    assert(bitwise_equal(a.motility.swim_speed, b.motility.swim_speed));
    assert(bitwise_equal(a.motility.run_timer, b.motility.run_timer));
    assert(a.motility.is_stopped == b.motility.is_stopped);
    assert(bitwise_equal(a.motility.stop_timer, b.motility.stop_timer));
    assert(bitwise_equal(a.motility.prev_carbon, b.motility.prev_carbon));
    assert(bitwise_equal(a.motility.prev_oxygen, b.motility.prev_oxygen));
    assert(bitwise_equal(a.motility.prev_ai2, b.motility.prev_ai2));
    assert(a.genome.bi_loci.size() == b.genome.bi_loci.size());
    for (size_t j = 0; j < a.genome.bi_loci.size(); ++j) {
      const auto& la = a.genome.bi_loci[j];
      const auto& lb = b.genome.bi_loci[j];
      assert(la.toxin_id == lb.toxin_id);
      assert(la.is_nuclease == lb.is_nuclease);
      assert(bitwise_equal(la.protease_half_life, lb.protease_half_life));
      assert(bitwise_equal(la.burst_size, lb.burst_size));
      assert(bitwise_equal(la.phage_induction_rate, lb.phage_induction_rate));
    }
  }
}

void assert_fields_bitwise_equal(const ChemicalField& lhs,
                                 const ChemicalField& rhs) {
  assert(lhs.num_species() == rhs.num_species());
  assert(lhs.ncells() == rhs.ncells());
  for (Int s = 0; s < lhs.num_species(); ++s) {
    for (Int c = 0; c < lhs.ncells(); ++c) {
      assert(bitwise_equal(lhs.conc(s, c), rhs.conc(s, c)));
    }
  }
}

SimulationConfig make_stochastic_config(const std::string& restart_dir) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.advection.mucus_thickness = 25e-6;
  cfg.advection.distal_length = 50e-6;
  cfg.advection.radial_turnover = 2.0e4;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  cfg.time.bio_dt = 60.0;
  cfg.time.total_time = 20.0 * 60.0;
  cfg.time.output_interval = 60.0;
  cfg.seed = 20260472;

  cfg.hdf5.enabled = false;
  cfg.restart.enabled = !restart_dir.empty();
  cfg.restart.directory = restart_dir;
  cfg.restart.interval_steps = 10;

  // Enable rich stochastic dynamics: SOS induction & lysis bursts, receptor
  // killing, conjugation, mutation, motility, and immigration.
  cfg.fixes.bacteriocin.sos_basal_rate = 2.0e-3;
  cfg.fixes.bacteriocin.sos_lysis_prob = 0.25;
  cfg.fixes.conjugation.base_transfer_rate = 5.0e-3;
  cfg.fixes.receptor.btuB_competitor_conc = 0.0;
  cfg.cell_bio.motility.enabled = true;
  cfg.immigration.enabled = true;
  cfg.immigration.rate = 0.05;
  cfg.immigration.strain_index = 1;

  cfg.initial_strains.clear();
  for (auto& chem : cfg.chemicals) {
    if (chem.name == "carbon") {
      chem.initial_conc = 5.0e-2;
      chem.boundary_conc = 5.0e-2;
      chem.z_gradient_enabled = false;
    }
  }
  SimulationConfig::InitialStrain producer;
  producer.type = 1;
  producer.count = 16;
  producer.mu_max = 1.2e-3;
  producer.plasmids = {"ColE1", "ColE2"};
  producer.conjugative = true;
  cfg.initial_strains.push_back(producer);

  SimulationConfig::InitialStrain sensitive;
  sensitive.type = 2;
  sensitive.count = 16;
  sensitive.mu_max = 1.2e-3;
  sensitive.plasmids = {};
  sensitive.conjugative = false;
  cfg.initial_strains.push_back(sensitive);

  return cfg;
}

void test_bitwise_checkpoint_replay() {
  const std::string base_path =
      resolve_test_h5_path("GUTIBM_TEST_RNG_REPLAY_H5", "rng_replay_base");
  const fs::path scratch =
      fs::path(base_path).parent_path() / "gutibm_rng_replay_artifacts";
  const fs::path restart_dir = scratch / "ckpt";
  fs::remove_all(scratch);
  fs::create_directories(restart_dir);

  constexpr Real dt = 60.0;
  constexpr int total_steps = 20;
  constexpr int split_step = 10;

  // 1. Uninterrupted 20-step reference run.
  SimulationConfig full_cfg = make_stochastic_config("");
  Simulation full_sim;
  full_sim.init(full_cfg);
  std::vector<StepEvents> full_tail_events;
  full_tail_events.reserve(static_cast<size_t>(total_steps - split_step));
  for (int s = 1; s <= total_steps; ++s) {
    full_sim.step(dt);
    if (s > split_step) {
      full_tail_events.push_back(full_sim.step_events());
    }
    full_sim.reset_step_events_after_summary(s, full_sim.time());
  }

  // Verify the stochastic scenario actually exercised stochastic events.
  const StepEvents& full_cum = full_sim.cumulative_events();
  assert(full_cum.divisions > 0);
  assert(full_cum.sos_inductions > 0 || full_cum.mortality_lysis > 0);
  assert(full_cum.immigrations > 0);

  // 2. First 10 steps + closed HDF5 checkpoint at step 10.
  SimulationConfig first_cfg = make_stochastic_config(restart_dir.string());
  first_cfg.time.total_time = static_cast<Real>(split_step) * dt;
  Simulation first_sim;
  first_sim.init(first_cfg);
  for (int s = 1; s <= split_step; ++s) {
    first_sim.step(dt);
    first_sim.reset_step_events_after_summary(s, first_sim.time());
  }
  const fs::path ckpt_path = restart_dir / "step_000010.h5";
  const bool wrote = HDF5Writer::write_closed_restart(
      first_sim, ckpt_path.string(), split_step, first_sim.time(), dt, true);
  assert(wrote);
  assert(fs::exists(ckpt_path));

  // Verify snapshot contains RNG state.
  HDF5CheckpointSnapshot snap = HDF5Reader::load_snapshot(ckpt_path.string(), "");
  assert(snap.metadata.has_rng_state);
  assert(!snap.metadata.simulation_rng_states.empty());
  assert(!snap.metadata.simulation_rng_states[0].empty());
  assert(!snap.metadata.immigration_rng_state.empty());

  // 3. Restart from step 10 checkpoint and run remaining 10 steps.
  SimulationConfig resume_cfg = make_stochastic_config("");
  Simulation resumed_sim;
  resumed_sim.init_from_checkpoint(resume_cfg, ckpt_path.string(), "");
  assert(resumed_sim.step_count() == split_step);
  assert(bitwise_equal(resumed_sim.time(), static_cast<Real>(split_step) * dt));

  std::vector<StepEvents> resumed_tail_events;
  resumed_tail_events.reserve(static_cast<size_t>(total_steps - split_step));
  for (int s = split_step + 1; s <= total_steps; ++s) {
    resumed_sim.step(dt);
    resumed_tail_events.push_back(resumed_sim.step_events());
    resumed_sim.reset_step_events_after_summary(s, resumed_sim.time());
  }

  // 4. Bitwise comparisons on step-by-step StepEvents, cumulative StepEvents,
  //    agent state/positions, and chemical concentrations.
  assert(full_tail_events.size() == resumed_tail_events.size());
  for (size_t i = 0; i < full_tail_events.size(); ++i) {
    assert_step_events_bitwise_equal(full_tail_events[i], resumed_tail_events[i]);
  }
  assert_step_events_bitwise_equal(full_sim.cumulative_events(),
                                   resumed_sim.cumulative_events());
  assert_agents_bitwise_equal(full_sim.agents(), resumed_sim.agents());
  assert_fields_bitwise_equal(full_sim.chemical_field(),
                              resumed_sim.chemical_field());

  std::cout << "  PASS: 20-step uninterrupted vs 10+ckpt+10 replay is bitwise identical"
            << " (agents=" << full_sim.agents().size()
            << ", divisions=" << full_cum.divisions
            << ", sos=" << full_cum.sos_inductions
            << ", lysis=" << full_cum.mortality_lysis
            << ", colicin_kills=" << full_cum.mortality_colicin
            << ", immigrations=" << full_cum.immigrations << ")\n";

  // 5. Legacy checkpoint backward compatibility: strip RNG, toxin_bursts, and
  //    optional extended datasets from a copy of the checkpoint and confirm
  //    init_from_checkpoint loads and steps cleanly.
#ifdef GUTIBM_HDF5
  const fs::path legacy_path = restart_dir / "legacy_step_000010.h5";
  fs::copy_file(ckpt_path, legacy_path, fs::copy_options::overwrite_existing);
  hid_t file = H5Fopen(legacy_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  assert(file >= 0);
  auto remove_if_exists = [file](const char* path) {
    if (H5Lexists(file, path, H5P_DEFAULT) > 0) {
      const herr_t rc = H5Ldelete(file, path, H5P_DEFAULT);
      assert(rc >= 0);
    }
  };
  remove_if_exists("/summary/step_000010/rng");
  remove_if_exists("/summary/step_000010/toxin_bursts");
  remove_if_exists("/agents/step_000010/mass");
  remove_if_exists("/agents/step_000010/outer_radius");
  remove_if_exists("/agents/step_000010/age");
  remove_if_exists("/agents/step_000010/sos_timer");
  remove_if_exists("/agents/step_000010/death_time");
  remove_if_exists("/agents/step_000010/pending_growth_chemistry_biomass");
  remove_if_exists("/agents/step_000010/receptor_expr");
  remove_if_exists("/agents/step_000010/receptor_expr_base");
  remove_if_exists("/agents/step_000010/motility_dir_x");
  remove_if_exists("/agents/step_000010/motility_dir_y");
  remove_if_exists("/agents/step_000010/motility_dir_z");
  remove_if_exists("/agents/step_000010/motility_swim_speed");
  remove_if_exists("/agents/step_000010/motility_run_timer");
  remove_if_exists("/agents/step_000010/motility_is_stopped");
  remove_if_exists("/agents/step_000010/motility_stop_timer");
  remove_if_exists("/agents/step_000010/motility_prev_carbon");
  remove_if_exists("/agents/step_000010/motility_prev_oxygen");
  remove_if_exists("/agents/step_000010/motility_prev_ai2");
  remove_if_exists("/agents/step_000010/sampled_nuclease_conc");
  remove_if_exists("/genome/step_000010/bi_protease_half_life");
  remove_if_exists("/genome/step_000010/bi_is_nuclease");
  remove_if_exists("/genome/step_000010/bi_burst_size");
  H5Fclose(file);

  HDF5CheckpointSnapshot legacy_snap =
      HDF5Reader::load_snapshot(legacy_path.string(), "");
  assert(!legacy_snap.metadata.has_rng_state);
  assert(!legacy_snap.metadata.has_toxin_bursts);
  assert(legacy_snap.agents.mass.empty());

  Simulation legacy_sim;
  legacy_sim.init_from_checkpoint(resume_cfg, legacy_path.string(), "");
  assert(legacy_sim.step_count() == split_step);
  assert(legacy_sim.agents().size() == first_sim.agents().size());
  for (int s = split_step + 1; s <= split_step + 5; ++s) {
    legacy_sim.step(dt);
  }
  assert(legacy_sim.step_count() == split_step + 5);
  std::cout << "  PASS: legacy checkpoint without RNG/extended datasets loads cleanly\n";
#endif

  fs::remove_all(scratch);
}


void test_checkpoint_fork_reseed_and_tag_repartition() {
#ifdef GUTIBM_HDF5
  const std::string base_path =
      resolve_test_h5_path("GUTIBM_TEST_RNG_FORK_H5", "rng_fork_base");
  const fs::path scratch =
      fs::path(base_path).parent_path() / "gutibm_rng_fork_artifacts";
  const fs::path restart_dir = scratch / "ckpt";
  fs::remove_all(scratch);
  fs::create_directories(restart_dir);

  constexpr Real dt = 60.0;
  constexpr int burn_steps = 6;
  constexpr int fork_steps = 6;

  // Burn-in run WITHOUT immigration enabled (seed = 1001).
  SimulationConfig burn_cfg = make_stochastic_config(restart_dir.string());
  burn_cfg.seed = 1001;
  burn_cfg.immigration.enabled = false;
  burn_cfg.time.total_time = static_cast<Real>(burn_steps) * dt;

  Simulation burn_sim;
  burn_sim.init(burn_cfg);
  for (int s = 1; s <= burn_steps; ++s) {
    burn_sim.step(dt);
    burn_sim.reset_step_events_after_summary(s, burn_sim.time());
  }
  const fs::path ckpt_path = restart_dir / "step_000006.h5";
  assert(HDF5Writer::write_closed_restart(
      burn_sim, ckpt_path.string(), burn_steps, burn_sim.time(), dt, true));

  const HDF5CheckpointSnapshot snap =
      HDF5Reader::load_snapshot(ckpt_path.string(), "");
  assert(snap.metadata.has_rng_seed);
  assert(snap.metadata.rng_seed == 1001ULL);
  assert(snap.metadata.has_immigration_enabled);
  assert(!snap.metadata.immigration_was_enabled);

  // Fork A (seed = 2002) and Fork B (seed = 3003) from the same burn-in
  // checkpoint, enabling a single-pulse immigration schedule (schedule="pulse",
  // step=0, count=3), which fires when relative_step == 0 at fork start.
  auto run_fork = [&](uint64_t fork_seed, Simulation& out_sim) {
    SimulationConfig fcfg = burn_cfg;
    fcfg.restart.enabled = false;
    fcfg.restart.directory.clear();
    fcfg.seed = fork_seed;
    fcfg.immigration.enabled = true;
    fcfg.immigration.schedule = "pulse";
    fcfg.immigration.step = 0;
    fcfg.immigration.count = 3;
    fcfg.immigration.rate = 0.0;
    fcfg.immigration.strain_index = 1;
    out_sim.init_from_checkpoint(fcfg, ckpt_path.string(), "");
    assert(out_sim.immigration().start_step() == burn_steps);
    for (int s = burn_steps + 1; s <= burn_steps + fork_steps; ++s) {
      out_sim.step(dt);
      out_sim.reset_step_events_after_summary(s, out_sim.time());
    }
  };

  Simulation fork_a;
  Simulation fork_b;
  run_fork(2002, fork_a);
  run_fork(3003, fork_b);

  // Single-pulse immigration at fork start must have injected 3 immigrants on
  // the first post-checkpoint step.
  assert(fork_a.cumulative_events().immigrations == 3);
  assert(fork_b.cumulative_events().immigrations == 3);

  // Distinct fork seeds must produce distinct stochastic trajectories instead
  // of replaying the burn-in RNG state identically.
  bool diverged = fork_a.agents().size() != fork_b.agents().size();
  for (Int i = 0; !diverged && i < fork_a.agents().size(); ++i) {
    if (!bitwise_equal_vec3(fork_a.agents()[i].x, fork_b.agents()[i].x) ||
        !bitwise_equal(fork_a.agents()[i].biomass,
                       fork_b.agents()[i].biomass)) {
      diverged = true;
    }
  }
  assert(diverged);

  // Simulate an MPI rank-repartitioned checkpoint (e.g., 2-rank checkpoint
  // loaded on a 1-rank process) where next_tag = {1, 2} would collide with live
  // founder tags if next_tags.size() != domain.nprocs() were not guarded, and
  // verify bounded h5::read_str_vec on a non-null-terminated string dataset.
  const fs::path repart_path = restart_dir / "repart_step_000006.h5";
  fs::copy_file(ckpt_path, repart_path, fs::copy_options::overwrite_existing);
  {
    hid_t fid = H5Fopen(repart_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    assert(fid >= 0);
    if (H5Lexists(fid, "/summary/step_000006/rng/next_tag", H5P_DEFAULT) > 0) {
      assert(H5Ldelete(fid, "/summary/step_000006/rng/next_tag",
                       H5P_DEFAULT) >= 0);
    }
    const std::vector<int64_t> two_rank_tags = {1, 2};
    h5::write_vec(fid, "/summary/step_000006/rng/next_tag", H5T_NATIVE_INT64,
                  two_rank_tags);

    hsize_t dims[1] = {2};
    hid_t sp = H5Screate_simple(1, dims, nullptr);
    hid_t tp = H5Tcopy(H5T_C_S1);
    H5Tset_size(tp, h5::kRngStrLen);
    hid_t ds = H5Dcreate2(fid, "/summary/step_000006/rng/non_nul_probe", tp, sp,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(ds >= 0);
    std::vector<char> raw(2 * h5::kRngStrLen, 'X');
    assert(H5Dwrite(ds, tp, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw.data()) >= 0);
    H5Dclose(ds);
    H5Tclose(tp);
    H5Sclose(sp);

    const auto bounded =
        h5::read_str_vec(fid, "/summary/step_000006/rng/non_nul_probe");
    assert(bounded.size() == 2);
    assert(bounded[0].size() == h5::kRngStrLen);
    assert(bounded[1].size() == h5::kRngStrLen);
    H5Fclose(fid);
  }

  Simulation repart_sim;
  SimulationConfig repart_cfg = burn_cfg;
  repart_cfg.restart.enabled = false;
  repart_cfg.immigration.enabled = true;
  repart_cfg.immigration.schedule = "pulse";
  repart_cfg.immigration.step = 0;
  repart_cfg.immigration.count = 4;
  repart_cfg.immigration.rate = 0.0;
  repart_sim.init_from_checkpoint(repart_cfg, repart_path.string(), "");
  for (int s = burn_steps + 1; s <= burn_steps + 4; ++s) {
    repart_sim.step(dt);
  }
  std::unordered_set<TagID> seen_tags;
  for (Int i = 0; i < repart_sim.agents().size(); ++i) {
    const TagID tag = repart_sim.agents()[i].identity.tag;
    assert(tag > 0);
    assert(seen_tags.insert(tag).second);
  }
  std::cout << "  PASS: seed-fork RNG divergence, fork pulse immigration, "
               "rank-repartition tag uniqueness, and bounded read_str_vec\n";
  fs::remove_all(scratch);
#endif
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GUTIBM_MPI
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  std::cout << "=== Layer-1 HDF5 RNG Bitwise Replay Tests ===\n";
  test_bitwise_checkpoint_replay();
  test_checkpoint_fork_reseed_and_tag_repartition();
  std::cout << "All HDF5 RNG bitwise replay tests passed.\n";
#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif
  return 0;
}
