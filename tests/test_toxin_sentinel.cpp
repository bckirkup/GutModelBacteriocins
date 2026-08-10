/* -----------------------------------------------------------------------
   GutIBM – Toxin-sensitive producer-colony challenge
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "path_utils.h"
#include "simulation.h"
#include <cassert>
#include <format>
#include <iostream>
#include <string>
#include <vector>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

using namespace gutibm;

namespace {

struct ChallengeResult {
  Int colicin_kills = 0;
  Int divisions = 0;
};

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
#endif

SimulationConfig challenge_config(Int producer_count) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {400e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.bio_dt = 60.0;
  cfg.time.total_time = 600.0;
  cfg.time.output_interval = 600.0;
  cfg.seed = 90210 + static_cast<uint64_t>(producer_count);
  cfg.hdf5.enabled = false;
  cfg.advection.radial_turnover = 1.0e12;
  cfg.advection.distal_transit_time = 1.0e12;
  cfg.cell_bio.motility.enabled = false;
  cfg.cell_bio.cdi.enabled = false;
  cfg.advection.crypts_enabled = false;
  cfg.fixes.mechanics.hertzian_enabled = false;
  cfg.fixes.bacteriocin.sos_basal_rate = 1.0;
  cfg.fixes.receptor.kill_rate_colicin = 1.0;

  cfg.initial_strains.clear();

  SimulationConfig::InitialStrain producers;
  producers.type = 1;
  producers.count = producer_count;
  producers.mu_max = 5.0e-4;
  producers.plasmids = {"ColE1"};
  producers.conjugative = false;
  cfg.initial_strains.push_back(producers);

  SimulationConfig::InitialStrain targets;
  targets.type = 2;
  targets.count = 20;
  targets.mu_max = 5.0e-4;
  targets.conjugative = false;
  cfg.initial_strains.push_back(targets);
  return cfg;
}

ChallengeResult run_challenge(Int producer_count) {
#ifndef GUTIBM_HDF5
  (void)producer_count;
  return {};
#else
  SimulationConfig cfg = challenge_config(producer_count);
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = resolve_test_h5_path("GUTIBM_TOXIN_SENTINEL_H5",
                                            "toxin_sentinel");
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 0;
  cfg.hdf5.schedule.grid = 0;
  cfg.hdf5.schedule.lineage = 0;
  cfg.hdf5.schedule.genome = 0;

  Simulation sim;
  sim.init(cfg);

  const Vec3 producer_center = {75e-6, 50e-6, 25e-6};
  const Vec3 target_center = {275e-6, 50e-6, 25e-6};
  for (Agent& agent : sim.agents()) {
    agent.x = agent.identity.type == 1 ? producer_center : target_center;
    agent.flags.in_crypt = false;
    if (agent.identity.type == 1) {
      agent.state = PhenoState::SOS_INDUCED;
      agent.timers.sos_timer = 300.0;
    }
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid(agent.x, ix, iy, iz);
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  }

  sim.run();

  hid_t file = H5Fopen(cfg.hdf5.filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  Int colicin_kills = 0;
  Int divisions = 0;
  for (Int step = 0; step <= 10; ++step) {
    const std::string prefix = std::format("summary/step_{:06}/events/", step);
    htri_t exists = 0;
    H5E_BEGIN_TRY {
      exists = H5Lexists(file, (prefix + "colicin_kills").c_str(), H5P_DEFAULT);
    } H5E_END_TRY;
    if (exists <= 0) {
      break;
    }
    colicin_kills += read_event(file, prefix + "colicin_kills");
    divisions += read_event(file, prefix + "divisions");
  }
  H5Fclose(file);

  ChallengeResult result;
  result.colicin_kills = colicin_kills;
  result.divisions = divisions;
  return result;
#endif
}

}  // namespace

int main() {
#ifndef GUTIBM_HDF5
  std::cout << "HDF5 disabled — skipping toxin sentinel challenge.\n";
  return 0;
#else
  const std::vector<Int> producer_counts = {1, 10, 100, 1000};
  std::vector<ChallengeResult> results;
  results.reserve(producer_counts.size());
  for (const Int count : producer_counts) {
    results.push_back(run_challenge(count));
  }

  std::cout << "=== Toxin Sentinel Colony Challenge ===\n";
  for (size_t i = 0; i < producer_counts.size(); ++i) {
    std::cout << "producers=" << producer_counts[i]
              << " colicin_kills=" << results[i].colicin_kills
              << " divisions=" << results[i].divisions << "\n";
  }

  assert(results[0].colicin_kills <= 1);
  assert(results[1].colicin_kills >= 5);
  assert(results[2].colicin_kills >= 10);
  assert(results[3].colicin_kills >= results[2].colicin_kills);
  assert(results[3].colicin_kills >= 10);
  for (const ChallengeResult& result : results) {
    assert(result.divisions == 0);
  }

  std::cout << "Toxin-sensitive colony challenge passed.\n";
  return 0;
#endif
}
