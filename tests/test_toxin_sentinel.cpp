/* -----------------------------------------------------------------------
   GutIBM – Toxin-sensitive producer-colony challenge
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "path_utils.h"
#include "simulation.h"
#include "species_names.h"
#include "greens_function.h"
#include <cassert>
#include <cmath>
#include <format>
#include <iomanip>
#include <iostream>
#include <numbers>
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
  Int target_count = 100;
  Real target_toxin = 0.0;
  Real greens_function_toxin = 0.0;
  Real source_to_cell_distance = 0.0;
  Real expected_kill_fraction = 0.0;
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

Real read_grid_value(hid_t file, const std::string& path,
                     const Domain& domain, Int cell) {
  hid_t dataset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dataset >= 0);
  std::vector<Real> values(static_cast<size_t>(domain.ncells()));
  assert(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, values.data()) >= 0);
  H5Dclose(dataset);
  return values[static_cast<size_t>(cell)];
}
#endif

SimulationConfig challenge_config(Int producer_count) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {100e-6, 50e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.periodic = {false, false, false};
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
  cfg.enabled_fixes = {"bacteriocin", "receptor"};
  cfg.qssa.toxin_cutoff = 60e-6;
  cfg.fixes.bacteriocin.sos_basal_rate = 1.0;

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
  targets.count = 100;
  targets.mu_max = 5.0e-4;
  targets.conjugative = false;
  cfg.initial_strains.push_back(targets);
  return cfg;
}

ChallengeResult run_challenge(Int producer_count, Real target_distance) {
#ifndef GUTIBM_HDF5
  return {};
#else
  SimulationConfig cfg = challenge_config(producer_count);
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = resolve_test_h5_path("GUTIBM_TOXIN_SENTINEL_H5",
                                            "toxin_sentinel");
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 0;
  cfg.hdf5.schedule.grid = 1;
  cfg.hdf5.schedule.grid_species = {species::BACTERIOCIN_BTUB};
  cfg.hdf5.schedule.lineage = 0;
  cfg.hdf5.schedule.genome = 0;

  Simulation sim;
  sim.init(cfg);

  const Vec3 producer_center = {25e-6, 25e-6, 25e-6};
  const Vec3 target_center = {25e-6 + target_distance, 25e-6, 25e-6};
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
  const Int target_ix = static_cast<Int>(
      std::lround(target_center[0] / cfg.domain.grid_dx));
  const Int target_iy = static_cast<Int>(
      std::lround(target_center[1] / cfg.domain.grid_dx));
  const Int target_iz = static_cast<Int>(
      std::lround(target_center[2] / cfg.domain.grid_dx));
  const Int target_cell = sim.domain().cell_index(target_ix, target_iy, target_iz);
  const std::string grid_path =
      std::string("grid/step_000006/") + species::BACTERIOCIN_BTUB;
  const Real target_toxin =
      read_grid_value(file, grid_path, sim.domain(), target_cell);
  GreensFunction greens_function;
  greens_function.init(sim.domain(), sim.advection());
  GreensFunctionParams params;
  params.diff_coeff = 4.0e-11;
  params.retardation = 50.0;
  params.source_rate = 1.0e5 / AVOGADRO / 300.0
      * std::exp(-60.0 / 300.0);
  params.decay_rate = std::numbers::ln2 / 1800.0;
  const Real greens_function_toxin = greens_function.concentration_bounded(
      producer_center, sim.domain().cell_center(target_ix, target_iy, target_iz),
      params);
  const Vec3 target_cell_center =
      sim.domain().cell_center(target_ix, target_iy, target_iz);
  const Real dx = target_cell_center[0] - producer_center[0];
  const Real dy = target_cell_center[1] - producer_center[1];
  const Real dz = target_cell_center[2] - producer_center[2];
  constexpr Real apparent_kd = 5.005e-4;
  constexpr Real kill_rate = 1.0e-3;
  constexpr Real bio_dt = 60.0;
  constexpr Real release_tau = 300.0;
  constexpr Real toxin_half_life = 1800.0;
  constexpr Real decay_rate = std::numbers::ln2 / toxin_half_life;
  Real cumulative_hazard = 0.0;
  for (Int step = 0; step < 5; ++step) {
    const Real age = static_cast<Real>(step) * bio_dt;
    const Real toxin = target_toxin
        * std::exp(-age / release_tau - decay_rate * age);
    const Real occupancy = toxin / (apparent_kd + toxin);
    cumulative_hazard += kill_rate * occupancy * bio_dt;
  }
  H5Fclose(file);

  ChallengeResult result;
  result.colicin_kills = colicin_kills;
  result.divisions = divisions;
  result.target_toxin = target_toxin;
  result.greens_function_toxin = greens_function_toxin;
  result.source_to_cell_distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  result.expected_kill_fraction = 1.0 - std::exp(-cumulative_hazard);
  return result;
#endif
}

}  // namespace

int main() {
#ifndef GUTIBM_HDF5
  std::cout << "HDF5 disabled — skipping toxin sentinel challenge.\n";
  return 0;
#else
  const std::vector<Int> producer_counts = {1, 10, 100, 1000, 10000};
  const std::vector<Real> distances = {10e-6, 50e-6};
  for (size_t distance_index = 0; distance_index < distances.size();
       ++distance_index) {
    const Real distance = distances[distance_index];
    std::vector<ChallengeResult> results;
    results.reserve(producer_counts.size());
    for (const Int count : producer_counts) {
      results.push_back(run_challenge(count, distance));
    }

    std::cout << "distance_um=" << distance * 1e6 << "\n";
    std::cout << std::setprecision(10);
    for (size_t i = 0; i < producer_counts.size(); ++i) {
      const Real simulated_fraction =
          static_cast<Real>(results[i].colicin_kills) / results[i].target_count;
      std::cout << "producers=" << producer_counts[i]
                << " killed_fraction=" << simulated_fraction
                << " expected_fraction=" << results[i].expected_kill_fraction
                << " target_toxin=" << results[i].target_toxin
                << " gf_toxin=" << results[i].greens_function_toxin
                << " source_cell_distance_um="
                << results[i].source_to_cell_distance * 1e6
                << " divisions=" << results[i].divisions << "\n";
    }

    if (distance_index == 1U) {
      const Real simulated_fraction =
          static_cast<Real>(results[2].colicin_kills)
          / results[2].target_count;
      const Real sampled_per_source =
          results[2].target_toxin / producer_counts[2];
      assert(results[2].colicin_kills <= 5);
      assert(results[3].colicin_kills >= 10);
      assert(std::abs(sampled_per_source / results[2].greens_function_toxin
                      - 1.0) <= 0.1);
      assert(std::abs(simulated_fraction
                      - results[2].expected_kill_fraction) <= 0.05);
    }
    for (const ChallengeResult& result : results) {
      assert(result.divisions == 0);
    }
  }
  std::cout << "Toxin-sensitive colony challenge passed.\n";
  return 0;
#endif
}
