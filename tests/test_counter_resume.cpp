#include "hdf5_reader.h"
#include "hdf5_writer.h"
#include "input_parser.h"
#include "path_utils.h"
#include "simulation.h"

#include <cassert>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

using namespace gutibm;

#ifdef GUTIBM_HDF5
namespace {

double scalar(hid_t file, const std::string& path) {
  hid_t dataset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dataset >= 0);
  std::array<double, 16> values{};
  assert(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, values.data()) >= 0);
  H5Dclose(dataset);
  return values[0];
}

SimulationConfig config() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {25e-6, 25e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.hdf5.enabled = false;
  cfg.restart.enabled = false;
  cfg.initial_strains.clear();
  cfg.initial_strains.push_back({1, 1, 5e-4, {}, false, 0, 0, {}});
  return cfg;
}

void populate_window(Simulation& sim, Int kills, Real boundary) {
  sim.step_events().mortality_colicin = kills;
  sim.chemical_field().flux_accounting().boundary_interval[0] = boundary;
  sim.chemical_field().flux_accounting().vbf_source_interval[0] = boundary + 1.0;
  sim.chemical_field().flux_accounting().vbf_sink_interval[0] = boundary + 2.0;
  sim.chemical_field().flux_accounting().agent_uptake_interval[0] =
      boundary + 3.0;
  sim.chemical_field().flux_accounting().maintenance_interval[0] =
      boundary + 4.0;
  sim.chemical_field().flux_accounting().maintenance_shortfall_interval[0] =
      boundary + 5.0;
  sim.chemical_field().flux_accounting().maintenance_limited_agents_interval[0] =
      boundary + 6.0;
}

void assert_window(const std::string& path, Int step, Int interval_kills,
                   Int cumulative_kills, double boundary_interval,
                   double boundary_cumulative, double start_time) {
  hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  const std::string prefix = "summary/step_" +
      (step == 1 ? std::string("000001") : std::string("000002"));
  assert(scalar(file, prefix + "/events/mortality_colicin") == interval_kills);
  assert(scalar(file, prefix + "/events/cumulative_mortality_colicin")
         == cumulative_kills);
  assert(scalar(file, prefix + "/events/interval_start_step")
         == static_cast<double>(step));
  assert(scalar(file, prefix + "/events/interval_end_step")
         == static_cast<double>(step));
  assert(std::abs(scalar(file, prefix + "/events/interval_start_time")
                  - start_time) < 1e-12);
  assert(std::abs(scalar(file, prefix + "/nutrient_flux/boundary_interval")
                  - boundary_interval) < 1e-12);
  assert(std::abs(scalar(file, prefix + "/nutrient_flux/boundary_cumulative")
                  - boundary_cumulative) < 1e-12);
  assert(std::abs(scalar(file, prefix + "/nutrient_flux/vbf_source_interval")
                  - boundary_interval - 1.0) < 1e-12);
  assert(std::abs(scalar(file, prefix + "/nutrient_flux/vbf_source_cumulative")
                  - boundary_cumulative - (step == 1 ? 1.0 : 2.0)) < 1e-12);
  assert(std::abs(scalar(file, prefix + "/nutrient_flux/vbf_sink_interval")
                  - boundary_interval - 2.0) < 1e-12);
  assert(std::abs(scalar(file, prefix + "/nutrient_flux/vbf_sink_cumulative")
                  - boundary_cumulative - (step == 1 ? 2.0 : 4.0)) < 1e-12);
  assert(std::abs(scalar(file, prefix + "/nutrient_flux/agent_uptake_interval")
                  - boundary_interval - 3.0) < 1e-12);
  assert(std::abs(scalar(file, prefix + "/nutrient_flux/agent_uptake_cumulative")
                  - boundary_cumulative - (step == 1 ? 3.0 : 6.0)) < 1e-12);
  assert(std::abs(scalar(file, prefix + "/nutrient_flux/maintenance_interval")
                  - boundary_interval - 4.0) < 1e-12);
  assert(std::abs(scalar(file, prefix + "/nutrient_flux/maintenance_cumulative")
                  - boundary_cumulative - (step == 1 ? 4.0 : 8.0)) < 1e-12);
  assert(std::abs(scalar(file,
                         prefix + "/nutrient_flux/maintenance_shortfall_interval")
                  - boundary_interval - 5.0) < 1e-12);
  assert(std::abs(scalar(file,
                         prefix + "/nutrient_flux/maintenance_shortfall_cumulative")
                  - boundary_cumulative - (step == 1 ? 5.0 : 10.0)) < 1e-12);
  assert(std::abs(scalar(
                     file,
                     prefix + "/nutrient_flux/maintenance_limited_agents_interval")
                  - boundary_interval - 6.0) < 1e-12);
  assert(std::abs(scalar(
                     file,
                     prefix + "/nutrient_flux/maintenance_limited_agents_cumulative")
                  - boundary_cumulative - (step == 1 ? 6.0 : 12.0)) < 1e-12);
  assert(scalar(file, prefix + "/nutrient_flux/interval_start_step")
         == static_cast<double>(step));
  assert(scalar(file, prefix + "/nutrient_flux/interval_end_step")
         == static_cast<double>(step));
  H5Fclose(file);
}

}  // namespace
#endif

int main() {
#ifndef GUTIBM_HDF5
  std::cout << "HDF5 disabled — skipping counter resume test.\n";
  return 0;
#else
  const std::filesystem::path root =
      std::filesystem::path(resolve_test_h5_path(
          "GUTIBM_COUNTER_RESUME_H5", "counter_resume")).parent_path();
  std::filesystem::create_directories(root);
  const std::filesystem::path first = root / "window_1.h5";
  const std::filesystem::path second = root / "window_2.h5";

  Simulation first_sim;
  first_sim.init(config());
  first_sim.set_event_window_start(1, 0.0);
  populate_window(first_sim, 2, 3.0);
  assert(HDF5Writer::write_closed_restart(
      first_sim, first.string(), 1, 60.0, 60.0));
  assert_window(first.string(), 1, 2, 2, 3.0, 3.0, 0.0);

  SimulationConfig resume_cfg = config();
  Simulation resumed;
  resumed.init_from_checkpoint(resume_cfg, first.string(), "");
  assert(resumed.cumulative_events().mortality_colicin == 2);
  populate_window(resumed, 5, 7.0);
  assert(HDF5Writer::write_closed_restart(
      resumed, second.string(), 2, 120.0, 60.0));
  assert_window(second.string(), 2, 5, 7, 7.0, 10.0, 60.0);
  std::cout << "PASS: event and nutrient counters survive a two-window resume\n";
  return 0;
#endif
}
