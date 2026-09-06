#include "hdf5_writer.h"
#include "hdf5_test_helpers.h"
#include "input_parser.h"
#include "path_utils.h"
#include "simulation.h"
#include "step_events.h"

#include <cassert>
#include <array>
#include <filesystem>
#include <iostream>
#include <string>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

using namespace gutibm;

#ifdef GUTIBM_HDF5
namespace {

int32_t read_i32(hid_t file, const std::string& path) {
  hid_t dataset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dataset >= 0);
  int32_t value = 0;
  assert(H5Dread(dataset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 &value) >= 0);
  H5Dclose(dataset);
  return value;
}

std::array<int32_t, MAX_AGENT_TYPES> read_i32_array(hid_t file,
                                                    const std::string& path) {
  hid_t dataset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dataset >= 0);
  std::array<int32_t, MAX_AGENT_TYPES> values{};
  assert(H5Dread(dataset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 values.data()) >= 0);
  H5Dclose(dataset);
  return values;
}

SimulationConfig make_cfg(const std::string& h5_path) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {30e-6, 30e-6, 30e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.bio_dt = 60.0;
  cfg.time.total_time = 60.0;
  cfg.time.output_interval = 60.0;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = h5_path;
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 0;
  cfg.hdf5.schedule.grid = 0;
  cfg.hdf5.schedule.lineage = 0;
  cfg.hdf5.schedule.genome = 0;
  cfg.hdf5.schedule.provenance = 0;
  cfg.restart.enabled = false;
  cfg.gpu.enabled = false;
  cfg.dysbiosis_threshold = 0.0;
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s1;
  s1.type = 1;
  s1.count = 2;
  s1.mu_max = 5e-4;
  SimulationConfig::InitialStrain s2;
  s2.type = 2;
  s2.count = 3;
  s2.mu_max = 5e-4;
  cfg.initial_strains.push_back(s1);
  cfg.initial_strains.push_back(s2);
  return cfg;
}

void assert_by_type_matches(const StepEvents& events,
                            const std::array<int32_t, MAX_AGENT_TYPES>& stored) {
  for (Int i = 0; i < MAX_AGENT_TYPES; ++i) {
    assert(stored[static_cast<size_t>(i)]
           == events.divisions_by_type[static_cast<size_t>(i)]);
  }
}

}  // namespace
#endif

int main() {
#ifndef GUTIBM_HDF5
  std::cout << "test_divisions_by_type: SKIPPED (HDF5 disabled)\n";
  return 0;
#else
  const std::filesystem::path path =
      resolve_test_h5_path("GUTIBM_DIVISIONS_BY_TYPE_H5", "divisions_by_type");
  std::filesystem::remove(path);

  // --- Metabolism path: supersize biomass so every live agent divides once ---
  {
    Simulation sim;
    sim.init(make_cfg(path.string()));
    const Real initial_mass =
        sphere_mass(CELL_RADIUS_DEFAULT, CELL_DENSITY_DEFAULT);
    Int type1 = 0;
    Int type2 = 0;
    for (Agent& agent : sim.agents()) {
      if (agent.state == PhenoState::DEAD || agent.flags.is_ghost) continue;
      agent.biomass = initial_mass * 3.0;
      if (agent.identity.type == 1) ++type1;
      if (agent.identity.type == 2) ++type2;
    }
    assert(type1 == 2);
    assert(type2 == 3);

    // write_hdf5_step (and cumulative commit) runs inside run(), not step().
    assert(sim.run() == 0);

    const StepEvents& events = sim.cumulative_events();
    assert(events.divisions_by_type[1] == type1);
    assert(events.divisions_by_type[2] == type2);
    Int by_type_sum = 0;
    for (Int v : events.divisions_by_type) by_type_sum += v;
    assert(by_type_sum == events.divisions);
    assert(events.divisions == type1 + type2);

    hid_t file = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    assert(file >= 0);
    const std::string prefix = "summary/step_000001/events/";
    assert(read_i32(file, prefix + "divisions") == events.divisions);
    assert(read_i32(file, prefix + "cumulative_divisions") == events.divisions);
    assert_by_type_matches(events,
                           read_i32_array(file, prefix + "divisions_by_type"));
    assert_by_type_matches(
        events, read_i32_array(file, prefix + "cumulative_divisions_by_type"));
    H5Fclose(file);
  }

  // --- Restart continuity of cumulative_divisions_by_type ---
  {
    const auto first = path.parent_path() / "divisions_by_type_window1.h5";
    const auto second = path.parent_path() / "divisions_by_type_window2.h5";
    std::filesystem::remove(first);
    std::filesystem::remove(second);

    SimulationConfig cfg = make_cfg(first.string());
    cfg.hdf5.enabled = false;
    Simulation first_sim;
    first_sim.init(cfg);
    first_sim.set_event_window_start(1, 0.0);
    first_sim.step_events().record_division(1);
    first_sim.step_events().record_division(1);
    first_sim.step_events().record_division(2);
    assert(HDF5Writer::write_closed_restart(
        first_sim, first.string(), 1, 60.0, 60.0));

    Simulation resumed;
    resumed.init_from_checkpoint(cfg, first.string(), "");
    assert(resumed.cumulative_events().divisions == 3);
    assert(resumed.cumulative_events().divisions_by_type[1] == 2);
    assert(resumed.cumulative_events().divisions_by_type[2] == 1);
    resumed.step_events().record_division(1);
    resumed.step_events().record_division(3);
    assert(HDF5Writer::write_closed_restart(
        resumed, second.string(), 2, 120.0, 60.0));

    hid_t file = H5Fopen(second.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    assert(file >= 0);
    const std::string prefix = "summary/step_000002/events/";
    assert(read_i32(file, prefix + "divisions") == 2);
    assert(read_i32(file, prefix + "cumulative_divisions") == 5);
    const auto interval = read_i32_array(file, prefix + "divisions_by_type");
    const auto cumulative =
        read_i32_array(file, prefix + "cumulative_divisions_by_type");
    assert(interval[1] == 1);
    assert(interval[3] == 1);
    assert(cumulative[1] == 3);
    assert(cumulative[2] == 1);
    assert(cumulative[3] == 1);
    H5Fclose(file);
  }

  std::cout << "test_divisions_by_type: PASSED\n";
  return 0;
#endif
}
