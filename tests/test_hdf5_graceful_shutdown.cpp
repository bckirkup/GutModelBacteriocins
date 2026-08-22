/* -----------------------------------------------------------------------
   GutIBM – Graceful shutdown / valid HDF5 finalize (Spec 4)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "path_utils.h"
#include "stop_signal.h"

#include <cassert>
#include <filesystem>
#include <iostream>
#include <string>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

using namespace gutibm;

int main() {
#ifndef GUTIBM_HDF5
  std::cout << "HDF5 disabled — skipping graceful shutdown tests.\n";
  return 0;
#else
  install_stop_signal_handlers();
  gutibm_reset_stop_request();

  const std::string filename = resolve_test_h5_path("GUTIBM_SHUTDOWN_H5", "graceful_shutdown");

  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {30e-6, 30e-6, 15e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.total_time = 1.0e9;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 60.0;
  cfg.seed = 555;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = filename;
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 1;
  cfg.initial_strains[0].count = 4;
  cfg.initial_strains.resize(1);
  cfg.restart.enabled = true;
  cfg.restart.directory = filename + ".restart";
  cfg.restart.interval_steps = 360;
  std::filesystem::create_directories(cfg.restart.directory);

  Simulation sim;
  sim.init(cfg);
  gutibm_request_stop();
  sim.run();
  assert(sim.termination_cause() == TerminationCause::StopRequested);

  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  assert(H5Lexists(file, "summary/step_000000/time", H5P_DEFAULT) > 0);
  H5Fclose(file);

  const std::string incomplete_filename =
      resolve_test_h5_path("GUTIBM_INCOMPLETE_H5", "incomplete_marker");
  cfg.hdf5.filename = incomplete_filename;
  {
    Simulation incomplete;
    incomplete.init(cfg);
    incomplete.step(60.0);
    assert(incomplete.termination_cause() ==
           TerminationCause::IncompleteUnknown);
  }
  hid_t incomplete_file =
      H5Fopen(incomplete_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(incomplete_file >= 0);
  hid_t incomplete_cause = H5Dopen2(
      incomplete_file, "run_provenance/termination_cause_code", H5P_DEFAULT);
  int32_t incomplete_code = -1;
  assert(incomplete_cause >= 0);
  assert(H5Dread(incomplete_cause, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, &incomplete_code) >= 0);
  H5Dclose(incomplete_cause);
  H5Fclose(incomplete_file);
  assert(incomplete_code == static_cast<int32_t>(
      to_underlying(TerminationCause::IncompleteUnknown)));

  const std::string restart_output = filename + ".restart_output";
  cfg.hdf5.filename = restart_output;
  gutibm_reset_stop_request();
  Simulation restart_sim;
  restart_sim.init(cfg);
  restart_sim.step(60.0);
  gutibm_request_stop();
  restart_sim.run();

  const std::string restart_path =
      cfg.restart.directory + "/step_000001.h5";
  assert(std::filesystem::is_regular_file(restart_path));
  hid_t restart = H5Fopen(restart_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(restart >= 0);
  assert(H5Lexists(restart, "agents/step_000001", H5P_DEFAULT) > 0);
  H5Fclose(restart);

  std::cout << "All graceful shutdown tests passed.\n";
  return 0;
#endif
}
