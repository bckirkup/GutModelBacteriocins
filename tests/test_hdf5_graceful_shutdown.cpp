/* -----------------------------------------------------------------------
   GutIBM – Graceful shutdown / valid HDF5 finalize (Spec 4)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "path_utils.h"
#include "stop_signal.h"

#include <cassert>
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

  Simulation sim;
  sim.init(cfg);
  gutibm_request_stop();
  sim.run();

  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  assert(H5Lexists(file, "summary/step_000000/time", H5P_DEFAULT) > 0);
  H5Fclose(file);

  std::cout << "All graceful shutdown tests passed.\n";
  return 0;
#endif
}
