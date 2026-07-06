/* -----------------------------------------------------------------------
   GutIBM – HDF5 schedule layer tests (Spec 4)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "path_utils.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

using namespace gutibm;

namespace {

bool dataset_exists(hid_t file, const std::string& path) {
#ifdef GUTIBM_HDF5
  return H5Lexists(file, path.c_str(), H5P_DEFAULT) > 0;
#else
  (void)file; (void)path;
  return false;
#endif
}

}  // namespace

int main() {
#ifndef GUTIBM_HDF5
  std::cout << "HDF5 disabled — skipping schedule tests.\n";
  return 0;
#else
  const std::string filename = resolve_test_h5_path("GUTIBM_SCHEDULE_H5", "schedule_test");

  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {40e-6, 40e-6, 20e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.total_time = 120.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 60.0;
  cfg.seed = 111;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = filename;
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 0;
  cfg.hdf5.schedule.grid = 0;
  cfg.hdf5.schedule.lineage = 0;
  cfg.hdf5.schedule.genome = 0;
  cfg.initial_strains[0].count = 4;
  cfg.initial_strains.resize(1);

  Simulation sim;
  sim.init(cfg);
  sim.run();

  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  assert(dataset_exists(file, "summary/step_000000/time"));
  assert(!dataset_exists(file, "agents/step_000000/id"));
  assert(!dataset_exists(file, "grid/step_000000/carbon"));
  H5Fclose(file);

  std::cout << "All HDF5 schedule tests passed.\n";
  return 0;
#endif
}
