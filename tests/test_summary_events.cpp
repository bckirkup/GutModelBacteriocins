/* -----------------------------------------------------------------------
   GutIBM – Step event counter tests (Spec 4 summary layer)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "path_utils.h"

#include <cassert>
#include <iostream>
#include <string>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

using namespace gutibm;

namespace {

#ifdef GUTIBM_HDF5
int32_t read_event(hid_t file, const std::string& path) {
  hid_t dset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dset >= 0);
  int32_t value = 0;
  H5Dread(dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
  H5Dclose(dset);
  return value;
}
#endif

}  // namespace

int main() {
#ifndef GUTIBM_HDF5
  std::cout << "HDF5 disabled — skipping summary event tests.\n";
  return 0;
#else
  const std::string filename = resolve_test_h5_path("GUTIBM_EVENTS_H5", "summary_events");

  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {40e-6, 40e-6, 20e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.total_time = 3600.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 60.0;
  cfg.seed = 333;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = filename;
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 0;
  cfg.hdf5.schedule.grid = 0;
  cfg.hdf5.schedule.lineage = 0;
  cfg.hdf5.schedule.genome = 0;
  cfg.initial_strains[0].count = 20;
  cfg.initial_strains.resize(1);

  Simulation sim;
  sim.init(cfg);
  sim.run();

  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);

  const int32_t divisions = read_event(file, "summary/step_000001/events/divisions");
  assert(divisions >= 0);

  H5Fclose(file);
  std::cout << "All summary event counter tests passed.\n";
  return 0;
#endif
}
