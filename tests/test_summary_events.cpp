/* -----------------------------------------------------------------------
   GutIBM – Step event counter tests (Spec 4 summary layer)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "path_utils.h"

#include <cassert>
#include <array>
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

#ifdef GUTIBM_HDF5
int32_t read_event(hid_t file, const std::string& path) {
  hid_t dset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dset >= 0);
  int32_t value = 0;
  H5Dread(dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
  H5Dclose(dset);
  return value;
}

bool dataset_exists(hid_t file, const std::string& path) {
  const htri_t exists = H5Lexists(file, path.c_str(), H5P_DEFAULT);
  return exists > 0;
}

std::vector<double> read_vector(hid_t file, const std::string& path) {
  hid_t dset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dset >= 0);
  hid_t space = H5Dget_space(dset);
  std::array<hsize_t, 1> dims = {0};
  assert(H5Sget_simple_extent_ndims(space) == 1);
  H5Sget_simple_extent_dims(space, dims.data(), nullptr);
  std::vector<double> values(static_cast<size_t>(dims[0]));
  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          values.data());
  H5Sclose(space);
  H5Dclose(dset);
  return values;
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
  const std::string flux_prefix = "summary/step_000001/nutrient_flux/";
  assert(dataset_exists(file, flux_prefix + "species_names"));
  assert(dataset_exists(file, flux_prefix + "boundary_area_flux_interval"));
  const std::vector<double> boundary =
      read_vector(file, flux_prefix + "boundary_interval");
  const std::vector<double> cumulative =
      read_vector(file, flux_prefix + "boundary_cumulative");
  assert(!boundary.empty());
  assert(boundary.size() == cumulative.size());
  assert(dataset_exists(file, flux_prefix + "interval_start_time"));
  assert(dataset_exists(file, flux_prefix + "interval_end_time"));

  H5Fclose(file);
  std::cout << "All summary event counter tests passed.\n";
  return 0;
#endif
}
