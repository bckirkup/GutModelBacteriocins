/* -----------------------------------------------------------------------
   GutIBM – Per-kill provenance output test
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "path_utils.h"
#include "simulation.h"
#include "species_names.h"

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

using namespace gutibm;

#ifdef GUTIBM_HDF5
namespace {

std::vector<double> read_double_vector(hid_t file, const std::string& path) {
  hid_t dataset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dataset >= 0);
  hid_t space = H5Dget_space(dataset);
  assert(space >= 0);
  hsize_t extent = 0;
  assert(H5Sget_simple_extent_dims(space, &extent, nullptr) == 1);
  std::vector<double> values(static_cast<size_t>(extent));
  if (extent > 0) {
    assert(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                   H5P_DEFAULT, values.data()) >= 0);
  }
  H5Sclose(space);
  H5Dclose(dataset);
  return values;
}

std::vector<int32_t> read_int_vector(hid_t file, const std::string& path) {
  hid_t dataset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dataset >= 0);
  hid_t space = H5Dget_space(dataset);
  assert(space >= 0);
  hsize_t extent = 0;
  assert(H5Sget_simple_extent_dims(space, &extent, nullptr) == 1);
  std::vector<int32_t> values(static_cast<size_t>(extent));
  if (extent > 0) {
    assert(H5Dread(dataset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL,
                   H5P_DEFAULT, values.data()) >= 0);
  }
  H5Sclose(space);
  H5Dclose(dataset);
  return values;
}

int32_t read_scalar(hid_t file, const std::string& path) {
  hid_t dataset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dataset >= 0);
  int32_t value = 0;
  assert(H5Dread(dataset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, &value) >= 0);
  H5Dclose(dataset);
  return value;
}

int count_cause(const std::vector<int32_t>& causes, ProvenanceCause cause) {
  int count = 0;
  for (const int32_t value : causes) {
    if (value == static_cast<int32_t>(cause)) ++count;
  }
  return count;
}

}  // namespace
#endif

int main() {
#ifndef GUTIBM_HDF5
  std::cout << "HDF5 disabled — skipping kill provenance tests.\n";
  return 0;
#else
  const std::string filename =
      resolve_test_h5_path("GUTIBM_KILL_PROVENANCE_H5", "kill_provenance");

  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {40e-6, 40e-6, 20e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.total_time = 60.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 60.0;
  cfg.seed = 12345;
  cfg.enabled_fixes = {"bacteriocin", "receptor"};
  cfg.qssa.microcin_secretion = 1.0e-12;
  cfg.fixes.receptor.kill_rate_microcin = 1.0;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = filename;
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 0;
  cfg.hdf5.schedule.grid = 0;
  cfg.hdf5.schedule.lineage = 0;
  cfg.hdf5.schedule.genome = 0;
  cfg.hdf5.schedule.provenance = 1;
  cfg.initial_strains.clear();
  cfg.initial_strains.push_back({1, 1, 5e-4, {"MccV"}});
  cfg.initial_strains.push_back({2, 1, 5e-4, {}});

  Simulation sim;
  sim.init(cfg);
  Agent& producer = sim.agents()[0];
  Agent& victim = sim.agents()[1];
  producer.receptor_expr[to_underlying(ReceptorType::CirA)] = 0.0;
  victim.receptor_expr[to_underlying(ReceptorType::BtuB)] = 0.0;
  victim.receptor_expr[to_underlying(ReceptorType::FepA)] = 0.0;
  victim.receptor_expr[to_underlying(ReceptorType::FhuA)] = 0.0;
  victim.x = producer.x;
  victim.grid_cell = producer.grid_cell;

  sim.run();

  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  const auto causes =
      read_int_vector(file, "provenance/step_000001/cause");
  const auto concentrations =
      read_double_vector(file, "provenance/step_000001/toxin_concentration");
  const auto occupancies =
      read_double_vector(file, "provenance/step_000001/toxin_occupancy");
  const auto hazards =
      read_double_vector(file, "provenance/step_000001/toxin_hazard");
  assert(read_scalar(file, "summary/step_000001/events/colicin_kills") == 1);
  assert(read_scalar(file, "summary/step_000001/events/cdi_kills") == 0);
  assert(read_scalar(file, "summary/step_000001/events/washout_deaths") == 0);
  assert(read_scalar(file, "summary/step_000001/events/boundary_deaths") == 0);
  assert(read_scalar(file, "summary/step_000001/events/starvation_deaths") == 0);
  assert(causes.size() == 1);
  assert(causes[0] == static_cast<int32_t>(ProvenanceCause::COLICIN));
  assert(count_cause(causes, ProvenanceCause::COLICIN) == 1);
  assert(count_cause(causes, ProvenanceCause::CDI) == 0);
  assert(count_cause(causes, ProvenanceCause::WASHOUT) == 0);
  assert(count_cause(causes, ProvenanceCause::BOUNDARY) == 0);
  assert(count_cause(causes, ProvenanceCause::STARVATION) == 0);
  assert(concentrations.size() == 4);
  assert(occupancies.size() == 4);
  assert(hazards.size() == 4);
  assert(concentrations[2] > 0.0);
  assert(occupancies[2] > 0.0);
  assert(hazards[2] > 0.0);
  H5Fclose(file);

  std::cout << "All kill provenance tests passed.\n";
  return 0;
#endif
}
