/* -----------------------------------------------------------------------
   GutIBM – Per-kill provenance output test
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "path_utils.h"
#include "simulation.h"
#include "species_names.h"

#include <cassert>
#include <filesystem>
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

int count_cause(const std::vector<int32_t>& causes, ProvenanceCause cause);

void test_lysis_provenance() {
  const std::string filename =
      resolve_test_h5_path("GUTIBM_LYSIS_PROVENANCE_H5", "lysis_provenance");
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {40e-6, 40e-6, 20e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.time.total_time = 360.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 60.0;
  cfg.seed = 54321;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = filename;
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.provenance = 1;
  cfg.initial_strains.clear();
  cfg.enabled_fixes = {"bacteriocin"};
  cfg.fixes.bacteriocin.sos_basal_rate = 1.0;
  cfg.advection.radial_turnover = 1.0e12;
  cfg.advection.distal_transit_time = 1.0e12;
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 1;
  strain.mu_max = 5e-4;
  strain.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(strain);
  SimulationConfig::InitialStrain survivor;
  survivor.type = 2;
  survivor.count = 1;
  survivor.mu_max = 5e-4;
  cfg.initial_strains.push_back(survivor);

  Simulation sim;
  sim.init(cfg);
  sim.run();

  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  const std::string provenance_path = std::format(
      "provenance/step_{:06}/cause", sim.step_count());
  const auto causes = read_int_vector(file, provenance_path);
  assert(count_cause(causes, ProvenanceCause::LYSIS) == 1);
  H5Fclose(file);
}

double read_double_scalar(hid_t file, const std::string& path) {
  hid_t dataset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dataset >= 0);
  double value = 0.0;
  assert(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, &value) >= 0);
  H5Dclose(dataset);
  return value;
}

int count_cause(const std::vector<int32_t>& causes, ProvenanceCause cause) {
  int count = 0;
  for (const int32_t value : causes) {
    if (value == to_underlying(cause)) ++count;
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
  cfg.restart.enabled = true;
  cfg.restart.directory = std::filesystem::path(filename).parent_path().string()
      + "/kill_provenance_restart";
  cfg.restart.interval_steps = 1;
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
  assert(read_scalar(file, "summary/step_000001/events/mortality_colicin") == 1);
  assert(read_scalar(file, "summary/step_000001/events/mortality_cdi") == 0);
  assert(read_scalar(file, "summary/step_000001/events/outflow_washout") == 0);
  assert(read_scalar(file, "summary/step_000001/events/outflow_boundary") == 0);
  assert(read_scalar(file, "summary/step_000001/events/interval_start_step") == 1);
  assert(read_scalar(file, "summary/step_000001/events/interval_end_step") == 1);
  assert(read_double_scalar(
             file, "summary/step_000001/events/interval_start_time") == 0.0);
  assert(read_double_scalar(
             file, "summary/step_000001/events/interval_end_time") == 60.0);
  assert(read_scalar(
             file, "summary/step_000001/events/cumulative_mortality_colicin") == 1);
  assert(causes.size() == 1);
  assert(causes[0] == static_cast<int32_t>(ProvenanceCause::COLICIN));
  assert(count_cause(causes, ProvenanceCause::COLICIN) == 1);
  assert(count_cause(causes, ProvenanceCause::CDI) == 0);
  assert(count_cause(causes, ProvenanceCause::WASHOUT) == 0);
  assert(count_cause(causes, ProvenanceCause::BOUNDARY) == 0);
  assert(count_cause(causes, ProvenanceCause::LYSIS) == 0);
  assert(concentrations.size() == 4);
  assert(occupancies.size() == 4);
  assert(hazards.size() == 4);
  assert(concentrations[2] > 0.0);
  assert(occupancies[2] > 0.0);
  assert(hazards[2] > 0.0);
  H5Fclose(file);

  const std::string restart_filename = cfg.restart.directory + "/step_000001.h5";
  hid_t restart_file = H5Fopen(
      restart_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(restart_file >= 0);
  const auto restart_causes =
      read_int_vector(restart_file, "provenance/step_000001/cause");
  assert(restart_causes.size() == 1);
  assert(read_scalar(
             restart_file,
             "summary/step_000001/events/interval_start_step") == 1);
  assert(read_scalar(
             restart_file,
             "summary/step_000001/events/interval_end_step") == 1);
  assert(read_scalar(
             restart_file,
             "summary/step_000001/events/cumulative_mortality_colicin") == 1);
  H5Fclose(restart_file);

  test_lysis_provenance();
  std::cout << "All kill provenance tests passed.\n";
  return 0;
#endif
}
