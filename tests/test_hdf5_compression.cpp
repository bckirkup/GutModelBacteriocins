/* -----------------------------------------------------------------------
   GutIBM – HDF5 gzip compression tests (Spec 4)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "path_utils.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sys/stat.h>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

using namespace gutibm;

namespace {

size_t file_bytes(const std::string& path) {
  struct stat st {};
  if (stat(path.c_str(), &st) != 0) return 0;
  return static_cast<size_t>(st.st_size);
}

}  // namespace

int main() {
#ifndef GUTIBM_HDF5
  std::cout << "HDF5 disabled — skipping compression tests.\n";
  return 0;
#else
  const std::string raw_file = resolve_test_h5_path("GUTIBM_COMPRESS_RAW_H5", "compress_raw");
  const std::string gzip_file = resolve_test_h5_path("GUTIBM_COMPRESS_GZIP_H5", "compress_gzip");

  auto make_cfg = [](const std::string& path, const std::string& compression) {
    SimulationConfig cfg = InputParser::default_config();
    cfg.domain.hi = {80e-6, 80e-6, 40e-6};
    cfg.domain.grid_dx = 2e-6;
    cfg.time.total_time = 60.0;
    cfg.time.bio_dt = 60.0;
    cfg.time.output_interval = 60.0;
    cfg.seed = 222;
    cfg.hdf5.enabled = true;
    cfg.hdf5.filename = path;
    cfg.hdf5.compression = compression;
    cfg.hdf5.compression_level = 4;
    cfg.hdf5.schedule.summary = 0;
    cfg.hdf5.schedule.agents = 0;
    cfg.hdf5.schedule.grid = 1;
    cfg.hdf5.schedule.lineage = 0;
    cfg.hdf5.schedule.genome = 0;
    cfg.hdf5.schedule.grid_species = {"all"};
    cfg.initial_strains[0].count = 6;
    cfg.initial_strains.resize(1);
    return cfg;
  };

  Simulation sim_raw;
  sim_raw.init(make_cfg(raw_file, "none"));
  sim_raw.run();

  Simulation sim_gzip;
  sim_gzip.init(make_cfg(gzip_file, "gzip"));
  sim_gzip.run();

  const size_t raw_size = file_bytes(raw_file);
  const size_t gzip_size = file_bytes(gzip_file);
  assert(raw_size > 0);
  assert(gzip_size > 0);
  assert(gzip_size < raw_size);

  hid_t gz = H5Fopen(gzip_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(gz >= 0);
  assert(H5Lexists(gz, "grid/step_000001/carbon", H5P_DEFAULT) > 0);
  H5Fclose(gz);

  std::cout << "All HDF5 compression tests passed.\n";
  return 0;
#endif
}
