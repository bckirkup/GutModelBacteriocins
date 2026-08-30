/* -----------------------------------------------------------------------
   GutIBM – requested HDF5 output fails closed (Spec 4)
   ----------------------------------------------------------------------- */

#include "domain.h"
#include "error.h"
#include "hdf5_writer.h"
#include "path_utils.h"

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

namespace {

Domain make_domain() {
  DomainConfig cfg;
  cfg.hi = {20e-6, 20e-6, 20e-6};
  cfg.grid_dx = 5e-6;
  Domain domain;
  domain.init(cfg);
  return domain;
}

}  // namespace

int main() {
#ifndef GUTIBM_HDF5
  std::cout << "HDF5 disabled — skipping fail-closed tests.\n";
  return 0;
#else
  const Domain domain = make_domain();

  HDF5Config disabled_cfg;
  disabled_cfg.enabled = false;
  disabled_cfg.filename = "/missing/hdf5/output.h5";
  HDF5Writer disabled;
  disabled.init(disabled_cfg, domain);
  assert(!disabled.is_enabled());

  const std::string valid_path =
      resolve_test_h5_path("GUTIBM_FAIL_CLOSED_H5", "fail_closed");
  HDF5Config valid_cfg;
  valid_cfg.filename = valid_path;
  valid_cfg.schedule.summary = 1;
  valid_cfg.schedule.agents = 0;
  valid_cfg.schedule.grid = 0;
  valid_cfg.schedule.lineage = 0;
  valid_cfg.schedule.genome = 0;
  valid_cfg.schedule.provenance = 0;
  HDF5Writer valid;
  valid.init(valid_cfg, domain);
  assert(valid.is_enabled());
  valid.finalize();

  hid_t file = H5Fopen(valid_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  H5Fclose(file);
  assert(std::filesystem::remove(valid_path));
  assert(!std::filesystem::exists(valid_path));

  const std::string marker = secure_temp_file("gutibm_fail_closed_");
  std::filesystem::remove(marker);
  const std::filesystem::path invalid_path =
      std::filesystem::path(marker + "_missing") / "output.h5";
  HDF5Config invalid_cfg = valid_cfg;
  invalid_cfg.filename = invalid_path.string();
  HDF5Writer invalid;
  bool threw = false;
  std::string message;
  try {
    invalid.init(invalid_cfg, domain);
  } catch (const IOError& error) {
    threw = true;
    message = error.what();
  }
  assert(threw);
  assert(message.find(invalid_cfg.filename) != std::string::npos);
  assert(message.find("output directory") != std::string::npos);

  std::cout << "All HDF5 fail-closed tests passed.\n";
  return 0;
#endif
}
