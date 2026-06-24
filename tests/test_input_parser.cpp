/* -----------------------------------------------------------------------
   GutIBM – Example input file parsing tests
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include <string>

#ifndef GUTIBM_SOURCE_DIR
#define GUTIBM_SOURCE_DIR "."
#endif

using namespace gutibm;

void test_single_colony_example() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/examples/single_colony/input.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(std::abs(cfg.total_time - 86400.0) < 1e-6);
  assert(std::abs(cfg.bio_dt - 60.0) < 1e-6);
  assert(cfg.seed == 12345);
  assert(std::abs(cfg.domain.hi[0] - 0.001) < 1e-12);
  assert(cfg.hdf5.filename == "single_colony_output.h5");
  assert(cfg.hdf5.dump_every == 60);
  std::cout << "  test_single_colony_example: PASSED\n";
}

void test_diversity_paradox_example() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/examples/diversity_paradox/input.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(std::abs(cfg.total_time - 604800.0) < 1e-6);
  assert(cfg.seed == 42);
  assert(std::abs(cfg.domain.hi[0] - 0.002) < 1e-12);
  assert(cfg.hdf5.filename == "diversity_paradox_output.h5");
  assert(cfg.qssa.use_fmm == true);
  assert(std::abs(cfg.qssa.fmm_theta - 0.5) < 1e-12);
  assert(std::abs(cfg.qssa.toxin_cutoff - 200e-6) < 1e-15);
  assert(std::abs(cfg.qssa.nutrient_cutoff - 50e-6) < 1e-15);
  std::cout << "  test_diversity_paradox_example: PASSED\n";
}

void test_single_colony_peristaltic_keys() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/examples/single_colony/input.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(cfg.advection.peristaltic_enabled == true);
  assert(std::abs(cfg.advection.peristaltic_period - 20.0) < 1e-12);
  assert(std::abs(cfg.advection.peristaltic_amplitude - 0.5) < 1e-12);
  assert(std::abs(cfg.advection.peristaltic_wavelength - 0.001) < 1e-15);
  std::cout << "  test_single_colony_peristaltic_keys: PASSED\n";
}

void test_fmm_peristaltic_fixture() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/parser_fmm_peristaltic.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(std::abs(cfg.total_time - 3600.0) < 1e-6);
  assert(cfg.seed == 99);
  assert(cfg.advection.peristaltic_enabled == true);
  assert(std::abs(cfg.advection.peristaltic_period - 15.0) < 1e-12);
  assert(std::abs(cfg.advection.peristaltic_amplitude - 0.3) < 1e-12);
  assert(std::abs(cfg.advection.peristaltic_wavelength - 0.0005) < 1e-15);
  assert(cfg.qssa.use_fmm == true);
  assert(std::abs(cfg.qssa.fmm_theta - 0.3) < 1e-12);
  assert(std::abs(cfg.qssa.toxin_cutoff - 100e-6) < 1e-15);
  assert(std::abs(cfg.qssa.nutrient_cutoff - 25e-6) < 1e-15);
  std::cout << "  test_fmm_peristaltic_fixture: PASSED\n";
}

int main() {
  std::cout << "=== Input Parser Example Tests ===\n";
  test_single_colony_example();
  test_single_colony_peristaltic_keys();
  test_diversity_paradox_example();
  test_fmm_peristaltic_fixture();
  std::cout << "All input parser example tests passed.\n";
  return 0;
}
