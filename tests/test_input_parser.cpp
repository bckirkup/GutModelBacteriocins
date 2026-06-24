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
  std::cout << "  test_diversity_paradox_example: PASSED\n";
}

int main() {
  std::cout << "=== Input Parser Example Tests ===\n";
  test_single_colony_example();
  test_diversity_paradox_example();
  std::cout << "All input parser example tests passed.\n";
  return 0;
}
