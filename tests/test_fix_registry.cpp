/* -----------------------------------------------------------------------
   GutIBM – Fix registry unit tests
   ----------------------------------------------------------------------- */

#include "fix_registry.h"
#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <iostream>

using namespace gutibm;

void test_default_fix_order() {
  auto names = FixRegistry::registered_names();
  assert(names.size() == 6);
  assert(names[0] == "metabolism");
  assert(names[1] == "bacteriocin");
  assert(names[2] == "receptor");
  assert(names[3] == "conjugation");
  assert(names[4] == "mutation");
  assert(names[5] == "mechanics");
  std::cout << "  test_default_fix_order: PASSED\n";
}

void test_create_all_fixes() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.hdf5.enabled = false;
  cfg.initial_strains.clear();

  Simulation sim;
  sim.init(cfg);

  auto names = sim.fix_names();
  assert(names.size() == 6);
  assert(names[0] == "metabolism");
  assert(names[5] == "mechanics");
  std::cout << "  test_create_all_fixes: PASSED\n";
}

void test_subset_fixes() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.hdf5.enabled = false;
  cfg.total_time = 120.0;
  cfg.initial_strains.clear();

  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 5;
  strain.mu_max = 5.0e-4;
  cfg.initial_strains.push_back(strain);

  cfg.enabled_fixes = {"metabolism", "mechanics"};

  Simulation sim;
  sim.init(cfg);
  sim.run();

  auto names = sim.fix_names();
  assert(names.size() == 2);
  assert(names[0] == "metabolism");
  assert(names[1] == "mechanics");
  std::cout << "  test_subset_fixes: PASSED\n";
}

int main() {
  std::cout << "=== Fix Registry Tests ===\n";
  test_default_fix_order();
  test_create_all_fixes();
  test_subset_fixes();
  std::cout << "All fix registry tests passed.\n";
  return 0;
}
