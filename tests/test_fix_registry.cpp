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

  // Six fixes registered and initialized without crash
  std::cout << "  test_create_all_fixes: PASSED\n";
}

int main() {
  std::cout << "=== Fix Registry Tests ===\n";
  test_default_fix_order();
  test_create_all_fixes();
  std::cout << "All fix registry tests passed.\n";
  return 0;
}
