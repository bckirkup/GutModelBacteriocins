/* -----------------------------------------------------------------------
   GutIBM – Host/device mechanics participation predicate tests
   ----------------------------------------------------------------------- */

#include "mechanics_participation.h"
#include "types.h"

#include <cassert>
#include <iostream>

using namespace gutibm;

static_assert(kDeadStateValue == to_underlying(PhenoState::DEAD));

int main() {
  constexpr double sim_time = 1000.0;
  constexpr double persistence = 300.0;

  assert(mechanics_participates(
      to_underlying(PhenoState::NORMAL), -1.0, sim_time, 0, 0.0));
  assert(mechanics_participates(
      to_underlying(PhenoState::NORMAL), 1000.0, sim_time, 0, 0.0));

  assert(mechanics_participates(
      kDeadStateValue, sim_time - 1.0, sim_time, 1, persistence));
  assert(!mechanics_participates(
      kDeadStateValue, sim_time - 1.0, sim_time, 0, persistence));
  assert(mechanics_participates(
      kDeadStateValue, sim_time - persistence + 1.0e-9, sim_time, 1,
      persistence));
  assert(!mechanics_participates(
      kDeadStateValue, sim_time - persistence - 1.0e-9, sim_time, 1,
      persistence));
  assert(!mechanics_participates(
      kDeadStateValue, -1.0, sim_time, 1, persistence));

  std::cout << "All mechanics participation tests passed.\n";
  return 0;
}
