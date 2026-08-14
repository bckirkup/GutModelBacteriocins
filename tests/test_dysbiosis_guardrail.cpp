/* -----------------------------------------------------------------------
   GutIBM – Operating-envelope dysbiosis guardrail tests
   ----------------------------------------------------------------------- */

#include "dysbiosis_guard.h"

#include <cassert>
#include <iostream>
#include <vector>

using namespace gutibm;

int main() {
  constexpr Real threshold = 1.0e8;
  constexpr Int sample_count = 7;
  constexpr Real interval = 300.0;

  const std::vector accelerating{
      1.01e8, 1.02e8, 1.04e8, 1.07e8, 1.11e8, 1.16e8, 1.22e8};
  DysbiosisGuard accelerating_guard(threshold, interval, sample_count);
  accelerating_guard.reset(0.0);
  for (Int i = 0; i < sample_count; ++i) {
    accelerating_guard.observe(i * interval, accelerating[static_cast<size_t>(i)]);
  }
  assert(accelerating_guard.halted());

  const std::vector plateau{
      1.01e8, 1.02e8, 1.03e8, 1.03e8, 1.03e8, 1.03e8, 1.03e8};
  DysbiosisGuard plateau_guard(threshold, interval, sample_count);
  plateau_guard.reset(0.0);
  for (Int i = 0; i < sample_count; ++i) {
    plateau_guard.observe(i * interval, plateau[static_cast<size_t>(i)]);
  }
  assert(!plateau_guard.halted());

  const std::vector decelerating{
      1.01e8, 1.04e8, 1.06e8, 1.07e8, 1.075e8, 1.078e8, 1.079e8};
  DysbiosisGuard decelerating_guard(threshold, interval, sample_count);
  decelerating_guard.reset(0.0);
  for (Int i = 0; i < sample_count; ++i) {
    decelerating_guard.observe(i * interval,
                               decelerating[static_cast<size_t>(i)]);
  }
  assert(!decelerating_guard.halted());

  const std::vector dip_then_accelerate{
      1.05e8, 1.02e8, 1.03e8, 1.05e8, 1.08e8, 1.12e8, 1.17e8};
  DysbiosisGuard dip_guard(threshold, interval, sample_count);
  dip_guard.reset(0.0);
  for (Int i = 0; i < sample_count; ++i) {
    dip_guard.observe(i * interval,
                      dip_then_accelerate[static_cast<size_t>(i)]);
  }
  assert(!dip_guard.halted());

  DysbiosisGuard short_interval_guard(threshold, 1.0, sample_count);
  short_interval_guard.reset(0.0);
  short_interval_guard.observe(0.0, 1.01e8);
  short_interval_guard.observe(60.0, 1.02e8);
  assert(short_interval_guard.density_history().size() == 2);
  short_interval_guard.reset(60.0);
  assert(short_interval_guard.density_history().empty());
  assert(!short_interval_guard.halted());

  std::cout << "Operating-envelope dysbiosis guardrail tests passed.\n";
  return 0;
}
