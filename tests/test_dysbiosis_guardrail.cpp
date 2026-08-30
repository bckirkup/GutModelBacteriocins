/* -----------------------------------------------------------------------
   GutIBM – Operating-envelope dysbiosis guardrail tests
   ----------------------------------------------------------------------- */

#include "dysbiosis_guard.h"

#include <cassert>
#include <iostream>
#include <vector>

using namespace gutibm;

static Real observe_series(DysbiosisGuard& guard,
                           const std::vector<Real>& densities,
                           Real interval) {
  for (size_t i = 0; i < densities.size(); ++i) {
    const Real time = static_cast<Real>(i) * interval;
    guard.observe(time, densities[i]);
    if (guard.halted()) return time;
  }
  return -1.0;
}

int main() {
  constexpr Real threshold = 1.0e8;
  constexpr Int sample_count = 7;
  constexpr Real interval = 300.0;

  const std::vector accelerating{
      1.01e8, 1.02e8, 1.04e8, 1.07e8, 1.11e8, 1.16e8, 1.22e8};
  DysbiosisGuard accelerating_guard(threshold, interval, sample_count);
  accelerating_guard.reset(0.0);
  observe_series(accelerating_guard, accelerating, interval);
  assert(accelerating_guard.halted());

  const std::vector plateau{
      1.01e8, 1.02e8, 1.03e8, 1.03e8, 1.03e8, 1.03e8, 1.03e8};
  DysbiosisGuard plateau_guard(threshold, interval, sample_count);
  plateau_guard.reset(0.0);
  observe_series(plateau_guard, plateau, interval);
  assert(!plateau_guard.halted());

  const std::vector decelerating{
      1.01e8, 1.04e8, 1.06e8, 1.07e8, 1.075e8, 1.078e8, 1.079e8};
  DysbiosisGuard decelerating_guard(threshold, interval, sample_count);
  decelerating_guard.reset(0.0);
  observe_series(decelerating_guard, decelerating, interval);
  assert(!decelerating_guard.halted());

  const std::vector noisy_rise{
      1.01e8, 1.05e8, 1.04e8, 1.08e8, 1.12e8, 1.15e8, 1.21e8};
  DysbiosisGuard dip_guard(threshold, interval, sample_count);
  dip_guard.reset(0.0);
  observe_series(dip_guard, noisy_rise, interval);
  assert(dip_guard.halted());
  assert(dip_guard.density_rate_cells_per_mL_per_s() > 0.0);

  const std::vector below_threshold{
      0.90e8, 0.93e8, 0.92e8, 0.96e8, 0.99e8, 0.98e8, 1.00e8, 1.03e8};
  DysbiosisGuard below_guard(threshold, interval, sample_count);
  below_guard.reset(0.0);
  observe_series(below_guard, below_threshold, interval);
  assert(!below_guard.halted());

  const std::vector declining{
      1.30e8, 1.27e8, 1.25e8, 1.20e8, 1.18e8, 1.12e8, 1.08e8};
  DysbiosisGuard declining_guard(threshold, interval, sample_count);
  declining_guard.reset(0.0);
  observe_series(declining_guard, declining, interval);
  assert(!declining_guard.halted());

  const std::vector flat{
      1.20e8, 1.20e8, 1.20e8, 1.20e8, 1.20e8, 1.20e8, 1.20e8};
  DysbiosisGuard flat_guard(threshold, interval, sample_count);
  flat_guard.reset(0.0);
  observe_series(flat_guard, flat, interval);
  assert(!flat_guard.halted());

  const std::vector growth_rates{1.0e6, 2.0e6, 4.0e6};
  Real previous_halt_time = 1.0e30;
  for (const Real growth_rate : growth_rates) {
    std::vector<Real> growth;
    for (Int i = 0; i < 20; ++i) {
      growth.push_back(0.90e8 + growth_rate * static_cast<Real>(i));
    }
    DysbiosisGuard growth_guard(threshold, interval, sample_count);
    growth_guard.reset(0.0);
    const Real halt_time = observe_series(growth_guard, growth, interval);
    assert(halt_time >= 0.0);
    assert(halt_time <= previous_halt_time);
    previous_halt_time = halt_time;
  }

  const std::vector zero_growth(20, 0.90e8);
  DysbiosisGuard zero_growth_guard(threshold, interval, sample_count);
  zero_growth_guard.reset(0.0);
  observe_series(zero_growth_guard, zero_growth, interval);
  assert(!zero_growth_guard.halted());

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
