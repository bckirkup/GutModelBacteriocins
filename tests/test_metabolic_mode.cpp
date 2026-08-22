#include "metabolic_mode.h"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
  using gutibm::metabolic_mode::acid_inhibition;
  using gutibm::metabolic_mode::fermentation_fraction;
  using gutibm::metabolic_mode::interpolate;
  using gutibm::metabolic_mode::relax;

  const double low_oxygen = fermentation_fraction(0.0, 1.0e-4, 3.0e-4);
  const double medium_oxygen = fermentation_fraction(0.5, 1.0e-4, 3.0e-4);
  const double high_oxygen = fermentation_fraction(1.0, 1.0e-4, 3.0e-4);
  assert(low_oxygen >= 0.0 && low_oxygen <= 1.0);
  assert(low_oxygen > medium_oxygen && medium_oxygen > high_oxygen);

  const double overflow = fermentation_fraction(1.0, 6.0e-4, 3.0e-4);
  assert(overflow > high_oxygen);
  assert(interpolate(1.0, 4.1, overflow)
         > interpolate(1.0, 4.1, high_oxygen));
  const double respiring_cost = interpolate(1.0, 4.1, high_oxygen);
  const double fermenting_cost = interpolate(1.0, 4.1, low_oxygen);
  assert(fermenting_cost > respiring_cost);
  assert(fermenting_cost <= 4.1);

  const double relaxed = relax(0.0, 1.0, 3600.0, 3600.0);
  assert(std::abs(relaxed - (1.0 - std::exp(-1.0))) < 1.0e-12);
  const double acid_low = acid_inhibition(20.0, 5.0, 4.76, 50.0, 0.8);
  const double acid_mid = acid_inhibition(80.0, 5.0, 4.76, 50.0, 0.8);
  const double acid_high = acid_inhibition(160.0, 5.0, 4.76, 50.0, 0.8);
  assert(acid_low < acid_mid && acid_mid < acid_high);
  assert(acid_high <= 0.8);
  assert(acid_inhibition(80.0, 5.0, 4.76, 50.0, 0.8)
         > acid_inhibition(80.0, 6.0, 4.76, 50.0, 0.8));
  assert(acid_inhibition(0.0, 6.0, 4.76, 50.0, 0.8) == 0.0);

  std::cout << "metabolic mode helper tests passed\n";
  return 0;
}
