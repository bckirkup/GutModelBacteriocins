/* -----------------------------------------------------------------------
   GutIBM – Operating-envelope dysbiosis guardrail tests
   ----------------------------------------------------------------------- */

#include "simulation.h"

#include <cassert>
#include <iostream>
#include <vector>

using namespace gutibm;

int main() {
  constexpr Real threshold = 1.0e8;
  constexpr Int sample_count = 7;

  const std::vector<Real> accelerating{
      1.01e8, 1.02e8, 1.04e8, 1.07e8, 1.11e8, 1.16e8, 1.22e8};
  assert(is_accelerating_density_window(accelerating, threshold,
                                        sample_count));

  const std::vector<Real> plateau{
      1.01e8, 1.02e8, 1.03e8, 1.03e8, 1.03e8, 1.03e8, 1.03e8};
  assert(!is_accelerating_density_window(plateau, threshold, sample_count));

  const std::vector<Real> decelerating{
      1.01e8, 1.04e8, 1.06e8, 1.07e8, 1.075e8, 1.078e8, 1.079e8};
  assert(!is_accelerating_density_window(decelerating, threshold,
                                         sample_count));

  std::cout << "Operating-envelope dysbiosis guardrail tests passed.\n";
  return 0;
}
