#ifndef GUTIBM_VBF_CARBON_SINK_H
#define GUTIBM_VBF_CARBON_SINK_H

#include "types.h"

#include <cmath>

#ifdef __CUDACC__
#define GUTIBM_VBF_HOST_DEVICE __host__ __device__
#else
#define GUTIBM_VBF_HOST_DEVICE
#endif

namespace gutibm::vbf {

GUTIBM_VBF_HOST_DEVICE inline Real implicit_carbon_sink(
    Real concentration, Real vmax, Real km, Real dt) {
  if (concentration <= 0.0 || vmax <= 0.0 || dt <= 0.0) {
    return 0.0;
  }
  if (km <= 0.0) {
    return concentration / dt;
  }
  const Real linear_term = km + dt * vmax - concentration;
  const Real discriminant = linear_term * linear_term
      + 4.0 * concentration * km;
  const Real root = std::sqrt(discriminant);
  const Real concentration_after_sink = linear_term >= 0.0
      ? 2.0 * concentration * km / (linear_term + root)
      : 0.5 * (-linear_term + root);
  return (concentration - concentration_after_sink) / dt;
}

}  // namespace gutibm::vbf

#undef GUTIBM_VBF_HOST_DEVICE

#endif  // GUTIBM_VBF_CARBON_SINK_H
