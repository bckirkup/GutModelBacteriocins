#ifndef GUTIBM_VBF_CARBON_SINK_H
#define GUTIBM_VBF_CARBON_SINK_H

#include <cmath>

#ifdef __CUDACC__
#define GUTIBM_VBF_HOST_DEVICE __host__ __device__
#else
#define GUTIBM_VBF_HOST_DEVICE
#endif

namespace gutibm::vbf {

GUTIBM_VBF_HOST_DEVICE inline double implicit_carbon_sink(
    double concentration, double vmax, double km, double dt) {
  if (concentration <= 0.0 || vmax <= 0.0 || dt <= 0.0) {
    return 0.0;
  }
  if (km <= 0.0) {
    return concentration / dt;
  }
  const double linear_term = km + dt * vmax - concentration;
  const double discriminant = linear_term * linear_term
      + 4.0 * concentration * km;
  const double root = std::sqrt(discriminant);
  const double concentration_after_sink = linear_term >= 0.0
      ? 2.0 * concentration * km / (linear_term + root)
      : 0.5 * (-linear_term + root);
  return (concentration - concentration_after_sink) / dt;
}

}  // namespace gutibm::vbf

#undef GUTIBM_VBF_HOST_DEVICE

#endif  // GUTIBM_VBF_CARBON_SINK_H
