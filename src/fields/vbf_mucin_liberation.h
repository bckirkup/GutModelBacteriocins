#ifndef GUTIBM_VBF_MUCIN_LIBERATION_H
#define GUTIBM_VBF_MUCIN_LIBERATION_H

#include <cmath>

#ifdef __CUDACC__
#define GUTIBM_VBF_HOST_DEVICE __host__ __device__
#else
#define GUTIBM_VBF_HOST_DEVICE
#endif

namespace gutibm::vbf {

// Exact backward-Euler Monod mucin liberation rate over timestep dt:
// solves M_new = M - dt * vmax * M_new / (km + M_new), where
// vmax = k_liberation * vbf_density * z_weight.
// Unconditionally guarantees 0 <= rate * dt <= max(0, M), preventing negative
// mucin concentrations and unphysical mass creation while recovering
// vmax * M / (km + M) + O(dt) when dt * vmax << M + km.
GUTIBM_VBF_HOST_DEVICE inline double implicit_mucin_liberation(
    double mucin_conc, double k_liberation, double vbf_density,
    double z_weight, double km, double dt) {
  if (mucin_conc <= 0.0 || k_liberation <= 0.0 || vbf_density <= 0.0
      || z_weight <= 0.0 || dt <= 0.0) {
    return 0.0;
  }
  const double vmax = k_liberation * vbf_density * z_weight;
  if (vmax <= 0.0) {
    return 0.0;
  }

  double rate = 0.0;
  if (km <= 0.0) {
    rate = (vmax * dt <= mucin_conc) ? vmax : (mucin_conc / dt);
  } else {
    const double linear_term = km + dt * vmax - mucin_conc;
    const double discriminant =
        linear_term * linear_term + 4.0 * mucin_conc * km;
    const double root = std::sqrt(discriminant);
    const double mucin_after = linear_term >= 0.0
        ? (2.0 * mucin_conc * km) / (linear_term + root)
        : 0.5 * (-linear_term + root);
    if (mucin_after > 0.5 * mucin_conc) {
      rate = vmax * mucin_after / (km + mucin_after);
    } else {
      rate = (mucin_conc - mucin_after) / dt;
    }
  }

  // Guard against 1-ULP IEEE-754 rounding overshoot when mucin_after << M.
  constexpr double kOneMinusEps = 1.0 - 2.220446049250313e-16;
  while (mucin_conc - rate * dt < 0.0 && rate > 0.0) {
    rate *= kOneMinusEps;
  }
  return rate;
}

}  // namespace gutibm::vbf

#undef GUTIBM_VBF_HOST_DEVICE

#endif  // GUTIBM_VBF_MUCIN_LIBERATION_H
