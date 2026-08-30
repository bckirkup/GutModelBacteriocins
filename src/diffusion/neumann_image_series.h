/* -----------------------------------------------------------------------
   GutIBM – Shared host/device two-wall Neumann image series
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_NEUMANN_IMAGE_SERIES_H
#define GUTIBM_NEUMANN_IMAGE_SERIES_H

#include <cmath>

#ifdef __CUDACC__
#define GUTIBM_NEUMANN_HOST_DEVICE __host__ __device__
#else
#define GUTIBM_NEUMANN_HOST_DEVICE
#endif

namespace gutibm::neumann {

constexpr int kMaxImageShells = 512;
constexpr int kLowScreeningShells = 8;
constexpr double kLowScreeningFloorThreshold = 0.0225;
constexpr double kRelativeTolerance = 1.0e-10;

struct ImageSeriesBudget {
  int max_shells;
  bool low_screening_floor;
};

GUTIBM_NEUMANN_HOST_DEVICE inline double series_screening(
    double d_eff, double decay_rate, double flow_magnitude, double flow_z) {
  const double screened_speed = sqrt(
      flow_magnitude * flow_magnitude + 4.0 * d_eff * decay_rate);
  return (screened_speed - fabs(flow_z)) / (2.0 * d_eff);
}

GUTIBM_NEUMANN_HOST_DEVICE inline ImageSeriesBudget image_series_budget(
    double d_eff, double decay_rate, double flow_magnitude, double flow_z,
    double height, double rel_tol = kRelativeTolerance) {
  const double k_series = series_screening(
      d_eff, decay_rate, flow_magnitude, flow_z);
  const double scaled_screening = k_series * height;
  if (!(scaled_screening >= kLowScreeningFloorThreshold)) {
    return {kLowScreeningShells, true};
  }
  if (!(rel_tol > 0.0)) {
    return {kMaxImageShells, false};
  }
  const double estimated_shells = ceil(
      log(1.0 / rel_tol) / (2.0 * scaled_screening));
  double shell_value = estimated_shells;
  if (shell_value < 1.0) {
    shell_value = 1.0;
  } else if (shell_value > static_cast<double>(kMaxImageShells)) {
    shell_value = static_cast<double>(kMaxImageShells);
  }
  const auto shell_count = static_cast<int>(shell_value);
  return {shell_count, false};
}

// Kernel evaluates one image at image_z. reflected is nonzero for family B,
// whose odd z-reflection reverses the z component of the flow.
template <typename Kernel>
GUTIBM_NEUMANN_HOST_DEVICE inline double sum_image_series(
    double source_z, double z_lo, double z_hi, Kernel kernel,
    double rel_tol = kRelativeTolerance, int max_shells = kMaxImageShells,
    int* shell_count = nullptr, int* cap_hit = nullptr) {
  const double height = z_hi - z_lo;
  const double reflected_source_z = 2.0 * z_lo - source_z;
  double total = kernel(source_z, 0) + kernel(reflected_source_z, 1);
  bool converged = false;
  int completed_shells = 0;

  // Family B at m=1 is 2*z_hi - source_z, the upper-wall reflection.
  for (int m = 1; m <= max_shells; ++m) {
    const double offset = 2.0 * static_cast<double>(m) * height;
    double shell = 0.0;
    shell += kernel(source_z + offset, 0);
    shell += kernel(source_z - offset, 0);
    shell += kernel(reflected_source_z + offset, 1);
    shell += kernel(reflected_source_z - offset, 1);
    total += shell;
    completed_shells = m;
    if (rel_tol > 0.0 && fabs(shell) < rel_tol * fabs(total)) {
      converged = true;
      break;
    }
  }

  if (shell_count != nullptr) {
    *shell_count = completed_shells;
  }
  if (cap_hit != nullptr) {
    *cap_hit = converged ? 0 : 1;
  }
  return total;
}

template <typename Kernel>
GUTIBM_NEUMANN_HOST_DEVICE inline double sum_image_series(
    double source_z, double z_lo, double z_hi, Kernel kernel,
    double d_eff, double decay_rate, double flow_magnitude, double flow_z,
    double rel_tol = kRelativeTolerance, int* shell_count = nullptr,
    int* cap_hit = nullptr, int* low_screening_floor = nullptr) {
  const ImageSeriesBudget budget = image_series_budget(
      d_eff, decay_rate, flow_magnitude, flow_z, z_hi - z_lo, rel_tol);
  const double total = sum_image_series(
      source_z, z_lo, z_hi, kernel, rel_tol, budget.max_shells,
      shell_count, cap_hit);
  if (low_screening_floor != nullptr) {
    *low_screening_floor = budget.low_screening_floor ? 1 : 0;
  }
  return total;
}

}  // namespace gutibm::neumann

#undef GUTIBM_NEUMANN_HOST_DEVICE

#endif  // GUTIBM_NEUMANN_IMAGE_SERIES_H
