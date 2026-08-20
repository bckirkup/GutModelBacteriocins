/* -----------------------------------------------------------------------
   GutIBM – Shared host/device agent uptake limitation closed form
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_UPTAKE_LIMIT_H
#define GUTIBM_UPTAKE_LIMIT_H

#ifdef __CUDACC__
#define GUTIBM_UPTAKE_HOST_DEVICE __host__ __device__
#else
#define GUTIBM_UPTAKE_HOST_DEVICE
#endif

namespace gutibm {

enum class UptakeLimitMode : int {
  None = 0,
  Sherwood = 1,
  Voxel = 2,
};

namespace uptake {

constexpr double kPi = 3.14159265358979323846;

// Returns the per-step uptake ceiling in mol, or a negative value when the
// selected model imposes no ceiling.
GUTIBM_UPTAKE_HOST_DEVICE inline double allowed_uptake_mol(
    int mode, double concentration, double effective_diffusivity,
    double agent_radius, double cell_volume, double dt) {
  const double available = concentration > 0.0 ? concentration : 0.0;
  if (mode == static_cast<int>(UptakeLimitMode::Sherwood)) {
    if (effective_diffusivity <= 0.0 || agent_radius <= 0.0 || dt <= 0.0) {
      return 0.0;
    }
    return 4.0 * kPi * effective_diffusivity * agent_radius * available * dt;
  }
  if (mode == static_cast<int>(UptakeLimitMode::Voxel)) {
    if (cell_volume <= 0.0) return 0.0;
    return available * cell_volume;
  }
  return -1.0;
}

GUTIBM_UPTAKE_HOST_DEVICE inline double limit_fraction(double allowed,
                                                       double demanded) {
  if (allowed < 0.0) return 1.0;
  if (demanded <= 0.0) return 1.0;
  if (allowed <= 0.0) return 0.0;
  const double fraction = allowed / demanded;
  return fraction < 1.0 ? fraction : 1.0;
}

GUTIBM_UPTAKE_HOST_DEVICE inline double effective_diffusivity(
    double diff_coeff, double retardation) {
  return retardation > 0.0 ? diff_coeff / retardation : diff_coeff;
}

}  // namespace uptake
}  // namespace gutibm

#undef GUTIBM_UPTAKE_HOST_DEVICE

#endif  // GUTIBM_UPTAKE_LIMIT_H
