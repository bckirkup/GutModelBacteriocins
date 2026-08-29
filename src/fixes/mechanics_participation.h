/* -----------------------------------------------------------------------
   GutIBM – Shared host/device mechanics participation predicate
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_MECHANICS_PARTICIPATION_H
#define GUTIBM_MECHANICS_PARTICIPATION_H

#ifdef __CUDACC__
#define GUTIBM_MECHANICS_HOST_DEVICE __host__ __device__
#else
#define GUTIBM_MECHANICS_HOST_DEVICE
#endif

namespace gutibm {

inline constexpr int kDeadStateValue = 3;

GUTIBM_MECHANICS_HOST_DEVICE inline bool mechanics_participates(
    int state, double death_time, double sim_time, int cdi_enabled,
    double corpse_persistence) {
  if (state != kDeadStateValue) return true;
  return cdi_enabled != 0 && death_time >= 0.0
      && (sim_time - death_time) < corpse_persistence;
}

}  // namespace gutibm

#undef GUTIBM_MECHANICS_HOST_DEVICE

#endif  // GUTIBM_MECHANICS_PARTICIPATION_H
