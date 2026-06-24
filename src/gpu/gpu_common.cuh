#ifndef GUTIBM_GPU_COMMON_CUH
#define GUTIBM_GPU_COMMON_CUH

#include "gpu_types.h"

namespace gutibm {
namespace gpu {

static constexpr int N_IMAGES = 3;
static constexpr double PI_GPU = 3.14159265358979323846;

__device__ inline double clampd(double v, double lo, double hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

__device__ inline double peristaltic_factor(const AdvectionParams& adv, double x) {
  if (!adv.peristaltic_enabled) return 1.0;
  double phase = 2.0 * PI_GPU * adv.current_time / adv.peristaltic_period;
  if (adv.peristaltic_wavelength > 0.0) {
    phase -= 2.0 * PI_GPU * x / adv.peristaltic_wavelength;
  }
  return 1.0 + adv.peristaltic_amplitude * sin(phase);
}

__device__ inline void flow_velocity(const AdvectionParams& adv, double x, double z,
                                     double& vx, double& vy, double& vz) {
  if (adv.crypts_enabled && z < adv.lo_z + adv.crypt_depth) {
    vx = vy = vz = 0.0;
    return;
  }
  double zn = clampd((z - adv.lo_z) / adv.h, 0.0, 1.0);
  double pf = peristaltic_factor(adv, x);
  double vr = adv.v_radial_max * pow(zn, adv.profile_alpha);
  double vd = adv.v_distal_max * pow(zn, adv.profile_alpha);
  vx = vd * pf;
  vy = 0.0;
  vz = vr * pf;
}

__device__ inline double single_kernel(const double src[3], const double tgt[3],
                                       double D_eff, double Q,
                                       const double flow[3]) {
  double dx = tgt[0] - src[0];
  double dy = tgt[1] - src[1];
  double dz = tgt[2] - src[2];
  double r = sqrt(dx * dx + dy * dy + dz * dz);
  if (r < 1.0e-9) return 0.0;

  double U_mag = sqrt(flow[0] * flow[0] + flow[1] * flow[1] + flow[2] * flow[2]);
  double U_dot_r = flow[0] * dx + flow[1] * dy + flow[2] * dz;
  double prefactor = Q / (4.0 * PI_GPU * D_eff * r);
  double exponent = (U_dot_r - U_mag * r) / (2.0 * D_eff);
  exponent = exponent < -500.0 ? -500.0 : exponent;
  return prefactor * exp(exponent);
}

__device__ inline double concentration_bounded(const double src[3], const double tgt[3],
                                               const GfSourceParams& p,
                                               const DomainParams& dom,
                                               const AdvectionParams& adv) {
  double D_eff = p.diff_coeff / p.retardation;
  double flow[3];
  flow_velocity(adv, src[0], src[2], flow[0], flow[1], flow[2]);
  double Q = p.source_rate;
  double total = single_kernel(src, tgt, D_eff, Q, flow);

  for (int n = 1; n <= N_IMAGES; ++n) {
    double img[3];
    double dz_span = dom.z_hi - dom.z_lo;

    img[0] = src[0]; img[1] = src[1];
    img[2] = 2.0 * dom.z_lo - src[2] - 2.0 * n * dz_span;
    total += single_kernel(img, tgt, D_eff, Q, flow);

    img[2] = 2.0 * dom.z_hi - src[2] + 2.0 * n * dz_span;
    total += single_kernel(img, tgt, D_eff, Q, flow);

    img[2] = 2.0 * dom.z_lo - src[2] + 2.0 * n * dz_span;
    total += single_kernel(img, tgt, D_eff, Q, flow);

    img[2] = 2.0 * dom.z_hi - src[2] - 2.0 * n * dz_span;
    total += single_kernel(img, tgt, D_eff, Q, flow);
  }

  return total > 0.0 ? total : 0.0;
}

__device__ inline int cell_index(const DomainParams& dom, int ix, int iy, int iz) {
  return iz * (dom.nx * dom.ny) + iy * dom.nx + ix;
}

__device__ inline void cell_center(const DomainParams& dom, int ix, int iy, int iz,
                                   double out[3]) {
  out[0] = dom.lo[0] + (ix + 0.5) * dom.dx;
  out[1] = dom.lo[1] + (iy + 0.5) * dom.dx;
  out[2] = dom.lo[2] + (iz + 0.5) * dom.dx;
}

}  // namespace gpu
}  // namespace gutibm

#endif  // GUTIBM_GPU_COMMON_CUH
