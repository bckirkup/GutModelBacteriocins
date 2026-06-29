#ifndef GUTIBM_GPU_TYPES_H
#define GUTIBM_GPU_TYPES_H

namespace gutibm {
namespace gpu {

struct DomainParams {
  int nx;
  int ny;
  int nz;
  double dx;
  double lo[3];       // NOLINT(modernize-avoid-c-arrays) — passed to CUDA kernels
  bool periodic[3];   // NOLINT(modernize-avoid-c-arrays) — passed to CUDA kernels
  double z_lo;
  double z_hi;
};

struct AdvectionParams {
  double v_radial_max;
  double v_distal_max;
  double h;
  double lo_z;
  double profile_alpha;
  bool crypts_enabled;
  double crypt_depth;
  bool peristaltic_enabled;
  double peristaltic_period;
  double peristaltic_amplitude;
  double peristaltic_wavelength;
  double current_time;
};

struct GfSourceParams {
  double diff_coeff;
  double source_rate;
  double retardation;
};

}  // namespace gpu
}  // namespace gutibm

#endif  // GUTIBM_GPU_TYPES_H
