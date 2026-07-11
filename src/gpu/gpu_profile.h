#ifndef GUTIBM_GPU_PROFILE_H
#define GUTIBM_GPU_PROFILE_H

namespace gutibm {

struct GpuTransferProfile {
  double h2d_s = 0.0;
  double d2h_s = 0.0;
};

bool gpu_transfer_profiling_enabled();
void gpu_transfer_profile_set_enabled(bool enabled);
void gpu_transfer_profile_reset();
void gpu_transfer_record_h2d(double seconds);
void gpu_transfer_record_d2h(double seconds);
GpuTransferProfile gpu_transfer_profile_snapshot();

}  // namespace gutibm

#endif  // GUTIBM_GPU_PROFILE_H
