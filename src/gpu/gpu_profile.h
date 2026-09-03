#ifndef GUTIBM_GPU_PROFILE_H
#define GUTIBM_GPU_PROFILE_H

#include <string>
#include <vector>

namespace gutibm {

struct GpuTransferProfile {
  double h2d_s = 0.0;
  double d2h_s = 0.0;
  unsigned long long h2d_bytes = 0;
  unsigned long long d2h_bytes = 0;
  unsigned long long h2d_calls = 0;
  unsigned long long d2h_calls = 0;
  double slab_x_roundtrip_s = 0.0;
};

struct GpuTransferSiteProfile {
  std::string label;
  double h2d_s = 0.0;
  double d2h_s = 0.0;
  unsigned long long h2d_bytes = 0;
  unsigned long long d2h_bytes = 0;
  unsigned long long h2d_calls = 0;
  unsigned long long d2h_calls = 0;
};

class GpuTransferSite {
 public:
  explicit GpuTransferSite(const char* label);
  ~GpuTransferSite();
  GpuTransferSite(const GpuTransferSite&) = delete;
  GpuTransferSite& operator=(const GpuTransferSite&) = delete;
  GpuTransferSite(GpuTransferSite&&) = delete;
  GpuTransferSite& operator=(GpuTransferSite&&) = delete;

 private:
  const char* previous_;
};

bool gpu_transfer_profiling_enabled();
void gpu_transfer_profile_set_enabled(bool enabled);
void gpu_transfer_profile_reset();
void gpu_transfer_record_h2d(double seconds, unsigned long long bytes);
void gpu_transfer_record_d2h(double seconds, unsigned long long bytes);
GpuTransferProfile gpu_transfer_profile_snapshot();
std::vector<GpuTransferSiteProfile> gpu_transfer_site_profiles();
void gpu_transfer_record_slab_x_roundtrip(double seconds);

}  // namespace gutibm

#endif  // GUTIBM_GPU_PROFILE_H
