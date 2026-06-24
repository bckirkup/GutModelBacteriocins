#include "dispatch.h"

namespace gutibm {

namespace {
GpuConfig g_cfg;
DeviceContext g_device;
bool g_initialized = false;
}  // namespace

void gpu_set_config(const GpuConfig& cfg) { g_cfg = cfg; }

const GpuConfig& gpu_config() { return g_cfg; }

const DeviceContext& gpu_device() { return g_device; }

bool gpu_init_for_rank(int mpi_rank, int mpi_nprocs) {
  g_initialized = false;
  if (!g_cfg.enabled) {
    return false;
  }

  int dev_id = g_cfg.device_id;
  if (dev_id < 0) {
    int count = DeviceContext::device_count();
    if (count <= 0) {
      return false;
    }
    dev_id = mpi_rank % count;
  }

  if (!g_device.init(dev_id)) {
    return false;
  }

  (void)mpi_nprocs;
  g_initialized = true;
  return true;
}

bool gpu_runtime_enabled() {
#ifdef GUTIBM_CUDA
  return g_cfg.enabled && g_initialized && g_device.active();
#else
  return false;
#endif
}

}  // namespace gutibm
