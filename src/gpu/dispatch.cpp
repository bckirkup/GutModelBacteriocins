#include "dispatch.h"

namespace gutibm {

namespace {
struct GpuDispatchState {
  GpuConfig cfg;
  DeviceContext device;
  bool initialized = false;
};

GpuDispatchState& gpu_dispatch_state() {
  static GpuDispatchState state;
  return state;
}
}  // namespace

void gpu_set_config(const GpuConfig& cfg) { gpu_dispatch_state().cfg = cfg; }

const GpuConfig& gpu_config() { return gpu_dispatch_state().cfg; }

const DeviceContext& gpu_device() { return gpu_dispatch_state().device; }

bool gpu_init_for_rank(int mpi_rank, int mpi_nprocs) {
  auto& state = gpu_dispatch_state();
  state.initialized = false;
  if (!state.cfg.enabled) {
    return false;
  }

  int dev_id = state.cfg.device_id;
  if (dev_id < 0) {
    int count = DeviceContext::device_count();
    if (count <= 0) {
      return false;
    }
    dev_id = mpi_rank % count;
  }

  if (!state.device.init(dev_id)) {
    return false;
  }

  (void)mpi_nprocs;
  state.initialized = true;
  return true;
}

bool gpu_runtime_enabled() {
#ifdef GUTIBM_CUDA
  const auto& state = gpu_dispatch_state();
  return state.cfg.enabled && state.initialized && state.device.active();
#else
  return false;
#endif
}

}  // namespace gutibm
