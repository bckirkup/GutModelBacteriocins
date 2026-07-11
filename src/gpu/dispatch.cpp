#include "dispatch.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

namespace gutibm {

namespace {
struct GpuDispatchState {
  GpuConfig cfg;
  DeviceContext device;
  bool initialized = false;
#ifdef GUTIBM_CUDA
  cudaStream_t compute_stream = nullptr;
#endif
};

GpuDispatchState& gpu_dispatch_state() {
  static GpuDispatchState state;
  return state;
}

#ifdef GUTIBM_CUDA
void destroy_compute_stream(GpuDispatchState& state) {
  if (state.compute_stream != nullptr) {
    cudaStreamDestroy(state.compute_stream);
    state.compute_stream = nullptr;
  }
}
#endif
}  // namespace

void gpu_set_config(const GpuConfig& cfg) { gpu_dispatch_state().cfg = cfg; }

const GpuConfig& gpu_config() { return gpu_dispatch_state().cfg; }

const DeviceContext& gpu_device() { return gpu_dispatch_state().device; }

bool gpu_init_for_rank(int mpi_rank, int mpi_nprocs) {
  auto& state = gpu_dispatch_state();
  state.initialized = false;
#ifdef GUTIBM_CUDA
  destroy_compute_stream(state);
#endif
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

#ifdef GUTIBM_CUDA
  if (cudaStreamCreate(&state.compute_stream) != cudaSuccess) {
    state.device = DeviceContext{};
    return false;
  }
#endif

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

#ifdef GUTIBM_CUDA
cudaStream_t gpu_compute_stream() {
  return gpu_dispatch_state().compute_stream;
}
#endif

void gpu_sync_compute() {
#ifdef GUTIBM_CUDA
  cudaStream_t stream = gpu_dispatch_state().compute_stream;
  if (stream != nullptr) {
    cudaStreamSynchronize(stream);
  }
#endif
}

}  // namespace gutibm
