#include "device.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

namespace gutibm {

namespace {
thread_local std::string g_last_error;
}  // namespace

std::string DeviceContext::last_error() { return g_last_error; }

int DeviceContext::device_count() {
#ifdef GUTIBM_CUDA
  int count = 0;
  if (cudaGetDeviceCount(&count) != cudaSuccess) {
    return 0;
  }
  return count;
#else
  return 0;
#endif
}

bool DeviceContext::init(int device_id) {
#ifdef GUTIBM_CUDA
  int count = device_count();
  if (count <= 0) {
    g_last_error = "No CUDA devices found";
    active_ = false;
    return false;
  }
  if (device_id < 0 || device_id >= count) {
    g_last_error = "Invalid CUDA device id";
    active_ = false;
    return false;
  }
  if (cudaError_t err = cudaSetDevice(device_id); err != cudaSuccess) {
    g_last_error = cudaGetErrorString(err);
    active_ = false;
    return false;
  }
  device_id_ = device_id;
  active_ = true;
  g_last_error.clear();
  return true;
#else
  (void)device_id;
  g_last_error = "CUDA support not compiled in";
  active_ = false;
  return false;
#endif
}

}  // namespace gutibm
