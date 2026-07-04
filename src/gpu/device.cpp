#include "device.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

namespace gutibm {

namespace {
std::string& last_error_storage() {
  thread_local std::string storage;
  return storage;
}
}  // namespace

std::string DeviceContext::last_error() { return last_error_storage(); }

int DeviceContext::device_count() {
#ifdef GUTIBM_CUDA
  int count = 0;
  if (cudaError_t err = cudaGetDeviceCount(&count); err != cudaSuccess) {
    last_error_storage() = cudaGetErrorString(err);
    return 0;
  }
  last_error_storage().clear();
  return count;
#else
  return 0;
#endif
}

bool DeviceContext::init(int device_id) {
#ifdef GUTIBM_CUDA
  int count = device_count();
  if (count <= 0) {
    if (last_error_storage().empty()) {
      last_error_storage() = "No CUDA devices found";
    }
    active_ = false;
    return false;
  }
  if (device_id < 0 || device_id >= count) {
    last_error_storage() = "Invalid CUDA device id";
    active_ = false;
    return false;
  }
  if (cudaError_t err = cudaSetDevice(device_id); err != cudaSuccess) {
    last_error_storage() = cudaGetErrorString(err);
    active_ = false;
    return false;
  }
  device_id_ = device_id;
  active_ = true;
  last_error_storage().clear();
  return true;
#else
  (void)device_id;
  last_error_storage() = "CUDA support not compiled in";
  active_ = false;
  return false;
#endif
}

}  // namespace gutibm
