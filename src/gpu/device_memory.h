#ifndef GUTIBM_DEVICE_MEMORY_H
#define GUTIBM_DEVICE_MEMORY_H

#include <chrono>
#include <cstddef>
#include "error.h"
#include "gpu_profile.h"
#include <vector>

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

namespace gutibm {

inline void gpu_check_error(const char* what) {
#ifdef GUTIBM_CUDA
  if (cudaError_t err = cudaGetLastError(); err != cudaSuccess) {
    throw Error(std::string(what) + ": " + cudaGetErrorString(err));
  }
#endif
  (void)what;
}

template <typename T>
class DeviceBuffer {
 public:
  DeviceBuffer() = default;
  explicit DeviceBuffer(size_t count) { allocate(count); }
  ~DeviceBuffer() { free(); }

  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;

  DeviceBuffer(DeviceBuffer&& other) noexcept { move_from(other); }
  DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
    if (this != &other) {
      free();
      move_from(other);
    }
    return *this;
  }

  void allocate(size_t count) {
    free();
    count_ = count;
    if (count_ == 0) return;
#ifdef GUTIBM_CUDA
    cudaMalloc(&ptr_, count_ * sizeof(T));
    gpu_check_error("cudaMalloc");
#else
    ptr_ = nullptr;
#endif
  }

  void free() {
#ifdef GUTIBM_CUDA
    if (host_ptr_) {
      cudaFreeHost(host_ptr_);
      host_ptr_ = nullptr;
    }
    if (ptr_) {
      cudaFree(ptr_);
      ptr_ = nullptr;
    }
#endif
    count_ = 0;
  }

  void upload(const T* host, size_t count) {
    if (count > count_) allocate(count);
#ifdef GUTIBM_CUDA
    if (count > 0) {
      const auto t0 = std::chrono::steady_clock::now();
      cudaStream_t stream = nullptr;
      cudaMemcpyAsync(ptr_, host, count * sizeof(T), cudaMemcpyHostToDevice, stream);
      cudaStreamSynchronize(stream);
      gpu_check_error("cudaMemcpy H2D");
      if (gpu_transfer_profiling_enabled()) {
        const auto t1 = std::chrono::steady_clock::now();
        gpu_transfer_record_h2d(
            std::chrono::duration<double>(t1 - t0).count());
      }
    }
#else
    (void)host;
#endif
  }

  void upload(const std::vector<T>& host) { upload(host.data(), host.size()); }

  void download(T* host, size_t count) const {
#ifdef GUTIBM_CUDA
    if (count > 0) {
      const auto t0 = std::chrono::steady_clock::now();
      cudaStream_t stream = nullptr;
      cudaMemcpyAsync(host, ptr_, count * sizeof(T), cudaMemcpyDeviceToHost, stream);
      cudaStreamSynchronize(stream);
      gpu_check_error("cudaMemcpy D2H");
      if (gpu_transfer_profiling_enabled()) {
        const auto t1 = std::chrono::steady_clock::now();
        gpu_transfer_record_d2h(
            std::chrono::duration<double>(t1 - t0).count());
      }
    }
#else
    (void)host;
    (void)count;
#endif
  }

  void download(std::vector<T>& host) const {
    host.resize(count_);
    download(host.data(), count_);
  }

  T*       data() { return ptr_; }
  const T* data() const { return ptr_; }
  size_t   size() const { return count_; }

 private:
  T*     ptr_      = nullptr;
  T*     host_ptr_ = nullptr;
  size_t count_    = 0;

  void move_from(DeviceBuffer& other) noexcept {
    ptr_      = other.ptr_;
    host_ptr_ = other.host_ptr_;
    count_    = other.count_;
    other.ptr_      = nullptr;
    other.host_ptr_ = nullptr;
    other.count_    = 0;
  }
};

}  // namespace gutibm

#endif  // GUTIBM_DEVICE_MEMORY_H
