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
class PinnedHostBuffer {
 public:
  PinnedHostBuffer() = default;
  explicit PinnedHostBuffer(size_t count) { allocate(count); }
  ~PinnedHostBuffer() { free(); }

  PinnedHostBuffer(const PinnedHostBuffer&) = delete;
  PinnedHostBuffer& operator=(const PinnedHostBuffer&) = delete;

  PinnedHostBuffer(PinnedHostBuffer&& other) noexcept { move_from(other); }
  PinnedHostBuffer& operator=(PinnedHostBuffer&& other) noexcept {
    if (this != &other) {
      free();
      move_from(other);
    }
    return *this;
  }

  void allocate(size_t count) {
    if (count == count_ && (count_ == 0 || data() != nullptr)) {
      return;
    }
    free();
    count_ = count;
    if (count_ == 0) return;
#ifdef GUTIBM_CUDA
    void* raw = nullptr;
    const cudaError_t status = cudaMallocHost(&raw, count_ * sizeof(T));
    if (status == cudaSuccess && raw != nullptr) {
      pinned_ptr_ = static_cast<T*>(raw);
      pinned_ = true;
      gpu_transfer_set_pinned_staging_active(true);
      return;
    }
    // Clear the failed cudaMallocHost error (e.g., RLIMIT_MEMLOCK / container
    // ulimit -l / no visible CUDA device) and fall back to pageable storage.
    (void)cudaGetLastError();
#endif
    fallback_.resize(count_);
    pinned_ = false;
  }

  void ensure_capacity(size_t count) {
    if (count > count_) {
      allocate(count);
    }
  }

  void free() noexcept {
#ifdef GUTIBM_CUDA
    if (pinned_ptr_ != nullptr) {
      cudaFreeHost(pinned_ptr_);
      pinned_ptr_ = nullptr;
    }
#endif
    fallback_.clear();
    fallback_.shrink_to_fit();
    count_ = 0;
    pinned_ = false;
  }

  T* data() noexcept {
    return pinned_ ? pinned_ptr_ : fallback_.data();
  }
  const T* data() const noexcept {
    return pinned_ ? pinned_ptr_ : fallback_.data();
  }
  size_t size() const noexcept { return count_; }
  bool empty() const noexcept { return count_ == 0; }
  bool is_pinned() const noexcept { return pinned_; }

 private:
  T* pinned_ptr_ = nullptr;
  std::vector<T> fallback_;
  size_t count_ = 0;
  bool pinned_ = false;

  void move_from(PinnedHostBuffer& other) noexcept {
    pinned_ptr_ = other.pinned_ptr_;
    fallback_ = std::move(other.fallback_);
    count_ = other.count_;
    pinned_ = other.pinned_;
    other.pinned_ptr_ = nullptr;
    other.count_ = 0;
    other.pinned_ = false;
  }
};

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
    if (ptr_) {
      cudaFree(ptr_);
      ptr_ = nullptr;
    }
#endif
    count_ = 0;
  }

  void zero() {
#ifdef GUTIBM_CUDA
    if (ptr_ != nullptr && count_ > 0) {
      if (const cudaError_t err = cudaMemset(ptr_, 0, count_ * sizeof(T));
          err != cudaSuccess) {
        throw Error(std::string("cudaMemset: ") + cudaGetErrorString(err));
      }
    }
#endif
  }

  void upload(const T* host, size_t count) {
    if (count > count_) allocate(count);
#ifdef GUTIBM_CUDA
    if (count > 0) {
      if (gpu_transfer_profiling_enabled()) {
        cudaStreamSynchronize(nullptr);
        gpu_check_error("cudaStreamSynchronize before H2D profiling");
      }
      const auto t0 = std::chrono::steady_clock::now();
      cudaStream_t stream = nullptr;
      cudaMemcpyAsync(ptr_, host, count * sizeof(T), cudaMemcpyHostToDevice, stream);
      cudaStreamSynchronize(stream);
      gpu_check_error("cudaMemcpy H2D");
      if (gpu_transfer_profiling_enabled()) {
        const auto t1 = std::chrono::steady_clock::now();
        gpu_transfer_record_h2d(
            std::chrono::duration<double>(t1 - t0).count(),
            static_cast<unsigned long long>(count * sizeof(T)));
      }
    }
#else
    (void)host;
#endif
  }

  void upload(const std::vector<T>& host) { upload(host.data(), host.size()); }

  void upload(const PinnedHostBuffer<T>& host, size_t count) {
    upload(host.data(), count);
  }

  void download(T* host, size_t count) const {
#ifdef GUTIBM_CUDA
    if (count > 0) {
      if (gpu_transfer_profiling_enabled()) {
        cudaStreamSynchronize(nullptr);
        gpu_check_error("cudaStreamSynchronize before D2H profiling");
      }
      const auto t0 = std::chrono::steady_clock::now();
      cudaStream_t stream = nullptr;
      cudaMemcpyAsync(host, ptr_, count * sizeof(T), cudaMemcpyDeviceToHost, stream);
      cudaStreamSynchronize(stream);
      gpu_check_error("cudaMemcpy D2H");
      if (gpu_transfer_profiling_enabled()) {
        const auto t1 = std::chrono::steady_clock::now();
        gpu_transfer_record_d2h(
            std::chrono::duration<double>(t1 - t0).count(),
            static_cast<unsigned long long>(count * sizeof(T)));
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

  void download(PinnedHostBuffer<T>& host, size_t count) const {
    host.ensure_capacity(count);
    download(host.data(), count);
  }

  T*       data() { return ptr_; }
  const T* data() const { return ptr_; }
  size_t   size() const { return count_; }

 private:
  T*     ptr_      = nullptr;
  size_t count_    = 0;

  void move_from(DeviceBuffer& other) noexcept {
    ptr_   = other.ptr_;
    count_ = other.count_;
    other.ptr_   = nullptr;
    other.count_ = 0;
  }
};

}  // namespace gutibm

#endif  // GUTIBM_DEVICE_MEMORY_H
