#ifndef GUTIBM_GPU_TEST_SUPPORT_H
#define GUTIBM_GPU_TEST_SUPPORT_H

#include "device.h"

#include <cstdlib>
#include <iostream>
#include <string_view>

namespace gutibm::test {

inline bool gpu_required() {
  const char* value = std::getenv("REQUIRE_GPU");
  return value != nullptr && std::string_view(value) == "1";
}

inline int require_gpu(std::string_view test_name) {
#ifndef GUTIBM_CUDA
  if (gpu_required()) {
    std::cerr << test_name
              << ": FAILED (REQUIRE_GPU=1 but CUDA was not compiled in)\n";
    return 1;
  }
  std::cout << test_name << ": SKIPPED (CUDA not compiled in)\n";
  return 77;
#else
  if (DeviceContext::device_count() <= 0) {
    if (gpu_required()) {
      std::cerr << test_name
                << ": FAILED (REQUIRE_GPU=1 but no CUDA device is available)\n";
      return 1;
    }
    std::cout << test_name << ": SKIPPED (no CUDA device)\n";
    return 77;
  }
  return 0;
#endif
}

}  // namespace gutibm::test

#endif  // GUTIBM_GPU_TEST_SUPPORT_H
