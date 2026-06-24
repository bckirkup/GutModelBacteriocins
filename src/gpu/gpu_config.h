#ifndef GUTIBM_GPU_CONFIG_H
#define GUTIBM_GPU_CONFIG_H

namespace gutibm {

struct GpuConfig {
  bool enabled     = false;  // runtime toggle (requires CUDA build)
  int  device_id   = -1;     // -1 = auto-select from MPI rank
};

}  // namespace gutibm

#endif  // GUTIBM_GPU_CONFIG_H
