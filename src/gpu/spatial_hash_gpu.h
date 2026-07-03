#ifndef GUTIBM_SPATIAL_HASH_GPU_H
#define GUTIBM_SPATIAL_HASH_GPU_H

#include "types.h"
#include "device_memory.h"
#include <vector>

namespace gutibm {

class AgentPoolGpu;

struct SpatialHashGpu {
  bool active_ = false;

  Int num_cells_x = 0;
  Int num_cells_y = 0;
  Int num_cells_z = 0;
  Real cell_size = 10.0e-6;
  Vec3 lo{};

  DeviceBuffer<int> cell_offsets;
  DeviceBuffer<int> cell_counts;
  DeviceBuffer<int> sorted_agent_indices;
  DeviceBuffer<int> cell_keys;

  void build(const AgentPoolGpu& agents, Int num_agents) const;
  void set_active(bool v) { active_ = v; }
  bool active() const { return active_; }
};

bool gpu_build_spatial_hash(const AgentPoolGpu& agents, Int num_agents,
                            Vec3 lo, Vec3 hi, Real cell_size,
                            SpatialHashGpu& out);

}  // namespace gutibm

#endif  // GUTIBM_SPATIAL_HASH_GPU_H
