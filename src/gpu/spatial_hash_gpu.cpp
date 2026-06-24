#include "spatial_hash_gpu.h"
#include "agent_pool_gpu.h"
#include "dispatch.h"
#include "gpu_kernels.h"
#include "device_memory.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

#include <algorithm>
#include <cmath>
#include <vector>

namespace gutibm {

bool gpu_build_spatial_hash(const AgentPoolGpu& agents, Int num_agents,
                            Vec3 lo, Vec3 hi, Real cell_size,
                            SpatialHashGpu& out) {
#ifndef GUTIBM_CUDA
  (void)agents;
  (void)num_agents;
  (void)lo;
  (void)hi;
  (void)cell_size;
  (void)out;
  return false;
#else
  if (!gpu_runtime_enabled() || num_agents <= 0) return false;

  out.lo = lo;
  out.cell_size = cell_size;
  out.num_cells_x = static_cast<Int>(std::ceil((hi[0] - lo[0]) / cell_size));
  out.num_cells_y = static_cast<Int>(std::ceil((hi[1] - lo[1]) / cell_size));
  out.num_cells_z = static_cast<Int>(std::ceil((hi[2] - lo[2]) / cell_size));

  out.cell_keys.allocate(static_cast<size_t>(num_agents));
  out.sorted_agent_indices.allocate(static_cast<size_t>(num_agents));

  gpu::launch_spatial_hash_build_kernel(
      agents.x(), agents.y(), agents.z(), agents.state(),
      out.cell_keys.data(), out.sorted_agent_indices.data(), num_agents,
      lo[0], lo[1], lo[2], cell_size,
      out.num_cells_x, out.num_cells_y, out.num_cells_z, 0);
  cudaDeviceSynchronize();
  gpu_check_error("spatial_hash_build_kernel");

  std::vector<int> keys(static_cast<size_t>(num_agents));
  std::vector<int> indices(static_cast<size_t>(num_agents));
  out.cell_keys.download(keys);
  out.sorted_agent_indices.download(indices);

  std::vector<int> order(static_cast<size_t>(num_agents));
  for (Int i = 0; i < num_agents; ++i) order[static_cast<size_t>(i)] = i;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
    return keys[static_cast<size_t>(a)] < keys[static_cast<size_t>(b)];
  });
  for (Int i = 0; i < num_agents; ++i) {
    indices[static_cast<size_t>(i)] = order[static_cast<size_t>(i)];
  }
  out.sorted_agent_indices.upload(indices);

  out.set_active(true);
  return true;
#endif
}

void SpatialHashGpu::build(const AgentPoolGpu& agents, Int num_agents) {
  (void)agents;
  (void)num_agents;
}

}  // namespace gutibm
