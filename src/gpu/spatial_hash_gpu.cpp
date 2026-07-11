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
#include <numeric>
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
  const Int num_cells = out.num_cells_x * out.num_cells_y * out.num_cells_z;

  out.cell_keys.allocate(static_cast<size_t>(num_agents));
  out.sorted_agent_indices.allocate(static_cast<size_t>(num_agents));

  gpu::launch_spatial_hash_build_kernel(
      agents.x(), agents.y(), agents.z(), agents.state(),
      out.cell_keys.data(), out.sorted_agent_indices.data(), num_agents,
      lo[0], lo[1], lo[2], cell_size,
      out.num_cells_x, out.num_cells_y, out.num_cells_z, gpu_compute_stream());
  gpu_sync_compute();
  gpu_check_error("spatial_hash_build_kernel");

  std::vector<int> keys(static_cast<size_t>(num_agents));
  out.cell_keys.download(keys);

  std::vector<int> counts(static_cast<size_t>(num_cells), 0);
  for (Int i = 0; i < num_agents; ++i) {
    const int key = keys[static_cast<size_t>(i)];
    if (key >= 0 && key < num_cells) {
      counts[static_cast<size_t>(key)]++;
    }
  }

  std::vector<int> offsets(static_cast<size_t>(num_cells + 1), 0);
  for (Int c = 0; c < num_cells; ++c) {
    offsets[static_cast<size_t>(c + 1)] =
        offsets[static_cast<size_t>(c)] + counts[static_cast<size_t>(c)];
  }

  std::vector<int> cursor = offsets;
  std::vector<int> sorted(static_cast<size_t>(num_agents), -1);
  for (Int i = 0; i < num_agents; ++i) {
    const int key = keys[static_cast<size_t>(i)];
    if (key < 0 || key >= num_cells) continue;
    const int pos = cursor[static_cast<size_t>(key)]++;
    sorted[static_cast<size_t>(pos)] = i;
  }

  out.cell_offsets.allocate(static_cast<size_t>(num_cells + 1));
  out.cell_counts.allocate(static_cast<size_t>(num_cells));
  out.cell_offsets.upload(offsets);
  out.cell_counts.upload(counts);
  out.sorted_agent_indices.upload(sorted);

  out.set_active(true);
  return true;
#endif
}

void SpatialHashGpu::build(const AgentPoolGpu& agents, Int num_agents) const {
  (void)agents;
  (void)num_agents;
}

}  // namespace gutibm
