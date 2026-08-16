#include "mechanics_gpu.h"
#include "agent_pool_gpu.h"
#include "domain.h"
#include "fix_mechanics.h"
#include "dispatch.h"
#include "device_memory.h"
#include "gpu_kernels.h"
#include "gpu_types.h"

#include <type_traits>

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

namespace gutibm {

static_assert(std::is_same_v<Int, int>,
              "GPU clamp counts must match the host Int type");

namespace {

#ifdef GUTIBM_CUDA
struct MechanicsDeviceScratch {
  DeviceBuffer<double> dx;
  DeviceBuffer<double> dy;
  DeviceBuffer<double> dz;
  DeviceBuffer<int> clamp_count;
};

MechanicsDeviceScratch& mechanics_device_scratch() {
  static MechanicsDeviceScratch scratch;
  return scratch;
}
#endif

}  // namespace

bool gpu_run_mechanics(AgentPoolGpu& agents, Int num_agents,
                       const SpatialHashGpu& hash, const Domain& domain,
                       const MechanicsConfig& cfg, Real dt, Real viscosity,
                       Int& clamp_count) {
  clamp_count = 0;
#ifndef GUTIBM_CUDA
  (void)agents;
  (void)num_agents;
  (void)hash;
  (void)domain;
  (void)cfg;
  (void)dt;
  (void)viscosity;
  (void)clamp_count;
  return false;
#else
  if (!gpu_runtime_enabled() || num_agents <= 0 || !hash.active()) return false;
  if (hash.cell_offsets.size() == 0 || hash.sorted_agent_indices.size() == 0) {
    return false;
  }

  auto& scratch = mechanics_device_scratch();
  scratch.dx.allocate(static_cast<size_t>(num_agents));
  scratch.dy.allocate(static_cast<size_t>(num_agents));
  scratch.dz.allocate(static_cast<size_t>(num_agents));
  scratch.clamp_count.allocate(1);

  gpu::MechanicsLaunchParams params{};
  params.hertz_k = cfg.hertz_k;
  params.hertzian_enabled = cfg.hertzian_enabled ? 1 : 0;
  params.adhesion_enabled = cfg.adhesion_enabled ? 1 : 0;
  params.adhesion_strength = cfg.adhesion_strength;
  params.adhesion_range = cfg.adhesion_range;
  params.viscosity = viscosity;
  params.max_displacement_fraction =
      kMechanicsMaxDisplacementRadiusFraction;
  params.dt = dt;
  params.lo0 = domain.lo()[0];
  params.lo1 = domain.lo()[1];
  params.lo2 = domain.lo()[2];
  params.hi0 = domain.hi()[0];
  params.hi1 = domain.hi()[1];
  params.hi2 = domain.hi()[2];
  params.periodic_x = domain.config().periodic[0] ? 1 : 0;
  params.periodic_y = domain.config().periodic[1] ? 1 : 0;
  params.periodic_z = domain.config().periodic[2] ? 1 : 0;
  params.cell_size = hash.cell_size;
  params.nx_cells = hash.num_cells_x;
  params.ny_cells = hash.num_cells_y;
  params.nz_cells = hash.num_cells_z;

  cudaStream_t stream = gpu_compute_stream();
  gpu::launch_mechanics_clear_kernel(
      scratch.dx.data(), scratch.dy.data(), scratch.dz.data(),
      scratch.clamp_count.data(), num_agents, stream);
  gpu::launch_mechanics_forces_kernel(
      agents.x(), agents.y(), agents.z(),
      agents.radius(), agents.state(),
      hash.cell_offsets.data(), hash.sorted_agent_indices.data(),
      scratch.dx.data(), scratch.dy.data(), scratch.dz.data(),
      num_agents, params, stream);
  gpu::launch_mechanics_apply_kernel(
      agents.x(), agents.y(), agents.z(),
      agents.radius(),
      scratch.dx.data(), scratch.dy.data(), scratch.dz.data(),
      scratch.clamp_count.data(), num_agents, params, stream);

  gpu_sync_compute();
  gpu_check_error("mechanics_kernel");
  scratch.clamp_count.download(&clamp_count, 1);
  return true;
#endif
}

}  // namespace gutibm
