#include "qssa_gpu.h"
#include "species_names.h"
#include "agent_pool_gpu.h"
#include "chemical_field_gpu.h"
#include "chemical_field.h"
#include "chem_environment_config.h"
#include "dispatch.h"
#include "domain.h"
#include "gpu_kernels.h"
#include "device_memory.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

namespace gutibm {

bool gpu_solve_nutrient_depletion(const AgentPoolGpu& agents,
                                  Int num_agents,
                                  ChemicalFieldGpu& chem_gpu,
                                  const ChemicalField& chem,
                                  const OxygenConfig& oxygen,
                                  const Domain& domain) {
#ifndef GUTIBM_CUDA
  (void)agents;
  (void)num_agents;
  (void)chem_gpu;
  (void)chem;
  (void)oxygen;
  (void)domain;
  return false;
#else
  if (!gpu_runtime_enabled() || !chem_gpu.active()) return false;
  if (!oxygen.enabled) return true;

  const Int i_oxygen = chem.find(species::OXYGEN);
  if (i_oxygen < 0) return true;

  const double cell_vol = domain.dx() * domain.dx() * domain.dx();
  if (cell_vol <= 0.0 || num_agents <= 0) return true;

  double* d_reac_oxygen = chem_gpu.reac_device(i_oxygen);
  if (d_reac_oxygen == nullptr) return false;

  gpu::launch_o2_depletion_kernel(
      d_reac_oxygen,
      agents.mu_realized(),
      agents.grid_cell(),
      agents.state(),
      num_agents,
      oxygen.q_consumption,
      oxygen.q_maintenance,
      1.0 / cell_vol,
      nullptr);

  cudaDeviceSynchronize();
  gpu_check_error("gpu_solve_nutrient_depletion");
  return true;
#endif
}

}  // namespace gutibm
