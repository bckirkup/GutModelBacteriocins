#include "qssa_gpu.h"
#include "agent_pool_gpu.h"
#include "chemical_field_gpu.h"
#include "chemical_field.h"
#include "dispatch.h"
#include "gpu_kernels.h"
#include "device_memory.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

namespace gutibm {

bool gpu_solve_nutrient_depletion(const AgentPoolGpu& agents, Int num_agents,
                                  ChemicalFieldGpu& chem_gpu, ChemicalField& chem) {
#ifndef GUTIBM_CUDA
  (void)agents;
  (void)num_agents;
  (void)chem_gpu;
  (void)chem;
  return false;
#else
  if (!gpu_runtime_enabled() || num_agents <= 0) return false;

  Int i_iron = chem.find("iron");
  Int i_b12 = chem.find("b12");
  Int i_carbon = chem.find("carbon");

  gpu::launch_nutrient_depletion_kernel(
      agents.grid_cell(), agents.mu_realized(), agents.biomass(), agents.state(),
      i_iron >= 0 ? chem_gpu.reac_device(i_iron) : nullptr,
      i_b12 >= 0 ? chem_gpu.reac_device(i_b12) : nullptr,
      i_carbon >= 0 ? chem_gpu.reac_device(i_carbon) : nullptr,
      num_agents,
      i_iron >= 0, i_b12 >= 0, i_carbon >= 0, nullptr);

  cudaDeviceSynchronize();
  gpu_check_error("nutrient_depletion_kernel");
  return true;
#endif
}

}  // namespace gutibm
