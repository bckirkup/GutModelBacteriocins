#include "qssa_gpu.h"
#include "species_names.h"
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

// Spec 6 — carbon/iron uptake is applied by the GPU metabolism kernel
// (mirroring FixMetabolism::grow_agent), and B12/corrinoid is no longer
// depleted. The chemistry-phase nutrient-depletion kernel that previously
// re-applied carbon/iron/B12 here duplicated that uptake (double-counting) and
// has been removed. Returning true signals "handled" so the CPU fallback is
// not run on the GPU path. (Agent O2 respiration is CPU-only; on the GPU path
// it is not applied, matching prior behavior.)
bool gpu_solve_nutrient_depletion(const AgentPoolGpu& agents, Int num_agents,
                                  ChemicalFieldGpu& chem_gpu, const ChemicalField& chem) {
  (void)agents;
  (void)num_agents;
  (void)chem_gpu;
  (void)chem;
#ifndef GUTIBM_CUDA
  return false;
#else
  return gpu_runtime_enabled();
#endif
}

}  // namespace gutibm
