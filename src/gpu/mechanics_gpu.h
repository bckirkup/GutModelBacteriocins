#ifndef GUTIBM_MECHANICS_GPU_H
#define GUTIBM_MECHANICS_GPU_H

#include "types.h"
#include "spatial_hash_gpu.h"

namespace gutibm {

class AgentPoolGpu;
class Domain;
struct MechanicsConfig;

bool gpu_run_mechanics(AgentPoolGpu& agents, Int num_agents,
                       const SpatialHashGpu& hash, const Domain& domain,
                       const MechanicsConfig& cfg, Real dt);

}  // namespace gutibm

#endif  // GUTIBM_MECHANICS_GPU_H
