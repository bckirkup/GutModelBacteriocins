#ifndef GUTIBM_RECEPTOR_GPU_H
#define GUTIBM_RECEPTOR_GPU_H

#include "types.h"
#include <vector>

namespace gutibm {

class AgentPool;
class AgentPoolGpu;
class ChemicalFieldGpu;
class ChemicalField;
struct ReceptorConfig;

bool gpu_compute_receptor_kill_probs_host_packed(
    const AgentPoolGpu& agents,
    const AgentPool& pool,
    const ChemicalFieldGpu& chem_gpu,
    const ChemicalField& chem,
    const ReceptorConfig& cfg,
    double dt,
    std::vector<double>& kill_probs_out);

}  // namespace gutibm

#endif  // GUTIBM_RECEPTOR_GPU_H
