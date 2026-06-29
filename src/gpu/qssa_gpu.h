#ifndef GUTIBM_QSSA_GPU_H
#define GUTIBM_QSSA_GPU_H

#include "types.h"

namespace gutibm {

class AgentPoolGpu;
class ChemicalFieldGpu;
class ChemicalField;

bool gpu_solve_nutrient_depletion(const AgentPoolGpu& agents, Int num_agents,
                                  ChemicalFieldGpu& chem_gpu, const ChemicalField& chem);

}  // namespace gutibm

#endif  // GUTIBM_QSSA_GPU_H
