#ifndef GUTIBM_VBF_GPU_H
#define GUTIBM_VBF_GPU_H

#include "types.h"

namespace gutibm {

class ChemicalField;
class ChemicalFieldGpu;
class Domain;
class VBF;
class AgentPoolGpu;

struct AcetateConfig;
struct MucinConfig;
struct OxygenConfig;
struct VbfFluxTotals;

bool gpu_apply_vbf_coupling(ChemicalFieldGpu& chem_gpu,
                            ChemicalField& chem,
                            const Domain& domain,
                            const VBF& vbf,
                            const OxygenConfig& oxygen,
                            const AcetateConfig& acetate,
                            const MucinConfig& mucin,
                            const AgentPoolGpu& agents,
                            VbfFluxTotals& totals, Real dt);

}  // namespace gutibm

#endif  // GUTIBM_VBF_GPU_H
