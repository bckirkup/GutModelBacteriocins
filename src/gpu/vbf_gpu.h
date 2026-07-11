#ifndef GUTIBM_VBF_GPU_H
#define GUTIBM_VBF_GPU_H

#include "types.h"

namespace gutibm {

class ChemicalField;
class ChemicalFieldGpu;
class Domain;
class VBF;

struct AcetateConfig;
struct MucinConfig;
struct OxygenConfig;

bool gpu_apply_vbf_coupling(ChemicalFieldGpu& chem_gpu,
                            const ChemicalField& chem,
                            const Domain& domain,
                            const VBF& vbf,
                            const OxygenConfig& oxygen,
                            const AcetateConfig& acetate,
                            const MucinConfig& mucin);

}  // namespace gutibm

#endif  // GUTIBM_VBF_GPU_H
