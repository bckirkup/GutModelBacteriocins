#ifndef GUTIBM_FMM_GPU_H
#define GUTIBM_FMM_GPU_H

#include "types.h"
#include <vector>

namespace gutibm {

class Domain;
class FMM;

bool gpu_accumulate_far_field_local(const FMM& fmm,
                                    const Domain& domain,
                                    int expansion_order,
                                    const std::vector<Real>& near_conc,
                                    std::vector<Real>& toxin_conc);

bool gpu_accumulate_far_field_local_device(const FMM& fmm,
                                           const Domain& domain,
                                           int expansion_order,
                                           double* d_near_conc,
                                           double* d_out_conc);

}  // namespace gutibm

#endif  // GUTIBM_FMM_GPU_H
