#ifndef GUTIBM_DIFFUSION_GPU_H
#define GUTIBM_DIFFUSION_GPU_H

#include "types.h"

namespace gutibm {

class ChemicalField;
class ChemicalSpec;
class Domain;

// Apply one backward-Euler directional-splitting diffusion step for a single
// species on the GPU. Mirrors ChemicalField::apply_diffusion for one spec.
// Returns false when CUDA is unavailable or line lengths exceed the limit.
bool gpu_apply_species_diffusion(const Domain& domain,
                                 const ChemicalSpec& spec,
                                 std::vector<Real>& concentration,
                                 Real dt);

}  // namespace gutibm

#endif  // GUTIBM_DIFFUSION_GPU_H
