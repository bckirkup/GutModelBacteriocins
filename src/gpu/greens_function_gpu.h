#ifndef GUTIBM_GREENS_FUNCTION_GPU_H
#define GUTIBM_GREENS_FUNCTION_GPU_H

#include "types.h"
#include "greens_function.h"
#include <vector>

namespace gutibm {

class Domain;
class AdvectionField;

// GPU-accelerated Green's function superposition. Returns false if GPU unavailable.
bool gpu_superpose_to_grid(
    const Domain& domain,
    const AdvectionField& adv,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    std::vector<Real>& grid_conc,
    Real cutoff_radius);

// Deposit near-field superposition directly into a device concentration buffer.
bool gpu_superpose_to_device(
    const Domain& domain,
    const AdvectionField& adv,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    double* d_grid_conc,
    Real cutoff_radius);

}  // namespace gutibm

#endif  // GUTIBM_GREENS_FUNCTION_GPU_H
