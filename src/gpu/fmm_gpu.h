#ifndef GUTIBM_FMM_GPU_H
#define GUTIBM_FMM_GPU_H

#include "types.h"
#include "greens_function.h"

#include "fmm.h"

namespace gutibm {

class AdvectionField;
class ChemicalFieldGpu;
class Domain;

// Device-resident GPU FMM far-field accumulation directly into chem_gpu's
// concentration buffer for toxin_species_idx. Supports both locals_ready()==true
// (Mode A: device octree lookup + local expansion + residual U-list) and
// locals_ready()==false (Mode B: stack-based GPU octree multipole traversal),
// across both replicated and MPI slab chemical layouts. Uploads only the
// compact O(N_nodes + N_sources) octree; performs zero full-grid H2D/D2H transfers.
bool gpu_accumulate_far_field_device(
    const FMM& fmm,
    const Domain& domain,
    const AdvectionField& adv,
    const GreensFunctionParams& avg_params,
    Real theta,
    Real near_cutoff,
    ChemicalFieldGpu* chem_gpu,
    Int toxin_species_idx,
    FmmNearCutoffShape cutoff_shape = FmmNearCutoffShape::Sphere);

struct FmmDeviceLayout {
  Int storage_nx = 0;
  Int owned_storage_x_begin = 0;
  Int owned_storage_x_end = 0;
  Int local_global_x_begin = 0;
};

// Low-level variant operating on a caller-owned device concentration buffer
// and explicit replicated or slab storage layout.
bool gpu_accumulate_far_field_device(
    const FMM& fmm,
    const Domain& domain,
    const AdvectionField& adv,
    const GreensFunctionParams& avg_params,
    Real theta,
    Real near_cutoff,
    double* d_conc,
    const FmmDeviceLayout& layout,
    FmmNearCutoffShape cutoff_shape = FmmNearCutoffShape::Sphere);

}  // namespace gutibm

#endif  // GUTIBM_FMM_GPU_H
