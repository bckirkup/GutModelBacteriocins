#ifndef GUTIBM_DIFFUSION_GPU_H
#define GUTIBM_DIFFUSION_GPU_H

#include "types.h"

namespace gutibm {

class ChemicalSpec;
class ChemicalField;
class Domain;

enum class DiffusionPhase {
  Replicated,
  SlabPreX,
  SlabPostX,
};

struct SlabDiffusionContext {
  ChemicalField& field;
  Int spec_index;
  Int storage_nx;
  Int owned_storage_x_begin;
  Int owned_storage_x_end;
};

// True when nx, ny, and nz-1 fit the PCR shared-memory line limit.
bool gpu_diffusion_line_lengths_supported(const Domain& domain);

// Apply one backward-Euler directional-splitting diffusion step for a single
// species directly on a device concentration buffer. Mirrors
// ChemicalField::apply_diffusion for one spec. Returns false when CUDA is
// unavailable, the species does not diffuse, or line lengths exceed the limit.
bool gpu_apply_species_diffusion_device(const Domain& domain,
                                        const ChemicalSpec& spec,
                                        double* d_conc,
                                        double* d_injected_amount,
                                        Real dt);

bool gpu_apply_species_diffusion_slab_device(
    const Domain& domain, const ChemicalSpec& spec, double* d_conc,
    double* d_injected_amount, Real dt, SlabDiffusionContext& context);

// Host-buffer variant: upload, diffuse on device, download.
bool gpu_apply_species_diffusion(const Domain& domain,
                                 const ChemicalSpec& spec,
                                 std::vector<Real>& concentration,
                                 Real dt);

}  // namespace gutibm

#endif  // GUTIBM_DIFFUSION_GPU_H
