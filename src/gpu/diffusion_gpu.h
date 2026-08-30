#ifndef GUTIBM_DIFFUSION_GPU_H
#define GUTIBM_DIFFUSION_GPU_H

#include "chemical_field.h"
#include "gpu_kernels.h"
#include "types.h"

namespace gutibm {

class Domain;

enum class DiffusionPhase {
  Replicated,
  SlabPreX,
  SlabPostX,
};

struct GpuSlabDiffusionContext {
  ChemicalField& field;
  Int spec_index;
  Int storage_nx;
  Int owned_storage_x_begin;
  Int owned_storage_x_end;
};

// Host-compilable checks shared by the CPU and GPU diffusion dispatch paths.
int diffusion_z_line_length(
    const Domain& domain, EpithelialBoundaryMode mode);
bool diffusion_line_lengths_within(
    const Domain& domain, EpithelialBoundaryMode mode, int max_line);
bool diffusion_all_species_within(
    const Domain& domain, const ChemicalField& field, int max_line);

// Route B is deliberately narrower than ordinary device diffusion.  Delivery
// uses a line-local Sherman--Morrison correction and therefore stays on the
// device only for replicated, single-rank fields whose lines fit 512 entries.
bool delivery_route_b_eligible(
    const Domain& domain, const ChemicalField& field,
    int max_line = gpu::kMaxDeliveryLineLength);
bool gpu_delivery_species_eligible(
    const Domain& domain, const ChemicalSpec& spec, Real dt);

// True when the selected z solve fits the PCR shared-memory line limit.
bool gpu_diffusion_line_lengths_supported(
    const Domain& domain,
    EpithelialBoundaryMode mode = EpithelialBoundaryMode::Dirichlet);

// Apply one backward-Euler directional-splitting diffusion step for a single
// species directly on a device concentration buffer. Mirrors
// ChemicalField::apply_diffusion for one spec. Returns false when CUDA is
// unavailable, the species does not diffuse, or line lengths exceed the limit.
bool gpu_apply_species_diffusion_device(const Domain& domain,
                                        const ChemicalSpec& spec,
                                        double* d_conc,
                                        double* d_injected_amount,
                                        Real dt);

bool gpu_apply_species_delivery_device(
    const Domain& domain, const ChemicalSpec& spec, double* d_conc,
    const double* d_sink, const double* d_prescribed, double* d_realized,
    double* d_boundary_injected, double* d_gradient_source,
    double* d_reaction_clip, Real dt, bool prescribed_active);

bool gpu_apply_species_diffusion_slab_device(
    const Domain& domain, const ChemicalSpec& spec, double* d_conc,
    double* d_injected_amount, Real dt, GpuSlabDiffusionContext& context);

// Host-buffer variant: upload, diffuse on device, download.
bool gpu_apply_species_diffusion(const Domain& domain,
                                 const ChemicalSpec& spec,
                                 std::vector<Real>& concentration,
                                 Real dt);

}  // namespace gutibm

#endif  // GUTIBM_DIFFUSION_GPU_H
