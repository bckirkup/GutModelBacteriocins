#ifndef GUTIBM_GPU_KERNELS_H
#define GUTIBM_GPU_KERNELS_H

#include "gpu_types.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#else
using cudaStream_t = void*;
#endif

namespace gutibm::gpu {

void launch_superpose_kernel(
    const double* src_x, const double* src_y, const double* src_z,
    const GfSourceParams* params, double* grid_conc,
    const DomainParams& dom, const AdvectionParams& adv,
    int num_sources, int span, cudaStream_t stream);

void launch_field_update_kernel(
    double* conc, const double* reac, int ncells, int num_species,
    double dt, cudaStream_t stream);

void launch_apply_boundaries_kernel(
    double* conc, int nx, int ny, int nz, int num_species,
    const double* boundary_conc, cudaStream_t stream);

void launch_grid_coupling_kernel(
    const double* x, const double* y, const double* z, int* grid_cell,
    const int* state, double lo0, double lo1, double lo2, double dx,
    int nx, int ny, int nz, int num_agents, cudaStream_t stream);

void launch_metabolism_kernel(
    const double* conc_carbon, const double* conc_iron, const double* conc_b12,
    const double* conc_acetate, const double* conc_eut,
    double* reac_carbon, double* reac_iron, double* reac_b12,
    double* mu_realized, double* biomass, double* radius, double* mass, double* age,
    const int* grid_cell, const int* state,
    const double* mu_max, const double* km_b12, const double* km_carbon,
    const double* receptor_expr, const double* ligand_affinity,
    const int* bi_loci_count, const double* plasmid_amelioration,
    int num_agents, double dt, double dx, double cell_density,
    double km_iron_primary, double km_iron_iroN, double km_iron_iutA, double km_iron_fiu,
    double maintenance_rate, double metE_penalty, double metE_acetate_max_factor,
    double metE_acetate_km, double eut_max_penalty, double eut_km,
    double yield_carbon, double yield_iron, double yield_b12,
    cudaStream_t stream);

void launch_spatial_hash_build_kernel(
    const double* x, const double* y, const double* z, const int* state,
    int* cell_keys, int* sorted_indices, int num_agents,
    double lo0, double lo1, double lo2, double cell_size,
    int nx_cells, int ny_cells, int nz_cells, cudaStream_t stream);

}  // namespace gutibm::gpu

#endif  // GUTIBM_GPU_KERNELS_H
