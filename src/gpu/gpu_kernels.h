#ifndef GUTIBM_GPU_KERNELS_H
#define GUTIBM_GPU_KERNELS_H

#include "gpu_types.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#else
using cudaStream_t = void*;
#endif

namespace gutibm::gpu {

#ifdef GUTIBM_CUDA
#ifdef __CUDACC__
#define GUTIBM_GPU_DEVICE __device__
#else
#define GUTIBM_GPU_DEVICE
#endif
GUTIBM_GPU_DEVICE inline int map_global_cell_to_storage(
    int global_cell, int global_nx, int global_ny, int storage_nx,
    int owned_global_x_begin, int owned_global_x_end,
    int owned_storage_x_begin) {
  if (global_cell < 0) return -1;
  const int global_plane = global_nx * global_ny;
  const int iz = global_cell / global_plane;
  const int rem = global_cell % global_plane;
  const int iy = rem / global_nx;
  const int ix = rem % global_nx;
  int local_x = owned_storage_x_begin + ix - owned_global_x_begin;
  if (local_x < 0 && owned_global_x_end == global_nx) {
    local_x += global_nx;
  } else if (local_x >= storage_nx && owned_global_x_begin == 0) {
    local_x -= global_nx;
  }
  if (local_x < 0 || local_x >= storage_nx) return -1;
  return iz * (storage_nx * global_ny) + iy * storage_nx + local_x;
}
#undef GUTIBM_GPU_DEVICE
#endif

void launch_superpose_kernel(
    const double* src_x, const double* src_y, const double* src_z,
    const GfSourceParams* params, double* grid_conc,
    const DomainParams& dom, const AdvectionParams& adv,
    int num_sources, int span_x, int span_y, int span_z,
    cudaStream_t stream);

void launch_field_update_kernel(
    double* conc, const double* reac, int ncells, int num_species,
    double dt, double* reaction_clip, double cell_volume,
    int storage_nx, int global_ny, int global_nz,
    int owned_x_begin, int owned_x_end,
    cudaStream_t stream);

void launch_apply_boundaries_kernel(
    double* conc, int nx, int ny, int nz, int num_species,
    const double* boundary_conc, cudaStream_t stream);

void launch_grid_coupling_kernel(
    const double* x, const double* y, const double* z, int* grid_cell,
    const int* state, double lo0, double lo1, double lo2,
    double dx_x, double dx_y, double dx_z,
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
    int num_agents, double dt, double cell_volume, double cell_density,
    double km_iron_primary, double km_iron_iroN, double km_iron_iutA, double km_iron_fiu,
    double maintenance_rate, double metE_penalty, double metE_acetate_max_factor,
    double metE_acetate_km, double eut_max_penalty, double eut_km,
    double yield_carbon, double yield_iron, double yield_b12,
    int o2_enabled, double o2_boost_max, double o2_Km,
    const double* conc_oxygen, double* agent_uptake_totals,
    int global_nx, int global_ny, int storage_nx,
    int owned_global_x_begin, int owned_global_x_end,
    int owned_storage_x_begin,
    cudaStream_t stream);

void launch_spatial_hash_build_kernel(
    const double* x, const double* y, const double* z, const int* state,
    int* cell_keys, int* sorted_indices, int num_agents,
    double lo0, double lo1, double lo2, double cell_size,
    int nx_cells, int ny_cells, int nz_cells, cudaStream_t stream);

void launch_diffuse_x_periodic(double* conc, int nx, int ny, int nz,
                               double alpha, double gamma, double corner,
                               double denominator, const double* correction,
                               cudaStream_t stream);
void launch_diffuse_y_periodic(double* conc, int storage_nx, int ny, int nz,
                               int owned_x_begin, int owned_x_end,
                               double alpha, double gamma, double corner,
                               double denominator, const double* correction,
                               cudaStream_t stream);
void launch_diffuse_z_bounded(double* conc, int storage_nx, int ny, int nz,
                              int owned_x_begin, int owned_x_end,
                              double alpha, double boundary_conc,
                              double cell_volume, double* face_exchange,
                              cudaStream_t stream);
void launch_set_epithelial_boundary(double* conc, int storage_nx, int ny,
                                    int owned_x_begin, int owned_x_end,
                                    double boundary_conc, double cell_volume,
                                    double* injected_amount,
                                    cudaStream_t stream);
void launch_set_luminal_neumann(double* conc, int storage_nx, int ny, int nz,
                                int owned_x_begin, int owned_x_end,
                                cudaStream_t stream);
void launch_shift_z_gradient(double* conc, int storage_nx, int ny, int nz,
                             int owned_x_begin, int owned_x_end, double dx_z,
                             double initial_conc, double lambda,
                             double boundary_conc, double scale,
                             cudaStream_t stream);
void launch_clamp_nonneg(double* conc, int storage_nx, int ny, int nz,
                         int owned_x_begin, int owned_x_end,
                         cudaStream_t stream);

void launch_o2_depletion_kernel(double* reac_oxygen,
                                const double* mu_realized,
                                const int* grid_cell,
                                const int* state,
                                int num_agents,
                                double q_consumption,
                                double q_maintenance,
                                double inv_cell_vol,
                                int global_nx, int global_ny, int storage_nx,
                                int owned_global_x_begin,
                                int owned_global_x_end,
                                int owned_storage_x_begin,
                                cudaStream_t stream);

void launch_vbf_coupling_kernel(int ncells,
                                const VbfLaunchParams& params,
                                double* reac_carbon,
                                const double* conc_carbon,
                                double* reac_iron,
                                const double* conc_iron,
                                double* reac_oxygen,
                                const double* conc_oxygen,
                                double* reac_acetate,
                                double* reac_mucin,
                                const double* conc_mucin,
                                double* vbf_totals,
                                double dt,
                                cudaStream_t stream);

void launch_fmm_far_local_kernel(
    const double* leaf_local,
    const double* leaf_center,
    const int* cell_leaf,
    const double* near_conc,
    double* out_conc,
    int ncells,
    int num_leaves,
    int coeffs_per_leaf,
    int expansion_order,
    double lo0, double lo1, double lo2, double dx_x, double dx_y, double dx_z,
    int nx, int ny,
    cudaStream_t stream);

void launch_mechanics_clear_kernel(double* dx, double* dy, double* dz,
                                   int* clamp_count, int num_agents,
                                   cudaStream_t stream);
void launch_mechanics_forces_kernel(
    const double* x, const double* y, const double* z,
    const double* radius, const int* state,
    const int* cell_offsets, const int* sorted_indices,
    double* dx, double* dy, double* dz,
    int num_agents, const MechanicsLaunchParams& params,
    cudaStream_t stream);
void launch_mechanics_apply_kernel(
    double* x, double* y, double* z,
    const double* radius,
    const double* dx, const double* dy, const double* dz,
    int* clamp_count, int num_agents,
    const MechanicsLaunchParams& params, cudaStream_t stream);

void launch_receptor_kill_prob_kernel(
    const int* grid_cell,
    const int* state,
    const double* receptor_expr,
    const double* ligand_affinity,
    const double* toxin_affinity,
    const double* immunity_eff,
    const double* conc_btuB,
    const double* conc_fepA,
    const double* conc_cirA,
    const double* conc_fhuA,
    const double* conc_b12,
    const double* conc_ferric_enterobactin,
    const double* conc_ferrichrome,
    double* kill_probs,
    int num_agents,
    double dt,
    double kd_b12_btuB,
    double kd_colicinE_btuB,
    double kd_enterobactin,
    double kd_colicinB_fepA,
    double kd_lin_enterobactin,
    double kd_colicinIa_cirA,
    double kd_colicinM_fhuA,
    double kd_ferrichrome,
    double cirA_linearized_fraction,
    double kill_rate_colicin,
    double kill_rate_microcin,
    int global_nx, int global_ny, int storage_nx,
    int owned_global_x_begin, int owned_global_x_end,
    int owned_storage_x_begin,
    cudaStream_t stream);

int diffusion_max_line_length();

}  // namespace gutibm::gpu

#endif  // GUTIBM_GPU_KERNELS_H
