#include "gpu_kernels.h"
#include <cuda_runtime.h>
#include <cmath>

namespace gutibm {
namespace gpu {

namespace {

__device__ double dynamic_mucin_liberation(double mucin_conc,
                                           const VbfLaunchParams& p) {
  const double substrate =
      mucin_conc / (p.mucin_Km_degradation + mucin_conc);
  return p.mucin_k_liberation * p.vbf_density * substrate;
}

__device__ void apply_vbf_at_cell(int cell,
                                  int iz,
                                  double static_liberation,
                                  double z_weight,
                                  const VbfLaunchParams& p,
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
                                  double cell_volume,
                                  double dt) {
  if (p.use_dynamic_mucin && p.mucin_enabled && reac_mucin && conc_mucin) {
    const double liberation =
        dynamic_mucin_liberation(conc_mucin[cell], p);
    reac_mucin[cell] -= liberation;
    if (reac_carbon) {
      reac_carbon[cell] += liberation;
      if (vbf_totals) atomicAdd(&vbf_totals[0], liberation * cell_volume * dt);
    }
  } else if (reac_carbon) {
    reac_carbon[cell] += static_liberation;
    if (vbf_totals) {
      atomicAdd(&vbf_totals[0], static_liberation * cell_volume * dt);
    }
  }

  if (reac_carbon && conc_carbon && p.carbon_sink_vmax > 0.0) {
    const double c = conc_carbon[cell];
    const double sink = p.carbon_sink_vmax * c / (p.carbon_sink_km + c);
    reac_carbon[cell] -= sink;
    if (vbf_totals) atomicAdd(&vbf_totals[1], sink * cell_volume * dt);
  }

  if (reac_iron && conc_iron) {
    const double sink = p.nutrient_sink * conc_iron[cell];
    reac_iron[cell] -= sink;
    if (vbf_totals) atomicAdd(&vbf_totals[2], sink * cell_volume * dt);
  }

  if (p.oxygen_enabled && reac_oxygen && conc_oxygen) {
    const double sink = p.oxygen_vbf_sink * conc_oxygen[cell];
    reac_oxygen[cell] -= sink;
    if (vbf_totals) atomicAdd(&vbf_totals[3], sink * cell_volume * dt);
  }

  if (p.acetate_enabled && reac_acetate) {
    reac_acetate[cell] += p.acetate_vbf_production * z_weight;
    reac_acetate[cell] -= p.acetate_vbf_consumption;
    if (iz == 0) {
      reac_acetate[cell] -= p.acetate_epithelial_uptake;
    }
  }

  if (p.mucin_enabled && reac_mucin && iz == 0) {
    reac_mucin[cell] += p.mucin_secretion_rate;
  }
}

__global__ void o2_depletion_kernel(double* reac_oxygen,
                                    const double* mu_realized,
                                    const int* grid_cell,
                                    const int* state,
                                    int num_agents,
                                    double q_consumption,
                                    double q_maintenance,
                                    double inv_cell_vol, int global_nx,
                                    int global_ny, int storage_nx,
                                    int owned_global_x_begin,
                                    int owned_global_x_end,
                                    int owned_storage_x_begin) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  if (state[i] == 3) return;

  const int cell = grid_cell[i];
  if (cell < 0 || reac_oxygen == nullptr) return;
  const int local_cell = map_global_cell_to_storage(
      cell, global_nx, global_ny, storage_nx, owned_global_x_begin,
      owned_global_x_end, owned_storage_x_begin);
  if (local_cell < 0) return;

  const double mu = mu_realized[i] > 0.0 ? mu_realized[i] : 0.0;
  const double o2_use = q_consumption * mu + q_maintenance;
  atomicAdd(reac_oxygen + local_cell, -o2_use * inv_cell_vol);
}

__global__ void vbf_coupling_kernel(int ncells,
                                    VbfLaunchParams params,
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
                                    double dt) {
  const int cell = blockIdx.x * blockDim.x + threadIdx.x;
  if (cell >= ncells) return;

  const int local_nx = params.owned_x_end - params.owned_x_begin;
  const int iz = cell / (local_nx * params.ny);
  const int rem = cell % (local_nx * params.ny);
  const int iy = rem / local_nx;
  const int ix = params.owned_x_begin + rem % local_nx;
  const int storage_cell = iz * (params.storage_nx * params.ny)
      + iy * params.storage_nx + ix;
  const double z_rel = (iz + 0.5) * params.dx_z;
  const double z_weight = params.mucin_z_gradient_enabled
      ? exp(-z_rel / params.mucin_z_gradient_lambda)
      : 1.0;
  double static_liberation = 0.0;
  if (!params.use_dynamic_mucin) {
    static_liberation = params.mucin_liberation;
    if (params.mucin_z_gradient_enabled) {
      static_liberation *= exp(-z_rel / params.mucin_z_gradient_lambda);
    }
  }

  apply_vbf_at_cell(storage_cell, iz, static_liberation, z_weight, params,
                    reac_carbon, conc_carbon,
                    reac_iron, conc_iron,
                    reac_oxygen, conc_oxygen,
                    reac_acetate,
                    reac_mucin, conc_mucin, vbf_totals,
                    params.dx_x * params.dx_y * params.dx_z, dt);
}

}  // namespace

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
                                cudaStream_t stream) {
  if (num_agents <= 0 || reac_oxygen == nullptr) return;
  const int block = 256;
  const int grid = (num_agents + block - 1) / block;
  o2_depletion_kernel<<<grid, block, 0, stream>>>(
      reac_oxygen, mu_realized, grid_cell, state, num_agents,
      q_consumption, q_maintenance, inv_cell_vol, global_nx, global_ny,
      storage_nx, owned_global_x_begin, owned_global_x_end,
      owned_storage_x_begin);
}

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
                                cudaStream_t stream) {
  if (ncells <= 0) return;
  const int block = 256;
  const int grid = (ncells + block - 1) / block;
  vbf_coupling_kernel<<<grid, block, 0, stream>>>(
      ncells, params,
      reac_carbon, conc_carbon,
      reac_iron, conc_iron,
      reac_oxygen, conc_oxygen,
      reac_acetate,
      reac_mucin, conc_mucin, vbf_totals, dt);
}

}  // namespace gpu
}  // namespace gutibm
