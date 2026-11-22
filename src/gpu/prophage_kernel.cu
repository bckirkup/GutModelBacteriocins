/* -----------------------------------------------------------------------
   GutIBM – CUDA kernels for Spec 14 temperate prophage induction & virions
   ----------------------------------------------------------------------- */

#include "gpu_kernels.h"
#include <cmath>
#include <cuda_runtime.h>

namespace gutibm {
namespace gpu {

namespace {

constexpr int kStateSosInduced = 2;
constexpr int kStateDead = 3;

__device__ inline double atomic_add_double(double* addr, double val) {
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
  auto* address_as_ull = reinterpret_cast<unsigned long long*>(addr);
  unsigned long long old = *address_as_ull;
  unsigned long long assumed = 0;
  do {
    assumed = old;
    const double sum = __longlong_as_double(assumed) + val;
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(sum));
  } while (assumed != old);
  return __longlong_as_double(old);
#else
  return atomicAdd(addr, val);
#endif
}

__device__ inline double prophage_receptor_expression(
    const double* receptor_expr, int receptor_index, int num_agents, int i) {
  if (receptor_index < 0 || receptor_expr == nullptr) return 1.0;
  const double expr = receptor_expr[receptor_index * num_agents + i];
  return expr > 0.0 ? expr : 0.0;
}

__global__ void prophage_induction_prob_kernel(
    const int* state,
    const int* is_ghost,
    const int* has_inducible_prophage,
    const int* prophage_induced,
    const double* damage_signal,
    double* induction_probs,
    int num_agents,
    double k_spont,
    double k_max_damage,
    double K_damage,
    double hill_n,
    double dt) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;

  if (is_ghost[i] != 0 || state[i] == kStateDead || state[i] == kStateSosInduced
      || has_inducible_prophage[i] == 0 || prophage_induced[i] != 0) {
    induction_probs[i] = 0.0;
    return;
  }

  double hazard = k_spont > 0.0 ? k_spont : 0.0;
  const double d = damage_signal[i];
  if (d > 0.0 && k_max_damage > 0.0 && hill_n > 0.0) {
    if (K_damage <= 0.0) {
      hazard += k_max_damage;
    } else {
      const double ratio = K_damage / d;
      const double denom = 1.0 + pow(ratio, hill_n);
      if (isfinite(denom) && denom > 0.0) {
        hazard += k_max_damage / denom;
      }
    }
  }
  const double p = hazard > 0.0 ? -expm1(-hazard * dt) : 0.0;
  induction_probs[i] = p > 1.0 ? 1.0 : (p < 0.0 ? 0.0 : p);
}

__global__ void prophage_adsorption_accum_kernel(
    const int* grid_cell,
    const int* state,
    const int* is_ghost,
    const double* receptor_expr,
    double* cell_adsorption_rate,
    int num_agents,
    int receptor_index,
    double k_adsorb,
    double cell_volume,
    int global_nx, int global_ny, int storage_nx,
    int owned_global_x_begin, int owned_global_x_end,
    int owned_storage_x_begin) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  if (is_ghost[i] != 0 || state[i] == kStateDead || state[i] == kStateSosInduced) {
    return;
  }

  const int cell = map_global_cell_to_storage(
      grid_cell[i], global_nx, global_ny, storage_nx,
      owned_global_x_begin, owned_global_x_end, owned_storage_x_begin);
  if (cell < 0) return;

  const double expr =
      prophage_receptor_expression(receptor_expr, receptor_index, num_agents, i);
  const double rate = k_adsorb * expr / cell_volume;
  if (rate > 0.0) {
    atomic_add_double(cell_adsorption_rate + cell, rate);
  }
}

__global__ void prophage_cell_removal_kernel(
    const double* conc_virion,
    const double* cell_adsorption_rate,
    const double* washout_by_z,
    double* reac_delta,
    double* cell_adsorb_factor,
    double* totals_out,
    int storage_nx, int global_ny, int global_nz,
    int owned_storage_x_begin, int owned_storage_x_end,
    double k_virion_decay,
    double cell_volume,
    double dt) {
  const int owned_nx = owned_storage_x_end - owned_storage_x_begin;
  const int total_owned = owned_nx * global_ny * global_nz;
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= total_owned) return;

  const int plane = owned_nx * global_ny;
  const int iz = idx / plane;
  const int rem = idx % plane;
  const int iy = rem / owned_nx;
  const int ix = owned_storage_x_begin + (rem % owned_nx);
  const int cell = iz * (storage_nx * global_ny) + iy * storage_nx + ix;

  const double conc = conc_virion[cell];
  if (conc <= 0.0) {
    reac_delta[cell] = 0.0;
    cell_adsorb_factor[cell] = 0.0;
    return;
  }

  const double k_decay = (k_virion_decay > 0.0 ? k_virion_decay : 0.0)
      + washout_by_z[iz];
  const double a_total = cell_adsorption_rate[cell];
  const double lambda = k_decay + a_total;
  if (lambda <= 0.0) {
    reac_delta[cell] = 0.0;
    cell_adsorb_factor[cell] = 0.0;
    return;
  }

  const double frac_removed = -expm1(-lambda * dt);
  const double removed_conc = conc * frac_removed;
  const double removed_virions = removed_conc * cell_volume;
  const double decayed = removed_virions * (k_decay / lambda);
  const double adsorbed = removed_virions * (a_total / lambda);

  reac_delta[cell] = -removed_conc / dt;
  cell_adsorb_factor[cell] = removed_virions / lambda;

  if (decayed > 0.0) {
    atomic_add_double(totals_out, decayed);
  }
  if (adsorbed > 0.0) {
    atomic_add_double(totals_out + 1, adsorbed);
  }
}

__global__ void prophage_infect_prob_kernel(
    const int* grid_cell,
    const int* state,
    const int* is_ghost,
    const int* has_prophage_immunity,
    const double* receptor_expr,
    const double* cell_adsorb_factor,
    double* infect_probs,
    int num_agents,
    int receptor_index,
    double k_adsorb,
    double cell_volume,
    int global_nx, int global_ny, int storage_nx,
    int owned_global_x_begin, int owned_global_x_end,
    int owned_storage_x_begin) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;

  // Existing lysogens (inducible OR cryptic) have repressor-mediated
  // superinfection immunity: zero infection / zero lytic kill.
  if (is_ghost[i] != 0 || state[i] == kStateDead || state[i] == kStateSosInduced
      || has_prophage_immunity[i] != 0) {
    infect_probs[i] = 0.0;
    return;
  }

  const int cell = map_global_cell_to_storage(
      grid_cell[i], global_nx, global_ny, storage_nx,
      owned_global_x_begin, owned_global_x_end, owned_storage_x_begin);
  if (cell < 0) {
    infect_probs[i] = 0.0;
    return;
  }

  const double factor = cell_adsorb_factor[cell];
  if (factor <= 0.0) {
    infect_probs[i] = 0.0;
    return;
  }

  const double expr =
      prophage_receptor_expression(receptor_expr, receptor_index, num_agents, i);
  const double v_adsorbed = (k_adsorb * expr / cell_volume) * factor;
  const double p = v_adsorbed > 0.0 ? -expm1(-v_adsorbed) : 0.0;
  infect_probs[i] = p > 1.0 ? 1.0 : (p < 0.0 ? 0.0 : p);
}

}  // namespace

void launch_prophage_induction_prob_kernel(
    const int* state,
    const int* is_ghost,
    const int* has_inducible_prophage,
    const int* prophage_induced,
    const double* damage_signal,
    double* induction_probs,
    int num_agents,
    double k_spont,
    double k_max_damage,
    double K_damage,
    double hill_n,
    double dt,
    cudaStream_t stream) {
  if (num_agents <= 0) return;
  const int block = 256;
  const int grid = (num_agents + block - 1) / block;
  prophage_induction_prob_kernel<<<grid, block, 0, stream>>>(
      state, is_ghost, has_inducible_prophage, prophage_induced,
      damage_signal, induction_probs, num_agents, k_spont, k_max_damage,
      K_damage, hill_n, dt);
}

void launch_prophage_adsorption_accum_kernel(
    const int* grid_cell,
    const int* state,
    const int* is_ghost,
    const double* receptor_expr,
    double* cell_adsorption_rate,
    int num_agents,
    int receptor_index,
    double k_adsorb,
    double cell_volume,
    int global_nx, int global_ny, int storage_nx,
    int owned_global_x_begin, int owned_global_x_end,
    int owned_storage_x_begin,
    cudaStream_t stream) {
  if (num_agents <= 0) return;
  const int block = 256;
  const int grid = (num_agents + block - 1) / block;
  prophage_adsorption_accum_kernel<<<grid, block, 0, stream>>>(
      grid_cell, state, is_ghost, receptor_expr, cell_adsorption_rate,
      num_agents, receptor_index, k_adsorb, cell_volume,
      global_nx, global_ny, storage_nx, owned_global_x_begin,
      owned_global_x_end, owned_storage_x_begin);
}

void launch_prophage_cell_removal_kernel(
    const double* conc_virion,
    const double* cell_adsorption_rate,
    const double* washout_by_z,
    double* reac_delta,
    double* cell_adsorb_factor,
    double* totals_out,
    int storage_nx, int global_ny, int global_nz,
    int owned_storage_x_begin, int owned_storage_x_end,
    double k_virion_decay,
    double cell_volume,
    double dt,
    cudaStream_t stream) {
  const int owned_nx = owned_storage_x_end - owned_storage_x_begin;
  const int total_owned = owned_nx * global_ny * global_nz;
  if (total_owned <= 0) return;
  const int block = 256;
  const int grid = (total_owned + block - 1) / block;
  prophage_cell_removal_kernel<<<grid, block, 0, stream>>>(
      conc_virion, cell_adsorption_rate, washout_by_z, reac_delta,
      cell_adsorb_factor, totals_out, storage_nx, global_ny, global_nz,
      owned_storage_x_begin, owned_storage_x_end, k_virion_decay,
      cell_volume, dt);
}

void launch_prophage_infect_prob_kernel(
    const int* grid_cell,
    const int* state,
    const int* is_ghost,
    const int* has_prophage_immunity,
    const double* receptor_expr,
    const double* cell_adsorb_factor,
    double* infect_probs,
    int num_agents,
    int receptor_index,
    double k_adsorb,
    double cell_volume,
    int global_nx, int global_ny, int storage_nx,
    int owned_global_x_begin, int owned_global_x_end,
    int owned_storage_x_begin,
    cudaStream_t stream) {
  if (num_agents <= 0) return;
  const int block = 256;
  const int grid = (num_agents + block - 1) / block;
  prophage_infect_prob_kernel<<<grid, block, 0, stream>>>(
      grid_cell, state, is_ghost, has_prophage_immunity, receptor_expr,
      cell_adsorb_factor, infect_probs, num_agents, receptor_index,
      k_adsorb, cell_volume, global_nx, global_ny, storage_nx,
      owned_global_x_begin, owned_global_x_end, owned_storage_x_begin);
}

}  // namespace gpu
}  // namespace gutibm
