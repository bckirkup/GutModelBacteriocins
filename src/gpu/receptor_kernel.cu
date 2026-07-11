#include "gpu_kernels.h"
#include "gpu_common.cuh"
#include <cuda_runtime.h>
#include <cmath>

namespace gutibm {
namespace gpu {

__device__ double toxin_occupancy_device(double tox_conc, double ligand_conc,
                                         double kd_tox, double kd_ligand,
                                         double receptor_expr,
                                         double toxin_aff, double ligand_aff) {
  const double eff_kd_tox = kd_tox / fmax(toxin_aff, 1.0e-6);
  const double eff_kd_ligand = kd_ligand / fmax(ligand_aff, 1.0e-6);
  const double competitive = 1.0 + ligand_conc / eff_kd_ligand;
  const double apparent_kd = eff_kd_tox * competitive;
  return receptor_expr * tox_conc / (apparent_kd + tox_conc);
}

__global__ void receptor_kill_prob_kernel(
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
    const double* conc_iron,
    const double* conc_siderophore,
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
    double kill_rate_microcin) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  if (state[i] == 3 || state[i] == 4) {
    kill_probs[i] = 0.0;
    return;
  }

  const int cell = grid_cell[i];
  if (cell < 0) {
    kill_probs[i] = 0.0;
    return;
  }

  const double tox_btuB = conc_btuB ? conc_btuB[cell] : 0.0;
  const double tox_fepA = conc_fepA ? conc_fepA[cell] : 0.0;
  const double tox_cirA = conc_cirA ? conc_cirA[cell] : 0.0;
  const double tox_fhuA = conc_fhuA ? conc_fhuA[cell] : 0.0;
  if (tox_btuB <= 0.0 && tox_fepA <= 0.0 && tox_cirA <= 0.0 && tox_fhuA <= 0.0) {
    kill_probs[i] = 0.0;
    return;
  }

  const double* rex_base = receptor_expr;
  const double* lig_base = ligand_affinity;
  const double* tox = toxin_affinity ? toxin_affinity + i * 8 : nullptr;
  const double* imm = immunity_eff + i * 4;

  double total_kill = 0.0;

  if (tox_btuB > 0.0) {
    const double ligand = conc_b12 ? conc_b12[cell] : 0.0;
    const double occ = toxin_occupancy_device(
        tox_btuB, ligand, kd_colicinE_btuB, kd_b12_btuB,
        rex_base[0 * num_agents + i], tox ? tox[0] : 1.0, lig_base[0 * num_agents + i]);
    total_kill += kill_rate_colicin * occ * imm[0] * dt;
  }
  if (tox_fepA > 0.0) {
    const double ligand = conc_iron ? conc_iron[cell] : 0.0;
    const double occ = toxin_occupancy_device(
        tox_fepA, ligand, kd_colicinB_fepA, kd_enterobactin,
        rex_base[1 * num_agents + i], tox ? tox[1] : 1.0, lig_base[1 * num_agents + i]);
    total_kill += kill_rate_colicin * occ * imm[1] * dt;
  }
  if (tox_cirA > 0.0) {
    const double ligand = conc_siderophore
        ? conc_siderophore[cell] * cirA_linearized_fraction : 0.0;
    const double occ = toxin_occupancy_device(
        tox_cirA, ligand, kd_colicinIa_cirA, kd_lin_enterobactin,
        rex_base[6 * num_agents + i], tox ? tox[6] : 1.0, lig_base[6 * num_agents + i]);
    total_kill += kill_rate_microcin * occ * imm[2] * dt;
  }
  if (tox_fhuA > 0.0) {
    const double ligand = conc_iron ? conc_iron[cell] : 0.0;
    const double occ = toxin_occupancy_device(
        tox_fhuA, ligand, kd_colicinM_fhuA, kd_ferrichrome,
        rex_base[3 * num_agents + i], tox ? tox[3] : 1.0, lig_base[3 * num_agents + i]);
    total_kill += kill_rate_colicin * occ * imm[3] * dt;
  }

  const double p = 1.0 - exp(-total_kill);
  kill_probs[i] = p > 1.0 ? 1.0 : p;
}

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
    const double* conc_iron,
    const double* conc_siderophore,
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
    cudaStream_t stream) {
  if (num_agents <= 0) return;
  int block = 256;
  int grid = (num_agents + block - 1) / block;
  receptor_kill_prob_kernel<<<grid, block, 0, stream>>>(
      grid_cell, state, receptor_expr, ligand_affinity, toxin_affinity,
      immunity_eff, conc_btuB, conc_fepA, conc_cirA, conc_fhuA,
      conc_b12, conc_iron, conc_siderophore, kill_probs, num_agents, dt,
      kd_b12_btuB, kd_colicinE_btuB, kd_enterobactin, kd_colicinB_fepA,
      kd_lin_enterobactin, kd_colicinIa_cirA, kd_colicinM_fhuA, kd_ferrichrome,
      cirA_linearized_fraction, kill_rate_colicin, kill_rate_microcin);
}

}  // namespace gpu
}  // namespace gutibm
