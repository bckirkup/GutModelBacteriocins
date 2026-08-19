#include "gpu_kernels.h"
#include <cuda_runtime.h>
#include <cmath>

namespace gutibm {
namespace gpu {

static constexpr double PI_M = 3.14159265358979323846;

__device__ inline double maxd(double a, double b) { return a > b ? a : b; }

__global__ void metabolism_kernel(
    const double* conc_carbon, const double* conc_iron, const double* conc_b12,
    const double* conc_acetate, const double* conc_eut,
    const double* conc_oxygen,
    double* reac_carbon, double* reac_iron, double* reac_b12,
    double* reac_acetate,
    double* mu_realized, double* biomass, double* radius, double* mass, double* age,
    const int* grid_cell, const int* state,
    const double* mu_max, const double* km_b12, const double* km_carbon,
    double* receptor_expr, const double* receptor_expr_base,
    const double* ligand_affinity, const int* iron_receptor,
    const int* bi_loci_count, const double* plasmid_amelioration,
    int num_agents, int local_agent_count, int agent_stride,
    int receptor_count, double dt, double cell_volume,
    double cell_density,
    double km_iron_primary, double km_iron_iroN, double km_iron_iutA, double km_iron_fiu,
    double maintenance_rate, double metE_penalty, double metE_acetate_max_factor,
    double metE_acetate_km, double eut_max_penalty, double eut_km,
    double yield_carbon, double yield_iron, double yield_b12,
    int iron_uptake_enabled, int b12_uptake_enabled, int eut_enabled,
    int fur_enabled, double fur_Km, double fur_upregulation_max,
    double fur_receptor_max,
    int acetate_enabled, double acetate_overflow_threshold,
    double acetate_overflow_rate, double acetate_scavenge_rate,
    double acetate_scavenge_Km,
    int o2_enabled, double o2_boost_max, double o2_Km,
    double* agent_uptake_totals,
    int global_nx, int global_ny, int storage_nx,
    int owned_global_x_begin, int owned_global_x_end,
    int owned_storage_x_begin) {

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  if (state[i] == 3) return;

  int cell = map_global_cell_to_storage(
      grid_cell[i], global_nx, global_ny, storage_nx,
      owned_global_x_begin, owned_global_x_end, owned_storage_x_begin);
  if (cell < 0) {
    mu_realized[i] = 0.0;
    return;
  }

  double S_carbon = conc_carbon ? conc_carbon[cell] : 1.0;
  double S_iron = iron_uptake_enabled && conc_iron
      ? conc_iron[cell] : 1.0;
  double S_b12 = b12_uptake_enabled && conc_b12
      ? conc_b12[cell] : 1.0;

  if (fur_enabled && iron_uptake_enabled) {
    const double fur_factor = 1.0 + fur_upregulation_max * fur_Km
        / (fur_Km + S_iron);
    for (int r = 0; r < receptor_count; ++r) {
      const double base_expression =
          receptor_expr_base[r * agent_stride + i];
      receptor_expr[r * agent_stride + i] = iron_receptor[r]
          ? fmin(base_expression * fur_factor, fur_receptor_max)
          : base_expression;
    }
  }

  // receptor_expr / ligand_affinity are SoA: index = receptor * agent_stride
  // + agent (matches AgentPoolGpu::sync_from_host and receptor_kernel).
  double expr_fepA = receptor_expr[1 * agent_stride + i];
  double expr_iroN = receptor_expr[4 * agent_stride + i];
  double expr_iutA = receptor_expr[7 * agent_stride + i];
  double expr_fiu  = receptor_expr[5 * agent_stride + i];
  double expr_btuB = maxd(receptor_expr[0 * agent_stride + i], 0.01);

  double lig_fepA = maxd(ligand_affinity[1 * agent_stride + i], 0.01);
  double lig_btuB = maxd(ligand_affinity[0 * agent_stride + i], 0.01);
  double lig_iroN = maxd(ligand_affinity[4 * agent_stride + i], 0.01);
  double lig_iutA = maxd(ligand_affinity[7 * agent_stride + i], 0.01);
  double lig_fiu  = maxd(ligand_affinity[5 * agent_stride + i], 0.01);

  double monod_iron = 1.0;
  if (iron_uptake_enabled) {
    double iron_uptake = 0.0;
    iron_uptake += expr_fepA * lig_fepA * S_iron
        / (km_iron_primary + S_iron);
    iron_uptake += expr_iroN * lig_iroN * S_iron
        / (km_iron_iroN + S_iron);
    iron_uptake += expr_iutA * lig_iutA * S_iron
        / (km_iron_iutA + S_iron);
    iron_uptake += expr_fiu * lig_fiu * S_iron
        / (km_iron_fiu + S_iron);
    monod_iron = iron_uptake
        / (1.0 + expr_iroN + expr_iutA + expr_fiu);
  }

  double Km_b12  = km_b12[i] / (expr_btuB * lig_btuB);
  double Km_carb = km_carbon[i];
  double monod_carbon = S_carbon / (Km_carb + S_carbon);
  double monod_b12 = b12_uptake_enabled
      ? S_b12 / (Km_b12 + S_b12) : 1.0;

  double mu = mu_max[i] * monod_carbon * monod_iron * monod_b12;

  if (o2_enabled && conc_oxygen) {
    const double s_o2 = conc_oxygen[cell];
    const double monod_o2_boost =
        1.0 + o2_boost_max * s_o2 / (o2_Km + s_o2);
    mu *= monod_o2_boost;
  }

  if (expr_btuB < 0.5) {
    double metE_eff = metE_penalty;
    if (conc_acetate) {
      double acetate_conc = conc_acetate[cell];
      double acetate_factor = 1.0 + (metE_acetate_max_factor - 1.0)
          * acetate_conc / (metE_acetate_km + acetate_conc);
      metE_eff *= acetate_factor;
    }
    double eut_conc = eut_enabled && conc_eut ? conc_eut[cell] : 0.0;
    double eut_effect = eut_enabled
        ? eut_max_penalty * eut_conc / (eut_km + eut_conc) : 0.0;
    mu *= (1.0 - metE_eff - eut_effect);
  }

  int n_loci = bi_loci_count[i];
  if (n_loci > 0) {
    double per_locus = fmax(0.0, 0.02 - plasmid_amelioration[i]);
    double plasmid_cost = per_locus * n_loci;
    plasmid_cost = plasmid_cost > 0.10 ? 0.10 : plasmid_cost;
    mu *= (1.0 - plasmid_cost);
  }

  mu -= maintenance_rate;
  mu_realized[i] = mu;

  if (i >= local_agent_count) return;

  double d_biomass = mu * biomass[i] * dt;
  biomass[i] += d_biomass;
  biomass[i] = biomass[i] > 1.0e-20 ? biomass[i] : 1.0e-20;

  double vol = biomass[i] / cell_density;
  radius[i] = cbrt(3.0 * vol / (4.0 * PI_M));
  mass[i] = biomass[i];
  age[i] += dt;

  if (d_biomass <= 0.0 || dt <= 0.0) return;
  if (cell_volume <= 0.0) return;

  if (reac_carbon) {
    const double uptake = d_biomass * yield_carbon;
    atomicAdd(&reac_carbon[cell], -uptake / (cell_volume * dt));
    if (agent_uptake_totals) atomicAdd(&agent_uptake_totals[0], uptake);
  }
  if (iron_uptake_enabled && reac_iron) {
    const double uptake = d_biomass * yield_iron;
    atomicAdd(&reac_iron[cell], -uptake / (cell_volume * dt));
    if (agent_uptake_totals) atomicAdd(&agent_uptake_totals[1], uptake);
  }
  // Spec 6 §3 — B12/corrinoid is not depleted (constant bioavailable pool).
  // reac_b12 / yield_b12 retained in the signature for ABI stability but unused.
  (void)reac_b12;
  (void)yield_b12;
  if (acetate_enabled && reac_acetate) {
    const double acetate_conc = conc_acetate ? conc_acetate[cell] : 0.0;
    if (mu_realized[i] > acetate_overflow_threshold) {
      atomicAdd(&reac_acetate[cell],
                acetate_overflow_rate * biomass[i] / cell_volume);
    }
    const double scavenge = acetate_scavenge_rate * acetate_conc
        / (acetate_scavenge_Km + acetate_conc)
        * biomass[i] / cell_volume;
    atomicAdd(&reac_acetate[cell], -scavenge);
  }
}

void launch_metabolism_kernel(
    const double* conc_carbon, const double* conc_iron, const double* conc_b12,
    const double* conc_acetate, const double* conc_eut,
    const double* conc_oxygen,
    double* reac_carbon, double* reac_iron, double* reac_b12,
    double* reac_acetate,
    double* mu_realized, double* biomass, double* radius, double* mass, double* age,
    const int* grid_cell, const int* state,
    const double* mu_max, const double* km_b12, const double* km_carbon,
    double* receptor_expr, const double* receptor_expr_base,
    const double* ligand_affinity, const int* iron_receptor,
    const int* bi_loci_count, const double* plasmid_amelioration,
    int num_agents, int local_agent_count, int agent_stride,
    int receptor_count, double dt, double cell_volume,
    double cell_density,
    double km_iron_primary, double km_iron_iroN, double km_iron_iutA, double km_iron_fiu,
    double maintenance_rate, double metE_penalty, double metE_acetate_max_factor,
    double metE_acetate_km, double eut_max_penalty, double eut_km,
    double yield_carbon, double yield_iron, double yield_b12,
    int iron_uptake_enabled, int b12_uptake_enabled, int eut_enabled,
    int fur_enabled, double fur_Km, double fur_upregulation_max,
    double fur_receptor_max,
    int acetate_enabled, double acetate_overflow_threshold,
    double acetate_overflow_rate, double acetate_scavenge_rate,
    double acetate_scavenge_Km,
    int o2_enabled, double o2_boost_max, double o2_Km,
    double* agent_uptake_totals,
    int global_nx, int global_ny, int storage_nx,
    int owned_global_x_begin, int owned_global_x_end,
    int owned_storage_x_begin,
    cudaStream_t stream) {
  if (num_agents <= 0) return;
  int block = 256;
  int grid = (num_agents + block - 1) / block;
  metabolism_kernel<<<grid, block, 0, stream>>>(
      conc_carbon, conc_iron, conc_b12, conc_acetate, conc_eut,
      conc_oxygen,
      reac_carbon, reac_iron, reac_b12, reac_acetate,
      mu_realized, biomass, radius, mass, age,
      grid_cell, state, mu_max, km_b12, km_carbon,
      receptor_expr, receptor_expr_base, ligand_affinity, iron_receptor,
      bi_loci_count, plasmid_amelioration,
      num_agents, local_agent_count, agent_stride, receptor_count, dt,
      cell_volume, cell_density,
      km_iron_primary, km_iron_iroN, km_iron_iutA, km_iron_fiu,
      maintenance_rate, metE_penalty, metE_acetate_max_factor,
      metE_acetate_km, eut_max_penalty, eut_km,
      yield_carbon, yield_iron, yield_b12,
      iron_uptake_enabled, b12_uptake_enabled, eut_enabled,
      fur_enabled, fur_Km, fur_upregulation_max, fur_receptor_max,
      acetate_enabled, acetate_overflow_threshold, acetate_overflow_rate,
      acetate_scavenge_rate, acetate_scavenge_Km,
      o2_enabled, o2_boost_max, o2_Km, agent_uptake_totals,
      global_nx, global_ny, storage_nx, owned_global_x_begin,
      owned_global_x_end, owned_storage_x_begin);
}

}  // namespace gpu
}  // namespace gutibm
