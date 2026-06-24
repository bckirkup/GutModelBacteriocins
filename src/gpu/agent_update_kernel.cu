#include "gpu_kernels.h"
#include <cuda_runtime.h>
#include <cmath>

namespace gutibm {
namespace gpu {

static constexpr int NUM_RECEPTORS_GPU = 8;
static constexpr double PI_M = 3.14159265358979323846;

__device__ inline double maxd(double a, double b) { return a > b ? a : b; }

__global__ void metabolism_kernel(
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
    double yield_carbon, double yield_iron, double yield_b12) {

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  if (state[i] == 3) return;

  int cell = grid_cell[i];
  if (cell < 0) {
    mu_realized[i] = 0.0;
    return;
  }

  double S_carbon = conc_carbon ? conc_carbon[cell] : 1.0;
  double S_iron   = conc_iron ? conc_iron[cell] : 1.0;
  double S_b12    = conc_b12 ? conc_b12[cell] : 1.0;

  const double* rex = receptor_expr + i * NUM_RECEPTORS_GPU;
  const double* lig = ligand_affinity + i * NUM_RECEPTORS_GPU;

  double expr_fepA = rex[1];
  double expr_iroN = rex[4];
  double expr_iutA = rex[7];
  double expr_fiu  = rex[5];
  double expr_btuB = maxd(rex[0], 0.01);

  double lig_fepA = maxd(lig[1], 0.01);
  double lig_btuB = maxd(lig[0], 0.01);
  double lig_iroN = maxd(lig[4], 0.01);
  double lig_iutA = maxd(lig[7], 0.01);
  double lig_fiu  = maxd(lig[5], 0.01);

  double iron_uptake = 0.0;
  iron_uptake += expr_fepA * lig_fepA * S_iron / (km_iron_primary + S_iron);
  iron_uptake += expr_iroN * lig_iroN * S_iron / (km_iron_iroN + S_iron);
  iron_uptake += expr_iutA * lig_iutA * S_iron / (km_iron_iutA + S_iron);
  iron_uptake += expr_fiu  * lig_fiu  * S_iron / (km_iron_fiu  + S_iron);
  double monod_iron = iron_uptake / (1.0 + expr_iroN + expr_iutA + expr_fiu);

  double Km_b12  = km_b12[i] / (expr_btuB * lig_btuB);
  double Km_carb = km_carbon[i];
  double monod_carbon = S_carbon / (Km_carb + S_carbon);
  double monod_b12    = S_b12 / (Km_b12 + S_b12);

  double mu = mu_max[i] * monod_carbon * monod_iron * monod_b12;

  if (expr_btuB < 0.5) {
    double metE_eff = metE_penalty;
    if (conc_acetate) {
      double acetate_conc = conc_acetate[cell];
      double acetate_factor = 1.0 + (metE_acetate_max_factor - 1.0)
          * acetate_conc / (metE_acetate_km + acetate_conc);
      metE_eff *= acetate_factor;
    }
    double eut_conc = conc_eut ? conc_eut[cell] : 0.0;
    double eut_effect = eut_max_penalty * eut_conc / (eut_km + eut_conc);
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

  double d_biomass = mu * biomass[i] * dt;
  biomass[i] += d_biomass;
  biomass[i] = biomass[i] > 1.0e-20 ? biomass[i] : 1.0e-20;

  double vol = biomass[i] / cell_density;
  radius[i] = cbrt(3.0 * vol / (4.0 * PI_M));
  mass[i] = biomass[i];
  age[i] += dt;

  if (d_biomass <= 0.0) return;
  double cell_vol = dx * dx * dx;
  if (cell_vol <= 0.0) return;

  if (reac_carbon) {
    atomicAdd(&reac_carbon[cell], -d_biomass * yield_carbon / cell_vol);
  }
  if (reac_iron) {
    atomicAdd(&reac_iron[cell], -d_biomass * yield_iron / cell_vol);
  }
  if (reac_b12) {
    atomicAdd(&reac_b12[cell], -d_biomass * yield_b12 / cell_vol);
  }
}

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
    cudaStream_t stream) {
  int block = 256;
  int grid = (num_agents + block - 1) / block;
  metabolism_kernel<<<grid, block, 0, stream>>>(
      conc_carbon, conc_iron, conc_b12, conc_acetate, conc_eut,
      reac_carbon, reac_iron, reac_b12,
      mu_realized, biomass, radius, mass, age,
      grid_cell, state, mu_max, km_b12, km_carbon,
      receptor_expr, ligand_affinity, bi_loci_count, plasmid_amelioration,
      num_agents, dt, dx, cell_density,
      km_iron_primary, km_iron_iroN, km_iron_iutA, km_iron_fiu,
      maintenance_rate, metE_penalty, metE_acetate_max_factor,
      metE_acetate_km, eut_max_penalty, eut_km,
      yield_carbon, yield_iron, yield_b12);
}

}  // namespace gpu
}  // namespace gutibm
