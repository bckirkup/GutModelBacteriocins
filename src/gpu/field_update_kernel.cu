#include "gpu_kernels.h"
#include <cuda_runtime.h>
#include <cmath>

namespace gutibm {
namespace gpu {

__device__ inline int iz_cell_index(int ix, int iy, int iz, int nx, int ny) {
  return iz * (nx * ny) + iy * nx + ix;
}

__global__ void field_update_kernel(double* conc, const double* reac,
                                    int ncells, int num_species, double dt) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int total = ncells * num_species;
  if (idx >= total) return;
  int spec = idx / ncells;
  int cell = idx % ncells;
  double c = conc[spec * ncells + cell] + reac[spec * ncells + cell] * dt;
  conc[spec * ncells + cell] = c > 0.0 ? c : 0.0;
}

__global__ void apply_boundaries_kernel(double* conc, int nx, int ny, int nz,
                                        int num_species, const double* boundary_conc) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int face_cells = nx * ny;
  int total = face_cells * 2 * num_species;
  if (idx >= total) return;

  int spec = idx / (face_cells * 2);
  int rem = idx % (face_cells * 2);
  int face = rem / face_cells;
  int cell2d = rem % face_cells;
  int ix = cell2d % nx;
  int iy = cell2d / nx;
  int ncells = nx * ny * nz;
  int base = spec * ncells;

  if (face == 0) {
    int cidx = iz_cell_index(ix, iy, 0, nx, ny);
    conc[base + cidx] = boundary_conc[spec];
  } else if (nz >= 2) {
    int top = iz_cell_index(ix, iy, nz - 1, nx, ny);
    int below = iz_cell_index(ix, iy, nz - 2, nx, ny);
    conc[base + top] = conc[base + below];
  }
}

__global__ void nutrient_depletion_kernel(
    const int* grid_cell, const double* mu_realized, const double* biomass,
    const int* state, double* reac_iron, double* reac_b12, double* reac_carbon,
    int num_agents, bool has_iron, bool has_b12, bool has_carbon) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  if (state[i] == 3) return;  // DEAD
  int cell = grid_cell[i];
  if (cell < 0) return;
  double consumption = mu_realized[i] * biomass[i];
  if (has_iron) atomicAdd(&reac_iron[cell], -consumption * 1.0e-6);
  if (has_b12) atomicAdd(&reac_b12[cell], -consumption * 1.0e-9);
  if (has_carbon) atomicAdd(&reac_carbon[cell], -consumption * 0.5);
}

__global__ void grid_coupling_kernel(
    const double* x, const double* y, const double* z, int* grid_cell,
    const int* state, double lo0, double lo1, double lo2, double dx,
    int nx, int ny, int nz, int num_agents) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  if (state[i] == 3) {
    grid_cell[i] = -1;
    return;
  }
  int ix = static_cast<int>((x[i] - lo0) / dx);
  int iy = static_cast<int>((y[i] - lo1) / dx);
  int iz = static_cast<int>((z[i] - lo2) / dx);
  ix = ix < 0 ? 0 : (ix >= nx ? nx - 1 : ix);
  iy = iy < 0 ? 0 : (iy >= ny ? ny - 1 : iy);
  iz = iz < 0 ? 0 : (iz >= nz ? nz - 1 : iz);
  grid_cell[i] = iz * (nx * ny) + iy * nx + ix;
}

void launch_field_update_kernel(double* conc, const double* reac, int ncells,
                                int num_species, double dt, cudaStream_t stream) {
  int total = ncells * num_species;
  int block = 256;
  int grid = (total + block - 1) / block;
  field_update_kernel<<<grid, block, 0, stream>>>(conc, reac, ncells, num_species, dt);
}

void launch_apply_boundaries_kernel(double* conc, int nx, int ny, int nz,
                                    int num_species, const double* boundary_conc,
                                    cudaStream_t stream) {
  int total = nx * ny * 2 * num_species;
  int block = 256;
  int grid = (total + block - 1) / block;
  apply_boundaries_kernel<<<grid, block, 0, stream>>>(
      conc, nx, ny, nz, num_species, boundary_conc);
}

void launch_nutrient_depletion_kernel(
    const int* grid_cell, const double* mu_realized, const double* biomass,
    const int* state, double* reac_iron, double* reac_b12, double* reac_carbon,
    int num_agents, bool has_iron, bool has_b12, bool has_carbon,
    cudaStream_t stream) {
  if (num_agents <= 0) return;
  int block = 256;
  int grid = (num_agents + block - 1) / block;
  nutrient_depletion_kernel<<<grid, block, 0, stream>>>(
      grid_cell, mu_realized, biomass, state,
      reac_iron, reac_b12, reac_carbon,
      num_agents, has_iron, has_b12, has_carbon);
}

void launch_grid_coupling_kernel(
    const double* x, const double* y, const double* z, int* grid_cell,
    const int* state, double lo0, double lo1, double lo2, double dx,
    int nx, int ny, int nz, int num_agents, cudaStream_t stream) {
  if (num_agents <= 0) return;
  int block = 256;
  int grid = (num_agents + block - 1) / block;
  grid_coupling_kernel<<<grid, block, 0, stream>>>(
      x, y, z, grid_cell, state, lo0, lo1, lo2, dx, nx, ny, nz, num_agents);
}

}  // namespace gpu
}  // namespace gutibm
