#include "gpu_kernels.h"
#include <cuda_runtime.h>
#include <cmath>

namespace gutibm {
namespace gpu {

__device__ inline int iz_cell_index(int ix, int iy, int iz, int storage_nx,
                                    int ny) {
  return iz * (storage_nx * ny) + iy * storage_nx + ix;
}

__global__ void field_update_kernel(double* conc, const double* reac,
                                    int ncells, int num_species, double dt,
                                    double* reaction_clip,
                                    double cell_volume, int storage_nx,
                                    int global_ny, int global_nz,
                                    int owned_x_begin, int owned_x_end) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int owned_nx = owned_x_end - owned_x_begin;
  int total = owned_nx * global_ny * global_nz * num_species;
  if (idx >= total) return;
  int spec = idx / (total / num_species);
  int local = idx % (total / num_species);
  int iz = local / (owned_nx * global_ny);
  int rem = local % (owned_nx * global_ny);
  int iy = rem / owned_nx;
  int ix = owned_x_begin + rem % owned_nx;
  int cell = iz_cell_index(ix, iy, iz, storage_nx, global_ny);
  double c = conc[spec * ncells + cell] + reac[spec * ncells + cell] * dt;
  if (c < 0.0 && reaction_clip != nullptr) {
    atomicAdd(reaction_clip, -c * cell_volume);
  }
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
                                int num_species, double dt,
                                double* reaction_clip, double cell_volume,
                                int storage_nx, int global_ny, int global_nz,
                                int owned_x_begin, int owned_x_end,
                                cudaStream_t stream) {
  int total = (owned_x_end - owned_x_begin) * global_ny * global_nz
      * num_species;
  int block = 256;
  int grid = (total + block - 1) / block;
  field_update_kernel<<<grid, block, 0, stream>>>(
      conc, reac, ncells, num_species, dt, reaction_clip, cell_volume,
      storage_nx, global_ny, global_nz, owned_x_begin, owned_x_end);
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
