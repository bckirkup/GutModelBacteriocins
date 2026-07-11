#include "gpu_kernels.h"
#include "gpu_types.h"

#include <cuda_runtime.h>
#include <cmath>

namespace gutibm {
namespace gpu {

namespace {

__device__ inline void min_image_delta_device(double dx, double dy, double dz,
                                              double& out_dx, double& out_dy,
                                              double& out_dz,
                                              const MechanicsLaunchParams& p) {
  out_dx = dx;
  out_dy = dy;
  out_dz = dz;
  if (p.periodic_x) {
    const double sx = p.hi0 - p.lo0;
    while (out_dx > 0.5 * sx) out_dx -= sx;
    while (out_dx < -0.5 * sx) out_dx += sx;
  }
  if (p.periodic_y) {
    const double sy = p.hi1 - p.lo1;
    while (out_dy > 0.5 * sy) out_dy -= sy;
    while (out_dy < -0.5 * sy) out_dy += sy;
  }
  if (p.periodic_z) {
    const double sz = p.hi2 - p.lo2;
    while (out_dz > 0.5 * sz) out_dz -= sz;
    while (out_dz < -0.5 * sz) out_dz += sz;
  }
}

__device__ inline void apply_pair_displacement_atomic(
    int i, int j, double nx, double ny, double nz, double force_mag, double dt,
    double mi, double mj, double* dx, double* dy, double* dz) {
  const double inv_mi = 1.0 / mi;
  const double inv_mj = 1.0 / mj;
  const double inv_sum = 1.0 / (inv_mi + inv_mj);
  const double push_i = force_mag * dt * inv_mi * inv_mi * inv_sum;
  const double push_j = force_mag * dt * inv_mj * inv_mj * inv_sum;

  atomicAdd(&dx[i], -nx * push_i);
  atomicAdd(&dy[i], -ny * push_i);
  atomicAdd(&dz[i], -nz * push_i);
  atomicAdd(&dx[j], nx * push_j);
  atomicAdd(&dy[j], ny * push_j);
  atomicAdd(&dz[j], nz * push_j);
}

__device__ inline void apply_adhesion_atomic(
    int i, int j, double nx, double ny, double nz, double adhesion_force,
    double dt, double mi, double mj, double* dx, double* dy, double* dz) {
  const double inv_mi = 1.0 / mi;
  const double inv_mj = 1.0 / mj;
  const double inv_sum = 1.0 / (inv_mi + inv_mj);
  const double pull_i = adhesion_force * dt * inv_mi * inv_mi * inv_sum;
  const double pull_j = adhesion_force * dt * inv_mj * inv_mj * inv_sum;

  atomicAdd(&dx[i], nx * pull_i);
  atomicAdd(&dy[i], ny * pull_i);
  atomicAdd(&dz[i], nz * pull_i);
  atomicAdd(&dx[j], -nx * pull_j);
  atomicAdd(&dy[j], -ny * pull_j);
  atomicAdd(&dz[j], -nz * pull_j);
}

}  // namespace

__global__ void mechanics_clear_kernel(double* dx, double* dy, double* dz,
                                       int num_agents) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  dx[i] = 0.0;
  dy[i] = 0.0;
  dz[i] = 0.0;
}

__global__ void mechanics_forces_kernel(
    const double* x, const double* y, const double* z,
    const double* radius, const double* mass, const int* state,
    const int* cell_offsets, const int* sorted_indices,
    double* dx, double* dy, double* dz,
    int num_agents, MechanicsLaunchParams params) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  if (state[i] == 3) return;

  const double mi = mass[i] < 1.0e-30 ? 1.0e-30 : mass[i];

  int ix = static_cast<int>((x[i] - params.lo0) / params.cell_size);
  int iy = static_cast<int>((y[i] - params.lo1) / params.cell_size);
  int iz = static_cast<int>((z[i] - params.lo2) / params.cell_size);
  ix = ix < 0 ? 0 : (ix >= params.nx_cells ? params.nx_cells - 1 : ix);
  iy = iy < 0 ? 0 : (iy >= params.ny_cells ? params.ny_cells - 1 : iy);
  iz = iz < 0 ? 0 : (iz >= params.nz_cells ? params.nz_cells - 1 : iz);

  for (int dz_cell = -1; dz_cell <= 1; ++dz_cell) {
    int nz = iz + dz_cell;
    if (nz < 0 || nz >= params.nz_cells) continue;
    for (int dy_cell = -1; dy_cell <= 1; ++dy_cell) {
      int ny = iy + dy_cell;
      if (ny < 0 || ny >= params.ny_cells) continue;
      for (int dx_cell = -1; dx_cell <= 1; ++dx_cell) {
        int nx = ix + dx_cell;
        if (nx < 0 || nx >= params.nx_cells) continue;
        const int cell = nz * (params.nx_cells * params.ny_cells)
                         + ny * params.nx_cells + nx;
        const int begin = cell_offsets[cell];
        const int end = cell_offsets[cell + 1];
        for (int pos = begin; pos < end; ++pos) {
          const int j = sorted_indices[pos];
          if (j <= i) continue;
          if (state[j] == 3) continue;

          double delta_x = 0.0;
          double delta_y = 0.0;
          double delta_z = 0.0;
          min_image_delta_device(x[j] - x[i], y[j] - y[i], z[j] - z[i],
                                 delta_x, delta_y, delta_z, params);
          const double d2 = delta_x * delta_x + delta_y * delta_y
                            + delta_z * delta_z;
          if (d2 <= 0.0) continue;
          const double d = sqrt(d2);
          const double sum_r = radius[i] + radius[j];
          const double nx_dir = delta_x / d;
          const double ny_dir = delta_y / d;
          const double nz_dir = delta_z / d;
          const double mj = mass[j] < 1.0e-30 ? 1.0e-30 : mass[j];

          const double overlap = sum_r - d;
          if (overlap > 0.0) {
            const double force_mag = params.hertzian_enabled
                ? params.hertz_k * pow(overlap, 1.5)
                : params.hertz_k * overlap;
            apply_pair_displacement_atomic(
                i, j, nx_dir, ny_dir, nz_dir, force_mag, params.dt,
                mi, mj, dx, dy, dz);
          }

          if (params.adhesion_enabled) {
            const double gap = d - sum_r;
            if (gap > 0.0 && gap < params.adhesion_range) {
              const double adhesion_frac = 1.0 - (gap / params.adhesion_range);
              const double adhesion_force =
                  params.adhesion_strength * adhesion_frac;
              apply_adhesion_atomic(
                  i, j, nx_dir, ny_dir, nz_dir, adhesion_force, params.dt,
                  mi, mj, dx, dy, dz);
            }
          }
        }
      }
    }
  }
}

__global__ void mechanics_apply_kernel(double* x, double* y, double* z,
                                       const double* dx, const double* dy,
                                       const double* dz, int num_agents) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  x[i] += dx[i];
  y[i] += dy[i];
  z[i] += dz[i];
}

void launch_mechanics_clear_kernel(double* dx, double* dy, double* dz,
                                   int num_agents, cudaStream_t stream) {
  if (num_agents <= 0) return;
  int block = 256;
  int grid = (num_agents + block - 1) / block;
  mechanics_clear_kernel<<<grid, block, 0, stream>>>(dx, dy, dz, num_agents);
}

void launch_mechanics_forces_kernel(
    const double* x, const double* y, const double* z,
    const double* radius, const double* mass, const int* state,
    const int* cell_offsets, const int* sorted_indices,
    double* dx, double* dy, double* dz,
    int num_agents, const MechanicsLaunchParams& params,
    cudaStream_t stream) {
  if (num_agents <= 0) return;
  int block = 256;
  int grid = (num_agents + block - 1) / block;
  mechanics_forces_kernel<<<grid, block, 0, stream>>>(
      x, y, z, radius, mass, state,
      cell_offsets, sorted_indices,
      dx, dy, dz, num_agents, params);
}

void launch_mechanics_apply_kernel(
    double* x, double* y, double* z,
    const double* dx, const double* dy, const double* dz,
    int num_agents, cudaStream_t stream) {
  if (num_agents <= 0) return;
  int block = 256;
  int grid = (num_agents + block - 1) / block;
  mechanics_apply_kernel<<<grid, block, 0, stream>>>(
      x, y, z, dx, dy, dz, num_agents);
}

}  // namespace gpu
}  // namespace gutibm
