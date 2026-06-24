#include "gpu_kernels.h"
#include "gpu_common.cuh"
#include <cuda_runtime.h>

namespace gutibm {
namespace gpu {

__global__ void compute_cell_keys_kernel(
    const double* x, const double* y, const double* z, const int* state,
    int* cell_keys, int* sorted_indices, int num_agents,
    double lo0, double lo1, double lo2, double cell_size,
    int nx_cells, int ny_cells, int nz_cells) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= num_agents) return;
  sorted_indices[i] = i;
  if (state[i] == 3) {
    cell_keys[i] = -1;
    return;
  }
  int ix = static_cast<int>((x[i] - lo0) / cell_size);
  int iy = static_cast<int>((y[i] - lo1) / cell_size);
  int iz = static_cast<int>((z[i] - lo2) / cell_size);
  ix = ix < 0 ? 0 : (ix >= nx_cells ? nx_cells - 1 : ix);
  iy = iy < 0 ? 0 : (iy >= ny_cells ? ny_cells - 1 : iy);
  iz = iz < 0 ? 0 : (iz >= nz_cells ? nz_cells - 1 : iz);
  cell_keys[i] = iz * (nx_cells * ny_cells) + iy * nx_cells + ix;
}

void launch_spatial_hash_build_kernel(
    const double* x, const double* y, const double* z, const int* state,
    int* cell_keys, int* sorted_indices, int num_agents,
    double lo0, double lo1, double lo2, double cell_size,
    int nx_cells, int ny_cells, int nz_cells, cudaStream_t stream) {
  int block = 256;
  int grid = (num_agents + block - 1) / block;
  compute_cell_keys_kernel<<<grid, block, 0, stream>>>(
      x, y, z, state, cell_keys, sorted_indices, num_agents,
      lo0, lo1, lo2, cell_size, nx_cells, ny_cells, nz_cells);
}

}  // namespace gpu
}  // namespace gutibm
