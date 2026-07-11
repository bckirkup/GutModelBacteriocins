#include "gpu_kernels.h"
#include "gpu_common.cuh"
#include <cuda_runtime.h>

namespace gutibm {
namespace gpu {

__device__ double eval_local_expansion(const double* local,
                                       int order,
                                       double cx, double cy, double cz,
                                       double tx, double ty, double tz) {
  const double drx = tx - cx;
  const double dry = ty - cy;
  const double drz = tz - cz;
  double sum = 0.0;
  int idx = 0;
  for (int total = 0; total <= order; ++total) {
    for (int iz = 0; iz <= total; ++iz) {
      for (int iy = 0; iy <= total - iz; ++iy) {
        const int ix = total - iy - iz;
        double term = local[idx++];
        for (int d = 0; d < ix; ++d) term *= drx;
        for (int d = 0; d < iy; ++d) term *= dry;
        for (int d = 0; d < iz; ++d) term *= drz;
        sum += term;
      }
    }
  }
  return sum;
}

__global__ void fmm_far_local_kernel(
    const double* leaf_local,
    const double* leaf_center,
    const int* cell_leaf,
    const double* near_conc,
    double* out_conc,
    int ncells,
    int num_leaves,
    int coeffs_per_leaf,
    int expansion_order,
    double lo0, double lo1, double lo2, double dx,
    int nx, int ny) {
  int cell = blockIdx.x * blockDim.x + threadIdx.x;
  if (cell >= ncells) return;

  const int leaf = cell_leaf[cell];
  if (leaf < 0 || leaf >= num_leaves) return;

  const int iz = cell / (nx * ny);
  const int rem = cell - iz * nx * ny;
  const int iy = rem / nx;
  const int ix = rem - iy * nx;

  const double tx = lo0 + (static_cast<double>(ix) + 0.5) * dx;
  const double ty = lo1 + (static_cast<double>(iy) + 0.5) * dx;
  const double tz = lo2 + (static_cast<double>(iz) + 0.5) * dx;

  const double* local = leaf_local + static_cast<size_t>(leaf) * coeffs_per_leaf;
  const double* center = leaf_center + static_cast<size_t>(leaf) * 3;
  const double far = eval_local_expansion(
      local, expansion_order, center[0], center[1], center[2], tx, ty, tz);
  const double near = near_conc ? near_conc[cell] : 0.0;
  const double contrib = far > near ? far - near : 0.0;
  out_conc[cell] += contrib;
}

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
    double lo0, double lo1, double lo2, double dx,
    int nx, int ny,
    cudaStream_t stream) {
  if (ncells <= 0 || num_leaves <= 0) return;
  int block = 256;
  int grid = (ncells + block - 1) / block;
  fmm_far_local_kernel<<<grid, block, 0, stream>>>(
      leaf_local, leaf_center, cell_leaf, near_conc, out_conc,
      ncells, num_leaves, coeffs_per_leaf, expansion_order,
      lo0, lo1, lo2, dx, nx, ny);
}

}  // namespace gpu
}  // namespace gutibm
