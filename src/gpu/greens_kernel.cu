#include "gpu_kernels.h"
#include "gpu_common.cuh"
#include <cuda_runtime.h>

namespace gutibm {
namespace gpu {

__global__ void superpose_kernel(
    const double* __restrict__ src_x,
    const double* __restrict__ src_y,
    const double* __restrict__ src_z,
    const GfSourceParams* __restrict__ params,
    double* __restrict__ grid_conc,
    const double* __restrict__ robin_tables,
    DomainParams dom,
    AdvectionParams adv,
    int num_sources,
    int span_x, int span_y, int span_z,
    unsigned long long* cap_hits,
    unsigned long long* kernel_evaluations) {

  int sid = blockIdx.x;
  if (sid >= num_sources) return;

  int local = threadIdx.x + blockIdx.y * blockDim.x;
  int stencil_vol = (2 * span_x + 1) * (2 * span_y + 1) * (2 * span_z + 1);
  if (local >= stencil_vol) return;

  int dz = local / ((2 * span_x + 1) * (2 * span_y + 1)) - span_z;
  int rem = local % ((2 * span_x + 1) * (2 * span_y + 1));
  int dy = rem / (2 * span_x + 1) - span_y;
  int dx = rem % (2 * span_x + 1) - span_x;

  double src[3] = {src_x[sid], src_y[sid], src_z[sid]};

  if (dom.periodic[0] && dx - dom.nx >= -span_x) return;
  if (dom.periodic[1] && dy - dom.ny >= -span_y) return;

  int src_ix = static_cast<int>((src[0] - dom.lo[0]) / dom.dx_x);
  int src_iy = static_cast<int>((src[1] - dom.lo[1]) / dom.dx_y);
  int src_iz = static_cast<int>((src[2] - dom.lo[2]) / dom.dx_z);

  int iz = src_iz + dz;
  if (iz < 0 || iz >= dom.nz) return;

  int iy = src_iy + dy;
  if (dom.periodic[1]) {
    iy = ((iy % dom.ny) + dom.ny) % dom.ny;
  } else if (iy < 0 || iy >= dom.ny) {
    return;
  }

  int ix = src_ix + dx;
  if (dom.periodic[0]) {
    ix = ((ix % dom.nx) + dom.nx) % dom.nx;
  } else if (ix < 0 || ix >= dom.nx) {
    return;
  }

  double tgt[3];
  cell_center(dom, ix, iy, iz, tgt);
  double c = concentration_bounded(
      src, tgt, params[sid], dom, adv, robin_tables, cap_hits,
      kernel_evaluations);
  int idx = cell_index(dom, ix, iy, iz);
  atomicAdd(&grid_conc[idx], c);
}

void launch_superpose_kernel(
    const double* src_x, const double* src_y, const double* src_z,
    const GfSourceParams* params, double* grid_conc,
    const double* robin_tables,
    const DomainParams& dom, const AdvectionParams& adv,
    int num_sources, int span_x, int span_y, int span_z,
    cudaStream_t stream, unsigned long long* cap_hits,
    unsigned long long* kernel_evaluations) {

  if (num_sources == 0) return;
  int stencil_vol = (2 * span_x + 1) * (2 * span_y + 1) * (2 * span_z + 1);
  dim3 block(256);
  dim3 grid(num_sources, (stencil_vol + block.x - 1) / block.x);
  superpose_kernel<<<grid, block, 0, stream>>>(
      src_x, src_y, src_z, params, grid_conc, robin_tables, dom, adv,
      num_sources,
      span_x, span_y, span_z, cap_hits, kernel_evaluations);
}

}  // namespace gpu
}  // namespace gutibm
