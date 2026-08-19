#include "gpu_kernels.h"
#include <cuda_runtime.h>
#include <cmath>

namespace gutibm {
namespace gpu {

namespace {

constexpr int kMaxLineLength = 1024;

__device__ inline int cell_index(int ix, int iy, int iz, int storage_nx,
                                 int ny) {
  return iz * (storage_nx * ny) + iy * storage_nx + ix;
}

__device__ void pcr_solve(double* s_lower,
                            double* s_diag,
                            double* s_upper,
                            double* s_rhs,
                            int n,
                            int tid) {
  const int n_padded = blockDim.x;

  for (int stride = 1; stride < n_padded; stride <<= 1) {
    double a_new = 0.0;
    double d_new = s_diag[tid];
    double c_new = 0.0;
    double r_new = s_rhs[tid];

    if (tid >= stride && tid < n) {
      const double ratio = s_lower[tid] / s_diag[tid - stride];
      a_new = -ratio * s_lower[tid - stride];
      d_new -= ratio * s_upper[tid - stride];
      r_new -= ratio * s_rhs[tid - stride];
    }
    if (tid + stride < n) {
      const double ratio = s_upper[tid] / s_diag[tid + stride];
      c_new = -ratio * s_upper[tid + stride];
      d_new -= ratio * s_lower[tid + stride];
      r_new -= ratio * s_rhs[tid + stride];
    }

    __syncthreads();

    if (tid < n) {
      s_lower[tid] = a_new;
      s_diag[tid] = d_new;
      s_upper[tid] = c_new;
      s_rhs[tid] = r_new;
    }

    __syncthreads();
  }

  if (tid < n) {
    s_rhs[tid] /= s_diag[tid];
  }
}

__device__ void solve_periodic_line_n2(double* line, double alpha) {
  const double diagonal = 1.0 + 2.0 * alpha;
  const double off = 2.0 * alpha;
  const double det = 1.0 + 4.0 * alpha;
  const double first = (diagonal * line[0] + off * line[1]) / det;
  const double second = (off * line[0] + diagonal * line[1]) / det;
  line[0] = first;
  line[1] = second;
}

__device__ void solve_periodic_line(double* line,
                                    double* pcr_base,
                                    int n,
                                    double alpha,
                                    double gamma,
                                    double corner,
                                    double denominator,
                                    const double* correction,
                                    int tid) {
  if (n <= 1) return;
  if (n == 2) {
    if (tid == 0) solve_periodic_line_n2(line, alpha);
    return;
  }

  double* s_lower = pcr_base;
  double* s_diag = pcr_base + blockDim.x;
  double* s_upper = s_diag + blockDim.x;
  double* s_rhs = s_upper + blockDim.x;

  if (tid < n) {
    s_lower[tid] = (tid == 0) ? 0.0 : -alpha;
    s_upper[tid] = (tid == n - 1) ? 0.0 : -alpha;
    const double diagonal_value = 1.0 + 2.0 * alpha;
    s_diag[tid] = diagonal_value;
    if (tid == 0) s_diag[tid] -= gamma;
    if (tid == n - 1) s_diag[tid] -= corner * corner / gamma;
    s_rhs[tid] = line[tid];
  }

  __syncthreads();
  pcr_solve(s_lower, s_diag, s_upper, s_rhs, n, tid);
  __syncthreads();

  if (tid < n) {
    const double numerator = s_rhs[0] + corner * s_rhs[n - 1] / gamma;
    const double adjustment = numerator / denominator;
    line[tid] = s_rhs[tid] - adjustment * correction[tid];
  }
}

__device__ void solve_bounded_line(double* line,
                                   double* pcr_base,
                                   int n,
                                   double alpha,
                                   double boundary_conc,
                                   int tid) {
  double* s_lower = pcr_base;
  double* s_diag = pcr_base + blockDim.x;
  double* s_upper = s_diag + blockDim.x;
  double* s_rhs = s_upper + blockDim.x;

  if (tid < n) {
    s_lower[tid] = (tid == 0) ? 0.0 : -alpha;
    s_upper[tid] = (tid == n - 1) ? 0.0 : -alpha;
    s_diag[tid] = 1.0 + 2.0 * alpha;
    if (tid == n - 1) s_diag[tid] = 1.0 + alpha;
    s_rhs[tid] = line[tid];
    if (tid == 0) s_rhs[tid] += alpha * boundary_conc;
  }

  __syncthreads();
  pcr_solve(s_lower, s_diag, s_upper, s_rhs, n, tid);
  __syncthreads();

  if (tid < n) line[tid] = s_rhs[tid];
}

__device__ void solve_delivery_line(double* line, double* pcr_base, int n,
                                    double alpha, int boundary_mode,
                                    double boundary_conc, double beta,
                                    double flux_source, int tid) {
  double* s_lower = pcr_base;
  double* s_diag = pcr_base + blockDim.x;
  double* s_upper = s_diag + blockDim.x;
  double* s_rhs = s_upper + blockDim.x;

  if (tid < n) {
    s_lower[tid] = (tid == 0) ? 0.0 : -alpha;
    s_upper[tid] = (tid == n - 1) ? 0.0 : -alpha;
    s_diag[tid] = 1.0 + 2.0 * alpha;
    if (tid == 0) {
      s_diag[tid] = boundary_mode == 1
          ? 1.0 + alpha + beta : 1.0 + alpha;
    } else if (tid == n - 1) {
      s_diag[tid] = 1.0 + alpha;
    }
    s_rhs[tid] = line[tid];
    if (tid == 0) {
      s_rhs[tid] += boundary_mode == 1
          ? beta * boundary_conc : flux_source;
    }
  }

  __syncthreads();
  pcr_solve(s_lower, s_diag, s_upper, s_rhs, n, tid);
  __syncthreads();

  if (tid < n) line[tid] = s_rhs[tid];
}

__global__ void diffuse_x_periodic_kernel(double* conc,
                                          int nx,
                                          int ny,
                                          int nz,
                                          double alpha,
                                          double gamma,
                                          double corner,
                                          double denominator,
                                          const double* correction) {
  const int line_id = blockIdx.x;
  if (line_id >= ny * nz) return;
  const int iy = line_id / nz;
  const int iz = line_id % nz;
  const int tid = threadIdx.x;

  extern __shared__ double smem[];
  double* line = smem;
  double* pcr_base = smem + blockDim.x;

  if (tid < nx) {
    line[tid] = conc[cell_index(tid, iy, iz, nx, ny)];
  }
  __syncthreads();

  solve_periodic_line(line, pcr_base, nx, alpha, gamma, corner, denominator,
                      correction, tid);
  __syncthreads();

  if (tid < nx) {
    conc[cell_index(tid, iy, iz, nx, ny)] = line[tid];
  }
}

__global__ void diffuse_y_periodic_kernel(double* conc,
                                          int storage_nx,
                                          int ny,
                                          int nz,
                                          int owned_x_begin,
                                          int owned_x_end,
                                          double alpha,
                                          double gamma,
                                          double corner,
                                          double denominator,
                                          const double* correction) {
  const int line_id = blockIdx.x;
  const int local_nx = owned_x_end - owned_x_begin;
  if (line_id >= local_nx * nz) return;
  const int ix = owned_x_begin + line_id / nz;
  const int iz = line_id % nz;
  const int tid = threadIdx.x;

  extern __shared__ double smem[];
  double* line = smem;
  double* pcr_base = smem + blockDim.x;

  if (tid < ny) {
    line[tid] = conc[cell_index(ix, tid, iz, storage_nx, ny)];
  }
  __syncthreads();

  solve_periodic_line(line, pcr_base, ny, alpha, gamma, corner, denominator,
                      correction, tid);
  __syncthreads();

  if (tid < ny) {
    conc[cell_index(ix, tid, iz, storage_nx, ny)] = line[tid];
  }
}

__global__ void diffuse_z_bounded_kernel(double* conc,
                                         int storage_nx,
                                         int ny,
                                         int nz,
                                         int owned_x_begin,
                                         int owned_x_end,
                                         double alpha,
                                         double boundary_conc,
                                         int boundary_mode,
                                         double beta,
                                         double flux_source,
                                         double cell_volume,
                                         double* face_exchange) {
  const int line_id = blockIdx.x;
  const int local_nx = owned_x_end - owned_x_begin;
  if (line_id >= local_nx * ny) return;
  const int ix = owned_x_begin + line_id / ny;
  const int iy = line_id % ny;
  const int n = boundary_mode == 0 ? nz - 1 : nz;
  if (n <= 0) return;
  const int tid = threadIdx.x;

  extern __shared__ double smem[];
  double* line = smem;
  double* pcr_base = smem + blockDim.x;

  if (tid < n) {
    const int iz = boundary_mode == 0 ? tid + 1 : tid;
    line[tid] = conc[cell_index(ix, iy, iz, storage_nx, ny)];
  }
  __syncthreads();

  if (boundary_mode == 0) {
    solve_bounded_line(line, pcr_base, n, alpha, boundary_conc, tid);
  } else {
    solve_delivery_line(line, pcr_base, n, alpha, boundary_mode,
                        boundary_conc, beta, flux_source, tid);
  }
  __syncthreads();

  if (tid == 0 && face_exchange != nullptr) {
    const double realized = boundary_mode == 1
        ? beta * (boundary_conc - line[0]) * cell_volume
        : flux_source * cell_volume;
    atomicAdd(face_exchange, boundary_mode == 0
        ? alpha * (boundary_conc - line[0]) * cell_volume : realized);
  }

  if (tid < n) {
    const int iz = boundary_mode == 0 ? tid + 1 : tid;
    conc[cell_index(ix, iy, iz, storage_nx, ny)] = line[tid];
  }
}

__global__ void set_epithelial_boundary_kernel(double* conc,
                                               int storage_nx,
                                               int ny,
                                               int owned_x_begin,
                                               int owned_x_end,
                                               double boundary_conc,
                                               double cell_volume,
                                               double* injected_amount) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int local_nx = owned_x_end - owned_x_begin;
  const int face = local_nx * ny;
  if (idx >= face) return;
  const int ix = owned_x_begin + idx % local_nx;
  const int iy = idx / local_nx;
  const int cell = cell_index(ix, iy, 0, storage_nx, ny);
  if (injected_amount != nullptr) {
    atomicAdd(injected_amount, (boundary_conc - conc[cell]) * cell_volume);
  }
  conc[cell] = boundary_conc;
}

__global__ void set_luminal_neumann_kernel(double* conc, int storage_nx,
                                           int ny, int nz,
                                           int owned_x_begin,
                                           int owned_x_end) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int local_nx = owned_x_end - owned_x_begin;
  const int face = local_nx * ny;
  if (idx >= face || nz < 2) return;
  const int ix = owned_x_begin + idx % local_nx;
  const int iy = idx / local_nx;
  const int top = cell_index(ix, iy, nz - 1, storage_nx, ny);
  const int below = cell_index(ix, iy, nz - 2, storage_nx, ny);
  conc[top] = conc[below];
}

__global__ void shift_z_gradient_kernel(double* conc,
                                        int storage_nx,
                                        int ny,
                                        int nz,
                                        int owned_x_begin,
                                        int owned_x_end,
                                        double dx_z,
                                        double initial_conc,
                                        double lambda,
                                        double boundary_conc,
                                        double scale) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int local_nx = owned_x_end - owned_x_begin;
  const int ncells = local_nx * ny * nz;
  if (idx >= ncells) return;

  const int iz = idx / (local_nx * ny);
  const int rem = idx % (local_nx * ny);
  const int iy = rem / local_nx;
  const int ix = owned_x_begin + rem % local_nx;
  double shift = 0.0;
  if (iz == 0) {
    shift = scale * boundary_conc;
  } else {
    const int profile_iz = (nz >= 2 && iz == nz - 1) ? nz - 2 : iz;
    const double z_rel = (profile_iz + 0.5) * dx_z;
    shift = scale * initial_conc * exp(-z_rel / lambda);
  }
  conc[cell_index(ix, iy, iz, storage_nx, ny)] += shift;
}

__global__ void clamp_nonneg_kernel(double* conc, int storage_nx, int ny,
                                     int nz, int owned_x_begin,
                                     int owned_x_end) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int local_nx = owned_x_end - owned_x_begin;
  const int ncells = local_nx * ny * nz;
  if (idx >= ncells) return;
  const int iz = idx / (local_nx * ny);
  const int rem = idx % (local_nx * ny);
  const int iy = rem / local_nx;
  const int ix = owned_x_begin + rem % local_nx;
  conc[cell_index(ix, iy, iz, storage_nx, ny)] =
      fmax(conc[cell_index(ix, iy, iz, storage_nx, ny)], 0.0);
}

int next_pow2(int n) {
  int p = 1;
  while (p < n) p <<= 1;
  return p;
}

size_t periodic_smem_bytes(int line_len) {
  const int padded = next_pow2(line_len);
  return static_cast<size_t>(5 * padded) * sizeof(double);
}

}  // namespace

void launch_diffuse_x_periodic(double* conc,
                               int nx,
                               int ny,
                               int nz,
                               double alpha,
                               double gamma,
                               double corner,
                               double denominator,
                               const double* correction,
                               cudaStream_t stream) {
  if (nx <= 0 || ny <= 0 || nz <= 0) return;
  const int block = next_pow2(nx);
  const int grid = ny * nz;
  const size_t smem = periodic_smem_bytes(nx);
  diffuse_x_periodic_kernel<<<grid, block, smem, stream>>>(
      conc, nx, ny, nz, alpha, gamma, corner, denominator, correction);
}

void launch_diffuse_y_periodic(double* conc,
                               int storage_nx,
                               int ny,
                               int nz,
                               int owned_x_begin,
                               int owned_x_end,
                               double alpha,
                               double gamma,
                               double corner,
                               double denominator,
                               const double* correction,
                               cudaStream_t stream) {
  const int local_nx = owned_x_end - owned_x_begin;
  if (local_nx <= 0 || ny <= 0 || nz <= 0) return;
  const int block = next_pow2(ny);
  const int grid = local_nx * nz;
  const size_t smem = periodic_smem_bytes(ny);
  diffuse_y_periodic_kernel<<<grid, block, smem, stream>>>(
      conc, storage_nx, ny, nz, owned_x_begin, owned_x_end, alpha, gamma,
      corner, denominator, correction);
}

void launch_diffuse_z_bounded(double* conc,
                              int storage_nx,
                              int ny,
                              int nz,
                              int owned_x_begin,
                              int owned_x_end,
                              double alpha,
                              double boundary_conc,
                              int boundary_mode,
                              double beta,
                              double flux_source,
                              double cell_volume,
                              double* face_exchange,
                              cudaStream_t stream) {
  const int local_nx = owned_x_end - owned_x_begin;
  if (local_nx <= 0 || ny <= 0 || nz <= 1) return;
  const int n = boundary_mode == 0 ? nz - 1 : nz;
  const int block = next_pow2(n);
  const int grid = local_nx * ny;
  const size_t smem = static_cast<size_t>(5 * block) * sizeof(double);
  diffuse_z_bounded_kernel<<<grid, block, smem, stream>>>(
      conc, storage_nx, ny, nz, owned_x_begin, owned_x_end, alpha,
      boundary_conc, boundary_mode, beta, flux_source, cell_volume,
      face_exchange);
}

void launch_set_epithelial_boundary(double* conc,
                                    int storage_nx,
                                    int ny,
                                    int owned_x_begin,
                                    int owned_x_end,
                                    double boundary_conc,
                                    double cell_volume,
                                    double* injected_amount,
                                    cudaStream_t stream) {
  const int face = (owned_x_end - owned_x_begin) * ny;
  if (face <= 0) return;
  const int block = 256;
  const int grid = (face + block - 1) / block;
  set_epithelial_boundary_kernel<<<grid, block, 0, stream>>>(
      conc, storage_nx, ny, owned_x_begin, owned_x_end, boundary_conc,
      cell_volume, injected_amount);
}

void launch_set_luminal_neumann(double* conc,
                                int storage_nx,
                                int ny,
                                int nz,
                                int owned_x_begin,
                                int owned_x_end,
                                cudaStream_t stream) {
  const int face = (owned_x_end - owned_x_begin) * ny;
  if (face <= 0 || nz < 2) return;
  const int block = 256;
  const int grid = (face + block - 1) / block;
  set_luminal_neumann_kernel<<<grid, block, 0, stream>>>(
      conc, storage_nx, ny, nz, owned_x_begin, owned_x_end);
}

void launch_shift_z_gradient(double* conc,
                             int storage_nx,
                             int ny,
                             int nz,
                             int owned_x_begin,
                             int owned_x_end,
                             double dx_z,
                             double initial_conc,
                             double lambda,
                             double boundary_conc,
                             double scale,
                             cudaStream_t stream) {
  const int ncells = (owned_x_end - owned_x_begin) * ny * nz;
  if (ncells <= 0) return;
  const int block = 256;
  const int grid = (ncells + block - 1) / block;
  shift_z_gradient_kernel<<<grid, block, 0, stream>>>(
      conc, storage_nx, ny, nz, owned_x_begin, owned_x_end, dx_z, initial_conc,
      lambda, boundary_conc, scale);
}

void launch_clamp_nonneg(double* conc, int storage_nx, int ny, int nz,
                         int owned_x_begin, int owned_x_end,
                         cudaStream_t stream) {
  const int ncells = (owned_x_end - owned_x_begin) * ny * nz;
  if (ncells <= 0) return;
  const int block = 256;
  const int grid = (ncells + block - 1) / block;
  clamp_nonneg_kernel<<<grid, block, 0, stream>>>(
      conc, storage_nx, ny, nz, owned_x_begin, owned_x_end);
}

int diffusion_max_line_length() { return kMaxLineLength; }

}  // namespace gpu
}  // namespace gutibm
