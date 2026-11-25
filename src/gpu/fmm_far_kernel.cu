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
    double lo0, double lo1, double lo2, double dx_x, double dx_y, double dx_z,
    int nx, int ny) {
  int cell = blockIdx.x * blockDim.x + threadIdx.x;
  if (cell >= ncells) return;

  const int leaf = cell_leaf[cell];
  if (leaf < 0 || leaf >= num_leaves) return;

  const int iz = cell / (nx * ny);
  const int rem = cell - iz * nx * ny;
  const int iy = rem / nx;
  const int ix = rem - iy * nx;

  const double tx = lo0 + (static_cast<double>(ix) + 0.5) * dx_x;
  const double ty = lo1 + (static_cast<double>(iy) + 0.5) * dx_y;
  const double tz = lo2 + (static_cast<double>(iz) + 0.5) * dx_z;

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
    double lo0, double lo1, double lo2, double dx_x, double dx_y, double dx_z,
    int nx, int ny,
    cudaStream_t stream) {
  if (ncells <= 0 || num_leaves <= 0) return;
  int block = 256;
  int grid = (ncells + block - 1) / block;
  fmm_far_local_kernel<<<grid, block, 0, stream>>>(
      leaf_local, leaf_center, cell_leaf, near_conc, out_conc,
      ncells, num_leaves, coeffs_per_leaf, expansion_order,
      lo0, lo1, lo2, dx_x, dx_y, dx_z, nx, ny);
}

namespace {

constexpr double kSqrt3Gpu = 1.732050807568877293527446341505872367;
constexpr int kFmmTraversalStackCap = 128;

__device__ inline bool decode_owned_voxel(int tid,
                                          const FmmGridLayout& grid,
                                          const DomainParams& dom,
                                          double tgt[3],
                                          int& storage_cell) {
  const int owned_nx = grid.owned_storage_x_end - grid.owned_storage_x_begin;
  const int total_owned = owned_nx * grid.global_ny * grid.global_nz;
  if (tid >= total_owned) return false;

  const int plane = owned_nx * grid.global_ny;
  const int iz = tid / plane;
  const int rem = tid - iz * plane;
  const int iy = rem / owned_nx;
  const int local_offset = rem - iy * owned_nx;
  const int ix_storage = grid.owned_storage_x_begin + local_offset;
  const int ix_global = grid.local_global_x_begin + local_offset;

  tgt[0] = dom.lo[0] + (static_cast<double>(ix_global) + 0.5) * dom.dx_x;
  tgt[1] = dom.lo[1] + (static_cast<double>(iy) + 0.5) * dom.dx_y;
  tgt[2] = dom.lo[2] + (static_cast<double>(iz) + 0.5) * dom.dx_z;
  storage_cell = iz * (grid.storage_nx * grid.global_ny)
               + iy * grid.storage_nx + ix_storage;
  return true;
}

__device__ inline int find_containing_node_device(const FmmDeviceTreeView& tree,
                                                  const double tgt[3]) {
  int node_idx = 0;
  while (tree.is_leaf[node_idx] == 0) {
    const double* c = tree.center + static_cast<size_t>(node_idx) * 3;
    int oct = 0;
    if (tgt[0] >= c[0]) oct |= 1;
    if (tgt[1] >= c[1]) oct |= 2;
    if (tgt[2] >= c[2]) oct |= 4;
    const int child = tree.children[static_cast<size_t>(node_idx) * 8 + oct];
    if (child < 0) break;
    node_idx = child;
  }
  return node_idx;
}

__device__ inline double eval_gf_unit_device(const double src[3],
                                             const double tgt[3],
                                             const FmmFarLaunchParams& p) {
  GfSourceParams unit = p.avg_params;
  unit.source_rate = 1.0;
  if (robin::transfer_enabled(unit.lumen_transfer_length)) {
    const double d_eff = unit.diff_coeff / unit.retardation;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    flow_velocity(p.adv, src[0], src[2], vx, vy, vz);
    const double threshold = unit.peclet_z_threshold >= 0.0
        ? unit.peclet_z_threshold : neumann::kDriftEnvelopePeZ;
    const bool drift_exceeded = neumann::drift_envelope_exceeded(
        vz, p.dom.z_hi - p.dom.z_lo, d_eff, threshold);
    const bool use_drift = (unit.drift_correction != 0)
        || (unit.drift_auto_dispatch != 0 && fabs(vz) > 0.0 && drift_exceeded);
    unit.robin_table_index = use_drift
        ? p.robin_table_index_physical : p.robin_table_index_legacy;
  }
  return concentration_bounded(
      src, tgt, unit, p.dom, p.adv, p.robin_tables, p.robin_table_count,
      nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
}

__device__ inline double stencil_margin_device(int order,
                                                 double dist,
                                                 const FmmFarLaunchParams& p) {
  const double length_scale = dist > 1.0e-9 ? dist : 1.0e-9;
  if (order >= 3) return 3.0 * length_scale * p.fd_rel_3;
  if (order == 2) return 2.0 * length_scale * p.fd_rel_2;
  return length_scale * p.fd_rel_1;
}

__device__ inline bool crosses_periodic_branch_cut_device(
    double dx, double dy, double dz, double pad, const DomainParams& dom) {
  if (dom.periodic[0] && (fabs(dx) + pad >= 0.5 * dom.extent[0])) return true;
  if (dom.periodic[1] && (fabs(dy) + pad >= 0.5 * dom.extent[1])) return true;
  if (dom.periodic[2] && (fabs(dz) + pad >= 0.5 * dom.extent[2])) return true;
  return false;
}

__device__ inline int pos_to_grid_axis_device(double pos,
                                              double lo,
                                              double dx,
                                              int n) {
  int idx = static_cast<int>(floor((pos - lo) / dx));
  if (idx < 0) idx = 0;
  if (idx >= n) idx = n - 1;
  return idx;
}

__device__ inline int periodic_cell_dist_device(int a, int b, int n, bool per) {
  const int d = a >= b ? (a - b) : (b - a);
  if (!per) return d;
  const int wrap = n - d;
  return wrap < d ? wrap : d;
}

__device__ inline bool is_source_in_near_region_device(
    const double src[3],
    const double tgt[3],
    const FmmFarLaunchParams& p) {
  if (p.near_cutoff <= 0.0) return false;
  if (p.cutoff_shape == 0) {
    const double dx = minimum_image_delta(
        tgt[0] - src[0], p.dom.extent[0], p.dom.periodic[0]);
    const double dy = minimum_image_delta(
        tgt[1] - src[1], p.dom.extent[1], p.dom.periodic[1]);
    const double dz = minimum_image_delta(
        tgt[2] - src[2], p.dom.extent[2], p.dom.periodic[2]);
    return (dx * dx + dy * dy + dz * dz) <= (p.near_cutoff * p.near_cutoff);
  }
  const int sx = pos_to_grid_axis_device(src[0], p.dom.lo[0], p.dom.dx_x, p.dom.nx);
  const int sy = pos_to_grid_axis_device(src[1], p.dom.lo[1], p.dom.dx_y, p.dom.ny);
  const int sz = pos_to_grid_axis_device(src[2], p.dom.lo[2], p.dom.dx_z, p.dom.nz);
  const int tx = pos_to_grid_axis_device(tgt[0], p.dom.lo[0], p.dom.dx_x, p.dom.nx);
  const int ty = pos_to_grid_axis_device(tgt[1], p.dom.lo[1], p.dom.dx_y, p.dom.ny);
  const int tz = pos_to_grid_axis_device(tgt[2], p.dom.lo[2], p.dom.dx_z, p.dom.nz);
  return periodic_cell_dist_device(tx, sx, p.dom.nx, p.dom.periodic[0]) <= p.span_x
      && periodic_cell_dist_device(ty, sy, p.dom.ny, p.dom.periodic[1]) <= p.span_y
      && periodic_cell_dist_device(tz, sz, p.dom.nz, p.dom.periodic[2]) <= p.span_z;
}

__device__ inline bool pair_completely_inside_near_device(
    double bdx, double bdy, double bdz,
    double half_size,
    const FmmFarLaunchParams& p) {
  if (p.near_cutoff <= 0.0) return false;
  if (p.cutoff_shape == 0) {
    const double box_dist = sqrt(bdx * bdx + bdy * bdy + bdz * bdz);
    return box_dist + half_size * kSqrt3Gpu <= p.near_cutoff;
  }
  if (!(p.dom.periodic[0] && (2 * p.span_x + 1 >= p.dom.nx))) {
    if (fabs(bdx) + half_size + p.dom.dx_x
        > static_cast<double>(p.span_x) * p.dom.dx_x) {
      return false;
    }
  }
  if (!(p.dom.periodic[1] && (2 * p.span_y + 1 >= p.dom.ny))) {
    if (fabs(bdy) + half_size + p.dom.dx_y
        > static_cast<double>(p.span_y) * p.dom.dx_y) {
      return false;
    }
  }
  if (!(p.dom.periodic[2] && (2 * p.span_z + 1 >= p.dom.nz))) {
    if (fabs(bdz) + half_size + p.dom.dx_z
        > static_cast<double>(p.span_z) * p.dom.dx_z) {
      return false;
    }
  }
  return true;
}

__device__ inline bool pair_completely_outside_near_device(
    double bdx, double bdy, double bdz,
    double half_size,
    const FmmFarLaunchParams& p) {
  if (p.near_cutoff <= 0.0) return true;
  if (p.cutoff_shape == 0) {
    const double box_dist = sqrt(bdx * bdx + bdy * bdy + bdz * bdz);
    return box_dist - half_size * kSqrt3Gpu > p.near_cutoff;
  }
  if (!(p.dom.periodic[0] && (2 * p.span_x + 1 >= p.dom.nx))) {
    if (fabs(bdx) - half_size
        > static_cast<double>(p.span_x + 1) * p.dom.dx_x) {
      return true;
    }
  }
  if (!(p.dom.periodic[1] && (2 * p.span_y + 1 >= p.dom.ny))) {
    if (fabs(bdy) - half_size
        > static_cast<double>(p.span_y + 1) * p.dom.dx_y) {
      return true;
    }
  }
  if (!(p.dom.periodic[2] && (2 * p.span_z + 1 >= p.dom.nz))) {
    if (fabs(bdz) - half_size
        > static_cast<double>(p.span_z + 1) * p.dom.dx_z) {
      return true;
    }
  }
  return false;
}

__device__ double eval_multipole_node_device(const double* moments,
                                             int order,
                                             const double src[3],
                                             const double tgt[3],
                                             double dist,
                                             const FmmFarLaunchParams& p) {
  const double length_scale = dist > 1.0e-9 ? dist : 1.0e-9;
  const double g0 = eval_gf_unit_device(src, tgt, p);
  double sum = moments[0] * g0;
  if (order < 1) return sum;

  // Order 1: (1,0,0)->1, (0,1,0)->2, (0,0,1)->3.
  // Moments centered at center_of_source have identically zero dipoles.
  const double h1 = length_scale * p.fd_rel_1;
  const double inv_2h1 = 1.0 / (2.0 * h1);
  for (int d = 0; d < 3; ++d) {
    const double m1 = moments[1 + d];
    if (fabs(m1) <= 0.0) continue;
    double plus[3] = {src[0], src[1], src[2]};
    double minus[3] = {src[0], src[1], src[2]};
    plus[d] += h1;
    minus[d] -= h1;
    const double deriv = (eval_gf_unit_device(plus, tgt, p)
                          - eval_gf_unit_device(minus, tgt, p)) * inv_2h1;
    sum += m1 * deriv;
  }
  if (order < 2) return sum;

  // Order 2: (2,0,0)->4, (1,1,0)->5, (0,2,0)->6, (1,0,1)->7, (0,1,1)->8, (0,0,2)->9.
  // Diagonal uses 5-point 4th-order central stencil with Taylor factor 1/2! = 0.5
  // (combined divisor 24 * h2^2); mixed uses 4-point stencil at h2 = L * eps^(1/6).
  const double h2 = length_scale * p.fd_rel_2;
  const double inv_h2_sq = 1.0 / (h2 * h2);
  const double inv_24h2_sq = inv_h2_sq / 24.0;
  const double inv_4h2_sq = 0.25 * inv_h2_sq;
  const int diag_idx[3] = {4, 6, 9};
  for (int d = 0; d < 3; ++d) {
    const double m2 = moments[diag_idx[d]];
    if (fabs(m2) <= 0.0) continue;
    double p2[3] = {src[0], src[1], src[2]};
    double p1[3] = {src[0], src[1], src[2]};
    double m1[3] = {src[0], src[1], src[2]};
    double m2p[3] = {src[0], src[1], src[2]};
    p2[d] += 2.0 * h2;
    p1[d] += h2;
    m1[d] -= h2;
    m2p[d] -= 2.0 * h2;
    const double deriv = (-eval_gf_unit_device(p2, tgt, p)
                          + 16.0 * eval_gf_unit_device(p1, tgt, p)
                          - 30.0 * g0
                          + 16.0 * eval_gf_unit_device(m1, tgt, p)
                          - eval_gf_unit_device(m2p, tgt, p)) * inv_24h2_sq;
    sum += m2 * deriv;
  }
  const int mixed_d1[3] = {0, 0, 1};
  const int mixed_d2[3] = {1, 2, 2};
  const int mixed_idx[3] = {5, 7, 8};
  for (int m = 0; m < 3; ++m) {
    const double m2 = moments[mixed_idx[m]];
    if (fabs(m2) <= 0.0) continue;
    const int d1 = mixed_d1[m];
    const int d2 = mixed_d2[m];
    double pp[3] = {src[0], src[1], src[2]};
    double pm[3] = {src[0], src[1], src[2]};
    double mp[3] = {src[0], src[1], src[2]};
    double mm[3] = {src[0], src[1], src[2]};
    pp[d1] += h2; pp[d2] += h2;
    pm[d1] += h2; pm[d2] -= h2;
    mp[d1] -= h2; mp[d2] += h2;
    mm[d1] -= h2; mm[d2] -= h2;
    const double deriv = (eval_gf_unit_device(pp, tgt, p)
                          - eval_gf_unit_device(pm, tgt, p)
                          - eval_gf_unit_device(mp, tgt, p)
                          + eval_gf_unit_device(mm, tgt, p)) * inv_4h2_sq;
    sum += m2 * deriv;
  }
  if (order < 3) return sum;

  // Order 3: pure single-axis 6-point 4th-order stencil (3,0,0)->10,
  // (0,3,0)->13, (0,0,3)->19, with Taylor factor 1/3! = 1/6 (divisor 48 * h3^3).
  const double h3 = length_scale * p.fd_rel_3;
  const double inv_48h3_cu = 1.0 / (48.0 * h3 * h3 * h3);
  const int third_idx[3] = {10, 13, 19};
  for (int d = 0; d < 3; ++d) {
    const double m3 = moments[third_idx[d]];
    if (fabs(m3) <= 0.0) continue;
    double p3[3] = {src[0], src[1], src[2]};
    double p2[3] = {src[0], src[1], src[2]};
    double p1[3] = {src[0], src[1], src[2]};
    double m1[3] = {src[0], src[1], src[2]};
    double m2[3] = {src[0], src[1], src[2]};
    double m3p[3] = {src[0], src[1], src[2]};
    p3[d] += 3.0 * h3;
    p2[d] += 2.0 * h3;
    p1[d] += h3;
    m1[d] -= h3;
    m2[d] -= 2.0 * h3;
    m3p[d] -= 3.0 * h3;
    const double deriv = (-eval_gf_unit_device(p3, tgt, p)
                          + 8.0 * eval_gf_unit_device(p2, tgt, p)
                          - 13.0 * eval_gf_unit_device(p1, tgt, p)
                          + 13.0 * eval_gf_unit_device(m1, tgt, p)
                          - 8.0 * eval_gf_unit_device(m2, tgt, p)
                          + eval_gf_unit_device(m3p, tgt, p)) * inv_48h3_cu;
    sum += m3 * deriv;
  }
  return sum;
}

__device__ void traverse_far_device(int root_node,
                                    const FmmDeviceTreeView& tree,
                                    const double tgt[3],
                                    const FmmFarLaunchParams& p,
                                    double& accumulator) {
  int stack[kFmmTraversalStackCap];
  int sp = 0;
  stack[sp++] = root_node;

  const double near_cutoff = p.near_cutoff;

  while (sp > 0) {
    const int node_idx = stack[--sp];
    if (tree.total_source_strength[node_idx] <= 0.0) continue;

    const double hs = tree.half_size[node_idx];
    const double* box_pos = tree.center + static_cast<size_t>(node_idx) * 3;
    const double bdx = minimum_image_delta(
        tgt[0] - box_pos[0], p.dom.extent[0], p.dom.periodic[0]);
    const double bdy = minimum_image_delta(
        tgt[1] - box_pos[1], p.dom.extent[1], p.dom.periodic[1]);
    const double bdz = minimum_image_delta(
        tgt[2] - box_pos[2], p.dom.extent[2], p.dom.periodic[2]);

    if (pair_completely_inside_near_device(bdx, bdy, bdz, hs, p)) continue;

    const double* cos_pos = tree.center_of_source + static_cast<size_t>(node_idx) * 3;
    const double dx = minimum_image_delta(
        tgt[0] - cos_pos[0], p.dom.extent[0], p.dom.periodic[0]);
    const double dy = minimum_image_delta(
        tgt[1] - cos_pos[1], p.dom.extent[1], p.dom.periodic[1]);
    const double dz = minimum_image_delta(
        tgt[2] - cos_pos[2], p.dom.extent[2], p.dom.periodic[2]);
    const double dist = sqrt(dx * dx + dy * dy + dz * dz);

    const bool straddles_near = (near_cutoff > 0.0)
        && !pair_completely_outside_near_device(bdx, bdy, bdz, hs, p);
    const double pad = stencil_margin_device(tree.expansion_order, dist, p);
    const bool is_leaf = (tree.is_leaf[node_idx] != 0);
    const bool straddles_branch =
        crosses_periodic_branch_cut_device(dx, dy, dz, pad, p.dom);

    if (!is_leaf) {
      const double cell_size = 2.0 * hs;
      const double denom = dist > 1.0e-30 ? dist : 1.0e-30;
      const bool fails_mac = (cell_size / denom) >= p.theta;
      if (fails_mac || straddles_near || straddles_branch) {
        const int* children = tree.children + static_cast<size_t>(node_idx) * 8;
        for (int c = 7; c >= 0; --c) {
          if (children[c] >= 0 && sp < kFmmTraversalStackCap) {
            stack[sp++] = children[c];
          }
        }
        continue;
      }
    } else if (straddles_near || straddles_branch) {
      if (tree.leaf_src_offsets != nullptr && tree.leaf_src_pos != nullptr
          && tree.leaf_src_strength != nullptr) {
        const int sb = tree.leaf_src_offsets[node_idx];
        const int se = tree.leaf_src_offsets[node_idx + 1];
        for (int s = sb; s < se; ++s) {
          const double* sp_pos = tree.leaf_src_pos + static_cast<size_t>(s) * 3;
          if (is_source_in_near_region_device(sp_pos, tgt, p)) continue;
          const double q = tree.leaf_src_strength[s];
          if (q <= 0.0) continue;
          accumulator += q * eval_gf_unit_device(sp_pos, tgt, p);
        }
      }
      continue;
    }

    const double* moments = tree.multipole
        + static_cast<size_t>(node_idx) * tree.coeffs_per_node;
    accumulator += eval_multipole_node_device(
        moments, tree.expansion_order, cos_pos, tgt, dist, p);
  }
}

__global__ void fmm_far_local_device_kernel(
    FmmDeviceTreeView tree,
    FmmFarLaunchParams params,
    double* d_conc) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  double tgt[3];
  int storage_cell = 0;
  if (!decode_owned_voxel(tid, params.grid, params.dom, tgt, storage_cell)) {
    return;
  }

  const int node_idx = find_containing_node_device(tree, tgt);
  const double* local = tree.local
      + static_cast<size_t>(node_idx) * tree.coeffs_per_node;
  const double* center = tree.center + static_cast<size_t>(node_idx) * 3;
  double far = eval_local_expansion(
      local, tree.expansion_order, center[0], center[1], center[2],
      tgt[0], tgt[1], tgt[2]);

  if (tree.near_offsets != nullptr && tree.near_nodes != nullptr) {
    const int nb = tree.near_offsets[node_idx];
    const int ne = tree.near_offsets[node_idx + 1];
    for (int k = nb; k < ne; ++k) {
      traverse_far_device(tree.near_nodes[k], tree, tgt, params, far);
    }
  }

  if (far > 0.0) {
    d_conc[storage_cell] += far;
  }
}

__global__ void fmm_far_tree_device_kernel(
    FmmDeviceTreeView tree,
    FmmFarLaunchParams params,
    double* d_conc) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  double tgt[3];
  int storage_cell = 0;
  if (!decode_owned_voxel(tid, params.grid, params.dom, tgt, storage_cell)) {
    return;
  }

  double far = 0.0;
  traverse_far_device(0, tree, tgt, params, far);
  if (far > 0.0) {
    d_conc[storage_cell] += far;
  }
}

}  // namespace

void launch_fmm_far_local_device_kernel(
    const FmmDeviceTreeView& tree,
    const FmmFarLaunchParams& params,
    double* d_conc,
    cudaStream_t stream) {
  const int owned_nx = params.grid.owned_storage_x_end
                     - params.grid.owned_storage_x_begin;
  const int total_owned = owned_nx * params.grid.global_ny * params.grid.global_nz;
  if (total_owned <= 0 || tree.num_nodes <= 0 || d_conc == nullptr) return;
  constexpr int kBlock = 128;
  const int grid = (total_owned + kBlock - 1) / kBlock;
  fmm_far_local_device_kernel<<<grid, kBlock, 0, stream>>>(tree, params, d_conc);
}

void launch_fmm_far_tree_device_kernel(
    const FmmDeviceTreeView& tree,
    const FmmFarLaunchParams& params,
    double* d_conc,
    cudaStream_t stream) {
  const int owned_nx = params.grid.owned_storage_x_end
                     - params.grid.owned_storage_x_begin;
  const int total_owned = owned_nx * params.grid.global_ny * params.grid.global_nz;
  if (total_owned <= 0 || tree.num_nodes <= 0 || d_conc == nullptr) return;
  constexpr int kBlock = 128;
  const int grid = (total_owned + kBlock - 1) / kBlock;
  fmm_far_tree_device_kernel<<<grid, kBlock, 0, stream>>>(tree, params, d_conc);
}

}  // namespace gpu
}  // namespace gutibm
