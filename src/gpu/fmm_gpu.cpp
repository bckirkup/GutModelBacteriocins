#include "fmm_gpu.h"
#include "fmm.h"
#include "fmm_kernel.h"
#include "domain.h"
#include "dispatch.h"
#include "device_memory.h"
#include "gpu_kernels.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

#include <vector>

namespace gutibm {

namespace {

int leaf_count(const FMM& fmm) {
  int count = 0;
  for (int i = 0; i < fmm.num_nodes(); ++i) {
    if (fmm.node(i).is_leaf) ++count;
  }
  return count;
}

}  // namespace

bool gpu_accumulate_far_field_local(const FMM& fmm,
                                    const Domain& domain,
                                    int expansion_order,
                                    const std::vector<Real>& near_conc,
                                    std::vector<Real>& toxin_conc) {
#ifndef GUTIBM_CUDA
  (void)fmm;
  (void)domain;
  (void)expansion_order;
  (void)near_conc;
  (void)toxin_conc;
  return false;
#else
  if (!gpu_runtime_enabled() || !fmm.locals_ready() || fmm.empty()) return false;

  const int num_leaves = leaf_count(fmm);
  if (num_leaves <= 0) return false;

  const int coeffs_per_leaf = FMM::num_coefficients(expansion_order);
  const Int ncells = domain.ncells();
  if (static_cast<Int>(near_conc.size()) != ncells
      || static_cast<Int>(toxin_conc.size()) != ncells) {
    return false;
  }

  std::vector<int> leaf_ids(fmm.num_nodes(), -1);
  int leaf_idx = 0;
  for (int i = 0; i < fmm.num_nodes(); ++i) {
    if (fmm.node(i).is_leaf) leaf_ids[static_cast<size_t>(i)] = leaf_idx++;
  }

  std::vector<int> cell_leaf(static_cast<size_t>(ncells), -1);
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.nx(); ++ix) {
        const Vec3 tgt = domain.cell_center(ix, iy, iz);
        const Int cell = domain.cell_index(ix, iy, iz);
        int node = 0;
        while (!fmm.node(node).is_leaf) {
          const auto& n = fmm.node(node);
          int oct = 0;
          if (tgt[0] > n.center[0]) oct |= 1;
          if (tgt[1] > n.center[1]) oct |= 2;
          if (tgt[2] > n.center[2]) oct |= 4;
          if (n.children[oct] < 0) break;
          node = n.children[oct];
        }
        cell_leaf[static_cast<size_t>(cell)] = leaf_ids[static_cast<size_t>(node)];
      }
    }
  }

  std::vector<double> packed_local(static_cast<size_t>(num_leaves) * coeffs_per_leaf, 0.0);
  std::vector<double> packed_center(static_cast<size_t>(num_leaves) * 3, 0.0);
  for (int i = 0; i < fmm.num_nodes(); ++i) {
    if (!fmm.node(i).is_leaf) continue;
    const int lid = leaf_ids[static_cast<size_t>(i)];
    const auto& node = fmm.node(i);
    packed_center[static_cast<size_t>(lid) * 3 + 0] = node.center[0];
    packed_center[static_cast<size_t>(lid) * 3 + 1] = node.center[1];
    packed_center[static_cast<size_t>(lid) * 3 + 2] = node.center[2];
    for (int c = 0; c < coeffs_per_leaf; ++c) {
      packed_local[static_cast<size_t>(lid) * coeffs_per_leaf + c] = node.local[c];
    }
  }

  DeviceBuffer<double> d_leaf_local;
  DeviceBuffer<double> d_leaf_center;
  DeviceBuffer<int> d_cell_leaf;
  DeviceBuffer<double> d_near;
  DeviceBuffer<double> d_out;
  d_leaf_local.upload(packed_local);
  d_leaf_center.upload(packed_center);
  d_cell_leaf.upload(cell_leaf);
  d_near.upload(near_conc);
  d_out.upload(toxin_conc);

  gpu::launch_fmm_far_local_kernel(
      d_leaf_local.data(), d_leaf_center.data(), d_cell_leaf.data(),
      d_near.data(), d_out.data(),
      static_cast<int>(ncells), num_leaves, coeffs_per_leaf, expansion_order,
      domain.lo()[0], domain.lo()[1], domain.lo()[2], domain.dx(),
      domain.nx(), domain.ny(), gpu_compute_stream());
  gpu_sync_compute();
  gpu_check_error("fmm_far_local_kernel");

  d_out.download(toxin_conc);
  return true;
#endif
}

bool gpu_accumulate_far_field_local_device(const FMM& fmm,
                                           const Domain& domain,
                                           int expansion_order,
                                           double* d_near_conc,
                                           double* d_out_conc) {
#ifndef GUTIBM_CUDA
  (void)fmm;
  (void)domain;
  (void)expansion_order;
  (void)d_near_conc;
  (void)d_out_conc;
  return false;
#else
  (void)fmm;
  (void)domain;
  (void)expansion_order;
  (void)d_near_conc;
  (void)d_out_conc;
  return false;
#endif
}

}  // namespace gutibm
