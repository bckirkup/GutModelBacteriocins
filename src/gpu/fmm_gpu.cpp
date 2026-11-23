#include "fmm_gpu.h"
#include "advection.h"
#include "chemical_field_gpu.h"
#include "device_memory.h"
#include "dispatch.h"
#include "domain.h"
#include "fmm.h"
#include "gpu_kernels.h"
#include "robin_correction_table.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

#include <cmath>
#include <limits>
#include <vector>

namespace gutibm {

namespace {

#ifdef GUTIBM_CUDA
gpu::DomainParams make_fmm_domain_params(const Domain& domain) {
  gpu::DomainParams p{};
  p.nx = domain.nx();
  p.ny = domain.ny();
  p.nz = domain.nz();
  p.dx_x = domain.dx_x();
  p.dx_y = domain.dx_y();
  p.dx_z = domain.dx_z();
  p.lo[0] = domain.lo()[0];
  p.lo[1] = domain.lo()[1];
  p.lo[2] = domain.lo()[2];
  p.extent[0] = domain.hi()[0] - domain.lo()[0];
  p.extent[1] = domain.hi()[1] - domain.lo()[1];
  p.extent[2] = domain.hi()[2] - domain.lo()[2];
  p.periodic[0] = domain.config().periodic[0];
  p.periodic[1] = domain.config().periodic[1];
  p.periodic[2] = domain.config().periodic[2];
  p.z_lo = domain.lo()[2];
  p.z_hi = domain.hi()[2];
  return p;
}

gpu::AdvectionParams make_fmm_advection_params(const AdvectionField& adv) {
  const auto& cfg = adv.config();
  gpu::AdvectionParams p{};
  p.h = cfg.mucus_thickness;
  p.lo_z = adv.lo_z();
  p.profile_alpha = cfg.profile_alpha;
  p.crypts_enabled = cfg.crypts_enabled;
  p.crypt_depth = cfg.crypt_depth;
  p.peristaltic_enabled = cfg.peristaltic_enabled;
  p.peristaltic_period = cfg.peristaltic_period;
  p.peristaltic_amplitude = cfg.peristaltic_amplitude;
  p.peristaltic_wavelength = cfg.peristaltic_wavelength;
  p.current_time = adv.current_time();
  p.v_radial_max = cfg.mucus_thickness / cfg.radial_turnover;
  p.v_distal_max = cfg.distal_length / cfg.distal_transit_time;
  return p;
}

gpu::GfSourceParams make_fmm_source_params(const GreensFunctionParams& avg) {
  gpu::GfSourceParams sp{};
  sp.diff_coeff = avg.diff_coeff;
  sp.source_rate = 1.0;
  sp.retardation = avg.retardation;
  sp.decay_rate = avg.decay_rate;
  sp.lumen_transfer_length = avg.lumen_transfer_length;
  sp.robin_cutoff = avg.robin_cutoff;
  sp.lumen_transfer_basis_free = avg.lumen_transfer_basis_free ? 1 : 0;
  sp.robin_table_index = -1;
  sp.image_series_relative_tolerance = avg.image_series_relative_tolerance;
  sp.image_series_max_shells =
      avg.image_series_legacy_reflections
              && !avg.image_series_max_shells_explicit
          ? kHistoricalLegacyImageSeriesShells
          : avg.image_series_max_shells;
  sp.image_series_legacy_reflections =
      avg.image_series_legacy_reflections ? 1 : 0;
  sp.drift_correction = avg.drift_correction ? 1 : 0;
  sp.drift_auto_dispatch = avg.drift_auto_dispatch ? 1 : 0;
  sp.peclet_z_threshold = avg.peclet_z_threshold;
  return sp;
}

void prepare_fmm_robin_tables(const Domain& domain,
                              const AdvectionField& adv,
                              const GreensFunctionParams& avg,
                              DeviceBuffer<double>& d_robin_tables,
                              gpu::FmmFarLaunchParams& launch_params) {
  if (!robin::transfer_enabled(avg.lumen_transfer_length)) return;

  const double d_eff = avg.diff_coeff / avg.retardation;
  const robin::TransferBasis basis = avg.lumen_transfer_basis_free
      ? robin::TransferBasis::Free
      : robin::TransferBasis::Effective;

  std::vector<double> table_values;
  const bool need_legacy = !avg.drift_correction;
  const bool need_physical = avg.drift_correction || avg.drift_auto_dispatch;
  int next_idx = 0;

  if (need_legacy) {
    const auto table = robin::global_table_cache().get(
        adv, domain.lo()[2], domain.hi()[2], avg.diff_coeff, d_eff,
        avg.decay_rate, avg.lumen_transfer_length, avg.robin_cutoff, basis,
        avg.image_series_relative_tolerance, avg.image_series_max_shells,
        avg.image_series_max_shells_explicit,
        avg.image_series_legacy_reflections, false);
    table_values.insert(table_values.end(),
                        table->legacy_values.begin(),
                        table->legacy_values.end());
    launch_params.robin_table_index_legacy = next_idx++;
  }
  if (need_physical) {
    const auto table = robin::global_table_cache().get(
        adv, domain.lo()[2], domain.hi()[2], avg.diff_coeff, d_eff,
        avg.decay_rate, avg.lumen_transfer_length, avg.robin_cutoff, basis,
        avg.image_series_relative_tolerance, avg.image_series_max_shells,
        avg.image_series_max_shells_explicit,
        avg.image_series_legacy_reflections, true);
    table_values.insert(table_values.end(),
                        table->physical_values.begin(),
                        table->physical_values.end());
    launch_params.robin_table_index_physical = next_idx++;
  }
  d_robin_tables.upload(table_values);
  launch_params.robin_tables = d_robin_tables.data();
  launch_params.robin_table_count = next_idx;
}

struct PackedTreeOffsets {
  size_t half_size = 0;
  size_t center_of_source = 0;
  size_t strength = 0;
  size_t multipole = 0;
  size_t local = 0;
  size_t is_leaf = 0;
  size_t near_offsets = 0;
  size_t near_nodes = 0;
};

PackedTreeOffsets pack_fmm_tree_buffers(const FMM& fmm,
                                        bool use_locals,
                                        std::vector<double>& real_tree,
                                        std::vector<int>& int_tree) {
  PackedTreeOffsets off{};
  const auto n = static_cast<size_t>(fmm.num_nodes());
  const auto ncoeff = static_cast<size_t>(
      FMM::num_coefficients(fmm.expansion_order()));

  // Pack all double arrays into a single contiguous H2D buffer.
  off.half_size = 3 * n;
  off.center_of_source = off.half_size + n;
  off.strength = off.center_of_source + 3 * n;
  off.multipole = off.strength + n;
  off.local = off.multipole + n * ncoeff;
  const size_t real_total = use_locals ? (off.local + n * ncoeff) : off.local;
  real_tree.assign(real_total, 0.0);

  size_t total_near = 0;
  if (use_locals) {
    for (int i = 0; i < fmm.num_nodes(); ++i) {
      total_near += fmm.node(i).near_list.size();
    }
  }

  // Pack all int arrays into a single contiguous H2D buffer.
  off.is_leaf = 8 * n;
  off.near_offsets = off.is_leaf + n;
  off.near_nodes = use_locals ? (off.near_offsets + n + 1) : off.near_offsets;
  const size_t int_total = use_locals ? (off.near_nodes + total_near) : off.near_offsets;
  int_tree.assign(int_total, 0);

  size_t near_cursor = 0;
  for (int i = 0; i < fmm.num_nodes(); ++i) {
    const auto idx = static_cast<size_t>(i);
    const FMMNode& node = fmm.node(i);

    real_tree[3 * idx + 0] = node.center[0];
    real_tree[3 * idx + 1] = node.center[1];
    real_tree[3 * idx + 2] = node.center[2];
    real_tree[off.half_size + idx] = node.half_size;
    real_tree[off.center_of_source + 3 * idx + 0] = node.center_of_source[0];
    real_tree[off.center_of_source + 3 * idx + 1] = node.center_of_source[1];
    real_tree[off.center_of_source + 3 * idx + 2] = node.center_of_source[2];
    real_tree[off.strength + idx] = node.total_source_strength;

    for (size_t k = 0; k < ncoeff && k < node.multipole.size(); ++k) {
      real_tree[off.multipole + idx * ncoeff + k] = node.multipole[k];
    }
    if (use_locals) {
      for (size_t k = 0; k < ncoeff && k < node.local.size(); ++k) {
        real_tree[off.local + idx * ncoeff + k] = node.local[k];
      }
    }

    for (size_t c = 0; c < 8; ++c) {
      int_tree[8 * idx + c] = node.children[c];
    }
    int_tree[off.is_leaf + idx] = node.is_leaf ? 1 : 0;

    if (use_locals) {
      int_tree[off.near_offsets + idx] = static_cast<int>(near_cursor);
      for (int src_idx : node.near_list) {
        int_tree[off.near_nodes + near_cursor++] = src_idx;
      }
    }
  }
  if (use_locals) {
    int_tree[off.near_offsets + n] = static_cast<int>(near_cursor);
  }
  return off;
}
#endif

}  // namespace

bool gpu_accumulate_far_field_device(const FMM& fmm,
                                     const Domain& domain,
                                     const AdvectionField& adv,
                                     const GreensFunctionParams& avg_params,
                                     Real theta,
                                     Real near_cutoff,
                                     double* d_conc,
                                     const FmmDeviceLayout& layout) {
#ifndef GUTIBM_CUDA
  (void)fmm;
  (void)domain;
  (void)adv;
  (void)avg_params;
  (void)theta;
  (void)near_cutoff;
  (void)d_conc;
  (void)layout;
  return false;
#else
  if (!gpu_runtime_enabled() || d_conc == nullptr) return false;
  if (fmm.empty()) return true;
  if (layout.owned_storage_x_end <= layout.owned_storage_x_begin) return true;

  GpuTransferSite site("qssa_fmm_tree");

  const bool use_locals = fmm.locals_match_query(theta, near_cutoff);
  std::vector<double> real_tree;
  std::vector<int> int_tree;
  const PackedTreeOffsets off = pack_fmm_tree_buffers(
      fmm, use_locals, real_tree, int_tree);

  DeviceBuffer<double> d_real_tree;
  DeviceBuffer<int> d_int_tree;
  d_real_tree.upload(real_tree);
  d_int_tree.upload(int_tree);

  gpu::FmmDeviceTreeView tree{};
  tree.center = d_real_tree.data();
  tree.half_size = d_real_tree.data() + off.half_size;
  tree.center_of_source = d_real_tree.data() + off.center_of_source;
  tree.total_source_strength = d_real_tree.data() + off.strength;
  tree.multipole = d_real_tree.data() + off.multipole;
  tree.children = d_int_tree.data();
  tree.is_leaf = d_int_tree.data() + off.is_leaf;
  tree.num_nodes = fmm.num_nodes();
  tree.coeffs_per_node = FMM::num_coefficients(fmm.expansion_order());
  tree.expansion_order = fmm.expansion_order();
  if (use_locals) {
    tree.local = d_real_tree.data() + off.local;
    tree.near_offsets = d_int_tree.data() + off.near_offsets;
    tree.near_nodes = (int_tree.size() > off.near_nodes)
        ? (d_int_tree.data() + off.near_nodes)
        : nullptr;
  }

  gpu::FmmFarLaunchParams launch{};
  launch.dom = make_fmm_domain_params(domain);
  launch.adv = make_fmm_advection_params(adv);
  launch.avg_params = make_fmm_source_params(avg_params);
  launch.grid.global_nx = domain.nx();
  launch.grid.global_ny = domain.ny();
  launch.grid.global_nz = domain.nz();
  launch.grid.storage_nx = layout.storage_nx;
  launch.grid.owned_storage_x_begin = layout.owned_storage_x_begin;
  launch.grid.owned_storage_x_end = layout.owned_storage_x_end;
  launch.grid.local_global_x_begin = layout.local_global_x_begin;
  launch.theta = theta;
  launch.near_cutoff = near_cutoff;

  const double eps = std::numeric_limits<double>::epsilon();
  launch.fd_rel_1 = std::cbrt(eps);
  launch.fd_rel_2 = std::pow(eps, 0.25);
  launch.fd_rel_3 = std::pow(eps, 0.2);

  DeviceBuffer<double> d_robin_tables;
  prepare_fmm_robin_tables(domain, adv, avg_params, d_robin_tables, launch);

  if (use_locals) {
    gpu::launch_fmm_far_local_device_kernel(
        tree, launch, d_conc, gpu_compute_stream());
    gpu_sync_compute();
    gpu_check_error("fmm_far_local_device_kernel");
  } else {
    gpu::launch_fmm_far_tree_device_kernel(
        tree, launch, d_conc, gpu_compute_stream());
    gpu_sync_compute();
    gpu_check_error("fmm_far_tree_device_kernel");
  }
  return true;
#endif
}

bool gpu_accumulate_far_field_device(const FMM& fmm,
                                     const Domain& domain,
                                     const AdvectionField& adv,
                                     const GreensFunctionParams& avg_params,
                                     Real theta,
                                     Real near_cutoff,
                                     ChemicalFieldGpu* chem_gpu,
                                     Int toxin_species_idx) {
  if (chem_gpu == nullptr || !chem_gpu->active()) return false;
  double* d_conc = chem_gpu->conc_device(toxin_species_idx);
  if (d_conc == nullptr) return false;
  FmmDeviceLayout layout{};
  layout.storage_nx = chem_gpu->storage_nx();
  layout.owned_storage_x_begin = chem_gpu->owned_storage_x_begin();
  layout.owned_storage_x_end = chem_gpu->owned_storage_x_end();
  layout.local_global_x_begin = chem_gpu->slab_mode()
      ? domain.local_grid_x_begin()
      : 0;
  return gpu_accumulate_far_field_device(
      fmm, domain, adv, avg_params, theta, near_cutoff, d_conc, layout);
}

}  // namespace gutibm
