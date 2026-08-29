#include "greens_function_gpu.h"
#include "dispatch.h"
#include "device_memory.h"
#include "gpu_types.h"
#include "gpu_kernels.h"
#include "robin_correction_table.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

#include "domain.h"
#include "advection.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <memory>
#include <vector>

namespace gutibm {
namespace {

#ifdef GUTIBM_CUDA
gpu::DomainParams make_domain_params(const Domain& domain) {
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

gpu::AdvectionParams make_advection_params(const AdvectionField& adv) {
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

  Real h = cfg.mucus_thickness;
  p.v_radial_max = h / cfg.radial_turnover;
  p.v_distal_max = cfg.distal_length / cfg.distal_transit_time;
  return p;
}

DeviceBuffer<double> robin_table_device_buffer;
uint64_t robin_table_device_generation = std::numeric_limits<uint64_t>::max();
std::vector<std::shared_ptr<const robin::Table>> robin_table_device_snapshot;
std::mutex robin_table_device_mutex;

bool launch_superpose(const Domain& domain,
                      const AdvectionField& adv,
                      const std::vector<Vec3>& sources,
                      const std::vector<GreensFunctionParams>& params,
                      const std::vector<Real>& strength_factors,
                      double* d_grid,
                      Real cutoff_radius, uint64_t* cap_hits,
                      uint64_t* kernel_evaluations) {
  if (!gpu_runtime_enabled() || sources.empty()) return false;
  if (d_grid == nullptr) return false;
  if (cap_hits != nullptr) *cap_hits = 0;
  if (kernel_evaluations != nullptr) *kernel_evaluations = 0;
  std::lock_guard<std::mutex> table_lock(robin_table_device_mutex);

  Int ncells = domain.ncells();
  cudaMemset(d_grid, 0, static_cast<size_t>(ncells) * sizeof(double));

  std::vector<double> sx(sources.size());
  std::vector<double> sy(sources.size());
  std::vector<double> sz(sources.size());
  std::vector<gpu::GfSourceParams> sp(params.size());
  std::vector<std::shared_ptr<const robin::Table>> robin_tables;
  for (size_t i = 0; i < sources.size(); ++i) {
    sx[i] = sources[i][0];
    sy[i] = sources[i][1];
    sz[i] = sources[i][2];
    Real scale = (i < strength_factors.size()) ? strength_factors[i] : 1.0;
    sp[i].diff_coeff = params[i].diff_coeff;
    sp[i].source_rate = params[i].source_rate * scale;
    sp[i].retardation = params[i].retardation;
    sp[i].decay_rate = params[i].decay_rate;
    sp[i].lumen_transfer_length = params[i].lumen_transfer_length;
    sp[i].robin_cutoff = params[i].robin_cutoff;
    sp[i].image_series_relative_tolerance =
        params[i].image_series_relative_tolerance;
    sp[i].image_series_max_shells = params[i].image_series_max_shells;
    sp[i].image_series_legacy_reflections =
        params[i].image_series_legacy_reflections ? 1 : 0;
    sp[i].lumen_transfer_basis_free =
        params[i].lumen_transfer_basis_free ? 1 : 0;
    sp[i].robin_table_index = -1;
    if (robin::transfer_enabled(params[i].lumen_transfer_length)) {
      const auto table = robin::global_table_cache().get(
          adv, domain.lo()[2], domain.hi()[2], params[i].diff_coeff,
          params[i].diff_coeff / params[i].retardation,
          params[i].decay_rate, params[i].lumen_transfer_length,
          params[i].robin_cutoff,
          params[i].lumen_transfer_basis_free
              ? robin::TransferBasis::Free
              : robin::TransferBasis::Effective);
      const auto existing = std::find_if(
          robin_tables.begin(), robin_tables.end(),
          [&table](const auto& existing_table) {
            return existing_table.get() == table.get();
          });
      if (existing == robin_tables.end()) {
        robin_tables.push_back(table);
        sp[i].robin_table_index =
            static_cast<int>(robin_tables.size() - 1);
      } else {
        sp[i].robin_table_index = static_cast<int>(
            std::distance(robin_tables.begin(), existing));
      }
    }
  }
  const auto cache_snapshot = robin::global_table_cache().snapshot();
  if (cache_snapshot.generation != robin_table_device_generation) {
    std::vector<double> robin_table_values;
    robin_table_values.reserve(
        cache_snapshot.tables.size() * robin::kTableValueCount);
    for (const auto& table : cache_snapshot.tables) {
      robin_table_values.insert(robin_table_values.end(),
                                table->values.begin(), table->values.end());
    }
    robin_table_device_buffer.upload(robin_table_values);
    robin_table_device_snapshot = cache_snapshot.tables;
    robin_table_device_generation = cache_snapshot.generation;
  }
  const auto& device_snapshot = robin_table_device_snapshot;
  for (auto& source : sp) {
    if (source.robin_table_index < 0) continue;
    const auto& table = robin_tables[static_cast<size_t>(
        source.robin_table_index)];
    const auto table_position = std::find_if(
        device_snapshot.begin(), device_snapshot.end(),
        [&table](const auto& candidate) {
          return candidate.get() == table.get();
        });
    source.robin_table_index = static_cast<int>(
        std::distance(device_snapshot.begin(), table_position));
  }

  DeviceBuffer<double> d_sx;
  DeviceBuffer<double> d_sy;
  DeviceBuffer<double> d_sz;
  DeviceBuffer<gpu::GfSourceParams> d_params;
  d_sx.upload(sx);
  d_sy.upload(sy);
  d_sz.upload(sz);
  d_params.upload(sp);

  const auto span_x = static_cast<int>(
      std::ceil(cutoff_radius / domain.dx_x()));
  const auto span_y = static_cast<int>(
      std::ceil(cutoff_radius / domain.dx_y()));
  const auto span_z = static_cast<int>(
      std::ceil(cutoff_radius / domain.dx_z()));
  const auto dom = make_domain_params(domain);
  const auto adv_p = make_advection_params(adv);
  DeviceBuffer<unsigned long long> d_cap_hits;
  DeviceBuffer<unsigned long long> d_kernel_evaluations;
  unsigned long long* cap_hits_device = nullptr;
  unsigned long long* kernel_evaluations_device = nullptr;
  if (cap_hits != nullptr) {
    d_cap_hits.allocate(1);
    cudaMemset(d_cap_hits.data(), 0, sizeof(unsigned long long));
    cap_hits_device = d_cap_hits.data();
  }
  if (kernel_evaluations != nullptr) {
    d_kernel_evaluations.allocate(1);
    cudaMemset(d_kernel_evaluations.data(), 0, sizeof(unsigned long long));
    kernel_evaluations_device = d_kernel_evaluations.data();
  }

  gpu::launch_superpose_kernel(
      d_sx.data(), d_sy.data(), d_sz.data(), d_params.data(), d_grid,
      robin_table_device_buffer.data(),
      dom, adv_p, static_cast<int>(sources.size()), span_x, span_y, span_z,
      gpu_compute_stream(), cap_hits_device, kernel_evaluations_device);

  gpu_sync_compute();
  gpu_check_error("superpose_kernel");
  if (cap_hits != nullptr) {
    unsigned long long count = 0;
    d_cap_hits.download(&count, 1);
    *cap_hits = static_cast<uint64_t>(count);
  }
  if (kernel_evaluations != nullptr) {
    unsigned long long count = 0;
    d_kernel_evaluations.download(&count, 1);
    *kernel_evaluations = static_cast<uint64_t>(count);
  }
  return true;
}
#endif

}  // namespace

bool gpu_superpose_to_grid(
    const Domain& domain,
    const AdvectionField& adv,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    std::vector<Real>& grid_conc,
    Real cutoff_radius, uint64_t* cap_hits,
    uint64_t* kernel_evaluations) {

#ifndef GUTIBM_CUDA
  (void)domain;
  (void)adv;
  (void)sources;
  (void)params;
  (void)strength_factors;
  (void)grid_conc;
  (void)cutoff_radius;
  (void)cap_hits;
  (void)kernel_evaluations;
  return false;
#else
  if (sources.empty()) {
    grid_conc.assign(domain.ncells(), 0.0);
    if (cap_hits != nullptr) *cap_hits = 0;
    if (kernel_evaluations != nullptr) *kernel_evaluations = 0;
    return true;
  }

  DeviceBuffer<double> d_grid;
  d_grid.allocate(static_cast<size_t>(domain.ncells()));
  if (!launch_superpose(domain, adv, sources, params, strength_factors,
                        d_grid.data(), cutoff_radius, cap_hits,
                        kernel_evaluations)) {
    return false;
  }
  d_grid.download(grid_conc);
  return true;
#endif
}

bool gpu_superpose_to_device(
    const Domain& domain,
    const AdvectionField& adv,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    double* d_grid_conc,
    Real cutoff_radius, uint64_t* cap_hits,
    uint64_t* kernel_evaluations) {

#ifndef GUTIBM_CUDA
  (void)domain;
  (void)adv;
  (void)sources;
  (void)params;
  (void)strength_factors;
  (void)d_grid_conc;
  (void)cutoff_radius;
  (void)cap_hits;
  (void)kernel_evaluations;
  return false;
#else
  if (sources.empty()) {
    if (d_grid_conc != nullptr) {
      cudaMemset(d_grid_conc, 0,
                 static_cast<size_t>(domain.ncells()) * sizeof(double));
    }
    if (cap_hits != nullptr) *cap_hits = 0;
    if (kernel_evaluations != nullptr) *kernel_evaluations = 0;
    return true;
  }
  return launch_superpose(domain, adv, sources, params, strength_factors,
                          d_grid_conc, cutoff_radius, cap_hits,
                          kernel_evaluations);
#endif
}

}  // namespace gutibm
