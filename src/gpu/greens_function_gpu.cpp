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
#include "error.h"
#include <algorithm>
#include <cmath>
#include <format>
#include <mutex>
#include <memory>
#include <ranges>
#include <string>
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
std::vector<std::shared_ptr<const robin::Table>> robin_table_device_set;
std::mutex robin_table_device_mutex;

bool launch_superpose(const Domain& domain,
                      const AdvectionField& adv,
                      const std::vector<Vec3>& sources,
                      const std::vector<GreensFunctionParams>& params,
                      const std::vector<Real>& strength_factors,
                      double* d_grid,
                      Real cutoff_radius, uint64_t* cap_hits,
                      uint64_t* kernel_evaluations,
                      uint64_t* low_screening_evaluations,
                      uint64_t* negative_field_count,
                      Real* most_negative_field) {
  if (!gpu_runtime_enabled() || sources.empty()) return false;
  if (d_grid == nullptr) return false;
  if (cap_hits != nullptr) *cap_hits = 0;
  if (kernel_evaluations != nullptr) *kernel_evaluations = 0;
  if (low_screening_evaluations != nullptr) *low_screening_evaluations = 0;
  if (negative_field_count != nullptr) *negative_field_count = 0;
  if (most_negative_field != nullptr) *most_negative_field = 0.0;
  std::scoped_lock table_lock(robin_table_device_mutex);

  Int ncells = domain.ncells();
  cudaMemset(d_grid, 0, static_cast<size_t>(ncells) * sizeof(double));

  std::vector<double> sx(sources.size());
  std::vector<double> sy(sources.size());
  std::vector<double> sz(sources.size());
  std::vector<gpu::GfSourceParams> sp(params.size());
  std::vector<std::shared_ptr<const robin::Table>> source_tables(params.size());
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
    sp[i].image_series_max_shells =
        params[i].image_series_legacy_reflections
                && !params[i].image_series_max_shells_explicit
            ? kHistoricalLegacyImageSeriesShells
            : params[i].image_series_max_shells;
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
      source_tables[i] = table;
    }
  }
  std::vector<std::shared_ptr<const robin::Table>> robin_tables;
  const std::vector<int> launch_indices = make_robin_launch_table_indices(
      source_tables, robin_tables);
  for (size_t i = 0; i < launch_indices.size(); ++i) {
    sp[i].robin_table_index = launch_indices[i];
  }
  if (const bool same_uploaded_set = robin_tables.size()
          == robin_table_device_set.size()
          && std::equal(
              robin_tables.begin(), robin_tables.end(),
              robin_table_device_set.begin(),
              [](const auto& current, const auto& uploaded) {
                return current.get() == uploaded.get();
              });
      !same_uploaded_set) {
    std::vector<double> robin_table_values;
    robin_table_values.reserve(
        robin_tables.size() * robin::kTableValueCount);
    for (const auto& table : robin_tables) {
      robin_table_values.insert(robin_table_values.end(),
                                table->values.begin(), table->values.end());
    }
    robin_table_device_buffer.upload(robin_table_values);
    robin_table_device_set = robin_tables;
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
  DeviceBuffer<unsigned int> d_robin_index_error;
  unsigned int* robin_index_error = nullptr;
  if (!robin_tables.empty()) {
    d_robin_index_error.allocate(1);
    cudaMemset(d_robin_index_error.data(), 0, sizeof(unsigned int));
    robin_index_error = d_robin_index_error.data();
  }
  DeviceBuffer<unsigned long long> d_cap_hits;
  DeviceBuffer<unsigned long long> d_low_screening_evaluations;
  DeviceBuffer<unsigned long long> d_negative_field_count;
  DeviceBuffer<double> d_most_negative_field;
  DeviceBuffer<unsigned long long> d_kernel_evaluations;
  unsigned long long* cap_hits_device = nullptr;
  unsigned long long* low_screening_evaluations_device = nullptr;
  unsigned long long* negative_field_count_device = nullptr;
  double* most_negative_field_device = nullptr;
  unsigned long long* kernel_evaluations_device = nullptr;
  if (cap_hits != nullptr) {
    d_cap_hits.allocate(1);
    cudaMemset(d_cap_hits.data(), 0, sizeof(unsigned long long));
    cap_hits_device = d_cap_hits.data();
  }
  if (low_screening_evaluations != nullptr) {
    d_low_screening_evaluations.allocate(1);
    cudaMemset(d_low_screening_evaluations.data(), 0,
               sizeof(unsigned long long));
    low_screening_evaluations_device = d_low_screening_evaluations.data();
  }
  if (negative_field_count != nullptr) {
    d_negative_field_count.allocate(1);
    cudaMemset(d_negative_field_count.data(), 0, sizeof(unsigned long long));
    negative_field_count_device = d_negative_field_count.data();
  }
  if (most_negative_field != nullptr) {
    d_most_negative_field.allocate(1);
    cudaMemset(d_most_negative_field.data(), 0, sizeof(double));
    most_negative_field_device = d_most_negative_field.data();
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
      static_cast<int>(robin_tables.size()), robin_index_error,
      gpu_compute_stream(), cap_hits_device, low_screening_evaluations_device,
      negative_field_count_device, most_negative_field_device,
      kernel_evaluations_device);

  gpu_sync_compute();
  gpu_check_error("superpose_kernel");
  if (robin_index_error != nullptr) {
    unsigned int error = 0;
    d_robin_index_error.download(&error, 1);
    if (error != 0) {
      throw SimulationError(std::format(
          "GPU Robin table index exceeded the launch-local table count "
          "(count={})",
          robin_tables.size()));
    }
  }
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
  if (low_screening_evaluations != nullptr) {
    unsigned long long count = 0;
    d_low_screening_evaluations.download(&count, 1);
    *low_screening_evaluations = static_cast<uint64_t>(count);
  }
  if (negative_field_count != nullptr) {
    unsigned long long count = 0;
    d_negative_field_count.download(&count, 1);
    *negative_field_count = static_cast<uint64_t>(count);
  }
  if (most_negative_field != nullptr) {
    d_most_negative_field.download(most_negative_field, 1);
  }
  return true;
}
#endif

}  // namespace

std::vector<int> make_robin_launch_table_indices(
    const std::vector<std::shared_ptr<const robin::Table>>& source_tables,
    std::vector<std::shared_ptr<const robin::Table>>& launch_tables) {
  std::vector indices(source_tables.size(), -1);
  for (size_t source = 0; source < source_tables.size(); ++source) {
    const auto& table = source_tables[source];
    if (table == nullptr) continue;
    const auto existing = std::ranges::find_if(
        launch_tables, [&table](const auto& candidate) {
          return candidate.get() == table.get();
        });
    if (existing == launch_tables.end()) {
      if (launch_tables.size() >= kMaximumRobinDeviceTables) {
        throw SimulationError(std::format(
            "Robin GPU launch references {} distinct tables; launch limit is "
            "{}",
            launch_tables.size() + 1, kMaximumRobinDeviceTables));
      }
      launch_tables.push_back(table);
      indices[source] = static_cast<int>(launch_tables.size() - 1);
    } else {
      indices[source] = static_cast<int>(
          std::distance(launch_tables.begin(), existing));
    }
  }
  return indices;
}

std::vector<size_t> robin_host_fallback_sources(
    const Domain& domain,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params) {
  const double cell_radius = std::min(
      {domain.dx_x(), domain.dx_y(), domain.dx_z()});
  std::vector<size_t> fallback;
  for (size_t source = 0; source < sources.size(); ++source) {
    if (!robin::transfer_enabled(params[source].lumen_transfer_length)) {
      continue;
    }
    if (robin::requires_direct_evaluation(
            sources[source][2], sources[source][2], 0.0,
            domain.lo()[2], domain.hi()[2], cell_radius)) {
      fallback.push_back(source);
    }
  }
  return fallback;
}

bool gpu_superpose_to_grid(
    const Domain& domain,
    const AdvectionField& adv,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    std::vector<Real>& grid_conc,
    Real cutoff_radius, uint64_t* cap_hits,
    uint64_t* kernel_evaluations,
    uint64_t* low_screening_evaluations,
    uint64_t* negative_field_count,
    Real* most_negative_field) {

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
  (void)low_screening_evaluations;
  (void)negative_field_count;
  (void)most_negative_field;
  return false;
#else
  if (sources.empty()) {
    grid_conc.assign(domain.ncells(), 0.0);
    if (cap_hits != nullptr) *cap_hits = 0;
    if (kernel_evaluations != nullptr) *kernel_evaluations = 0;
    if (low_screening_evaluations != nullptr) *low_screening_evaluations = 0;
    if (negative_field_count != nullptr) *negative_field_count = 0;
    if (most_negative_field != nullptr) *most_negative_field = 0.0;
    return true;
  }

  DeviceBuffer<double> d_grid;
  d_grid.allocate(static_cast<size_t>(domain.ncells()));
  if (!launch_superpose(domain, adv, sources, params, strength_factors,
                        d_grid.data(), cutoff_radius, cap_hits,
                        kernel_evaluations, low_screening_evaluations,
                        negative_field_count, most_negative_field)) {
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
    uint64_t* kernel_evaluations,
    uint64_t* low_screening_evaluations,
    uint64_t* negative_field_count,
    Real* most_negative_field) {

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
  (void)low_screening_evaluations;
  (void)negative_field_count;
  (void)most_negative_field;
  return false;
#else
  if (sources.empty()) {
    if (d_grid_conc != nullptr) {
      cudaMemset(d_grid_conc, 0,
                 static_cast<size_t>(domain.ncells()) * sizeof(double));
    }
    if (cap_hits != nullptr) *cap_hits = 0;
    if (kernel_evaluations != nullptr) *kernel_evaluations = 0;
    if (low_screening_evaluations != nullptr) *low_screening_evaluations = 0;
    if (negative_field_count != nullptr) *negative_field_count = 0;
    if (most_negative_field != nullptr) *most_negative_field = 0.0;
    return true;
  }
  return launch_superpose(domain, adv, sources, params, strength_factors,
                          d_grid_conc, cutoff_radius, cap_hits,
                          kernel_evaluations, low_screening_evaluations,
                          negative_field_count, most_negative_field);
#endif
}

}  // namespace gutibm
