#include "diffusion_gpu.h"
#include "chemical_field.h"
#include "domain.h"
#include "dispatch.h"
#include "device_memory.h"
#include "gpu_kernels.h"
#include "gpu_profile.h"
#include "tridiagonal_factorization.h"
#include <cassert>
#include <chrono>
#include <string>
#include <vector>

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

namespace gutibm {

namespace {

#ifdef GUTIBM_CUDA
struct PeriodicPcrCoeffs {
  double gamma = 0.0;
  double corner = 0.0;
  double denominator = 1.0;
  std::vector<double> correction;
};

PeriodicPcrCoeffs build_periodic_coeffs(int n, double alpha) {
  PeriodicPcrCoeffs out;
  if (n < 3) return out;

  const double diagonal_value = 1.0 + 2.0 * alpha;
  out.gamma = -diagonal_value;
  out.corner = -alpha;

  std::vector lower(static_cast<size_t>(n - 1), -alpha);
  std::vector upper(static_cast<size_t>(n - 1), -alpha);
  std::vector diagonal(static_cast<size_t>(n), diagonal_value);
  diagonal.front() -= out.gamma;
  diagonal.back() -= out.corner * out.corner / out.gamma;

  TridiagonalFactorization factorization;
  factorization.factorize(lower, diagonal, upper);

  out.correction.assign(static_cast<size_t>(n), 0.0);
  out.correction.front() = out.gamma;
  out.correction.back() = out.corner;
  factorization.solve_in_place(out.correction);
  out.denominator = 1.0 + out.correction.front()
      + out.corner * out.correction.back() / out.gamma;
  return out;
}
#endif

#ifdef GUTIBM_CUDA
bool species_diffusion_eligible(const ChemicalSpec& spec, Real dt,
                                const Domain& domain) {
  return dt > 0.0 && domain.dx_x() > 0.0 && domain.dx_y() > 0.0
      && domain.dx_z() > 0.0 && spec.diffuses();
}

bool apply_species_diffusion_on_device(const Domain& domain,
                                       const ChemicalSpec& spec,
                                       double* d_conc,
                                       double* d_injected_amount,
                                       Real dt, int storage_nx,
                                       int owned_x_begin, int owned_x_end,
                                       DiffusionPhase phase) {
  const int nx = storage_nx;
  const int ny = domain.ny();
  const int nz = domain.nz();
  if (const int ncells = storage_nx * domain.ny() * domain.nz();
      ncells <= 0 || d_conc == nullptr) return false;

  const Real effective_diffusion = spec.diff_coeff / spec.retardation;
  const Real alpha_x = effective_diffusion * dt
      / (domain.dx_x() * domain.dx_x());
  const Real alpha_y = effective_diffusion * dt
      / (domain.dx_y() * domain.dx_y());
  const Real alpha_z = effective_diffusion * dt
      / (domain.dx_z() * domain.dx_z());
  const bool preserve_gradient =
      spec.z_gradient_enabled && spec.z_gradient_lambda > 0.0;
  double diffusion_boundary = spec.boundary_conc;
  const auto boundary_mode = spec.epithelial_boundary_mode;
  const int boundary_mode_value = to_underlying(boundary_mode);
  const Real beta = boundary_mode == EpithelialBoundaryMode::Robin
      ? spec.epithelial_transfer_coeff * dt / domain.dx_z() : 0.0;
  const Real flux_source = boundary_mode == EpithelialBoundaryMode::Flux
      ? spec.epithelial_flux * dt / domain.dx_z() : 0.0;

  if (phase != DiffusionPhase::SlabPostX
      && boundary_mode == EpithelialBoundaryMode::Dirichlet) {
    gpu::launch_set_epithelial_boundary(
        d_conc, nx, ny, owned_x_begin, owned_x_end, spec.boundary_conc,
        domain.cell_volume(), d_injected_amount,
        gpu_compute_stream());
    if (preserve_gradient) {
      gpu::launch_set_luminal_neumann(
          d_conc, nx, ny, nz, owned_x_begin, owned_x_end,
          gpu_compute_stream());
      gpu::launch_shift_z_gradient(
          d_conc, nx, ny, nz, owned_x_begin, owned_x_end, domain.dx_z(),
          spec.initial_conc, spec.z_gradient_lambda, spec.boundary_conc, -1.0,
          gpu_compute_stream());
    }
  }
  if (phase == DiffusionPhase::SlabPreX) {
    return true;
  }
  if (preserve_gradient) diffusion_boundary = 0.0;

  const PeriodicPcrCoeffs y_coeffs = build_periodic_coeffs(ny, alpha_y);
  DeviceBuffer<double> d_corr_y;
  d_corr_y.upload(y_coeffs.correction);

  if (phase == DiffusionPhase::Replicated) {
    const PeriodicPcrCoeffs x_coeffs =
        build_periodic_coeffs(domain.nx(), alpha_x);
    DeviceBuffer<double> d_corr_x;
    d_corr_x.upload(x_coeffs.correction);
    gpu::launch_diffuse_x_periodic(
        d_conc, nx, ny, nz, alpha_x, x_coeffs.gamma, x_coeffs.corner,
        x_coeffs.denominator, d_corr_x.data(), gpu_compute_stream());
  }
  gpu::launch_diffuse_y_periodic(
      d_conc, nx, ny, nz, owned_x_begin, owned_x_end, alpha_y,
      y_coeffs.gamma, y_coeffs.corner, y_coeffs.denominator,
      d_corr_y.data(), gpu_compute_stream());
  gpu::launch_diffuse_z_bounded(
      d_conc, nx, ny, nz, owned_x_begin, owned_x_end, alpha_z,
      diffusion_boundary, boundary_mode_value, beta, flux_source,
      domain.cell_volume(), d_injected_amount,
      gpu_compute_stream());

  if (preserve_gradient) {
    gpu::launch_shift_z_gradient(
        d_conc, nx, ny, nz, owned_x_begin, owned_x_end, domain.dx_z(),
        spec.initial_conc,
        spec.z_gradient_lambda, spec.boundary_conc, 1.0, gpu_compute_stream());
    gpu::launch_set_luminal_neumann(
        d_conc, nx, ny, nz, owned_x_begin, owned_x_end, gpu_compute_stream());
  }

  gpu::launch_clamp_nonneg(
      d_conc, nx, ny, nz, owned_x_begin, owned_x_end, gpu_compute_stream());
  if (boundary_mode == EpithelialBoundaryMode::Dirichlet) {
    gpu::launch_set_epithelial_boundary(
        d_conc, nx, ny, owned_x_begin, owned_x_end, spec.boundary_conc, 0.0,
        nullptr, gpu_compute_stream());
  }
  return true;
}
#endif

}  // namespace

int diffusion_z_line_length(
    const Domain& domain, EpithelialBoundaryMode mode) {
  return mode == EpithelialBoundaryMode::Dirichlet
      ? domain.nz() - 1 : domain.nz();
}

bool diffusion_line_lengths_within(
    const Domain& domain, EpithelialBoundaryMode mode, int max_line) {
  if (max_line <= 0) return false;
  const int z_line = diffusion_z_line_length(domain, mode);
  return domain.nx() <= max_line && domain.ny() <= max_line
      && z_line <= max_line;
}

bool diffusion_all_species_within(
    const Domain& domain, const ChemicalField& field, int max_line) {
  for (Int s = 0; s < field.num_species(); ++s) {
    const ChemicalSpec& spec = field.spec(s);
    if (spec.diffuses() && !diffusion_line_lengths_within(
                                domain, spec.epithelial_boundary_mode, max_line)) {
      return false;
    }
  }
  return true;
}

bool delivery_route_b_eligible(
    const Domain& domain, const ChemicalField& field, int max_line) {
  if (max_line <= 0 || domain.nprocs() != 1 || field.slab_mode()) return false;
  for (Int s = 0; s < field.num_species(); ++s) {
    const ChemicalSpec& spec = field.spec(s);
    if (spec.diffuses() && spec.delivery_enabled
        && !diffusion_line_lengths_within(
            domain, spec.epithelial_boundary_mode, max_line)) {
      return false;
    }
  }
  return true;
}

bool gpu_diffusion_line_lengths_supported(
    const Domain& domain, EpithelialBoundaryMode mode) {
#ifdef GUTIBM_CUDA
  const int max_line = gpu::diffusion_max_line_length();
  return diffusion_line_lengths_within(domain, mode, max_line);
#else
  (void)domain;
  (void)mode;
  return false;
#endif
}

bool gpu_apply_species_diffusion_device(const Domain& domain,
                                        const ChemicalSpec& spec,
                                        double* d_conc,
                                        double* d_injected_amount,
                                        Real dt) {
#ifndef GUTIBM_CUDA
  (void)domain;
  (void)spec;
  (void)d_conc;
  (void)d_injected_amount;
  (void)dt;
  return false;
#else
  if (!gpu_runtime_enabled()) return false;
  if (!species_diffusion_eligible(spec, dt, domain)) return false;
  if (!gpu_diffusion_line_lengths_supported(
          domain, spec.epithelial_boundary_mode)) return false;
  return apply_species_diffusion_on_device(
      domain, spec, d_conc, d_injected_amount, dt, domain.nx(), 0,
      domain.nx(), DiffusionPhase::Replicated);
#endif
}

bool gpu_apply_species_delivery_device(
    const Domain& domain, const ChemicalSpec& spec, double* d_conc,
    const double* d_sink, const double* d_prescribed, double* d_realized,
    double* d_boundary_injected, double* d_gradient_source, Real dt,
    bool prescribed_active) {
#ifndef GUTIBM_CUDA
  (void)domain;
  (void)spec;
  (void)d_conc;
  (void)d_sink;
  (void)d_prescribed;
  (void)d_realized;
  (void)d_boundary_injected;
  (void)d_gradient_source;
  (void)dt;
  (void)prescribed_active;
  return false;
#else
  // Device rationing has no MPI reduction; Route B is single-rank only.
  assert(domain.nprocs() == 1);
  const int max_line = gpu::kMaxDeliveryLineLength;
  if (!gpu_runtime_enabled() || !spec.delivery_enabled
      || !species_diffusion_eligible(spec, dt, domain)
      || domain.nprocs() != 1 || domain.nx() > max_line
      || domain.ny() > max_line
      || !diffusion_line_lengths_within(
          domain, spec.epithelial_boundary_mode, max_line)) {
    return false;
  }
  const Real effective_diffusion = spec.diff_coeff / spec.retardation;
  const Real alpha_x = effective_diffusion * dt
      / (domain.dx_x() * domain.dx_x());
  const Real alpha_y = effective_diffusion * dt
      / (domain.dx_y() * domain.dx_y());
  const Real alpha_z = effective_diffusion * dt
      / (domain.dx_z() * domain.dx_z());
  const bool preserve_gradient =
      spec.z_gradient_enabled && spec.z_gradient_lambda > 0.0;
  const Real diffusion_boundary =
      preserve_gradient ? 0.0 : spec.boundary_conc;
  const int mode = to_underlying(spec.epithelial_boundary_mode);
  const Real beta = spec.epithelial_boundary_mode
      == EpithelialBoundaryMode::Robin
      ? spec.epithelial_transfer_coeff * dt / domain.dx_z() : 0.0;
  const Real flux_source = spec.epithelial_boundary_mode
      == EpithelialBoundaryMode::Flux
      ? spec.epithelial_flux * dt / domain.dx_z() : 0.0;
  const auto stream = gpu_compute_stream();
  if (spec.epithelial_boundary_mode == EpithelialBoundaryMode::Dirichlet) {
    gpu::launch_set_epithelial_boundary(
        d_conc, domain.nx(), domain.ny(), 0, domain.nx(),
        spec.boundary_conc, domain.cell_volume(), d_boundary_injected, stream);
  }
  if (preserve_gradient) {
    if (d_gradient_source != nullptr) {
      gpu::launch_set_luminal_neumann_accounted(
          d_conc, domain.nx(), domain.ny(), domain.nz(), 0, domain.nx(),
          domain.cell_volume(), d_gradient_source, stream);
      gpu::launch_shift_z_gradient_accounted(
          d_conc, domain.nx(), domain.ny(), domain.nz(), 0, domain.nx(),
          domain.dx_z(), spec.initial_conc, spec.z_gradient_lambda,
          spec.boundary_conc, -1.0, domain.cell_volume(), d_gradient_source,
          stream);
    } else {
      gpu::launch_set_luminal_neumann(
          d_conc, domain.nx(), domain.ny(), domain.nz(), 0, domain.nx(),
          stream);
      gpu::launch_shift_z_gradient(
          d_conc, domain.nx(), domain.ny(), domain.nz(), 0, domain.nx(),
          domain.dx_z(), spec.initial_conc, spec.z_gradient_lambda,
          spec.boundary_conc, -1.0, stream);
    }
  }
  if (!gpu::launch_diffuse_x_periodic_delivery(
          d_conc, d_sink, d_prescribed, d_realized, domain.nx(), domain.ny(),
          domain.nz(), alpha_x, dt / 3.0, domain.cell_volume(),
          domain.dx_z(), spec.initial_conc, spec.z_gradient_lambda,
          spec.boundary_conc, preserve_gradient, stream)) {
    return false;
  }
  if (!gpu::launch_diffuse_y_periodic_delivery(
          d_conc, d_sink, d_prescribed, d_realized, domain.nx(), domain.ny(),
          domain.nz(), 0, domain.nx(), alpha_y, dt / 3.0,
          domain.cell_volume(), domain.dx_z(), spec.initial_conc,
          spec.z_gradient_lambda, spec.boundary_conc, preserve_gradient,
          stream)) {
    return false;
  }
  if (!gpu::launch_diffuse_z_bounded_delivery(
          d_conc, d_sink, d_prescribed, d_realized, domain.nx(), domain.ny(),
          domain.nz(), 0, domain.nx(), alpha_z, mode, diffusion_boundary,
          spec.boundary_conc, beta, flux_source, dt / 3.0,
          domain.cell_volume(), domain.dx_z(), spec.initial_conc,
          spec.z_gradient_lambda, preserve_gradient, d_boundary_injected,
          stream)) {
    return false;
  }
  if (preserve_gradient) {
    if (d_gradient_source != nullptr) {
      gpu::launch_shift_z_gradient_accounted(
          d_conc, domain.nx(), domain.ny(), domain.nz(), 0, domain.nx(),
          domain.dx_z(), spec.initial_conc, spec.z_gradient_lambda,
          spec.boundary_conc, 1.0, domain.cell_volume(), d_gradient_source,
          stream);
      gpu::launch_set_luminal_neumann_accounted(
          d_conc, domain.nx(), domain.ny(), domain.nz(), 0, domain.nx(),
          domain.cell_volume(), d_gradient_source, stream);
    } else {
      gpu::launch_shift_z_gradient(
          d_conc, domain.nx(), domain.ny(), domain.nz(), 0, domain.nx(),
          domain.dx_z(), spec.initial_conc, spec.z_gradient_lambda,
          spec.boundary_conc, 1.0, stream);
      gpu::launch_set_luminal_neumann(
          d_conc, domain.nx(), domain.ny(), domain.nz(), 0, domain.nx(),
          stream);
    }
  }
  if (!prescribed_active) {
    gpu::launch_clamp_nonneg(
        d_conc, domain.nx(), domain.ny(), domain.nz(), 0, domain.nx(), stream);
  }
  if (spec.epithelial_boundary_mode == EpithelialBoundaryMode::Dirichlet) {
    gpu::launch_set_epithelial_boundary(
        d_conc, domain.nx(), domain.ny(), 0, domain.nx(),
        spec.boundary_conc, domain.cell_volume(), d_boundary_injected, stream);
  }
  gpu_sync_compute();
  gpu_check_error("gpu_apply_species_delivery_device");
  return true;
#endif
}

bool gpu_apply_species_diffusion_slab_device(
    const Domain& domain, const ChemicalSpec& spec, double* d_conc,
    double* d_injected_amount, Real dt, GpuSlabDiffusionContext& context) {
#ifndef GUTIBM_CUDA
  (void)domain;
  (void)spec;
  (void)d_conc;
  (void)d_injected_amount;
  (void)dt;
  (void)context;
  return false;
#else
  if (!gpu_runtime_enabled() || !species_diffusion_eligible(spec, dt, domain)
      || !gpu_diffusion_line_lengths_supported(
             domain, spec.epithelial_boundary_mode)) {
    return false;
  }
  if (!apply_species_diffusion_on_device(
          domain, spec, d_conc, d_injected_amount, dt, context.storage_nx,
          context.owned_storage_x_begin, context.owned_storage_x_end,
          DiffusionPhase::SlabPreX)) {
    return false;
  }
  gpu_sync_compute();
  gpu_check_error("slab diffusion boundary");
  const size_t count = static_cast<size_t>(context.storage_nx) * domain.ny()
      * domain.nz();
  const auto start = std::chrono::steady_clock::now();
  auto& host_concentration =
      context.field.mutable_species_concentration(context.spec_index);
  if (const cudaError_t d2h_status = cudaMemcpy(
          host_concentration.data(), d_conc, count * sizeof(double),
          cudaMemcpyDeviceToHost);
      d2h_status != cudaSuccess) {
    throw Error(std::string("slab diffusion D2H: ")
                + cudaGetErrorString(d2h_status));
  }
  gpu_transfer_record_d2h(
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
          .count());
  context.field.apply_periodic_x_diffusion(
      domain, dt, context.spec_index);
  const auto host_done = std::chrono::steady_clock::now();
  if (const cudaError_t h2d_status = cudaMemcpy(
          d_conc, host_concentration.data(), count * sizeof(double),
          cudaMemcpyHostToDevice);
      h2d_status != cudaSuccess) {
    throw Error(std::string("slab diffusion H2D: ")
                + cudaGetErrorString(h2d_status));
  }
  gpu_transfer_record_h2d(
      std::chrono::duration<double>(std::chrono::steady_clock::now() - host_done)
          .count());
  gpu_transfer_record_slab_x_roundtrip(
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
          .count());
  return apply_species_diffusion_on_device(
      domain, spec, d_conc, d_injected_amount, dt, context.storage_nx,
      context.owned_storage_x_begin, context.owned_storage_x_end,
      DiffusionPhase::SlabPostX);
#endif
}

bool gpu_apply_species_diffusion(const Domain& domain,
                                 const ChemicalSpec& spec,
                                 std::vector<Real>& concentration,
                                 Real dt) {
#ifndef GUTIBM_CUDA
  (void)domain;
  (void)spec;
  (void)concentration;
  (void)dt;
  return false;
#else
  if (!gpu_runtime_enabled()) return false;
  if (!species_diffusion_eligible(spec, dt, domain)) return false;
  if (!gpu_diffusion_line_lengths_supported(
          domain, spec.epithelial_boundary_mode)) return false;

  const int ncells = domain.ncells();
  if (ncells <= 0 || static_cast<int>(concentration.size()) < ncells) return false;

  DeviceBuffer<double> d_conc;
  d_conc.upload(concentration);
  if (!apply_species_diffusion_on_device(
          domain, spec, d_conc.data(), nullptr, dt, domain.nx(), 0, domain.nx(),
          DiffusionPhase::Replicated)) {
    return false;
  }

  gpu_sync_compute();
  gpu_check_error("gpu_apply_species_diffusion");

  d_conc.download(concentration.data(), static_cast<size_t>(ncells));
  return true;
#endif
}

}  // namespace gutibm
