#include "diffusion_gpu.h"
#include "chemical_field.h"
#include "domain.h"
#include "dispatch.h"
#include "device_memory.h"
#include "gpu_kernels.h"
#include <vector>

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

namespace gutibm {

namespace {

class TridiagonalFactorization {
 public:
  void factorize(const std::vector<Real>& lower,
                 const std::vector<Real>& diagonal,
                 const std::vector<Real>& upper) {
    diagonal_ = diagonal;
    upper_ = upper;
    multipliers_.assign(lower.size(), 0.0);
    for (size_t i = 1; i < diagonal_.size(); ++i) {
      const Real multiplier = lower[i - 1] / diagonal_[i - 1];
      multipliers_[i - 1] = multiplier;
      diagonal_[i] -= multiplier * upper_[i - 1];
    }
  }

  void solve_in_place(std::vector<Real>& values) const {
    for (size_t i = 1; i < values.size(); ++i) {
      values[i] -= multipliers_[i - 1] * values[i - 1];
    }
    values.back() /= diagonal_.back();
    for (size_t i = values.size() - 1; i > 0; --i) {
      values[i - 1] =
          (values[i - 1] - upper_[i - 1] * values[i]) / diagonal_[i - 1];
    }
  }

 private:
  std::vector<Real> diagonal_;
  std::vector<Real> upper_;
  std::vector<Real> multipliers_;
};

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

  std::vector<Real> lower(static_cast<size_t>(n - 1), -alpha);
  std::vector<Real> upper(static_cast<size_t>(n - 1), -alpha);
  std::vector<Real> diagonal(static_cast<size_t>(n), diagonal_value);
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

bool line_lengths_supported(const Domain& domain) {
  const int max_line = gpu::diffusion_max_line_length();
  return domain.nx() <= max_line && domain.ny() <= max_line
      && (domain.nz() - 1) <= max_line;
}

}  // namespace

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
  if (dt <= 0.0 || domain.dx() <= 0.0) return false;
  if (!spec.diffusion_enabled || spec.diff_coeff <= 0.0
      || spec.retardation <= 0.0) {
    return false;
  }
  if (!line_lengths_supported(domain)) return false;

  const int nx = domain.nx();
  const int ny = domain.ny();
  const int nz = domain.nz();
  const int ncells = domain.ncells();
  if (ncells <= 0 || static_cast<int>(concentration.size()) < ncells) return false;

  const Real dx2 = domain.dx() * domain.dx();
  const Real alpha = (spec.diff_coeff / spec.retardation) * dt / dx2;
  const bool preserve_gradient =
      spec.z_gradient_enabled && spec.z_gradient_lambda > 0.0;
  double diffusion_boundary = spec.boundary_conc;

  DeviceBuffer<double> d_conc;
  d_conc.upload(concentration);

  gpu::launch_set_epithelial_boundary(
      d_conc.data(), nx, ny, spec.boundary_conc, nullptr);

  if (preserve_gradient) {
    gpu::launch_set_luminal_neumann(d_conc.data(), nx, ny, nz, nullptr);
    gpu::launch_shift_z_gradient(
        d_conc.data(), nx, ny, nz, domain.dx(), spec.initial_conc,
        spec.z_gradient_lambda, spec.boundary_conc, -1.0, nullptr);
    diffusion_boundary = 0.0;
  }

  const PeriodicPcrCoeffs x_coeffs = build_periodic_coeffs(nx, alpha);
  const PeriodicPcrCoeffs y_coeffs = build_periodic_coeffs(ny, alpha);
  DeviceBuffer<double> d_corr_x;
  DeviceBuffer<double> d_corr_y;
  d_corr_x.upload(x_coeffs.correction);
  d_corr_y.upload(y_coeffs.correction);

  gpu::launch_diffuse_x_periodic(
      d_conc.data(), nx, ny, nz, alpha,
      x_coeffs.gamma, x_coeffs.corner, x_coeffs.denominator,
      d_corr_x.data(), nullptr);
  gpu::launch_diffuse_y_periodic(
      d_conc.data(), nx, ny, nz, alpha,
      y_coeffs.gamma, y_coeffs.corner, y_coeffs.denominator,
      d_corr_y.data(), nullptr);
  gpu::launch_diffuse_z_bounded(
      d_conc.data(), nx, ny, nz, alpha, diffusion_boundary, nullptr);

  if (preserve_gradient) {
    gpu::launch_shift_z_gradient(
        d_conc.data(), nx, ny, nz, domain.dx(), spec.initial_conc,
        spec.z_gradient_lambda, spec.boundary_conc, 1.0, nullptr);
    gpu::launch_set_luminal_neumann(d_conc.data(), nx, ny, nz, nullptr);
  }

  gpu::launch_clamp_nonneg(d_conc.data(), ncells, nullptr);
  gpu::launch_set_epithelial_boundary(
      d_conc.data(), nx, ny, spec.boundary_conc, nullptr);

  cudaDeviceSynchronize();
  gpu_check_error("gpu_apply_species_diffusion");

  d_conc.download(concentration.data(), static_cast<size_t>(ncells));
  return true;
#endif
}

}  // namespace gutibm
