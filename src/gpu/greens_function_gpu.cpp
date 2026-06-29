#include "greens_function_gpu.h"
#include "dispatch.h"
#include "device_memory.h"
#include "gpu_types.h"
#include "gpu_kernels.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

#include "domain.h"
#include "advection.h"
#include <vector>
#include <cmath>

namespace gutibm {
namespace {

#ifdef GUTIBM_CUDA
gpu::DomainParams make_domain_params(const Domain& domain) {
  gpu::DomainParams p{};
  p.nx = domain.nx();
  p.ny = domain.ny();
  p.nz = domain.nz();
  p.dx = domain.dx();
  p.lo[0] = domain.lo()[0];
  p.lo[1] = domain.lo()[1];
  p.lo[2] = domain.lo()[2];
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
#endif

}  // namespace

bool gpu_superpose_to_grid(
    const Domain& domain,
    const AdvectionField& adv,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    std::vector<Real>& grid_conc,
    Real cutoff_radius) {

#ifndef GUTIBM_CUDA
  (void)domain;
  (void)adv;
  (void)sources;
  (void)params;
  (void)grid_conc;
  (void)cutoff_radius;
  return false;
#else
  if (!gpu_runtime_enabled()) return false;
  if (sources.empty()) return true;

  Int ncells = domain.ncells();
  grid_conc.assign(ncells, 0.0);

  std::vector<double> sx(sources.size());
  std::vector<double> sy(sources.size());
  std::vector<double> sz(sources.size());
  std::vector<gpu::GfSourceParams> sp(params.size());
  for (size_t i = 0; i < sources.size(); ++i) {
    sx[i] = sources[i][0];
    sy[i] = sources[i][1];
    sz[i] = sources[i][2];
    sp[i].diff_coeff = params[i].diff_coeff;
    sp[i].source_rate = params[i].source_rate;
    sp[i].retardation = params[i].retardation;
  }

  DeviceBuffer<double> d_sx;
  DeviceBuffer<double> d_sy;
  DeviceBuffer<double> d_sz;
  DeviceBuffer<gpu::GfSourceParams> d_params;
  DeviceBuffer<double> d_grid;
  d_sx.upload(sx);
  d_sy.upload(sy);
  d_sz.upload(sz);
  d_params.upload(sp);
  d_grid.allocate(static_cast<size_t>(ncells));
  cudaMemset(d_grid.data(), 0, static_cast<size_t>(ncells) * sizeof(double));

  auto span = static_cast<int>(std::ceil(cutoff_radius / domain.dx()));
  auto dom = make_domain_params(domain);
  auto adv_p = make_advection_params(adv);

  gpu::launch_superpose_kernel(
      d_sx.data(), d_sy.data(), d_sz.data(), d_params.data(), d_grid.data(),
      dom, adv_p, static_cast<int>(sources.size()), span, nullptr);

  cudaDeviceSynchronize();
  gpu_check_error("superpose_kernel");

  d_grid.download(grid_conc);
  return true;
#endif
}

}  // namespace gutibm
