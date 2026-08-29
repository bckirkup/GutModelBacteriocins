/* -----------------------------------------------------------------------
   GutIBM – CPU/GPU parity for the shared Neumann image series
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "dispatch.h"
#include "domain.h"
#include "gpu_test_support.h"
#include "greens_function.h"
#include "greens_function_gpu.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

int main() {
  std::cout << "=== Neumann Image-Series GPU Parity Test ===\n";
  const int gpu_status = test::require_gpu("neumann_image_series_gpu");
  if (gpu_status != 0) return gpu_status;

  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {200.0e-6, 200.0e-6, 100.0e-6};
  dcfg.grid_dx = 10.0e-6;
  Domain domain;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 1.0e20;
  acfg.distal_transit_time = 1.0e20;
  acfg.mucus_thickness = 100.0e-6;
  acfg.distal_length = 200.0e-6;
  acfg.taylor_aris_enabled = false;
  AdvectionField adv;
  adv.init(acfg, domain);

  GreensFunction gf;
  gf.init(domain, adv);
  const std::vector<Vec3> sources = {
      {80.0e-6, 90.0e-6, 23.0e-6},
      {140.0e-6, 110.0e-6, 71.0e-6}};
  std::vector<GreensFunctionParams> params(2);
  for (auto& param : params) {
    param.diff_coeff = 4.0e-11;
    param.retardation = 1.0;
    param.source_rate = 1.0e-18;
    param.decay_rate = 1.0e-4;
  }

#ifndef GUTIBM_CUDA
  std::cout << "  test_neumann_image_series_gpu: SKIPPED (CUDA not compiled in)\n";
  return 77;
#else
  GpuConfig gpu_config;
  gpu_config.enabled = true;
  gpu_config.device_id = 0;
  gpu_set_config(gpu_config);
  if (!gpu_init_for_rank(0, 1)) return 77;

  std::vector<Real> cpu_grid;
  std::vector<Real> gpu_grid;
  gpu_config.enabled = false;
  gpu_set_config(gpu_config);
  gf.superpose_to_grid(sources, params, cpu_grid, 80.0e-6);
  gpu_config.enabled = true;
  gpu_set_config(gpu_config);
  if (!gpu_superpose_to_grid(
          domain, adv, sources, params,
          std::vector<Real>(sources.size(), 1.0), gpu_grid, 80.0e-6)) {
    std::cerr << "GPU Green's function evaluation failed\n";
    return 1;
  }

  Real max_relative_error = 0.0;
  for (size_t i = 0; i < cpu_grid.size(); ++i) {
    const Real denominator = std::max(std::abs(cpu_grid[i]), 1.0e-30);
    max_relative_error = std::max(
        max_relative_error, std::abs(cpu_grid[i] - gpu_grid[i]) / denominator);
  }
  if (!(max_relative_error < 1.0e-4)) {
    std::cerr << "CPU/GPU mismatch: " << max_relative_error << "\n";
    return 1;
  }
  std::cout << "  test_neumann_image_series_gpu: PASSED (max relative error="
            << max_relative_error << ")\n";

  for (auto& param : params) {
    param.decay_rate = 0.0;
  }
  gf.reset_image_series_cap_hits();
  gpu_config.enabled = false;
  gpu_set_config(gpu_config);
  std::vector<Real> unscreened_cpu_grid;
  gf.superpose_to_grid(
      sources, params, unscreened_cpu_grid, 80.0e-6);
  const uint64_t host_cap_hits = gf.image_series_cap_hits();

  gpu_config.enabled = true;
  gpu_set_config(gpu_config);
  std::vector<Real> unscreened_gpu_grid;
  uint64_t device_cap_hits = 0;
  if (!gpu_superpose_to_grid(
          domain, adv, sources, params,
          std::vector<Real>(sources.size(), 1.0), unscreened_gpu_grid,
          80.0e-6, &device_cap_hits)) {
    std::cerr << "GPU unscreened Green's function evaluation failed\n";
    return 1;
  }
  if (host_cap_hits == 0 || device_cap_hits == 0 ||
      host_cap_hits != device_cap_hits) {
    std::cerr << "CPU/GPU cap-hit mismatch: host=" << host_cap_hits
              << " device=" << device_cap_hits << "\n";
    return 1;
  }

  Real unscreened_max_relative_error = 0.0;
  for (size_t i = 0; i < unscreened_cpu_grid.size(); ++i) {
    const Real denominator =
        std::max(std::abs(unscreened_cpu_grid[i]), 1.0e-30);
    unscreened_max_relative_error = std::max(
        unscreened_max_relative_error,
        std::abs(unscreened_cpu_grid[i] - unscreened_gpu_grid[i]) /
            denominator);
  }
  if (!(unscreened_max_relative_error < 1.0e-4)) {
    std::cerr << "CPU/GPU unscreened mismatch: "
              << unscreened_max_relative_error << "\n";
    return 1;
  }
  std::cout << "  forced-cap parity: PASSED (cap hits=" << host_cap_hits
            << ", max relative error=" << unscreened_max_relative_error
            << ")\n";
  return 0;
#endif
}
