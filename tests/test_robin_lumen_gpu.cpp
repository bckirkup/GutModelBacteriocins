/* -----------------------------------------------------------------------
   GutIBM – CPU/GPU parity for the shared Robin correction table
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "domain.h"
#include "gpu_test_support.h"
#include "greens_function.h"
#include "greens_function_gpu.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using namespace gutibm;

int main() {
  std::cout << "=== Robin Lumen-Boundary GPU Parity Test ===\n";
  const int gpu_status = test::require_gpu("robin_lumen_gpu");
  if (gpu_status != 0) return gpu_status;

  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {200.0e-6, 200.0e-6, 100.0e-6};
  dcfg.grid_dx = 10.0e-6;
  Domain domain;
  domain.init(dcfg);
  AdvectionConfig acfg;
  acfg.radial_turnover = 5400.0;
  acfg.distal_transit_time = 43200.0;
  acfg.mucus_thickness = 100.0e-6;
  acfg.distal_length = 0.05;
  acfg.taylor_aris_enabled = false;
  AdvectionField adv;
  adv.init(acfg, domain);
  GreensFunction gf;
  gf.init(domain, adv);

#ifndef GUTIBM_CUDA
  std::cout << "  test_robin_lumen_gpu: SKIPPED (CUDA not compiled in)\n";
  return 77;
#else
  const std::vector<Vec3> sources = {
      {80.0e-6, 90.0e-6, 23.0e-6},
      {140.0e-6, 110.0e-6, 71.0e-6}};
  std::vector<GreensFunctionParams> params(2);
  for (auto& param : params) {
    param.diff_coeff = 4.0e-11;
    param.retardation = 1.0;
    param.source_rate = 1.0e-18;
    param.decay_rate = 1.0e-4;
    param.lumen_transfer_length = 100.0e-6;
    param.robin_cutoff = 80.0e-6;
  }
  std::vector<Real> cpu_grid;
  std::vector<Real> gpu_grid;
  GpuConfig gpu_config;
  gpu_config.enabled = true;
  gpu_set_config(gpu_config);
  if (!gpu_init_for_rank(0, 1)) return 77;
  gpu_config.enabled = false;
  gpu_set_config(gpu_config);
  gf.superpose_to_grid(sources, params, cpu_grid, 80.0e-6);
  gpu_config.enabled = true;
  gpu_set_config(gpu_config);
  if (!gpu_superpose_to_grid(
          domain, adv, sources, params,
          std::vector<Real>(sources.size(), 1.0), gpu_grid, 80.0e-6)) {
    std::cerr << "GPU Robin Green's-function evaluation failed\n";
    return 1;
  }
  Real maximum = 0.0;
  for (size_t i = 0; i < cpu_grid.size(); ++i) {
    const Real denominator = std::max(std::abs(cpu_grid[i]), 1.0e-30);
    maximum = std::max(
        maximum, std::abs(cpu_grid[i] - gpu_grid[i]) / denominator);
  }
  if (!(maximum < 1.0e-4)) {
    std::cerr << "Robin CPU/GPU mismatch: " << maximum << "\n";
    return 1;
  }

  for (auto& param : params) {
    param.lumen_transfer_length =
        std::numeric_limits<Real>::infinity();
  }
  gpu_config.enabled = false;
  gpu_set_config(gpu_config);
  gf.superpose_to_grid(sources, params, cpu_grid, 80.0e-6);
  gpu_config.enabled = true;
  gpu_set_config(gpu_config);
  if (!gpu_superpose_to_grid(
          domain, adv, sources, params,
          std::vector<Real>(sources.size(), 1.0), gpu_grid, 80.0e-6)) {
    std::cerr << "GPU sealed Green's-function evaluation failed\n";
    return 1;
  }
  maximum = 0.0;
  for (size_t i = 0; i < cpu_grid.size(); ++i) {
    const Real denominator = std::max(std::abs(cpu_grid[i]), 1.0e-30);
    maximum = std::max(
        maximum, std::abs(cpu_grid[i] - gpu_grid[i]) / denominator);
  }
  if (!(maximum < 1.0e-4)) {
    std::cerr << "Sealed CPU/GPU mismatch: " << maximum << "\n";
    return 1;
  }
  std::cout << "  test_robin_lumen_gpu: PASSED (max relative error="
            << maximum << ")\n";
  return 0;
#endif
}
