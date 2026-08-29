/* -----------------------------------------------------------------------
   GutIBM – CPU/GPU parity for the shared Robin correction table
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "dispatch.h"
#include "domain.h"
#include "gpu_test_support.h"
#include "greens_function.h"
#include "greens_function_gpu.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
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
  const auto run_large_working_set = [&](int count, Real cutoff_base) {
    std::vector<Vec3> many_sources;
    std::vector<GreensFunctionParams> many_params;
    many_sources.reserve(static_cast<size_t>(count));
    many_params.reserve(static_cast<size_t>(count));
    for (int i = 0; i < count; ++i) {
      many_sources.push_back(
          {80.0e-6 + (i % 10) * 2.0e-6,
           90.0e-6 + (i / 10) * 2.0e-6, 50.0e-6});
      auto param = params.front();
      param.robin_cutoff = cutoff_base + i * 0.1e-6;
      many_params.push_back(param);
    }
    std::vector<Real> many_grid;
    if (!gpu_superpose_to_grid(
            domain, adv, many_sources, many_params,
            std::vector<Real>(many_sources.size(), 1.0), many_grid,
            cutoff_base + count * 0.1e-6)) {
      std::cerr << "GPU Robin large working-set evaluation failed\n";
      std::exit(1);
    }
    for (const Real value : many_grid) {
      if (!std::isfinite(value)) {
        std::cerr << "GPU Robin large working-set result was not finite\n";
        std::exit(1);
      }
    }
  };

  run_large_working_set(65, 40.0e-6);
  run_large_working_set(128, 60.0e-6);

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

  const std::vector<Vec3> near_wall_sources = {
      {85.0e-6, 95.0e-6, 5.0e-6},
      {145.0e-6, 105.0e-6, 95.0e-6},
      {85.0e-6, 95.0e-6, 10.0e-6 - 1.0e-9},
      {145.0e-6, 105.0e-6, 10.0e-6 + 1.0e-9}};
  gpu_config.enabled = false;
  gpu_set_config(gpu_config);
  const uint64_t fallback_before = gf.robin_host_fallback_sources();
  gf.superpose_to_grid(
      near_wall_sources, params, cpu_grid, 80.0e-6);
  gpu_config.enabled = true;
  gpu_set_config(gpu_config);
  gf.superpose_to_grid(
      near_wall_sources, params, gpu_grid, 80.0e-6);
  if (gf.robin_host_fallback_sources() <= fallback_before) {
    std::cerr << "Robin near-wall host fallback was not exercised\n";
    return 1;
  }
  maximum = 0.0;
  for (size_t i = 0; i < cpu_grid.size(); ++i) {
    const Real denominator = std::max(std::abs(cpu_grid[i]), 1.0e-30);
    maximum = std::max(
        maximum, std::abs(cpu_grid[i] - gpu_grid[i]) / denominator);
  }
  if (!(maximum < 1.0e-4)) {
    std::cerr << "Near-wall Robin CPU/GPU mismatch: " << maximum << "\n";
    return 1;
  }
  std::cout << "  near-wall host fallback count="
            << (gf.robin_host_fallback_sources() - fallback_before)
            << " max relative error=" << maximum << "\n";

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
