/* -----------------------------------------------------------------------
   GutIBM – GPU Green's function parity test
   ----------------------------------------------------------------------- */

#include "greens_function.h"
#include "greens_function_gpu.h"
#include "dispatch.h"
#include "domain.h"
#include "gpu_test_support.h"
#include "advection.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

#ifdef GUTIBM_CUDA
static void cpu_superpose(const GreensFunction& gf,
                          const std::vector<Vec3>& sources,
                          const std::vector<GreensFunctionParams>& params,
                          std::vector<Real>& grid_conc,
                          Real cutoff) {
  gf.superpose_to_grid(sources, params, grid_conc, cutoff);
}
#endif

int main() {
  std::cout << "=== GPU Green's Function Tests ===\n";
  if (const int gpu_status = test::require_gpu("greens_function_gpu");
      gpu_status != 0) {
    return gpu_status;
  }

  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {200e-6, 200e-6, 50e-6};
  dcfg.grid_dx = 5e-6;

  Domain domain;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 5400.0;
  acfg.distal_transit_time = 43200.0;
  acfg.mucus_thickness = 50e-6;
  acfg.distal_length = 200e-6;

  AdvectionField adv;
  adv.init(acfg, domain);

  GreensFunction gf;
  gf.init(domain, adv);

  std::vector<Vec3> sources = {
      {50e-6, 50e-6, 25e-6},
      {120e-6, 80e-6, 25e-6},
  };
  std::vector<GreensFunctionParams> params(2);
  for (auto& p : params) {
    p.diff_coeff = 4e-11;
    p.source_rate = 1e-18;
    p.pI = 7.0;
    p.retardation = 2.0;
  }

  std::vector<Real> cpu_grid;
  std::vector<Real> gpu_grid;

#ifndef GUTIBM_CUDA
  std::cout << "  test_gpu_gf_parity: SKIPPED (CUDA not compiled in)\n";
  std::cout << "All GPU Green's function tests passed.\n";
  return 0;
#else
  Real cutoff = 80e-6;
  GpuConfig gcfg;
  gcfg.enabled = true;
  gcfg.device_id = 0;
  gpu_set_config(gcfg);

  if (!gpu_init_for_rank(0, 1)) return 1;

  // CPU reference (disable GPU dispatch temporarily)
  gcfg.enabled = false;
  gpu_set_config(gcfg);
  cpu_superpose(gf, sources, params, cpu_grid, cutoff);

  gcfg.enabled = true;
  gpu_set_config(gcfg);
  gpu_init_for_rank(0, 1);

  if (const auto strengths = std::vector<Real>(sources.size(), 1.0);
      !gpu_superpose_to_grid(domain, adv, sources, params, strengths, gpu_grid,
                             cutoff)) {
    std::cerr << "  test_gpu_gf_parity: FAILED (GPU kernel returned false)\n";
    return 1;
  }

  Real max_rel = 0.0;
  for (size_t i = 0; i < cpu_grid.size(); ++i) {
    Real denom = std::max(std::abs(cpu_grid[i]), 1e-30);
    Real rel = std::abs(cpu_grid[i] - gpu_grid[i]) / denom;
    max_rel = std::max(max_rel, rel);
  }

  if (max_rel > 1e-4) {
    std::cerr << "  test_gpu_gf_parity: FAILED (max_rel=" << max_rel << ")\n";
    return 1;
  }
  std::cout << "  test_gpu_gf_parity: PASSED (max_rel=" << max_rel << ")\n";

  // High-Pe_z ColE1 (retardation = 50.2) native GPU drift-corrected parity.
  DomainConfig hcfg;
  hcfg.lo = {0.0, 0.0, 0.0};
  hcfg.hi = {200.0e-6, 200.0e-6, 100.0e-6};
  hcfg.grid_dx = 10.0e-6;
  Domain high_domain;
  high_domain.init(hcfg);

  AdvectionConfig hacfg;
  hacfg.radial_turnover = 5400.0;
  hacfg.distal_transit_time = 43200.0;
  hacfg.mucus_thickness = 100.0e-6;
  hacfg.distal_length = 0.05;
  hacfg.profile_alpha = 1.5;
  AdvectionField high_adv;
  high_adv.init(hacfg, high_domain);

  GreensFunction high_gf;
  high_gf.init(high_domain, high_adv);

  const std::vector<Vec3> cole1_sources = {
      {80.0e-6, 80.0e-6, 10.0e-6},
      {100.0e-6, 100.0e-6, 50.0e-6},
      {120.0e-6, 120.0e-6, 90.0e-6},
  };
  std::vector<GreensFunctionParams> cole1_params(cole1_sources.size());
  for (auto& p : cole1_params) {
    p.diff_coeff = 4.0e-11;
    p.retardation = 50.2;
    p.source_rate = 1.0e-18;
    p.decay_rate = 5.0e-5;
    p.pI = 9.0;
    p.drift_correction = true;
  }

  std::vector<Real> cpu_drift_grid;
  std::vector<Real> gpu_drift_grid;
  gcfg.enabled = false;
  gpu_set_config(gcfg);
  high_gf.reset_low_screening_diagnostics();
  high_gf.superpose_to_grid(cole1_sources, cole1_params, cpu_drift_grid, 80.0e-6);

  gcfg.enabled = true;
  gpu_set_config(gcfg);
  high_gf.reset_low_screening_diagnostics();
  const uint64_t fallback_before = high_gf.robin_host_fallback_sources();
  high_gf.superpose_to_grid(cole1_sources, cole1_params, gpu_drift_grid, 80.0e-6);

  if (high_gf.robin_host_fallback_sources() != fallback_before) {
    std::cerr << "  test_gpu_colE1_drift_parity: FAILED (unexpected host fallback="
              << (high_gf.robin_host_fallback_sources() - fallback_before) << ")\n";
    return 1;
  }
  if (high_gf.drift_modal_evaluations() == 0) {
    std::cerr << "  test_gpu_colE1_drift_parity: FAILED (drift_modal_evaluations==0)\n";
    return 1;
  }
  Real num_sq = 0.0;
  Real den_sq = 0.0;
  for (size_t i = 0; i < cpu_drift_grid.size(); ++i) {
    const Real diff = cpu_drift_grid[i] - gpu_drift_grid[i];
    num_sq += diff * diff;
    den_sq += cpu_drift_grid[i] * cpu_drift_grid[i];
  }
  const Real drift_rel_l2 = std::sqrt(num_sq / std::max(den_sq, 1.0e-300));
  if (!(den_sq > 0.0) || drift_rel_l2 > 1.0e-6) {
    std::cerr << "  test_gpu_colE1_drift_parity: FAILED (rel_l2="
              << drift_rel_l2 << ")\n";
    return 1;
  }
  std::cout << "  test_gpu_colE1_drift_parity: PASSED (rel_l2="
            << drift_rel_l2 << ", modal_evals="
            << high_gf.drift_modal_evaluations() << ")\n";

  // Mixed ColE1 (Pe_z = 0.822 > 0.05) + ColE2 (R=1.2, Pe_z = 0.022 <= 0.05)
  // under drift_auto_dispatch = true.
  const std::vector<Vec3> mixed_sources = {
      {90.0e-6, 90.0e-6, 50.0e-6},
      {130.0e-6, 110.0e-6, 50.0e-6},
  };
  std::vector<GreensFunctionParams> mixed_params(2);
  mixed_params[0].diff_coeff = 4.0e-11;
  mixed_params[0].retardation = 50.2;
  mixed_params[0].source_rate = 1.0e-18;
  mixed_params[0].decay_rate = 5.0e-5;
  mixed_params[0].drift_auto_dispatch = true;
  mixed_params[0].peclet_z_threshold = 0.05;
  mixed_params[1].diff_coeff = 3.5e-11;
  mixed_params[1].retardation = 1.2;
  mixed_params[1].source_rate = 1.0e-18;
  mixed_params[1].decay_rate = 5.0e-5;
  mixed_params[1].drift_auto_dispatch = true;
  mixed_params[1].peclet_z_threshold = 0.05;

  std::vector<Real> cpu_auto_grid;
  std::vector<Real> gpu_auto_grid;
  gcfg.enabled = false;
  gpu_set_config(gcfg);
  high_gf.reset_low_screening_diagnostics();
  high_gf.superpose_to_grid(mixed_sources, mixed_params, cpu_auto_grid, 80.0e-6);

  gcfg.enabled = true;
  gpu_set_config(gcfg);
  high_gf.reset_low_screening_diagnostics();
  high_gf.superpose_to_grid(mixed_sources, mixed_params, gpu_auto_grid, 80.0e-6);

  if (high_gf.peclet_auto_dispatched_sources() != 1) {
    std::cerr << "  test_gpu_auto_dispatch_parity: FAILED (expected 1 dispatched, got "
              << high_gf.peclet_auto_dispatched_sources() << ")\n";
    return 1;
  }
  num_sq = 0.0;
  den_sq = 0.0;
  for (size_t i = 0; i < cpu_auto_grid.size(); ++i) {
    const Real diff = cpu_auto_grid[i] - gpu_auto_grid[i];
    num_sq += diff * diff;
    den_sq += cpu_auto_grid[i] * cpu_auto_grid[i];
  }
  const Real auto_rel_l2 = std::sqrt(num_sq / std::max(den_sq, 1.0e-300));
  if (!(den_sq > 0.0) || auto_rel_l2 > 1.0e-6) {
    std::cerr << "  test_gpu_auto_dispatch_parity: FAILED (rel_l2="
              << auto_rel_l2 << ")\n";
    return 1;
  }
  std::cout << "  test_gpu_auto_dispatch_parity: PASSED (rel_l2="
            << auto_rel_l2 << ", dispatched="
            << high_gf.peclet_auto_dispatched_sources() << ")\n";

  std::cout << "All GPU Green's function tests passed.\n";
  return 0;
#endif
}
