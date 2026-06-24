/* -----------------------------------------------------------------------
   GutIBM – GPU Green's function parity test
   ----------------------------------------------------------------------- */

#include "greens_function.h"
#include "greens_function_gpu.h"
#include "dispatch.h"
#include "domain.h"
#include "advection.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

static void cpu_superpose(const GreensFunction& gf,
                          const std::vector<Vec3>& sources,
                          const std::vector<GreensFunctionParams>& params,
                          std::vector<Real>& grid_conc,
                          Real cutoff) {
  gf.superpose_to_grid(sources, params, grid_conc, cutoff);
}

int main() {
  std::cout << "=== GPU Green's Function Tests ===\n";

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

  Real cutoff = 80e-6;
  std::vector<Real> cpu_grid, gpu_grid;

#ifndef GUTIBM_CUDA
  std::cout << "  test_gpu_gf_parity: SKIPPED (CUDA not compiled in)\n";
  std::cout << "All GPU Green's function tests passed.\n";
  return 0;
#else
  GpuConfig gcfg;
  gcfg.enabled = true;
  gcfg.device_id = 0;
  gpu_set_config(gcfg);

  if (!gpu_init_for_rank(0, 1)) {
    std::cout << "  test_gpu_gf_parity: SKIPPED (no CUDA device)\n";
    std::cout << "All GPU Green's function tests passed.\n";
    return 0;
  }

  // CPU reference (disable GPU dispatch temporarily)
  gcfg.enabled = false;
  gpu_set_config(gcfg);
  cpu_superpose(gf, sources, params, cpu_grid, cutoff);

  gcfg.enabled = true;
  gpu_set_config(gcfg);
  gpu_init_for_rank(0, 1);

  if (!gpu_superpose_to_grid(domain, adv, sources, params, gpu_grid, cutoff)) {
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
  std::cout << "All GPU Green's function tests passed.\n";
  return 0;
#endif
}
