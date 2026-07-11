/* -----------------------------------------------------------------------
   QSSA near-field CPU vs GPU parity test
   ----------------------------------------------------------------------- */

#include "qssa_solver.h"
#include "greens_function_gpu.h"
#include "dispatch.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

static Simulation make_qssa_sim(bool gpu_enabled, uint64_t seed = 77) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.hdf5.enabled = false;
  cfg.gpu.enabled = gpu_enabled;
  cfg.domain.hi = {120.0e-6, 120.0e-6, 60.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.qssa.toxin_cutoff = 60.0e-6;
  cfg.qssa.use_fmm = false;
  cfg.seed = seed;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 12;
  strain.mu_max = 5.5e-4;
  strain.plasmids = {"ColE1", "ColB"};
  strain.conjugative = true;
  cfg.initial_strains.push_back(strain);

  Simulation sim;
  sim.init(cfg);
  return sim;
}

int main() {
  std::cout << "=== QSSA GPU Near-Field Parity ===\n";

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  GpuConfig gcfg;
  gcfg.enabled = true;
  gcfg.device_id = 0;
  gpu_set_config(gcfg);

  if (!gpu_init_for_rank(0, 1)) {
    std::cout << "  SKIPPED (no CUDA device)\n";
    return 0;
  }

  Simulation sim_gpu = make_qssa_sim(true, 77);
  Simulation sim_cpu = make_qssa_sim(false, 77);

  const Int idx = sim_cpu.chemical_field().find(species::BACTERIOCIN_BTUB);
  if (idx < 0) {
    std::cout << "  FAIL: bacteriocin species missing\n";
    return 1;
  }

  sim_cpu.step(sim_cpu.config().time.bio_dt);
  sim_gpu.step(sim_gpu.config().time.bio_dt);

  const auto& chem_cpu = sim_cpu.chemical_field();
  const auto& chem_gpu = sim_gpu.chemical_field();
  const Int ncells = chem_cpu.ncells();

  double max_abs = 0.0;
  double max_rel = 0.0;
  int compared = 0;
  for (Int c = 0; c < ncells; ++c) {
    const double cpu_v = chem_cpu.conc(idx, c);
    const double gpu_v = chem_gpu.conc(idx, c);
    if (cpu_v <= 0.0 && gpu_v <= 0.0) continue;
    ++compared;
    const double abs_err = std::abs(cpu_v - gpu_v);
    const double denom = std::max(std::abs(cpu_v), 1.0e-30);
    max_abs = std::max(max_abs, abs_err);
    max_rel = std::max(max_rel, abs_err / denom);
  }

  if (compared < 10) {
    std::cout << "  FAIL: insufficient nonzero cells (" << compared << ")\n";
    return 1;
  }

  constexpr double kAbsTol = 1.0e-12;
  constexpr double kRelTol = 1.0e-4;
  if (max_abs > kAbsTol && max_rel > kRelTol) {
    std::cout << "  FAIL: max_abs=" << max_abs << " max_rel=" << max_rel << "\n";
    return 1;
  }

  std::cout << "  test_qssa_gpu_near_field_parity: PASSED"
            << " (cells=" << compared << " max_rel=" << max_rel << ")\n";
  return 0;
#endif
}
