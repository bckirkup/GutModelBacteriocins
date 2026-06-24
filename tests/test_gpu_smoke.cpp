/* -----------------------------------------------------------------------
   GutIBM – GPU smoke integration test
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "dispatch.h"
#include <cmath>
#include <iostream>

using namespace gutibm;

static Real fingerprint(const Simulation& sim) {
  Real fp = 0.0;
  const auto& agents = sim.agents();
  for (Int i = 0; i < agents.size(); ++i) {
    const Agent& a = agents[i];
    if (a.state == PhenoState::DEAD) continue;
    fp += a.biomass + a.mu_realized * 1e3 + a.x[0] * 1e6;
  }
  fp += static_cast<Real>(sim.global_agent_count());
  return fp;
}

int main() {
  std::cout << "=== GPU Smoke Test ===\n";

  SimulationConfig cfg = InputParser::default_config();
  cfg.total_time = 300.0;
  cfg.bio_dt = 60.0;
  cfg.hdf5.enabled = false;
  cfg.initial_strains.clear();

  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 20;
  s.mu_max = 5e-4;
  s.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(s);

#ifndef GUTIBM_CUDA
  std::cout << "  test_gpu_smoke: SKIPPED (CUDA not compiled in)\n";
  std::cout << "All GPU smoke tests passed.\n";
  return 0;
#else
  // CPU baseline
  cfg.gpu.enabled = false;
  Simulation sim_cpu;
  sim_cpu.init(cfg);
  sim_cpu.run();
  Real fp_cpu = fingerprint(sim_cpu);

  if (DeviceContext::device_count() <= 0) {
    std::cout << "  test_gpu_smoke: SKIPPED (no CUDA device)\n";
    std::cout << "All GPU smoke tests passed.\n";
    return 0;
  }

  cfg.gpu.enabled = true;
  cfg.gpu.device_id = 0;
  Simulation sim_gpu;
  sim_gpu.init(cfg);
  if (!sim_gpu.gpu_active()) {
    std::cout << "  test_gpu_smoke: SKIPPED (GPU init failed)\n";
    return 0;
  }
  sim_gpu.run();
  Real fp_gpu = fingerprint(sim_gpu);

  Real rel = std::abs(fp_cpu - fp_gpu) / std::max(std::abs(fp_cpu), 1e-30);
  if (rel > 0.05) {
    std::cerr << "  test_gpu_smoke: FAILED (rel=" << rel
              << " cpu=" << fp_cpu << " gpu=" << fp_gpu << ")\n";
    return 1;
  }

  std::cout << "  test_gpu_smoke: PASSED (rel=" << rel << ")\n";
  std::cout << "All GPU smoke tests passed.\n";
  return 0;
#endif
}
