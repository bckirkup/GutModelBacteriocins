#include "device.h"
#include "input_parser.h"
#include "simulation.h"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace gutibm;

namespace {

#ifdef GUTIBM_CUDA
SimulationConfig make_reproducibility_config(bool siderophore_enabled) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = 3004;
  cfg.time.total_time = 300.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 300.0;
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;
  cfg.chem_env.siderophore.enabled = siderophore_enabled;
  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {80e-6, 80e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 20;
  resident.mu_max = 5.5e-4;
  resident.plasmids = {"ColE1", "ColB"};
  resident.conjugative = true;
  cfg.initial_strains.push_back(resident);

  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = 10;
  immigrant.mu_max = 5.0e-4;
  immigrant.plasmids = {};
  immigrant.conjugative = false;
  cfg.initial_strains.push_back(immigrant);
  return cfg;
}

Simulation run_case(const SimulationConfig& cfg, bool gpu_enabled) {
  SimulationConfig run_cfg = cfg;
  run_cfg.gpu.enabled = gpu_enabled;
  run_cfg.gpu.device_id = 0;
  Simulation sim;
  sim.init(run_cfg);
  if (gpu_enabled) {
    assert(sim.gpu_active());
  }
  sim.run();
  return sim;
}

struct BiomassSummary {
  Int live = 0;
  Real total_biomass = 0.0;
};

BiomassSummary summarize(const Simulation& sim) {
  BiomassSummary summary;
  for (const Agent& agent : sim.agents()) {
    if (agent.state == PhenoState::DEAD) continue;
    ++summary.live;
    assert(std::isfinite(agent.biomass));
    assert(agent.biomass > 0.0);
    summary.total_biomass += agent.biomass;
  }
  assert(summary.live > 0);
  assert(summary.total_biomass > 0.0);
  return summary;
}

void print_difference(const char* label,
                      const BiomassSummary& first,
                      const BiomassSummary& second) {
  const Real absolute_difference =
      std::abs(first.total_biomass - second.total_biomass);
  const Real relative_difference =
      absolute_difference / first.total_biomass;
  std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
            << "  " << label
            << ": live_a=" << first.live
            << " live_b=" << second.live
            << " biomass_a=" << first.total_biomass
            << " biomass_b=" << second.total_biomass
            << " abs_diff=" << absolute_difference
            << " rel_diff=" << relative_difference << "\n";
}
#endif

}  // namespace

int main() {
  std::cout << "=== GPU Reproducibility Diagnostic ===\n";

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  if (DeviceContext::device_count() <= 0) {
    std::cout << "  SKIPPED (no CUDA device)\n";
    return 0;
  }

  const SimulationConfig fallback_cfg = make_reproducibility_config(true);
  const BiomassSummary fallback_gpu_a =
      summarize(run_case(fallback_cfg, true));
  const BiomassSummary fallback_gpu_b =
      summarize(run_case(fallback_cfg, true));
  print_difference("siderophore_on_gpu_repeat",
                   fallback_gpu_a, fallback_gpu_b);

  const BiomassSummary fallback_cpu_a =
      summarize(run_case(fallback_cfg, false));
  const BiomassSummary fallback_cpu_b =
      summarize(run_case(fallback_cfg, false));
  print_difference("siderophore_on_cpu_repeat",
                   fallback_cpu_a, fallback_cpu_b);
  assert(fallback_cpu_a.live == fallback_cpu_b.live);
  assert(fallback_cpu_a.total_biomass == fallback_cpu_b.total_biomass);

  const SimulationConfig metabolism_cfg =
      make_reproducibility_config(false);
  const BiomassSummary metabolism_cpu =
      summarize(run_case(metabolism_cfg, false));
  const BiomassSummary metabolism_gpu =
      summarize(run_case(metabolism_cfg, true));
  print_difference("siderophore_off_cpu_vs_gpu",
                   metabolism_cpu, metabolism_gpu);

  std::cout << "GPU reproducibility diagnostic complete.\n";
  return 0;
#endif
}
