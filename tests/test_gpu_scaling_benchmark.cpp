/* -----------------------------------------------------------------------
   GPU scaling benchmark smoke test.
   Mirrors scripts/run_gpu_scaling_benchmark.sh at CI-friendly scale.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "dispatch.h"

#include <cassert>
#include <iostream>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace {

using gutibm::InputParser;
using gutibm::Simulation;
using gutibm::SimulationConfig;

int mpi_rank() {
#ifdef GUTIBM_MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
#else
  return 0;
#endif
}

SimulationConfig bench_config(int agent_count, bool gpu_enabled) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.hdf5.enabled = false;
  cfg.time.total_time = 120.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 1.0e9;
  cfg.adaptive_dt.enabled = false;
  cfg.profile_steps = true;
  cfg.gpu.enabled = gpu_enabled;
  cfg.domain.hi = {200.0e-6, 200.0e-6, 100.0e-6};
  cfg.domain.grid_dx = 4.0e-6;
  cfg.qssa.toxin_cutoff = 80.0e-6;
  cfg.qssa.nutrient_cutoff = 40.0e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = agent_count;
  strain.mu_max = 5.5e-4;
  strain.plasmids = {"ColE1", "ColB"};
  strain.conjugative = true;
  cfg.initial_strains.push_back(strain);
  return cfg;
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GUTIBM_MPI
  MPI_Init(&argc, &argv);
#endif

  if (mpi_rank() != 0) {
#ifdef GUTIBM_MPI
    MPI_Finalize();
#endif
    return 0;
  }

  std::cout << "=== GPU Scaling Benchmark Smoke ===\n";

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  constexpr int kAgents = 200;
  constexpr int kSteps = 2;

  SimulationConfig cfg_cpu = bench_config(kAgents, false);
  Simulation sim_cpu;
  sim_cpu.init(cfg_cpu);
  for (int s = 0; s < kSteps; ++s) {
    sim_cpu.step(cfg_cpu.time.bio_dt);
  }
  const auto& prof_cpu = sim_cpu.step_profile();
  assert(prof_cpu.step_count == kSteps);
  std::cout << "  cpu_chemistry_s=" << prof_cpu.chemistry_s / kSteps << "\n";

  gutibm::GpuConfig gcfg;
  gcfg.enabled = true;
  gcfg.device_id = 0;
  gutibm::gpu_set_config(gcfg);
  if (!gutibm::gpu_init_for_rank(0, 1)) {
    std::cout << "  gpu: SKIPPED (no CUDA device)\n";
    return 0;
  }

  SimulationConfig cfg_gpu = bench_config(kAgents, true);
  Simulation sim_gpu;
  sim_gpu.init(cfg_gpu);
  for (int s = 0; s < kSteps; ++s) {
    sim_gpu.step(cfg_gpu.time.bio_dt);
  }
  const auto& prof_gpu = sim_gpu.step_profile();
  assert(prof_gpu.step_count == kSteps);
  std::cout << "  gpu_chemistry_s=" << prof_gpu.chemistry_s / kSteps
            << " gpu_h2d_s=" << prof_gpu.gpu_h2d_s / kSteps
            << " gpu_d2h_s=" << prof_gpu.gpu_d2h_s / kSteps << "\n";

  std::cout << "  test_gpu_scaling_benchmark: PASSED\n";
#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif
  return 0;
#endif
}
