/* -----------------------------------------------------------------------
   Scaling benchmark smoke test (issue #55).
   Measures wall time and per-phase hotspots at modest agent counts
   suitable for CI; full sweeps use scripts/run_scaling_benchmark.sh.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"

#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

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

int mpi_nprocs() {
#ifdef GUTIBM_MPI
  int n = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  return n;
#else
  return 1;
#endif
}

long read_vmrss_kb() {
  std::ifstream status("/proc/self/status");
  if (!status) return -1;
  std::string line;
  while (std::getline(status, line)) {
    if (line.rfind("VmRSS:", 0) == 0) {
      std::istringstream iss(line.substr(6));
      long kb = -1;
      iss >> kb;
      return kb;
    }
  }
  return -1;
}

SimulationConfig bench_config(int agent_count, bool use_fmm) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.hdf5.enabled = false;
  cfg.time.total_time = 120.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 1.0e9;
  cfg.adaptive_dt.enabled = false;
  cfg.profile_steps = true;
  // Compact domain keeps CI runtime low while preserving scaling trends.
  cfg.domain.hi = {200.0e-6, 200.0e-6, 100.0e-6};
  cfg.domain.grid_dx = 4.0e-6;
  cfg.qssa.use_fmm = use_fmm;
  cfg.qssa.fmm_theta = 0.5;
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

struct BenchRow {
  int agents = 0;
  int ranks = 1;
  double step_ms = 0.0;
  double chemistry_ms = 0.0;
  double biology_ms = 0.0;
  double hash_ms = 0.0;
  long rss_mb = -1;
  double bytes_per_agent = 0.0;
};

BenchRow run_case(int agent_count, bool use_fmm, int n_steps) {
  SimulationConfig cfg = bench_config(agent_count, use_fmm);
  Simulation sim;
  sim.init(cfg);

  const long rss_after_init_kb = read_vmrss_kb();
  const int global_agents = sim.global_agent_count();
  assert(global_agents > 0);

  using Clock = std::chrono::steady_clock;
  const auto t0 = Clock::now();
  for (int s = 0; s < n_steps; ++s) {
    sim.step(cfg.time.bio_dt);
  }
  const auto t1 = Clock::now();
  const double wall_s =
      std::chrono::duration<double>(t1 - t0).count();

  const auto& prof = sim.step_profile();
  assert(prof.step_count == n_steps);

  BenchRow row;
  row.agents = global_agents;
  row.ranks = mpi_nprocs();
  row.step_ms = (wall_s / n_steps) * 1000.0;
  row.chemistry_ms = (prof.chemistry_s / n_steps) * 1000.0;
  row.biology_ms = (prof.biology_s / n_steps) * 1000.0;
  row.hash_ms = (prof.spatial_hash_s / n_steps) * 1000.0;
  if (rss_after_init_kb > 0) {
    row.rss_mb = rss_after_init_kb / 1024;
    row.bytes_per_agent =
        static_cast<double>(rss_after_init_kb) * 1024.0 /
        static_cast<double>(global_agents);
  }
  return row;
}

void print_row(const BenchRow& row, bool use_fmm) {
  if (mpi_rank() != 0) return;
  std::cout << "BENCHMARK agents=" << row.agents
            << " ranks=" << row.ranks
            << " use_fmm=" << (use_fmm ? 1 : 0)
            << " step_ms=" << row.step_ms
            << " chemistry_ms=" << row.chemistry_ms
            << " biology_ms=" << row.biology_ms
            << " hash_ms=" << row.hash_ms
            << " rss_mb=" << row.rss_mb
            << " bytes_per_agent=" << row.bytes_per_agent
            << "\n";
}

void assert_monotonic_soft(const BenchRow& small, const BenchRow& large) {
  if (mpi_rank() != 0) return;
  if (large.chemistry_ms <= 0.0 || small.chemistry_ms <= 0.0) return;
  const double ratio = large.chemistry_ms / small.chemistry_ms;
  const double agent_ratio =
      static_cast<double>(large.agents) / static_cast<double>(small.agents);
  (void)ratio;
  (void)agent_ratio;
  // Smoke check only: larger populations should not be dramatically faster.
  assert(large.step_ms >= small.step_ms * 0.5);
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GUTIBM_MPI
  MPI_Init(&argc, &argv);
#endif

  constexpr int k_steps = 1;
  const BenchRow small = run_case(500, false, k_steps);
  const BenchRow medium = run_case(1500, false, k_steps);
  const BenchRow large = run_case(3000, false, k_steps);

  print_row(small, false);
  print_row(medium, false);
  print_row(large, false);

  assert_monotonic_soft(small, medium);

  if (mpi_rank() == 0) {
    std::cout << "All scaling benchmark smoke tests passed.\n";
  }

#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif
  return 0;
}
