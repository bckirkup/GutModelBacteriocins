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
#include <string_view>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace {

using gutibm::InputParser;
using gutibm::Simulation;
using gutibm::SimulationConfig;
using gutibm::ToxinBurstSource;

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

SimulationConfig bench_config(int agent_count, bool use_fmm,
                              std::string_view toxin_evaluation) {
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
  cfg.qssa.toxin_evaluation = toxin_evaluation;
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

SimulationConfig toxin_timing_config(std::string_view mode, int producer_count,
                                     int consumer_count, double domain_size,
                                     double grid_dx, double toxin_cutoff) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {domain_size, domain_size, domain_size};
  cfg.domain.grid_dx = grid_dx;
  cfg.qssa.toxin_cutoff = toxin_cutoff;
  cfg.qssa.toxin_evaluation = mode;
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain producers;
  producers.type = 1;
  producers.count = producer_count;
  producers.mu_max = 5.5e-4;
  producers.plasmids = {"ColB"};
  producers.conjugative = false;
  cfg.initial_strains.push_back(producers);
  SimulationConfig::InitialStrain consumers;
  consumers.type = 2;
  consumers.count = consumer_count;
  consumers.mu_max = 5.5e-4;
  consumers.conjugative = false;
  cfg.initial_strains.push_back(consumers);
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

BenchRow run_case(int agent_count, bool use_fmm, int n_steps,
                  const std::string& toxin_evaluation) {
  SimulationConfig cfg =
      bench_config(agent_count, use_fmm, toxin_evaluation);
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

void print_row(const BenchRow& row, bool use_fmm,
               const std::string& toxin_evaluation) {
  if (mpi_rank() != 0) return;
  std::cout << "BENCHMARK agents=" << row.agents
            << " ranks=" << row.ranks
            << " use_fmm=" << (use_fmm ? 1 : 0)
            << " toxin_evaluation=" << toxin_evaluation
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

ToxinBurstSource source_at(gutibm::Vec3 position) {
  ToxinBurstSource source;
  source.pos = position;
  source.creation_time = 0.0;
  source.release_tau = 300.0;
  source.target = gutibm::ReceptorType::BtuB;
  source.params.diff_coeff = 4e-11;
  source.params.retardation = 1.0;
  source.params.source_rate = 1e-18;
  source.params.decay_rate = 0.0;
  return source;
}

void run_toxin_timing_case() {
  constexpr int producer_count = 300;
  constexpr int consumer_count = 700;
  constexpr double domain_size = 200e-6;
  constexpr double grid_dx = 2e-6;
  constexpr double toxin_cutoff = 200e-6;

  std::vector<ToxinBurstSource> sources;
  sources.reserve(producer_count);
  for (int i = 0; i < producer_count; ++i) {
    const double x = (static_cast<double>(i % 10) + 0.5) * 20e-6;
    const double y = (static_cast<double>((i / 10) % 10) + 0.5) * 20e-6;
    const double z = (static_cast<double>(i / 100) + 0.5) * 20e-6;
    sources.push_back(source_at({x, y, z}));
  }

  Simulation grid_sim;
  Simulation agent_sim;
  const SimulationConfig grid_cfg = toxin_timing_config(
      "grid", producer_count, consumer_count, domain_size, grid_dx,
      toxin_cutoff);
  const SimulationConfig agent_cfg = toxin_timing_config(
      "agents", producer_count, consumer_count, domain_size, grid_dx,
      toxin_cutoff);
  grid_sim.init(grid_cfg);
  agent_sim.init(agent_cfg);

  using Clock = std::chrono::steady_clock;
  const auto grid_start = Clock::now();
  grid_sim.qssa().solve_all_bacteriocin_fields(
      grid_sim.agents(), sources, 0.0, grid_sim.config().chem_env.protease,
      grid_sim.advection(), grid_sim.chemical_field(), nullptr, true);
  const auto grid_end = Clock::now();
  const auto agent_start = Clock::now();
  agent_sim.qssa().solve_all_bacteriocin_fields(
      agent_sim.agents(), sources, 0.0, agent_sim.config().chem_env.protease,
      agent_sim.advection(), agent_sim.chemical_field(), nullptr, false);
  const auto agent_end = Clock::now();

  const double grid_ms =
      std::chrono::duration<double, std::milli>(grid_end - grid_start).count();
  const double agent_ms =
      std::chrono::duration<double, std::milli>(agent_end - agent_start)
          .count();
  if (mpi_rank() == 0) {
    std::cout << "TOXIN_TIMING label=producer_consumer"
              << " grid_ms=" << grid_ms
              << " agents_ms=" << agent_ms
              << " producers=" << producer_count
              << " consumers=" << consumer_count
              << " grid_dx_m=" << grid_dx
              << " toxin_cutoff_m=" << toxin_cutoff
              << " grid_cells="
              << grid_sim.chemical_field().global_ncells() << "\n";
  }
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GUTIBM_MPI
  MPI_Init(&argc, &argv);
#endif

  if (argc > 1 && std::string(argv[1]) == "--toxin-timing") {
    run_toxin_timing_case();
#ifdef GUTIBM_MPI
    MPI_Finalize();
#endif
    return 0;
  }

  constexpr int k_steps = 1;
  const BenchRow small = run_case(500, false, k_steps, "grid");
  const BenchRow medium = run_case(1500, false, k_steps, "grid");
  const BenchRow large = run_case(3000, false, k_steps, "grid");

  print_row(small, false, "grid");
  print_row(medium, false, "grid");
  print_row(large, false, "grid");

  assert_monotonic_soft(small, medium);

  if (mpi_rank() == 0) {
    std::cout << "All scaling benchmark smoke tests passed.\n";
  }

#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif
  return 0;
}
