/* -----------------------------------------------------------------------
   GutIBM – Progress and wall-clock heartbeat lines for Batch observability
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "types.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

using namespace gutibm;

SimulationConfig progress_config(Real output_interval, uint64_t seed) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {50e-6, 50e-6, 30e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.time.total_time = 180.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = output_interval;
  cfg.seed = seed;
  cfg.hdf5.enabled = false;
  cfg.advection.mucus_thickness = 30e-6;
  cfg.advection.distal_length = 50e-6;
  cfg.advection.radial_turnover = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;
  cfg.initial_strains.clear();

  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 8;
  resident.mu_max = 5.0e-4;
  resident.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(resident);
  return cfg;
}

void test_progress_line_fields() {
  SimulationConfig cfg = progress_config(60.0, 7);

  Simulation sim;
  sim.init(cfg);

  std::stringstream captured;
  std::streambuf* old_out = std::cout.rdbuf(captured.rdbuf());
  sim.run();
  std::cout.rdbuf(old_out);

  const std::string out = captured.str();
  assert(out.find("pct=") != std::string::npos);
  assert(out.find("rate=") != std::string::npos);
  assert(out.find("eta_s=") != std::string::npos);
  assert(out.find("global_agents=") != std::string::npos);
}

void test_heartbeat_when_progress_is_unscheduled() {
  SimulationConfig cfg = progress_config(3600.0, 8);

  Simulation sim;
  sim.init(cfg);

  std::stringstream captured;
  std::streambuf* old_out = std::cout.rdbuf(captured.rdbuf());
  sim.run();
  std::cout.rdbuf(old_out);

  const std::string out = captured.str();
  assert(out.find("Heartbeat step=") != std::string::npos);
  assert(out.find("wall_elapsed_s=") != std::string::npos);
  assert(out.find("Step 60 ") == std::string::npos);
}

int main() {
  test_progress_line_fields();
  test_heartbeat_when_progress_is_unscheduled();
  std::cout << "All progress-report tests passed\n";
  return 0;
}
