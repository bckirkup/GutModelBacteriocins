/* -----------------------------------------------------------------------
   GutIBM – VBF velocity-channel characterization test
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"

#include <cassert>
#include <iostream>

using namespace gutibm;

void test_vbf_velocity_remains_zero() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {60e-6, 60e-6, 40e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.advection.mucus_thickness = 40e-6;
  cfg.advection.distal_length = 60e-6;
  cfg.seed = 4242;
  cfg.hdf5.enabled = false;
  cfg.time.output_interval = 3600.0;
  cfg.time.total_time = 3600.0;
  cfg.initial_strains.clear();

  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 200;
  strain.mu_max = 5.5e-4;
  cfg.initial_strains.push_back(strain);
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);

  const Real dt = cfg.time.bio_dt;
  const int steps = static_cast<int>(cfg.time.total_time / dt);
  for (int step = 0; step < steps; ++step) {
    sim.step(dt);
    for (const Agent& agent : sim.agents()) {
      for (int k = 0; k < 3; ++k) {
        assert(agent.v[k] == 0.0);
      }
    }
  }

  std::cout << "  test_vbf_velocity_remains_zero: PASSED\n";
}

int main() {
  std::cout << "=== VBF Velocity Channel Test ===\n";
  // Characterization/invariant test: this holds before and after removal
  // because the measured velocity channel is numerically inert.
  test_vbf_velocity_remains_zero();
  std::cout << "All VBF velocity channel tests passed.\n";
  return 0;
}
