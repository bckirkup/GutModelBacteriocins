/* -----------------------------------------------------------------------
   GutIBM -- Adaptive timestep tests
   Verifies CFL-like adaptive dt selection under various activity levels.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "plasmid.h"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace gutibm;

// Helper: build a minimal config with adaptive dt enabled
static SimulationConfig make_adaptive_cfg() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;

  cfg.total_time      = 120.0;
  cfg.bio_dt          = 60.0;
  cfg.output_interval = 120.0;
  cfg.seed            = 7777;
  cfg.hdf5.enabled    = false;

  cfg.advection.mucus_thickness      = 50e-6;
  cfg.advection.distal_length        = 100e-6;
  cfg.advection.radial_turnover      = 5400.0;
  cfg.advection.distal_transit_time  = 43200.0;
  cfg.qssa.toxin_cutoff    = 50e-6;
  cfg.qssa.nutrient_cutoff = 25e-6;

  cfg.adaptive_dt_enabled = true;
  cfg.dt_min           = 1.0;
  cfg.dt_max           = 300.0;
  cfg.dt_safety        = 0.8;
  cfg.dt_growth_limit  = 0.1;

  cfg.initial_strains.clear();
  return cfg;
}

void test_disabled_returns_fixed_dt() {
  SimulationConfig cfg = make_adaptive_cfg();
  cfg.adaptive_dt_enabled = false;

  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 5; s.mu_max = 5e-4;
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);

  Real dt = sim.compute_adaptive_dt();
  assert(std::abs(dt - cfg.bio_dt) < 1e-12);

  std::cout << "  test_disabled_returns_fixed_dt: PASSED (dt=" << dt << ")\n";
}

void test_bounds_respected() {
  SimulationConfig cfg = make_adaptive_cfg();
  cfg.dt_min = 5.0;
  cfg.dt_max = 100.0;

  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 10; s.mu_max = 5e-4;
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);

  Real dt = sim.compute_adaptive_dt();
  assert(dt >= cfg.dt_min);
  assert(dt <= cfg.dt_max);

  std::cout << "  test_bounds_respected: PASSED (dt=" << dt
            << " in [" << cfg.dt_min << "," << cfg.dt_max << "])\n";
}

void test_high_growth_reduces_dt() {
  // High mu_max should reduce dt via growth_limit / mu constraint
  SimulationConfig cfg = make_adaptive_cfg();

  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 10; s.mu_max = 0.1;  // very high growth rate
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);

  // After init, mu_realized may be less than mu_max, but compute_adaptive_dt
  // uses mu_realized.  Step once so metabolism sets mu_realized.
  sim.step(1.0);

  Real dt = sim.compute_adaptive_dt();
  // dt should be well below dt_max
  assert(dt < cfg.dt_max);
  assert(dt >= cfg.dt_min);

  // Also verify it's substantially smaller than a quiescent run
  SimulationConfig cfg2 = make_adaptive_cfg();
  SimulationConfig::InitialStrain s2;
  s2.type = 1; s2.count = 2; s2.mu_max = 1e-6;  // very slow growth
  s2.plasmids = {}; s2.conjugative = false;
  cfg2.initial_strains.push_back(s2);

  Simulation sim2;
  sim2.init(cfg2);
  sim2.step(1.0);
  Real dt2 = sim2.compute_adaptive_dt();

  assert(dt < dt2);

  std::cout << "  test_high_growth_reduces_dt: PASSED"
            << " (fast=" << dt << " slow=" << dt2 << ")\n";
}

void test_quiescent_uses_large_dt() {
  // Very slow growth, few agents -> dt should be near dt_max
  SimulationConfig cfg = make_adaptive_cfg();
  cfg.dt_safety = 1.0;  // no safety reduction for cleaner test

  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 2; s.mu_max = 1e-8;  // extremely slow
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);

  Real dt = sim.compute_adaptive_dt();
  // growth_limit / mu = 0.1 / 1e-8 = 1e7, clamped to dt_max=300
  assert(std::abs(dt - cfg.dt_max) < 1e-6);

  std::cout << "  test_quiescent_uses_large_dt: PASSED (dt=" << dt << ")\n";
}

void test_safety_factor_applied() {
  SimulationConfig cfg = make_adaptive_cfg();
  cfg.dt_safety = 0.5;

  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 2; s.mu_max = 1e-8;
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);

  Real dt = sim.compute_adaptive_dt();
  // Without safety: dt_max=300; with safety 0.5: 150
  assert(std::abs(dt - 150.0) < 1e-6);

  std::cout << "  test_safety_factor_applied: PASSED (dt=" << dt << ")\n";
}

void test_adaptive_run_completes() {
  SimulationConfig cfg = make_adaptive_cfg();
  cfg.total_time = 300.0;

  SimulationConfig::InitialStrain resident;
  resident.type = 1; resident.count = 20; resident.mu_max = 5e-4;
  resident.plasmids = {"colicin_E1"}; resident.conjugative = false;
  cfg.initial_strains.push_back(resident);

  SimulationConfig::InitialStrain target;
  target.type = 2; target.count = 10; target.mu_max = 5e-4;
  target.plasmids = {}; target.conjugative = false;
  cfg.initial_strains.push_back(target);

  Simulation sim;
  sim.init(cfg);
  sim.run();

  // Should reach total_time
  assert(std::abs(sim.time() - cfg.total_time) < 1e-6);
  assert(sim.step_count() > 0);

  // With adaptive dt the step count won't be ceil(total_time/bio_dt)
  // It can be more or fewer depending on activity
  std::cout << "  test_adaptive_run_completes: PASSED"
            << " (steps=" << sim.step_count()
            << " time=" << sim.time() << ")\n";
}

void test_fixed_dt_run_still_works() {
  // Verify the refactored run() loop still works with fixed dt
  SimulationConfig cfg = make_adaptive_cfg();
  cfg.adaptive_dt_enabled = false;
  cfg.total_time = 180.0;
  cfg.bio_dt = 60.0;

  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 10; s.mu_max = 5e-4;
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);
  sim.run();

  assert(std::abs(sim.time() - 180.0) < 1e-6);
  // 180 / 60 = 3 steps
  assert(sim.step_count() == 3);

  std::cout << "  test_fixed_dt_run_still_works: PASSED"
            << " (steps=" << sim.step_count() << ")\n";
}

int main() {
  std::cout << "=== Adaptive Timestep Tests ===\n";
  test_disabled_returns_fixed_dt();
  test_bounds_respected();
  test_high_growth_reduces_dt();
  test_quiescent_uses_large_dt();
  test_safety_factor_applied();
  test_adaptive_run_completes();
  test_fixed_dt_run_still_works();
  std::cout << "All adaptive timestep tests passed.\n";
  return 0;
}
