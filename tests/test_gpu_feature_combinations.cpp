/* -----------------------------------------------------------------------
   GutIBM – GPU feature-combination smoke (Spec 9 PR4)
   Runs selected Spec 8 combo scenarios with gpu_enabled and asserts
   finite chemistry + population outcomes.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "dispatch.h"
#include "device.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

namespace {

SimulationConfig make_combo_config(uint64_t seed) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = seed;
  cfg.time.total_time = 300.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 300.0;
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;
  cfg.gpu.enabled = true;
  cfg.gpu.device_id = 0;

  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {80e-6, 80e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;

  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length = 80e-6;
  cfg.qssa.toxin_cutoff = 40e-6;
  cfg.qssa.nutrient_cutoff = 20e-6;

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

void assert_chemistry_sane(const Simulation& sim) {
  const auto& chem = sim.chemical_field();
  for (Int s = 0; s < chem.num_species(); ++s) {
    for (Int c = 0; c < chem.ncells(); ++c) {
      const Real val = chem.conc(s, c);
      assert(std::isfinite(val));
      assert(val >= 0.0);
      assert(val < 1e8);
    }
  }
}

void assert_population_sane(const Simulation& sim, Int min_expected = 1) {
  Int live = 0;
  for (const Agent& a : sim.agents()) {
    if (a.state == PhenoState::DEAD) continue;
    assert(std::isfinite(a.biomass));
    assert(a.biomass > 0.0);
    ++live;
  }
  assert(live >= min_expected);
}

Simulation run_gpu_combo(const SimulationConfig& cfg) {
  Simulation sim;
  sim.init(cfg);
  assert(sim.gpu_active());
  sim.run();
  return sim;
}

}  // namespace

void test_gpu_full_chemical_environment() {
  SimulationConfig cfg = make_combo_config(3001);
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.acetate.enabled = true;
  cfg.chem_env.mucin.enabled = true;
  cfg.chem_env.protease.enabled = true;

  Simulation sim = run_gpu_combo(cfg);
  assert_chemistry_sane(sim);
  assert_population_sane(sim);

  std::cout << "  test_gpu_full_chemical_environment: PASSED\n";
}

void test_gpu_adaptive_dt_with_crypts() {
  SimulationConfig cfg = make_combo_config(3002);
  cfg.chem_env.oxygen.enabled = true;
  cfg.adaptive_dt.enabled = true;
  cfg.adaptive_dt.min = 1.0;
  cfg.adaptive_dt.max = 120.0;
  cfg.advection.crypts_enabled = true;

  Simulation sim = run_gpu_combo(cfg);
  assert_chemistry_sane(sim);
  assert_population_sane(sim);

  std::cout << "  test_gpu_adaptive_dt_with_crypts: PASSED\n";
}

void test_gpu_kitchen_sink_light() {
  SimulationConfig cfg = make_combo_config(3003);
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.acetate.enabled = true;
  cfg.chem_env.mucin.enabled = true;
  cfg.adaptive_dt.enabled = true;
  cfg.advection.crypts_enabled = true;
  cfg.advection.peristaltic_enabled = true;
  cfg.time.total_time = 120.0;

  Simulation sim = run_gpu_combo(cfg);
  assert_chemistry_sane(sim);
  assert_population_sane(sim);

  std::cout << "  test_gpu_kitchen_sink_light: PASSED\n";
}

int main() {
  std::cout << "=== GPU Feature Combination Smoke Tests ===\n";

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  std::cout << "All GPU feature-combination tests passed.\n";
  return 0;
#else
  if (DeviceContext::device_count() <= 0) {
    std::cout << "  SKIPPED (no CUDA device)\n";
    std::cout << "All GPU feature-combination tests passed.\n";
    return 0;
  }

  test_gpu_full_chemical_environment();
  test_gpu_adaptive_dt_with_crypts();
  test_gpu_kitchen_sink_light();

  std::cout << "All GPU feature-combination tests passed.\n";
  return 0;
#endif
}
