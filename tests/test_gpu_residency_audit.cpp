#include "gpu_test_support.h"
#include "input_parser.h"
#include "simulation.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace gutibm;

namespace {

#ifdef GUTIBM_CUDA
SimulationConfig make_audit_config() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = 4242;
  cfg.time.total_time = 180.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 180.0;
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;
  cfg.gpu.enabled = true;
  cfg.gpu.device_id = 0;
  cfg.chem_env.siderophore.enabled = false;
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.acetate.enabled = true;
  cfg.chem_env.mucin.enabled = true;
  cfg.chem_env.protease.enabled = true;
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {80.0e-6, 40.0e-6, 40.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.qssa.use_fmm = true;
  cfg.qssa.fmm_theta = 0.5;
  cfg.qssa.toxin_cutoff = 80.0e-6;
  cfg.qssa.nutrient_cutoff = 40.0e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 8;
  resident.mu_max = 5.5e-4;
  resident.plasmids = {"ColE1", "ColB"};
  resident.conjugative = true;
  cfg.initial_strains.push_back(resident);

  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = 4;
  immigrant.mu_max = 5.0e-4;
  cfg.initial_strains.push_back(immigrant);
  return cfg;
}

void assert_chemical_field_sane(const Simulation& simulation) {
  const auto& field = simulation.chemical_field();
  for (Int spec = 0; spec < field.num_species(); ++spec) {
    for (Int cell = 0; cell < field.ncells(); ++cell) {
      const Real value = field.conc(spec, cell);
      assert(std::isfinite(value));
      assert(value >= 0.0);
    }
  }
}
#endif

}  // namespace

int main() {
  std::cout << "=== GPU Concentration Residency Audit ===\n";
  if (const int gpu_status = test::require_gpu("gpu_residency_audit");
      gpu_status != 0) {
    return gpu_status;
  }

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  setenv("GUTIBM_GPU_RESIDENCY_AUDIT", "1", 1);
  Simulation simulation;
  simulation.init(make_audit_config());
  assert(simulation.gpu_active());
  simulation.run();
  assert_chemical_field_sane(simulation);
  std::cout << "GPU concentration residency audit completed three steps.\n";
  return 0;
#endif
}
