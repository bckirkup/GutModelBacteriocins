/* -----------------------------------------------------------------------
   GutIBM – GPU production-path smoke (Spec 9 PR5 / issue #33)
   Exercises FMM hybrid (CPU far-field + GPU chemistry), production-scale
   domain sizing, and full chemical environment on the GPU path.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "dispatch.h"
#include "device.h"
#include "diffusion_gpu.h"
#include "domain.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

namespace {

SimulationConfig make_production_config(uint64_t seed) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = seed;
  cfg.time.total_time = 240.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 240.0;
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;
  cfg.gpu.enabled = true;
  cfg.gpu.device_id = 0;

  // diversity_paradox-scale domain (2 mm slab, 2 um resolution → nx=1000)
  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {2.0e-3, 2.0e-3, 100.0e-6};
  cfg.domain.grid_dx = 2.0e-6;
  cfg.domain.hash_cell_size = 10e-6;

  cfg.qssa.use_fmm = true;
  cfg.qssa.fmm_theta = 0.5;
  cfg.qssa.toxin_cutoff = 80e-6;
  cfg.qssa.nutrient_cutoff = 40e-6;

  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.acetate.enabled = true;
  cfg.chem_env.mucin.enabled = true;
  cfg.chem_env.protease.enabled = true;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 25;
  resident.mu_max = 5.5e-4;
  resident.plasmids = {"ColE1", "ColB"};
  resident.conjugative = true;
  cfg.initial_strains.push_back(resident);

  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = 12;
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

void assert_population_sane(const Simulation& sim) {
  Int live = 0;
  for (const Agent& a : sim.agents()) {
    if (a.state == PhenoState::DEAD) continue;
    assert(std::isfinite(a.biomass));
    assert(a.biomass > 0.0);
    ++live;
  }
  assert(live >= 1);
}

}  // namespace

void test_gpu_fmm_hybrid_production_domain() {
  SimulationConfig cfg = make_production_config(4001);

  Domain domain;
  domain.init(cfg.domain);
  assert(domain.nx() == 1000);
  assert(gpu_diffusion_line_lengths_supported(domain));

  Simulation sim;
  sim.init(cfg);
  assert(sim.gpu_active());
  sim.run();

  assert_chemistry_sane(sim);
  assert_population_sane(sim);

  std::cout << "  test_gpu_fmm_hybrid_production_domain: PASSED\n";
}

int main() {
  std::cout << "=== GPU Production Path Smoke Tests ===\n";

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  std::cout << "All GPU production-path tests passed.\n";
  return 0;
#else
  if (DeviceContext::device_count() <= 0) {
    std::cout << "  SKIPPED (no CUDA device)\n";
    std::cout << "All GPU production-path tests passed.\n";
    return 0;
  }

  test_gpu_fmm_hybrid_production_domain();

  std::cout << "All GPU production-path tests passed.\n";
  return 0;
#endif
}
