#include "device.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"
#include "gpu_test_support.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

namespace {

#ifdef GUTIBM_CUDA
SimulationConfig make_config(Real iron, uint64_t seed) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = seed;
  cfg.gpu.enabled = true;
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;
  cfg.time.total_time = 180.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 180.0;
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {40e-6, 40e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.cell_bio.fur.enabled = true;
  cfg.chem_env.siderophore.enabled = true;
  cfg.chem_env.acetate.enabled = true;
  cfg.chem_env.oxygen.enabled = true;
  cfg.fixes.metabolism.eut_enabled = true;
  cfg.fixes.metabolism.division_threshold = 100.0;
  cfg.quorum_sensing.enabled = true;
  cfg.enabled_fixes = {"metabolism", "quorum_sensing"};

  for (auto& species : cfg.chemicals) {
    if (species.name == gutibm::species::IRON) {
      species.initial_conc = iron;
      species.boundary_conc = iron;
    }
    if (species.name == gutibm::species::ACETATE) {
      species.initial_conc = 80.0;
      species.boundary_conc = 80.0;
    }
    if (species.name == gutibm::species::ETHANOLAMINE) {
      species.initial_conc = 1.0e-3;
      species.boundary_conc = 1.0e-3;
    }
  }

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 4;
  strain.mu_max = 5.0e-4;
  strain.plasmids = {};
  strain.conjugative = false;
  cfg.initial_strains.push_back(strain);
  return cfg;
}

Simulation run(const SimulationConfig& config, bool gpu) {
  SimulationConfig run_config = config;
  run_config.gpu.enabled = gpu;
  Simulation simulation;
  simulation.init(run_config);
  if (gpu) assert(simulation.gpu_active());
  simulation.run();
  return simulation;
}

Real receptor_expression(const Simulation& simulation) {
  const int receptor = to_underlying(ReceptorType::FepA);
  Real expression = 0.0;
  Int count = 0;
  for (const Agent& agent : simulation.agents()) {
    if (agent.state == PhenoState::DEAD || agent.flags.is_ghost) continue;
    expression += agent.receptor_expr[receptor];
    ++count;
  }
  assert(count > 0);
  return expression / static_cast<Real>(count);
}

Real mean_mu(const Simulation& simulation) {
  Real value = 0.0;
  Int count = 0;
  for (const Agent& agent : simulation.agents()) {
    if (agent.state == PhenoState::DEAD || agent.flags.is_ghost) continue;
    value += agent.mu_realized;
    ++count;
  }
  assert(count > 0);
  return value / static_cast<Real>(count);
}

Real maximum_ai2(const Simulation& simulation) {
  const Int ai2 = simulation.chemical_field().find(species::AI2);
  assert(ai2 >= 0);
  Real maximum = 0.0;
  for (Int cell = 0; cell < simulation.chemical_field().global_ncells();
       ++cell) {
    maximum = std::max(maximum,
                       simulation.chemical_field().conc_global(ai2, cell));
  }
  return maximum;
}

void assert_parity(const Simulation& cpu, const Simulation& gpu) {
  constexpr Real tolerance = 5.0e-8;
  assert(cpu.agents().size() == gpu.agents().size());
  for (Int i = 0; i < cpu.agents().size(); ++i) {
    assert(std::abs(cpu.agents()[i].mu_realized
                    - gpu.agents()[i].mu_realized) <= tolerance);
    assert(std::abs(cpu.agents()[i].biomass
                    - gpu.agents()[i].biomass) <= tolerance
           * std::max({1.0, std::abs(cpu.agents()[i].biomass),
                       std::abs(gpu.agents()[i].biomass)}));
  }

  const auto& cpu_chem = cpu.chemical_field();
  const auto& gpu_chem = gpu.chemical_field();
  assert(cpu_chem.num_species() == gpu_chem.num_species());
  for (Int species = 0; species < cpu_chem.num_species(); ++species) {
    for (Int cell = 0; cell < cpu_chem.global_ncells(); ++cell) {
      const Real expected = cpu_chem.conc_global(species, cell);
      const Real actual = gpu_chem.conc_global(species, cell);
      const Real scale = std::max({1.0, std::abs(expected), std::abs(actual)});
      assert(std::abs(expected - actual) <= tolerance * scale);
    }
  }
}

void test_fur_siderophore_metabolism_parity() {
  const SimulationConfig config = make_config(1.0e-8, 5511);
  Simulation cpu = run(config, false);
  Simulation gpu = run(config, true);
  assert(gpu.agents_gpu().metabolism_gpu_steps() > 0);
  assert(maximum_ai2(cpu) > 0.0);
  assert(maximum_ai2(gpu) > 0.0);
  assert_parity(cpu, gpu);
  std::cout << "  test_fur_siderophore_metabolism_parity: PASSED\n";
}

void test_fur_sensitivity() {
  const std::vector<Real> iron_values{1.0e-10, 1.0e-7, 1.0e-3};
  std::vector<Real> cpu_expression;
  std::vector<Real> gpu_expression;
  std::vector<Real> cpu_mu;
  std::vector<Real> gpu_mu;
  for (const Real iron : iron_values) {
    const SimulationConfig config = make_config(iron, 5512);
    Simulation cpu = run(config, false);
    Simulation gpu = run(config, true);
    assert(gpu.agents_gpu().metabolism_gpu_steps() > 0);
    cpu_expression.push_back(receptor_expression(cpu));
    gpu_expression.push_back(receptor_expression(gpu));
    cpu_mu.push_back(mean_mu(cpu));
    gpu_mu.push_back(mean_mu(gpu));
  }

  assert(cpu_expression[0] > cpu_expression[1]);
  assert(cpu_expression[1] > cpu_expression[2]);
  assert(gpu_expression[0] > gpu_expression[1]);
  assert(gpu_expression[1] > gpu_expression[2]);
  assert(cpu_mu[0] < cpu_mu[1]);
  assert(cpu_mu[1] < cpu_mu[2]);
  assert(gpu_mu[0] < gpu_mu[1]);
  assert(gpu_mu[1] < gpu_mu[2]);
  assert(cpu_expression[0] >= 5.0 - 1.0e-12);
  assert(gpu_expression[0] >= 5.0 - 1.0e-12);
  assert(std::abs((gpu_expression[0] - gpu_expression[2])
                  - (cpu_expression[0] - cpu_expression[2]))
         <= 0.2 * (cpu_expression[0] - cpu_expression[2]));
  assert(std::abs((gpu_mu[0] - gpu_mu[2]) - (cpu_mu[0] - cpu_mu[2]))
         <= 0.2 * std::abs(cpu_mu[0] - cpu_mu[2]));
  assert(!cpu_mu.empty());
  std::cout << "  test_fur_sensitivity: PASSED\n";
}
#endif

}  // namespace

int main() {
  std::cout << "=== GPU Fur/Siderophore Metabolism Parity ===\n";
  const int gpu_status = gutibm::test::require_gpu("gpu_metabolism_fur");
  if (gpu_status != 0) return gpu_status;
#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  test_fur_siderophore_metabolism_parity();
  test_fur_sensitivity();
  std::cout << "GPU Fur/siderophore metabolism tests passed.\n";
  return 0;
#endif
}
