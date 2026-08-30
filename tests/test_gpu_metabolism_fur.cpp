#include "device.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"
#include "gpu_test_support.h"
#include "gpu_diagnostic_format.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;
using gutibm::gpu_diagnostic::format_real;

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
    const Real cpu_mu = cpu.agents()[i].mu_realized;
    const Real gpu_mu = gpu.agents()[i].mu_realized;
    const Real mu_absolute_difference = std::abs(cpu_mu - gpu_mu);
    const Real mu_scale = std::max({1.0, std::abs(cpu_mu), std::abs(gpu_mu)});
    const Real mu_relative_difference = mu_absolute_difference / mu_scale;
    std::cerr << "[gpu_diag][gpu_metabolism_fur][agent]"
              << " index=" << i
              << " field=mu_realized"
              << " cpu=" << format_real(cpu_mu)
              << " gpu=" << format_real(gpu_mu)
              << " abs_diff=" << format_real(mu_absolute_difference)
              << " rel_diff=" << format_real(mu_relative_difference) << "\n";
    assert(mu_absolute_difference <= tolerance);

    const Real cpu_biomass = cpu.agents()[i].biomass;
    const Real gpu_biomass = gpu.agents()[i].biomass;
    const Real biomass_absolute_difference =
        std::abs(cpu_biomass - gpu_biomass);
    const Real biomass_scale =
        std::max({1.0, std::abs(cpu_biomass), std::abs(gpu_biomass)});
    const Real biomass_relative_difference =
        biomass_absolute_difference / biomass_scale;
    std::cerr << "[gpu_diag][gpu_metabolism_fur][agent]"
              << " index=" << i
              << " field=biomass"
              << " cpu=" << format_real(cpu_biomass)
              << " gpu=" << format_real(gpu_biomass)
              << " abs_diff=" << format_real(biomass_absolute_difference)
              << " rel_diff=" << format_real(biomass_relative_difference)
              << "\n";
    assert(biomass_absolute_difference <= tolerance * biomass_scale);
  }

  const auto& cpu_chem = cpu.chemical_field();
  const auto& gpu_chem = gpu.chemical_field();
  assert(cpu_chem.num_species() == gpu_chem.num_species());
  for (Int species = 0; species < cpu_chem.num_species(); ++species) {
    Real maximum_absolute_difference = 0.0;
    Real maximum_relative_difference = 0.0;
    Int maximum_absolute_cell = -1;
    Int maximum_relative_cell = -1;
    for (Int cell = 0; cell < cpu_chem.global_ncells(); ++cell) {
      const Real expected = cpu_chem.conc_global(species, cell);
      const Real actual = gpu_chem.conc_global(species, cell);
      const Real scale = std::max({1.0, std::abs(expected), std::abs(actual)});
      const Real absolute_difference = std::abs(expected - actual);
      const Real relative_difference = absolute_difference / scale;
      if (absolute_difference > maximum_absolute_difference) {
        maximum_absolute_difference = absolute_difference;
        maximum_absolute_cell = cell;
      }
      if (relative_difference > maximum_relative_difference) {
        maximum_relative_difference = relative_difference;
        maximum_relative_cell = cell;
      }
      if (absolute_difference > tolerance * scale) {
        std::cerr << "[gpu_diag][gpu_metabolism_fur][chemical]"
                  << " species=" << species
                  << " name=" << cpu_chem.spec(species).name
                  << " cell=" << cell
                  << " cpu=" << format_real(expected)
                  << " gpu=" << format_real(actual)
                  << " abs_diff=" << format_real(absolute_difference)
                  << " rel_diff=" << format_real(relative_difference) << "\n";
      }
      assert(absolute_difference <= tolerance * scale);
    }
    std::cerr << "[gpu_diag][gpu_metabolism_fur][chemical_summary]"
              << " species=" << species
              << " name=" << cpu_chem.spec(species).name
              << " max_abs_diff=" << format_real(maximum_absolute_difference)
              << " max_abs_cell=" << maximum_absolute_cell
              << " max_rel_diff=" << format_real(maximum_relative_difference)
              << " max_rel_cell=" << maximum_relative_cell << "\n";
  }
}

void assert_expression_parity(Real cpu_expression, Real gpu_expression,
                              Real iron) {
  constexpr Real tolerance = 5.0e-8;
  const Real absolute_difference = std::abs(cpu_expression - gpu_expression);
  const Real scale =
      std::max({1.0, std::abs(cpu_expression), std::abs(gpu_expression)});
  const Real relative_difference = absolute_difference / scale;
  std::cerr << "[gpu_diag][gpu_metabolism_fur][sensitivity]"
            << " iron=" << format_real(iron)
            << " field=fepA_expression"
            << " cpu=" << format_real(cpu_expression)
            << " gpu=" << format_real(gpu_expression)
            << " abs_diff=" << format_real(absolute_difference)
            << " rel_diff=" << format_real(relative_difference) << "\n";
  assert(absolute_difference <= tolerance);
}

void assert_mu_parity(Real cpu_mu, Real gpu_mu, Real iron) {
  constexpr Real tolerance = 5.0e-8;
  const Real absolute_difference = std::abs(cpu_mu - gpu_mu);
  const Real scale = std::max({1.0, std::abs(cpu_mu), std::abs(gpu_mu)});
  const Real relative_difference = absolute_difference / scale;
  std::cerr << "[gpu_diag][gpu_metabolism_fur][sensitivity]"
            << " iron=" << format_real(iron)
            << " field=mu_realized"
            << " cpu=" << format_real(cpu_mu)
            << " gpu=" << format_real(gpu_mu)
            << " abs_diff=" << format_real(absolute_difference)
            << " rel_diff=" << format_real(relative_difference) << "\n";
  assert(absolute_difference <= tolerance);
}

void assert_ordering_match(const std::vector<Real>& cpu_values,
                           const std::vector<Real>& gpu_values,
                           const std::vector<Real>& iron_values,
                           const char* field, bool require_decrease) {
  assert(cpu_values.size() == gpu_values.size());
  assert(cpu_values.size() == iron_values.size());
  for (size_t i = 1; i < cpu_values.size(); ++i) {
    const Real cpu_difference = cpu_values[i - 1] - cpu_values[i];
    const Real gpu_difference = gpu_values[i - 1] - gpu_values[i];
    std::cerr << "[gpu_diag][gpu_metabolism_fur][ordering]"
              << " field=" << field
              << " iron_high=" << format_real(iron_values[i - 1])
              << " iron_low=" << format_real(iron_values[i])
              << " cpu_high=" << format_real(cpu_values[i - 1])
              << " cpu_low=" << format_real(cpu_values[i])
              << " gpu_high=" << format_real(gpu_values[i - 1])
              << " gpu_low=" << format_real(gpu_values[i])
              << " cpu_delta=" << format_real(cpu_difference)
              << " gpu_delta=" << format_real(gpu_difference) << "\n";
    if (require_decrease) {
      assert(cpu_values[i - 1] > cpu_values[i]);
      assert(gpu_values[i - 1] > gpu_values[i]);
      continue;
    }
    if (cpu_values[i - 1] > cpu_values[i]) {
      assert(gpu_values[i - 1] > gpu_values[i]);
    } else if (cpu_values[i - 1] < cpu_values[i]) {
      assert(gpu_values[i - 1] < gpu_values[i]);
    } else {
      assert(gpu_values[i - 1] == gpu_values[i]);
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
  const std::vector iron_values{3.0e-5, 1.0e-4, 1.0e-3};
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
    assert_expression_parity(cpu_expression.back(), gpu_expression.back(), iron);
    assert_mu_parity(cpu_mu.back(), gpu_mu.back(), iron);
  }

  assert_ordering_match(cpu_expression, gpu_expression, iron_values,
                        "fepA_expression", true);
  assert_ordering_match(cpu_mu, gpu_mu, iron_values, "mu_realized", false);

  const std::vector cap_iron_values{1.0e-8, 1.0e-10};
  std::vector<Real> cap_cpu_expression;
  std::vector<Real> cap_gpu_expression;
  for (const Real iron : cap_iron_values) {
    const SimulationConfig config = make_config(iron, 5513);
    Simulation cpu = run(config, false);
    Simulation gpu = run(config, true);
    assert(gpu.agents_gpu().metabolism_gpu_steps() > 0);
    cap_cpu_expression.push_back(receptor_expression(cpu));
    cap_gpu_expression.push_back(receptor_expression(gpu));
    assert_expression_parity(cap_cpu_expression.back(), cap_gpu_expression.back(),
                             iron);
    std::cerr << "[gpu_diag][gpu_metabolism_fur][cap]"
              << " iron=" << format_real(iron)
              << " cpu=" << format_real(cap_cpu_expression.back())
              << " gpu=" << format_real(cap_gpu_expression.back())
              << " receptor_max="
              << format_real(config.cell_bio.fur.receptor_max) << "\n";
    assert(std::abs(cap_cpu_expression.back()
                    - config.cell_bio.fur.receptor_max)
           <= 5.0e-8);
    assert(std::abs(cap_gpu_expression.back()
                    - config.cell_bio.fur.receptor_max)
           <= 5.0e-8);
  }
  assert(std::abs(cap_cpu_expression[0] - cap_cpu_expression[1]) <= 5.0e-8);
  assert(std::abs(cap_gpu_expression[0] - cap_gpu_expression[1]) <= 5.0e-8);
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
