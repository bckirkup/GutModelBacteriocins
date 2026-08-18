/* -----------------------------------------------------------------------
   GutIBM – GPU feature-combination smoke (Spec 9 PR4)
   Runs selected Spec 8 combo scenarios with gpu_enabled and asserts
   finite chemistry + population outcomes.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "dispatch.h"
#include "device.h"
#include "gpu_diagnostic_format.h"
#include "species_names.h"
#include "gpu_test_support.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

using namespace gutibm;
using gutibm::gpu_diagnostic::format_real;

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
  // These combinations measure GPU metabolism; siderophore chemistry is
  // CPU-authoritative and has a dedicated fallback test below.
  cfg.chem_env.siderophore.enabled = false;

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

#ifdef GUTIBM_CUDA
void assert_flux_parity(const ChemicalField& reference,
                        const ChemicalField& candidate) {
  const auto& expected = reference.flux_accounting();
  const auto& actual = candidate.flux_accounting();
  const std::array<const std::vector<Real> NutrientFluxAccounting::*,
                   5> fields = {
      &NutrientFluxAccounting::boundary_interval,
      &NutrientFluxAccounting::vbf_source_interval,
      &NutrientFluxAccounting::vbf_sink_interval,
      &NutrientFluxAccounting::agent_uptake_interval,
      &NutrientFluxAccounting::reaction_clip_interval,
  };
  const std::array<const std::vector<Real> NutrientFluxAccounting::*,
                   5> cumulative_fields = {
      &NutrientFluxAccounting::boundary_cumulative,
      &NutrientFluxAccounting::vbf_source_cumulative,
      &NutrientFluxAccounting::vbf_sink_cumulative,
      &NutrientFluxAccounting::agent_uptake_cumulative,
      &NutrientFluxAccounting::reaction_clip_cumulative,
  };
  for (Int s = 0; s < reference.num_species(); ++s) {
    for (size_t i = 0; i < fields.size(); ++i) {
      const Real expected_value =
          (expected.*fields[i])[static_cast<size_t>(s)]
          + (expected.*cumulative_fields[i])[static_cast<size_t>(s)];
      const Real actual_value =
          (actual.*fields[i])[static_cast<size_t>(s)]
          + (actual.*cumulative_fields[i])[static_cast<size_t>(s)];
      const Real scale = std::max({1.0, std::abs(expected_value),
                                   std::abs(actual_value)});
      assert(std::abs(actual_value - expected_value) <= 1.0e-7 * scale);
    }
  }
}

void assert_chemistry_parity(const ChemicalField& reference,
                             const ChemicalField& candidate,
                             Real tolerance) {
  assert(reference.num_species() == candidate.num_species());
  for (Int s = 0; s < reference.num_species(); ++s) {
    for (Int cell = 0; cell < reference.global_ncells(); ++cell) {
      const Real expected = reference.conc(s, cell);
      const Real actual = candidate.conc_global(s, cell);
      const Real scale = std::max({1.0, std::abs(expected), std::abs(actual)});
      assert(std::abs(actual - expected) <= tolerance * scale);
    }
  }
}

void test_gpu_slab_equivalence_and_accounting() {
  SimulationConfig base = make_combo_config(3010);
  base.time.total_time = 120.0;
  base.chem_env.oxygen.enabled = true;
  base.chem_env.acetate.enabled = true;
  base.chem_env.mucin.enabled = true;
  base.domain.chemistry_stride = {2, 2, 1};

  SimulationConfig cpu_cfg = base;
  cpu_cfg.gpu.enabled = false;
  Simulation cpu;
  cpu.init(cpu_cfg);
  cpu.run();

  SimulationConfig replicated_cfg = base;
  Simulation replicated = run_gpu_combo(replicated_cfg);

  SimulationConfig slab_cfg = base;
  slab_cfg.chemistry_decomposition = "slab";
  slab_cfg.domain.chemistry_stride = {2, 2, 1};
  slab_cfg.domain.grid_halo_width = static_cast<Int>(
      std::ceil(slab_cfg.domain.ghost_width / (
          slab_cfg.domain.grid_dx * slab_cfg.domain.chemistry_stride[0])));
  Simulation slab = run_gpu_combo(slab_cfg);

  // Device and CPU directional diffusion use the same owned-cell arithmetic;
  // atomic uptake/VBF reductions can differ only in summation order.
  assert_chemistry_parity(replicated.chemical_field(),
                          slab.chemical_field(), 1.0e-9);
  assert_chemistry_parity(cpu.chemical_field(), replicated.chemical_field(),
                          1.0e-5);
  assert_flux_parity(replicated.chemical_field(), slab.chemical_field());

  const Int carbon = slab.chemical_field().find(species::CARBON);
  assert(carbon >= 0);
  std::array<Real, 3> sink_values = {0.0, 1.0e-7, 1.0e-5};
  std::array<Real, 3> carbon_means{};
  for (size_t i = 0; i < sink_values.size(); ++i) {
    SimulationConfig sensitivity_cfg = slab_cfg;
    sensitivity_cfg.seed = 3020;
    sensitivity_cfg.vbf.carbon_sink_vmax = sink_values[i];
    Simulation sensitivity = run_gpu_combo(sensitivity_cfg);
    Real sum = 0.0;
    for (Int cell = 0; cell < sensitivity.chemical_field().global_ncells();
         ++cell) {
      sum += sensitivity.chemical_field().conc_global(carbon, cell);
    }
    carbon_means[i] = sum
        / static_cast<Real>(sensitivity.chemical_field().global_ncells());
  }
  assert(carbon_means[0] > carbon_means[1]);
  assert(carbon_means[1] > carbon_means[2]);

  std::cout << "  test_gpu_slab_equivalence_and_accounting: PASSED\n";
}
#endif

}  // namespace

void test_gpu_siderophore_cpu_fallback() {
  SimulationConfig cfg = make_combo_config(3004);
  cfg.chem_env.siderophore.enabled = true;

  SimulationConfig cpu_cfg = cfg;
  cpu_cfg.gpu.enabled = false;
  Simulation cpu;
  cpu.init(cpu_cfg);
  cpu.run();

  Simulation gpu_cfg_sim;
  gpu_cfg_sim.init(cfg);
  assert(gpu_cfg_sim.gpu_active());
  gpu_cfg_sim.run();

  Int cpu_live = 0;
  Int gpu_live = 0;
  Real cpu_biomass = 0.0;
  Real gpu_biomass = 0.0;
  for (const Agent& agent : cpu.agents()) {
    if (agent.state == PhenoState::DEAD) continue;
    ++cpu_live;
    cpu_biomass += agent.biomass;
  }
  for (const Agent& agent : gpu_cfg_sim.agents()) {
    if (agent.state == PhenoState::DEAD) continue;
    ++gpu_live;
    gpu_biomass += agent.biomass;
  }
  assert(cpu_live > 0);
  assert(cpu_biomass > 0.0);
  assert(gpu_live == cpu_live);
  const Real biomass_rel_diff =
      std::abs(gpu_biomass - cpu_biomass) / cpu_biomass;
  // The deterministic CPU-vs-GPU numerical offset measured 4e-7--7e-7
  // across configurations. This 1e-5 bound leaves an order of headroom while
  // remaining two orders tighter than gpu_smoke's 5% bound. The offset comes
  // from the different diffusion solves, not nondeterminism; the
  // gpu_reproducibility target checks exact repeated GPU biomass for one
  // fixture/device, so that evidence is not a general reproducibility
  // guarantee for larger populations or more atomic contention.
  constexpr Real kBiomassRelativeTolerance = 1.0e-5;
  std::cout << "  test_gpu_siderophore_cpu_fallback: biomass_rel_diff="
            << biomass_rel_diff << "\n";
  if (!(biomass_rel_diff <= kBiomassRelativeTolerance)) {
    const Real biomass_absolute_difference =
        std::abs(gpu_biomass - cpu_biomass);
    std::cerr << "[gpu_diag][gpu_feature_combinations]"
              << " cpu_live=" << cpu_live << " gpu_live=" << gpu_live
              << " cpu_biomass=" << format_real(cpu_biomass)
              << " gpu_biomass=" << format_real(gpu_biomass)
              << " abs_diff=" << format_real(biomass_absolute_difference)
              << " rel_diff=" << format_real(biomass_rel_diff)
              << " tolerance=" << format_real(kBiomassRelativeTolerance)
              << "\n";
    gutibm::gpu_diagnostic::print_concentration_diagnostics(
        "gpu_feature_combinations", nullptr, cpu.chemical_field(),
        gpu_cfg_sim.chemical_field());
  }
  assert(biomass_rel_diff <= kBiomassRelativeTolerance);

  std::cout << "  test_gpu_siderophore_cpu_fallback: PASSED\n";
}

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

void test_gpu_slab_single_rank() {
  SimulationConfig cfg = make_combo_config(3005);
  cfg.chemistry_decomposition = "slab";
  cfg.domain.grid_halo_width = static_cast<Int>(
      std::ceil(cfg.domain.ghost_width / cfg.domain.grid_dx));
  Simulation sim = run_gpu_combo(cfg);
  assert_chemistry_sane(sim);
  assert_population_sane(sim);
  std::cout << "  test_gpu_slab_single_rank: PASSED\n";
}

int main() {
  std::cout << "=== GPU Feature Combination Smoke Tests ===\n";
  const int gpu_status = test::require_gpu("gpu_feature_combinations");
  if (gpu_status != 0) return gpu_status;

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  std::cout << "All GPU feature-combination tests passed.\n";
  return 0;
#else
  test_gpu_full_chemical_environment();
  test_gpu_adaptive_dt_with_crypts();
  test_gpu_kitchen_sink_light();
  test_gpu_siderophore_cpu_fallback();
  test_gpu_slab_single_rank();
  test_gpu_slab_equivalence_and_accounting();

  std::cout << "All GPU feature-combination tests passed.\n";
  return 0;
#endif
}
