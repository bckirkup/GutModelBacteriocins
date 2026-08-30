#include "chemical_field.h"
#include "error.h"
#include "input_parser.h"
#include "simulation.h"
#include "sim_fingerprint.h"
#include "species_names.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <set>
#include <string>

using namespace gutibm;

namespace {

std::set<std::string, std::less<>> names(const SimulationConfig& cfg) {
  std::set<std::string, std::less<>> result;
  for (const auto& spec : cfg.chemicals) result.insert(spec.name);
  return result;
}

void test_exact_species_sets() {
  const auto full = InputParser::default_config();
  const std::set<std::string, std::less<>> expected_full = {
      species::CARBON, species::IRON, species::B12,
      species::BACTERIOCIN_BTUB, species::BACTERIOCIN_FEPA,
      species::BACTERIOCIN_CIRA, species::BACTERIOCIN_FHUA,
      species::ACETATE, species::ETHANOLAMINE, species::SIDEROPHORE,
      species::FERRIC_ENTEROBACTIN};
  assert(names(full) == expected_full);

  auto nutrient = full;
  nutrient.species_subset = "nutrient_only";
  InputParser::finalize_config(nutrient);
  const std::set<std::string, std::less<>> expected_nutrient = {
      species::CARBON, species::IRON, species::B12, species::ACETATE,
      species::ETHANOLAMINE, species::SIDEROPHORE,
      species::FERRIC_ENTEROBACTIN};
  assert(names(nutrient) == expected_nutrient);

  auto carbon = full;
  carbon.species_subset = "carbon_only";
  InputParser::finalize_config(carbon);
  assert((names(carbon)
          == std::set<std::string, std::less<>>{species::CARBON}));
  std::cout << "  test_exact_species_sets: PASSED\n";
}

void test_required_species_rejection() {
  auto cfg = InputParser::default_config();
  cfg.chemicals.erase(
      std::remove_if(cfg.chemicals.begin(), cfg.chemicals.end(),
                     [](const ChemicalSpec& spec) {
                       return spec.name == species::CARBON;
                     }),
      cfg.chemicals.end());
  cfg.hdf5.enabled = false;
  bool threw = false;
  try {
    Simulation sim;
    sim.init(cfg);
  } catch (const ConfigError& ex) {
    const std::string message = ex.what();
    threw = message.find("fix_metabolism") != std::string::npos
        && message.find("carbon") != std::string::npos
        && message.find("chemistry.species_subset") != std::string::npos;
  }
  assert(threw);
  std::cout << "  test_required_species_rejection: PASSED\n";
}

void test_gpu_subset_rejection() {
  auto cfg = InputParser::default_config();
  cfg.species_subset = "nutrient_only";
  cfg.gpu.enabled = true;
  cfg.hdf5.enabled = false;
  bool rejected = false;
  try {
    Simulation sim;
    sim.init(cfg);
  } catch (const ConfigError& ex) {
    const std::string message = ex.what();
    rejected = message.find("gpu_enabled") != std::string::npos
        && message.find("species_subset") != std::string::npos;
  }
  assert(rejected);
  std::cout << "  test_gpu_subset_rejection: PASSED\n";
}

void test_subset_composition() {
  for (const std::string mode : {"slab", "agents"}) {
    auto cfg = InputParser::default_config();
    cfg.species_subset = "nutrient_only";
    cfg.hdf5.enabled = false;
    cfg.domain.hi = {50e-6, 50e-6, 25e-6};
    cfg.domain.grid_dx = 5e-6;
    if (mode == "slab") {
      cfg.chemistry_decomposition = "slab";
      // Slab cannot represent regularized support; opt into grid-dependent
      // single-voxel delivery explicitly.
      cfg.fixes.metabolism.delivery_far_field_radius = 0.0;
      cfg.domain.grid_halo_width = 2;
    } else {
      cfg.qssa.toxin_evaluation = "agents";
    }
    Simulation sim;
    sim.init(cfg);
    sim.step(1.0);
    assert(sim.agents().size() > 0);
    for (Int species_index = 0;
         species_index < sim.chemical_field().num_species(); ++species_index) {
      for (Int cell = 0; cell < sim.chemical_field().ncells(); ++cell) {
        const Real value = sim.chemical_field().conc(species_index, cell);
        assert(std::isfinite(value));
        assert(value >= 0.0);
      }
    }
  }
  std::cout << "  test_subset_composition: PASSED\n";
}

void test_subset_liveness() {
  std::array<uint64_t, 3> fingerprints{};
  size_t fingerprint_index = 0;
  for (const std::string subset : {"full", "nutrient_only", "carbon_only"}) {
    auto cfg = InputParser::default_config();
    cfg.species_subset = subset;
    cfg.hdf5.enabled = false;
    cfg.domain.hi = {50e-6, 50e-6, 25e-6};
    cfg.domain.grid_dx = 5e-6;
    Simulation sim;
    sim.init(cfg);
    sim.step(1.0);
    assert(sim.agents().size() > 0);
    fingerprints[fingerprint_index++] =
        test_util::simulation_fingerprint(sim);
    for (Int species_index = 0;
         species_index < sim.chemical_field().num_species(); ++species_index) {
      for (Int cell = 0; cell < sim.chemical_field().ncells(); ++cell) {
        const Real value = sim.chemical_field().conc(species_index, cell);
        assert(std::isfinite(value));
        assert(value >= 0.0);
      }
    }
  }
  assert(fingerprints[0] != fingerprints[1]);
  assert(fingerprints[0] != fingerprints[2]);
  std::cout << "  test_subset_liveness: PASSED\n";
}

}  // namespace

int main() {
  std::cout << "=== Species Subset Tests ===\n";
  test_exact_species_sets();
  test_required_species_rejection();
  test_gpu_subset_rejection();
  test_subset_composition();
  test_subset_liveness();
  return 0;
}
