#include "chemical_field.h"
#include "config_json.h"
#include "domain.h"
#include "input_parser.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace gutibm;

namespace {

int failures = 0;

void expect(bool condition, const std::string& message) {
  if (condition) return;
  std::cerr << "FAIL: " << message << "\n";
  ++failures;
}

bool near(Real actual, Real expected) {
  return std::abs(actual - expected) <= 1.0e-15;
}

void check_strain_defaults(const SimulationConfig::InitialStrain& strain) {
  expect(strain.type == 0, "default strain type must be zero");
  expect(strain.count == 0, "default strain count must be zero");
  expect(near(strain.mu_max, 5.0e-4), "default strain growth rate changed");
  expect(strain.plasmids.empty(), "default strain must have no plasmids");
  expect(!strain.conjugative, "default strain must not be conjugative");
  expect(strain.cdi_type == 0, "default strain CDI type must be zero");
  expect(strain.cdi_immunity == 0, "default strain CDI immunity must be zero");
}

void test_initial_strain_defaults_and_overrides() {
  check_strain_defaults(SimulationConfig::InitialStrain{});

  const auto defaults =
      ConfigJson::parse_initial_strains(R"({"initial_strains":[{}]})");
  expect(defaults.found, "empty strain object must parse");
  expect(defaults.strains.size() == 1, "empty strain object must produce one strain");
  if (defaults.strains.size() == 1) {
    check_strain_defaults(defaults.strains.front());
  }

  const auto overridden = ConfigJson::parse_initial_strains(
      R"json({"initial_strains":[{"type":7,"count":3,"mu_max":0.002,
      "plasmids":["ColE1"],"conjugative":true,"cdi_type":4,"cdi_immunity":5}]})json");
  expect(overridden.found, "overridden strain object must parse");
  expect(overridden.strains.size() == 1,
         "overridden strain object must produce one strain");
  if (overridden.strains.size() == 1) {
    const auto& strain = overridden.strains.front();
    expect(strain.type == 7, "strain type override was ignored");
    expect(strain.count == 3, "strain count override was ignored");
    expect(near(strain.mu_max, 0.002), "strain growth-rate override was ignored");
    expect(strain.plasmids.size() == 1 && strain.plasmids.front() == "ColE1",
           "strain plasmid override was ignored");
    expect(strain.conjugative, "strain conjugative override was ignored");
    expect(strain.cdi_type == 4, "strain CDI type override was ignored");
    expect(strain.cdi_immunity == 5, "strain CDI immunity override was ignored");
  }
}

Domain make_domain() {
  DomainConfig config;
  config.lo = {0.0, 0.0, 0.0};
  config.hi = {10.0e-6, 10.0e-6, 10.0e-6};
  config.grid_dx = 5.0e-6;
  config.hash_cell_size = 10.0e-6;
  Domain domain;
  domain.init(config);
  return domain;
}

void test_chemical_defaults_and_initial_concentration_sensitivity() {
  ChemicalSpec spec;
  expect(spec.name.empty(), "default chemical name must be empty");
  expect(near(spec.diff_coeff, 0.0), "default diffusion coefficient must be zero");
  expect(near(spec.retardation, 1.0), "default retardation must be one");
  expect(near(spec.initial_conc, 0.0), "default initial concentration must be zero");
  expect(near(spec.boundary_conc, 0.0), "default boundary concentration must be zero");
  expect(near(spec.decay_rate, 0.0), "default decay rate must be zero");
  expect(!spec.z_gradient_enabled, "default z gradient must be disabled");
  expect(near(spec.z_gradient_lambda, 25.0e-6),
         "default z-gradient length changed");
  expect(!spec.diffusion_enabled, "default diffusion must be disabled");

  const Domain domain = make_domain();
  ChemicalField default_field;
  default_field.init(domain, {spec});
  for (Int cell = 0; cell < default_field.ncells(); ++cell) {
    expect(near(default_field.conc(0, cell), 0.0),
           "default chemical field must initialize to zero");
  }

  spec.initial_conc = 0.25;
  ChemicalField configured_field;
  configured_field.init(domain, {spec});
  for (Int cell = 0; cell < configured_field.ncells(); ++cell) {
    expect(near(configured_field.conc(0, cell), 0.25),
           "initial concentration must change the initialized field");
  }
}

}  // namespace

int main() {
  test_initial_strain_defaults_and_overrides();
  test_chemical_defaults_and_initial_concentration_sensitivity();
  if (failures != 0) {
    std::cerr << failures << " data default check(s) failed\n";
    return 1;
  }
  std::cout << "All data default tests passed.\n";
  return 0;
}
