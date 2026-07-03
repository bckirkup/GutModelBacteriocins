/* -----------------------------------------------------------------------
   GutIBM – Oxygen gradient tests (Spec 1)
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "domain.h"
#include "input_parser.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

void test_oxygen_species_registered() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.oxygen.enabled = true;
  InputParser::finalize_config(cfg);

  bool found = false;
  for (const auto& spec : cfg.chemicals) {
    if (spec.name == "oxygen") {
      found = true;
      assert(spec.z_gradient_enabled);
      assert(std::abs(spec.boundary_conc - cfg.oxygen.epithelial_conc) < 1e-15);
    }
  }
  assert(found);
  std::cout << "  test_oxygen_species_registered: PASSED\n";
}

void test_oxygen_z_gradient_init() {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 100e-6};
  dcfg.grid_dx = 5e-6;

  Domain domain;
  domain.init(dcfg);

  const Real lambda = 25.0e-6;
  const Real c0 = 55.0e-6;

  ChemicalSpec oxygen;
  oxygen.name = "oxygen";
  oxygen.diff_coeff = 2.1e-9;
  oxygen.retardation = 1.0;
  oxygen.initial_conc = c0;
  oxygen.boundary_conc = c0;
  oxygen.z_gradient_enabled = true;
  oxygen.z_gradient_lambda = lambda;

  ChemicalField chem;
  chem.init(domain, {oxygen});

  for (Int iz = 0; iz < domain.nz(); ++iz) {
    Real z_rel = (iz + 0.5) * domain.dx();
    Real expected = c0 * std::exp(-z_rel / lambda);
    Int cell = domain.cell_index(0, 0, iz);
    Real actual = chem.conc(0, cell);
    assert(std::abs(actual - expected) / expected < 1e-12);
  }

  std::cout << "  test_oxygen_z_gradient_init: PASSED\n";
}

void test_oxygen_config_sensitivity() {
  SimulationConfig cfg_off = InputParser::default_config();
  cfg_off.oxygen.enabled = false;
  InputParser::finalize_config(cfg_off);

  SimulationConfig cfg_on = InputParser::default_config();
  cfg_on.oxygen.enabled = true;
  InputParser::finalize_config(cfg_on);

  Int count_off = 0;
  Int count_on = 0;
  for (const auto& s : cfg_off.chemicals) if (s.name == "oxygen") count_off++;
  for (const auto& s : cfg_on.chemicals) if (s.name == "oxygen") count_on++;
  assert(count_off == 0);
  assert(count_on == 1);

  std::cout << "  test_oxygen_config_sensitivity: PASSED\n";
}

int main() {
  std::cout << "=== Oxygen Gradient Tests ===\n";
  test_oxygen_species_registered();
  test_oxygen_z_gradient_init();
  test_oxygen_config_sensitivity();
  std::cout << "All oxygen gradient tests passed.\n";
  return 0;
}
