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
  cfg.chem_env.oxygen.enabled = true;
  InputParser::finalize_config(cfg);

  bool found = false;
  for (const auto& spec : cfg.chemicals) {
    if (spec.name == "oxygen") {
      found = true;
      assert(spec.z_gradient_enabled);
      assert(std::abs(spec.boundary_conc - cfg.chem_env.oxygen.epithelial_conc) < 1e-15);
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
    Real z_rel = (iz + 0.5) * domain.dx_z();
    Real expected = c0 * std::exp(-z_rel / lambda);
    Int cell = domain.cell_index(0, 0, iz);
    Real actual = chem.conc(0, cell);
    assert(std::abs(actual - expected) / expected < 1e-12);
  }

  std::cout << "  test_oxygen_z_gradient_init: PASSED\n";
}

void test_oxygen_config_sensitivity() {
  SimulationConfig cfg_off = InputParser::default_config();
  cfg_off.chem_env.oxygen.enabled = false;
  InputParser::finalize_config(cfg_off);

  SimulationConfig cfg_on = InputParser::default_config();
  cfg_on.chem_env.oxygen.enabled = true;
  InputParser::finalize_config(cfg_on);

  Int count_off = 0;
  Int count_on = 0;
  for (const auto& s : cfg_off.chemicals) if (s.name == "oxygen") count_off++;
  for (const auto& s : cfg_on.chemicals) if (s.name == "oxygen") count_on++;
  assert(count_off == 0);
  assert(count_on == 1);

  std::cout << "  test_oxygen_config_sensitivity: PASSED\n";
}

void test_oxygen_delivery_non_gradient_clip_accounting() {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 20e-6};
  dcfg.grid_dx = 5e-6;
  Domain domain;
  domain.init(dcfg);

  ChemicalSpec oxygen;
  oxygen.name = "oxygen";
  oxygen.diff_coeff = 2.1e-9;
  oxygen.initial_conc = 1.0e-3;
  oxygen.boundary_conc = 1.0e-3;
  oxygen.diffusion_enabled = true;
  oxygen.delivery_enabled = true;

  ChemicalField chem;
  chem.init(domain, {oxygen});
  const Int cell = domain.cell_index(1, 1, 1);
  chem.add_sink_rate_global(0, cell, 0.2);
  chem.apply_diffusion(domain, 60.0);

  assert(chem.sink_realized_global(0, cell) > 0.0);
  const auto& flux = chem.flux_accounting();
  assert(flux.reaction_clip_interval[0] == 0.0);
  std::cout << "  test_oxygen_delivery_non_gradient_clip_accounting: PASSED\n";
}

void test_oxygen_delivery_gradient_accounting() {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 100e-6};
  dcfg.grid_dx = 5e-6;
  Domain domain;
  domain.init(dcfg);

  ChemicalSpec oxygen;
  oxygen.name = "oxygen";
  oxygen.diff_coeff = 2.1e-9;
  oxygen.initial_conc = 55.0e-6;
  oxygen.boundary_conc = 55.0e-6;
  oxygen.z_gradient_enabled = true;
  oxygen.z_gradient_lambda = 25.0e-6;
  oxygen.diffusion_enabled = true;
  oxygen.delivery_enabled = true;

  ChemicalField chem;
  chem.init(domain, {oxygen});
  const Int cell = domain.cell_index(1, 1, 1);
  chem.add_sink_rate_global(0, cell, 0.2);
  chem.apply_diffusion(domain, 60.0);

  assert(chem.sink_realized_global(0, cell) > 0.0);
  assert(chem.flux_accounting().reaction_clip_interval[0] == 0.0);
  std::cout << "  test_oxygen_delivery_gradient_accounting: PASSED\n";
}

int main() {
  std::cout << "=== Oxygen Gradient Tests ===\n";
  test_oxygen_species_registered();
  test_oxygen_z_gradient_init();
  test_oxygen_config_sensitivity();
  test_oxygen_delivery_non_gradient_clip_accounting();
  test_oxygen_delivery_gradient_accounting();
  std::cout << "All oxygen gradient tests passed.\n";
  return 0;
}
