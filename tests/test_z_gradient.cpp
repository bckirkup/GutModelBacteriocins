/* -----------------------------------------------------------------------
   GutIBM – Tests for z-dependent nutrient gradient (Issue #14)
   Verifies:
     1. Carbon concentration decays exponentially from z=0
     2. VBF mucin liberation rate decays with z
     3. Agents near epithelium (z=0) have higher carbon availability
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "domain.h"
#include "vbf.h"
#include "input_parser.h"
#include "simulation.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

void test_carbon_z_gradient_init() {
  // Set up a small domain and verify the exponential profile
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 100e-6};
  dcfg.grid_dx = 5e-6;
  dcfg.hash_cell_size = 10e-6;

  Domain domain;
  domain.init(dcfg);

  Real lambda = 25.0e-6;
  Real C0 = 5.0e-3;

  ChemicalSpec carbon;
  carbon.name = "carbon";
  carbon.diff_coeff = 1.0e-9;
  carbon.retardation = 1.0;
  carbon.initial_conc = C0;
  carbon.boundary_conc = C0;
  carbon.decay_rate = 0.0;
  carbon.z_gradient_enabled = true;
  carbon.z_gradient_lambda = lambda;

  ChemicalField chem;
  chem.init(domain, {carbon});

  Int nz = domain.nz();

  // Verify concentration at each z-layer matches C0 * exp(-z_rel / lambda)
  for (Int iz = 0; iz < nz; ++iz) {
    Real z_rel = (iz + 0.5) * domain.dx();
    Real expected = C0 * std::exp(-z_rel / lambda);
    Int cell = domain.cell_index(0, 0, iz);
    Real actual = chem.conc(0, cell);
    Real rel_err = std::abs(actual - expected) / expected;
    assert(rel_err < 1e-12);
  }

  // Verify concentration at z=0 layer > concentration at z=nz-1 layer
  Int bottom = domain.cell_index(0, 0, 0);
  Int top = domain.cell_index(0, 0, nz - 1);
  assert(chem.conc(0, bottom) > chem.conc(0, top));

  // Verify monotonic decay
  for (Int iz = 1; iz < nz; ++iz) {
    Int c_prev = domain.cell_index(0, 0, iz - 1);
    Int c_curr = domain.cell_index(0, 0, iz);
    assert(chem.conc(0, c_prev) > chem.conc(0, c_curr));
  }

  std::cout << "  test_carbon_z_gradient_init: PASSED\n";
}

void test_uniform_when_gradient_disabled() {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 100e-6};
  dcfg.grid_dx = 5e-6;
  dcfg.hash_cell_size = 10e-6;

  Domain domain;
  domain.init(dcfg);

  Real C0 = 5.0e-3;
  ChemicalSpec carbon;
  carbon.name = "carbon";
  carbon.diff_coeff = 1.0e-9;
  carbon.retardation = 1.0;
  carbon.initial_conc = C0;
  carbon.boundary_conc = C0;
  carbon.decay_rate = 0.0;
  carbon.z_gradient_enabled = false;
  carbon.z_gradient_lambda = 25.0e-6;

  ChemicalField chem;
  chem.init(domain, {carbon});

  // All cells should have uniform concentration
  for (Int c = 0; c < chem.ncells(); ++c) {
    assert(std::abs(chem.conc(0, c) - C0) < 1e-15);
  }

  std::cout << "  test_uniform_when_gradient_disabled: PASSED\n";
}

void test_vbf_mucin_rate_z_dependent() {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 100e-6};
  dcfg.grid_dx = 5e-6;
  dcfg.hash_cell_size = 10e-6;

  Domain domain;
  domain.init(dcfg);

  VBFConfig vcfg;
  vcfg.mucin_liberation = 5.0e-5;
  vcfg.mucin_z_gradient_enabled = true;
  vcfg.mucin_z_gradient_lambda = 25.0e-6;

  VBF vbf;
  vbf.init(vcfg, domain);

  // At z=0, rate should be maximum
  Real rate_0 = vbf.mucin_rate(0.0);
  assert(std::abs(rate_0 - vcfg.mucin_liberation) < 1e-15);

  // At z=lambda, rate should be exp(-1) * base
  Real rate_lambda = vbf.mucin_rate(25.0e-6);
  Real expected = vcfg.mucin_liberation * std::exp(-1.0);
  assert(std::abs(rate_lambda - expected) / expected < 1e-12);

  // Rate must decrease with z
  Real prev = vbf.mucin_rate(0.0);
  for (int i = 1; i <= 10; ++i) {
    Real z = i * 10.0e-6;
    Real curr = vbf.mucin_rate(z);
    assert(curr < prev);
    prev = curr;
  }

  std::cout << "  test_vbf_mucin_rate_z_dependent: PASSED\n";
}

void test_vbf_mucin_rate_uniform_when_disabled() {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 100e-6};
  dcfg.grid_dx = 5e-6;
  dcfg.hash_cell_size = 10e-6;

  Domain domain;
  domain.init(dcfg);

  VBFConfig vcfg;
  vcfg.mucin_liberation = 5.0e-5;
  vcfg.mucin_z_gradient_enabled = false;

  VBF vbf;
  vbf.init(vcfg, domain);

  // All z positions should return the same rate
  for (int i = 0; i <= 10; ++i) {
    Real z = i * 10.0e-6;
    assert(std::abs(vbf.mucin_rate(z) - vcfg.mucin_liberation) < 1e-15);
  }

  std::cout << "  test_vbf_mucin_rate_uniform_when_disabled: PASSED\n";
}

void test_vbf_coupling_z_gradient() {
  // Verify that after VBF coupling, carbon reaction rates near z=0 are
  // higher than those further from the epithelium
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 100e-6};
  dcfg.grid_dx = 5e-6;
  dcfg.hash_cell_size = 10e-6;

  Domain domain;
  domain.init(dcfg);

  ChemicalSpec carbon;
  carbon.name = "carbon";
  carbon.diff_coeff = 1.0e-9;
  carbon.retardation = 1.0;
  carbon.initial_conc = 5.0e-3;
  carbon.boundary_conc = 5.0e-3;
  carbon.decay_rate = 0.0;
  carbon.z_gradient_enabled = true;
  carbon.z_gradient_lambda = 25.0e-6;

  ChemicalField chem;
  chem.init(domain, {carbon});

  VBFConfig vcfg;
  vcfg.mucin_liberation = 5.0e-5;
  vcfg.mucin_z_gradient_enabled = true;
  vcfg.mucin_z_gradient_lambda = 25.0e-6;

  VBF vbf;
  vbf.init(vcfg, domain);

  chem.zero_reactions();
  OxygenConfig oxygen;
  AcetateConfig acetate;
  MucinConfig mucin;
  vbf.apply_nutrient_coupling(chem, domain, 60.0, oxygen, acetate, mucin);

  // Carbon reaction rate near z=0 should be > rate near z=max
  Int bottom_cell = domain.cell_index(0, 0, 0);
  Int top_cell = domain.cell_index(0, 0, domain.nz() - 1);
  assert(chem.reac(0, bottom_cell) > chem.reac(0, top_cell));

  std::cout << "  test_vbf_coupling_z_gradient: PASSED\n";
}

void test_default_config_has_z_gradient() {
  SimulationConfig cfg = InputParser::default_config();

  // Carbon species should have z_gradient_enabled
  bool found_carbon = false;
  for (const auto& spec : cfg.chemicals) {
    if (spec.name == "carbon") {
      found_carbon = true;
      assert(spec.z_gradient_enabled);
      assert(std::abs(spec.z_gradient_lambda - 25.0e-6) < 1e-15);
    }
  }
  assert(found_carbon);

  // VBF should have mucin z-gradient enabled
  assert(cfg.vbf.mucin_z_gradient_enabled);
  assert(std::abs(cfg.vbf.mucin_z_gradient_lambda - 25.0e-6) < 1e-15);

  std::cout << "  test_default_config_has_z_gradient: PASSED\n";
}

void test_smoke_with_z_gradient() {
  // Run a mini simulation with z-gradient enabled (default) to ensure
  // it doesn't crash and produces reasonable results
  SimulationConfig cfg = InputParser::default_config();

  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {50e-6, 50e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;

  cfg.total_time = 120.0;
  cfg.bio_dt = 60.0;
  cfg.output_interval = 120.0;
  cfg.seed = 314;
  cfg.hdf5.enabled = false;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length = 50e-6;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 10; s.mu_max = 5e-4;
  s.plasmids = {}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);
  sim.run();

  assert(sim.step_count() > 0);

  // Carbon near epithelium (z=0) should still be > carbon far from it
  const auto& chem = sim.chemical_field();
  Int i_carbon = chem.find("carbon");
  assert(i_carbon >= 0);

  const auto& dom = sim.domain();
  Int bottom = dom.cell_index(0, 0, 0);
  Int top = dom.cell_index(0, 0, dom.nz() - 1);
  // z=0 is the Dirichlet boundary so its concentration is boundary_conc
  assert(chem.conc(i_carbon, bottom) >= chem.conc(i_carbon, top));

  std::cout << "  test_smoke_with_z_gradient: PASSED\n";
}

int main() {
  std::cout << "=== Z-Gradient Tests ===\n";
  test_carbon_z_gradient_init();
  test_uniform_when_gradient_disabled();
  test_vbf_mucin_rate_z_dependent();
  test_vbf_mucin_rate_uniform_when_disabled();
  test_vbf_coupling_z_gradient();
  test_default_config_has_z_gradient();
  test_smoke_with_z_gradient();
  std::cout << "All z-gradient tests passed.\n";
  return 0;
}
