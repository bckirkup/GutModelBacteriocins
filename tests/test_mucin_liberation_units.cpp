/* -----------------------------------------------------------------------
   GutIBM – Dynamic mucin liberation units test (Spec 1)
   ----------------------------------------------------------------------- */

#include "vbf.h"
#include "chemical_field.h"
#include "domain.h"
#include "input_parser.h"
#include "chem_environment_config.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

void test_dynamic_liberation_matches_static_scale() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.chem_env.mucin.enabled = true;
  InputParser::finalize_config(cfg);

  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 40e-6};
  dcfg.grid_dx = 5e-6;

  Domain domain;
  domain.init(dcfg);

  ChemicalSpec mucin_spec;
  mucin_spec.name = "mucin";
  mucin_spec.initial_conc = cfg.chem_env.mucin.initial_conc;
  mucin_spec.boundary_conc = cfg.chem_env.mucin.initial_conc;

  ChemicalSpec carbon_spec;
  carbon_spec.name = "carbon";
  carbon_spec.initial_conc = 0.0;

  ChemicalField chem;
  chem.init(domain, {carbon_spec, mucin_spec});

  VBF vbf;
  vbf.init(cfg.vbf, domain);
  OxygenConfig oxygen;
  AcetateConfig acetate;
  chem.zero_reactions();
  vbf.apply_nutrient_coupling(
      chem, domain, cfg.time.bio_dt, oxygen, acetate, cfg.chem_env.mucin);

  const Real dynamic_rate = chem.reac(0, domain.cell_index(0, 0, 1));
  const Real static_rate = cfg.vbf.mucin_liberation;
  assert(std::isfinite(dynamic_rate));
  assert(static_rate > 0.0);
  assert(std::abs(dynamic_rate) > static_rate / 10.0);
  assert(std::abs(dynamic_rate) < static_rate * 10.0);

  std::cout << "  test_dynamic_liberation_matches_static_scale: PASSED"
            << " (dynamic=" << dynamic_rate
            << " static=" << static_rate << ")\n";
}

int main() {
  std::cout << "=== Mucin Liberation Units Test ===\n";
  test_dynamic_liberation_matches_static_scale();
  std::cout << "All mucin liberation units tests passed.\n";
  return 0;
}
