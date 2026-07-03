/* -----------------------------------------------------------------------
   GutIBM – Dynamic mucin liberation tests (Spec 1)
   ----------------------------------------------------------------------- */

#include "vbf.h"
#include "chemical_field.h"
#include "domain.h"
#include "input_parser.h"
#include "chem_environment_config.h"
#include <cassert>
#include <iostream>

using namespace gutibm;

void test_mucin_liberation_scales_with_density() {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 40e-6};
  dcfg.grid_dx = 5e-6;

  Domain domain;
  domain.init(dcfg);

  ChemicalSpec mucin_spec;
  mucin_spec.name = "mucin";
  mucin_spec.initial_conc = 1.0e-2;
  mucin_spec.boundary_conc = 1.0e-2;

  ChemicalSpec carbon_spec;
  carbon_spec.name = "carbon";
  carbon_spec.initial_conc = 0.0;

  ChemicalField chem;
  chem.init(domain, {carbon_spec, mucin_spec});

  MucinConfig mucin_cfg;
  mucin_cfg.enabled = true;
  OxygenConfig oxygen;
  AcetateConfig acetate;

  VBFConfig vcfg_low;
  vcfg_low.density = 1.0e10;
  vcfg_low.use_dynamic_mucin = true;

  VBFConfig vcfg_high = vcfg_low;
  vcfg_high.density = 1.0e11;

  VBF vbf_low;
  vbf_low.init(vcfg_low, domain);
  chem.zero_reactions();
  vbf_low.apply_nutrient_coupling(chem, domain, 60.0, oxygen, acetate, mucin_cfg);
  const Real reac_low = chem.reac(0, domain.cell_index(0, 0, 1));

  ChemicalField chem2;
  chem2.init(domain, {carbon_spec, mucin_spec});
  VBF vbf_high;
  vbf_high.init(vcfg_high, domain);
  chem2.zero_reactions();
  vbf_high.apply_nutrient_coupling(chem2, domain, 60.0, oxygen, acetate, mucin_cfg);
  const Real reac_high = chem2.reac(0, domain.cell_index(0, 0, 1));

  assert(reac_high > reac_low);
  assert(reac_low > 0.0);

  std::cout << "  test_mucin_liberation_scales_with_density: PASSED"
            << " (low=" << reac_low << " high=" << reac_high << ")\n";
}

void test_mucin_species_registered() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.chem_env.mucin.enabled = true;
  InputParser::finalize_config(cfg);

  bool found_mucin = false;
  for (const auto& spec : cfg.chemicals) {
    if (spec.name == "mucin") found_mucin = true;
  }
  assert(found_mucin);
  assert(cfg.vbf.use_dynamic_mucin);

  std::cout << "  test_mucin_species_registered: PASSED\n";
}

int main() {
  std::cout << "=== Mucin Liberation Tests ===\n";
  test_mucin_liberation_scales_with_density();
  test_mucin_species_registered();
  std::cout << "All mucin liberation tests passed.\n";
  return 0;
}
