/* -----------------------------------------------------------------------
   GutIBM – QSSA nutrient stoichiometry config tests (Spec 0 §6)
   Verifies the formerly-hardcoded nutrient consumption stoichiometries are
   now driven by QSSAConfig and reach the reaction field (golden value +
   configuration sensitivity).
   ----------------------------------------------------------------------- */

#include "qssa_solver.h"
#include "chemical_field.h"
#include "domain.h"
#include "advection.h"
#include "agent.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

namespace {

Domain make_domain() {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 20e-6};
  dcfg.grid_dx = 5e-6;
  dcfg.hash_cell_size = 10e-6;
  Domain d;
  d.init(dcfg);
  return d;
}

ChemicalField make_chem(const Domain& domain) {
  ChemicalSpec carbon; carbon.name = "carbon"; carbon.diff_coeff = 1e-9;
  carbon.retardation = 1.0; carbon.initial_conc = 5e-3; carbon.boundary_conc = 5e-3;
  ChemicalSpec iron = carbon; iron.name = "iron";
  ChemicalSpec b12 = carbon; b12.name = "b12";
  ChemicalField chem;
  chem.init(domain, {carbon, iron, b12});
  return chem;
}

// Run one nutrient-depletion pass for a single agent and return the carbon
// reaction rate deposited at the agent's cell.
Real carbon_reac_for_stoich(const Domain& domain, Real carbon_stoich,
                            Real& out_consumption, Int& out_cell) {
  AdvectionField adv;
  AdvectionConfig acfg;
  acfg.radial_turnover = 1e20;
  acfg.distal_transit_time = 1e20;
  adv.init(acfg, domain);

  QSSAConfig cfg;
  cfg.carbon_stoichiometry = carbon_stoich;
  QSSASolver solver;
  solver.init(cfg, domain, adv);

  ChemicalField chem = make_chem(domain);
  chem.zero_reactions();

  Agent a = Agent::create_default(1, 1, {10e-6, 10e-6, 10e-6}, 5e-4);
  a.state = PhenoState::NORMAL;
  a.mu_realized = 3e-4;
  a.biomass = 2e-16;
  Int ix = 0, iy = 0, iz = 0;
  domain.pos_to_grid(a.x, ix, iy, iz);
  a.grid_cell = domain.cell_index(ix, iy, iz);
  out_cell = a.grid_cell;
  out_consumption = a.mu_realized * a.biomass;

  AgentPool agents;
  agents.push_back(a);

  OxygenConfig oxygen;
  solver.solve_nutrient_depletion(agents, chem, oxygen);

  const Int i_carbon = chem.find("carbon");
  assert(i_carbon >= 0);
  return chem.reac(i_carbon, out_cell);
}

}  // namespace

void test_carbon_stoichiometry_golden_and_sensitivity() {
  const Domain domain = make_domain();

  // Golden: default carbon stoichiometry (0.5) → reac = -consumption * 0.5.
  Real consumption = 0.0;
  Int cell = -1;
  const Real reac_default = carbon_reac_for_stoich(domain, 0.5, consumption, cell);
  assert(consumption > 0.0);
  assert(std::abs(reac_default - (-consumption * 0.5)) <= 1e-12 * std::abs(consumption));

  // Sensitivity: doubling the stoichiometry doubles the deposited consumption.
  Real c2 = 0.0; Int cell2 = -1;
  const Real reac_double = carbon_reac_for_stoich(domain, 1.0, c2, cell2);
  assert(std::abs(reac_double - (-consumption * 1.0)) <= 1e-12 * std::abs(consumption));
  assert(std::abs(reac_double - 2.0 * reac_default) <= 1e-12 * std::abs(reac_default));

  // Zeroing the stoichiometry removes carbon consumption entirely.
  Real c3 = 0.0; Int cell3 = -1;
  const Real reac_zero = carbon_reac_for_stoich(domain, 0.0, c3, cell3);
  assert(std::abs(reac_zero) <= 1e-30);

  std::cout << "  test_carbon_stoichiometry_golden_and_sensitivity: PASSED\n";
}

void test_qssa_fallback_defaults() {
  QSSAConfig cfg;
  assert(std::abs(cfg.iron_stoichiometry - 1e-6) < 1e-18);
  assert(std::abs(cfg.b12_stoichiometry - 1e-9) < 1e-21);
  assert(std::abs(cfg.carbon_stoichiometry - 0.5) < 1e-15);
  assert(std::abs(cfg.fallback_diff_coeff - 4e-11) < 1e-23);
  assert(std::abs(cfg.fallback_pI - 7.0) < 1e-12);
  assert(std::abs(cfg.fallback_retardation - 5.0) < 1e-12);
  std::cout << "  test_qssa_fallback_defaults: PASSED\n";
}

int main() {
  std::cout << "=== QSSA Stoichiometry Tests ===\n";
  test_carbon_stoichiometry_golden_and_sensitivity();
  test_qssa_fallback_defaults();
  std::cout << "All QSSA stoichiometry tests passed.\n";
  return 0;
}
