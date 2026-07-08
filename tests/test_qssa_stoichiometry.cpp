/* -----------------------------------------------------------------------
   GutIBM – QSSA nutrient-depletion (O2 respiration) tests

   Spec 6 removed the carbon/iron/B12 terms from QSSASolver::solve_nutrient_
   depletion — per-agent uptake of those nutrients is now owned solely by the
   metabolism Fix (yield-based), eliminating the double-count. The only
   per-agent term this function still applies is aerobic O2 respiration
   (q_consumption * mu_realized / cell_vol), which has no counterpart in the
   metabolism Fix. These tests pin that term with a golden value and a
   configuration-sensitivity check, and confirm carbon/iron/B12 are NOT touched
   here anymore.
   ----------------------------------------------------------------------- */

#include "qssa_solver.h"
#include "chemical_field.h"
#include "domain.h"
#include "advection.h"
#include "agent.h"
#include "chem_environment_config.h"
#include "species_names.h"
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
  ChemicalSpec carbon; carbon.name = species::CARBON; carbon.diff_coeff = 1e-9;
  carbon.retardation = 1.0; carbon.initial_conc = 5e-3; carbon.boundary_conc = 5e-3;
  ChemicalSpec iron = carbon; iron.name = species::IRON;
  ChemicalSpec b12 = carbon; b12.name = species::B12;
  ChemicalSpec oxygen = carbon; oxygen.name = species::OXYGEN;
  oxygen.initial_conc = 55e-6; oxygen.boundary_conc = 55e-6;
  ChemicalField chem;
  chem.init(domain, {carbon, iron, b12, oxygen});
  return chem;
}

// Run one nutrient-depletion pass for a single agent. Returns the O2 reaction
// rate deposited at the agent's cell; also reports carbon/iron/B12 reac so the
// caller can assert they are untouched.
struct DepletionResult {
  Real reac_o2 = 0.0;
  Real reac_carbon = 0.0;
  Real reac_iron = 0.0;
  Real reac_b12 = 0.0;
  Real consumption = 0.0;  // mu_realized (O2 term scales with this)
  Real cell_vol = 0.0;
};

DepletionResult run_depletion(const Domain& domain, Real q_consumption,
                              bool oxygen_enabled) {
  AdvectionField adv;
  AdvectionConfig acfg;
  acfg.radial_turnover = 1e20;
  acfg.distal_transit_time = 1e20;
  adv.init(acfg, domain);

  QSSAConfig cfg;
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

  AgentPool agents;
  agents.push_back(a);

  OxygenConfig oxygen;
  oxygen.enabled = oxygen_enabled;
  oxygen.q_consumption = q_consumption;
  solver.solve_nutrient_depletion(agents, chem, oxygen);

  DepletionResult r;
  r.consumption = std::max(a.mu_realized, 0.0);
  r.cell_vol = domain.dx() * domain.dx() * domain.dx();
  r.reac_o2 = chem.reac(chem.find(species::OXYGEN), a.grid_cell);
  r.reac_carbon = chem.reac(chem.find(species::CARBON), a.grid_cell);
  r.reac_iron = chem.reac(chem.find(species::IRON), a.grid_cell);
  r.reac_b12 = chem.reac(chem.find(species::B12), a.grid_cell);
  return r;
}

}  // namespace

void test_o2_respiration_golden_and_sensitivity() {
  const Domain domain = make_domain();

  // Golden: reac_o2 = -q_consumption * mu_realized / cell_vol.
  const Real q1 = 1.0e-14;
  const DepletionResult r1 = run_depletion(domain, q1, /*oxygen_enabled=*/true);
  const Real expected1 = -q1 * r1.consumption / r1.cell_vol;
  assert(r1.consumption > 0.0);
  assert(std::abs(r1.reac_o2 - expected1) <= 1e-9 * std::abs(expected1));
  assert(r1.reac_o2 < 0.0);

  // Sensitivity: doubling q_consumption doubles the O2 draw.
  const DepletionResult r2 = run_depletion(domain, 2.0 * q1, true);
  const Real expected2 = -2.0 * q1 * r2.consumption / r2.cell_vol;
  assert(std::abs(r2.reac_o2 - expected2) <= 1e-9 * std::abs(expected2));
  assert(std::abs(r2.reac_o2 - 2.0 * r1.reac_o2) <= 1e-9 * std::abs(r1.reac_o2));

  // Disabled oxygen removes the term entirely.
  const DepletionResult r0 = run_depletion(domain, q1, /*oxygen_enabled=*/false);
  assert(std::abs(r0.reac_o2) <= 1e-30);

  std::cout << "  test_o2_respiration_golden_and_sensitivity: PASSED\n";
}

void test_carbon_iron_b12_not_depleted_by_qssa() {
  const Domain domain = make_domain();
  const DepletionResult r = run_depletion(domain, 1.0e-14, /*oxygen_enabled=*/true);

  // Spec 6: solve_nutrient_depletion must not touch carbon/iron/B12 — those are
  // consumed only by the metabolism Fix. Only O2 is deposited here.
  assert(std::abs(r.reac_carbon) <= 1e-30);
  assert(std::abs(r.reac_iron) <= 1e-30);
  assert(std::abs(r.reac_b12) <= 1e-30);
  assert(r.reac_o2 < 0.0);

  std::cout << "  test_carbon_iron_b12_not_depleted_by_qssa: PASSED\n";
}

int main() {
  std::cout << "=== QSSA Nutrient-Depletion (O2) Tests ===\n";
  test_o2_respiration_golden_and_sensitivity();
  test_carbon_iron_b12_not_depleted_by_qssa();
  std::cout << "All QSSA nutrient-depletion tests passed.\n";
  return 0;
}
