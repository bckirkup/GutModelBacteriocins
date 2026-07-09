/* -----------------------------------------------------------------------
   GutIBM – QSSA nutrient-depletion (O2 respiration) tests

   Spec 6 removed the carbon/iron/B12 terms from QSSASolver::solve_nutrient_
   depletion — per-agent uptake of those nutrients is now owned solely by the
   metabolism Fix (yield-based), eliminating the double-count. The only
   per-agent term this function still applies is aerobic O2 respiration, which
   has no counterpart in the metabolism Fix. Respiration follows a Pirt model:
     reac_o2 = -(q_consumption * mu_realized + q_maintenance) / cell_vol
   The growth-associated term (q_consumption * mu) plus a basal maintenance
   term (q_maintenance) that is applied per living cell regardless of growth,
   so the O2 field tracks agent *density* (a non-growing cell still respires).
   These tests pin the term with golden values, config-sensitivity checks, and
   a density-coupling check, and confirm carbon/iron/B12 are NOT touched here.
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
                              bool oxygen_enabled, Real q_maintenance = 0.0,
                              Real mu_realized = 3e-4, int n_agents = 1) {
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

  Agent proto = Agent::create_default(1, 1, {10e-6, 10e-6, 10e-6}, 5e-4);
  proto.state = PhenoState::NORMAL;
  proto.mu_realized = mu_realized;
  proto.biomass = 2e-16;
  Int ix = 0, iy = 0, iz = 0;
  domain.pos_to_grid(proto.x, ix, iy, iz);
  proto.grid_cell = domain.cell_index(ix, iy, iz);

  AgentPool agents;
  for (int k = 0; k < n_agents; ++k) {
    Agent a = proto;
    a.identity.tag = static_cast<TagID>(k + 1);
    agents.push_back(a);
  }

  OxygenConfig oxygen;
  oxygen.enabled = oxygen_enabled;
  oxygen.q_consumption = q_consumption;
  oxygen.q_maintenance = q_maintenance;
  solver.solve_nutrient_depletion(agents, chem, oxygen);

  DepletionResult r;
  r.consumption = std::max(mu_realized, 0.0);
  r.cell_vol = domain.dx() * domain.dx() * domain.dx();
  r.reac_o2 = chem.reac(chem.find(species::OXYGEN), proto.grid_cell);
  r.reac_carbon = chem.reac(chem.find(species::CARBON), proto.grid_cell);
  r.reac_iron = chem.reac(chem.find(species::IRON), proto.grid_cell);
  r.reac_b12 = chem.reac(chem.find(species::B12), proto.grid_cell);
  return r;
}

}  // namespace

void test_o2_respiration_golden_and_sensitivity() {
  const Domain domain = make_domain();

  // Golden (growth-associated term only): reac_o2 = -q_consumption*mu/cell_vol.
  const Real q1 = 1.0e-14;
  const DepletionResult r1 = run_depletion(domain, q1, /*oxygen_enabled=*/true,
                                           /*q_maintenance=*/0.0);
  const Real expected1 = -q1 * r1.consumption / r1.cell_vol;
  assert(r1.consumption > 0.0);
  assert(std::abs(r1.reac_o2 - expected1) <= 1e-9 * std::abs(expected1));
  assert(r1.reac_o2 < 0.0);

  // Sensitivity: doubling q_consumption doubles the growth-associated O2 draw.
  const DepletionResult r2 = run_depletion(domain, 2.0 * q1, true, 0.0);
  const Real expected2 = -2.0 * q1 * r2.consumption / r2.cell_vol;
  assert(std::abs(r2.reac_o2 - expected2) <= 1e-9 * std::abs(expected2));
  assert(std::abs(r2.reac_o2 - 2.0 * r1.reac_o2) <= 1e-9 * std::abs(r1.reac_o2));

  // Combined Pirt model: reac_o2 = -(q_consumption*mu + q_maintenance)/cell_vol.
  const Real qm = 5.0e-18;
  const DepletionResult rc = run_depletion(domain, q1, true, qm);
  const Real expectedc = -(q1 * rc.consumption + qm) / rc.cell_vol;
  assert(std::abs(rc.reac_o2 - expectedc) <= 1e-9 * std::abs(expectedc));
  // Maintenance strictly increases the draw beyond the growth-only term.
  assert(rc.reac_o2 < r1.reac_o2);

  // Disabled oxygen removes the term entirely.
  const DepletionResult r0 = run_depletion(domain, q1, /*oxygen_enabled=*/false,
                                           qm);
  assert(std::abs(r0.reac_o2) <= 1e-30);

  std::cout << "  test_o2_respiration_golden_and_sensitivity: PASSED\n";
}

// Regression for the bug Edison flagged: respiration was purely growth-coupled
// (reac_o2 = -q_consumption*mu/cell_vol), so a present-but-non-growing cell
// (mu_realized <= 0) drew zero O2 and the field stopped tracking density. The
// maintenance term must keep O2 consumption proportional to cell *count* even
// when growth is zero.
void test_o2_maintenance_tracks_density() {
  const Domain domain = make_domain();
  const Real qm = 1.0e-18;

  // Non-growing cell: growth term is zero, but maintenance still respires.
  const DepletionResult stalled =
      run_depletion(domain, 1.0e-14, /*oxygen_enabled=*/true, qm,
                    /*mu_realized=*/0.0, /*n_agents=*/1);
  assert(stalled.consumption == 0.0);
  assert(stalled.reac_o2 < 0.0);  // still consumes O2 despite mu==0
  const Real expected1 = -qm / stalled.cell_vol;
  assert(std::abs(stalled.reac_o2 - expected1) <= 1e-9 * std::abs(expected1));

  // Density coupling: N non-growing cells in a cell draw N * the maintenance
  // rate. This is what makes the O2 field track agent density.
  const DepletionResult four =
      run_depletion(domain, 1.0e-14, true, qm, /*mu_realized=*/0.0,
                    /*n_agents=*/4);
  assert(std::abs(four.reac_o2 - 4.0 * stalled.reac_o2)
         <= 1e-9 * std::abs(4.0 * stalled.reac_o2));

  std::cout << "  test_o2_maintenance_tracks_density: PASSED\n";
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
  test_o2_maintenance_tracks_density();
  test_carbon_iron_b12_not_depleted_by_qssa();
  std::cout << "All QSSA nutrient-depletion tests passed.\n";
  return 0;
}
