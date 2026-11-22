/* -----------------------------------------------------------------------
   GutIBM – Metabolism wet-mass convention and delivery-mode carbon closure
   tests (GAP-07 / #469).

   Verifies:
     1. Agent::biomass is standardized as wet cell mass
        (volume * CELL_DENSITY_DEFAULT = 1,100 kg/m^3), matching Agent::mass
        and radius updates during growth.
     2. Exact ChemicalField carbon molar balance (residual ~ 0) and zero
        reaction clipping (reaction_clip == 0) over a multi-step
        delivery-mode growth window with closure.enforce_reaction_clip = true.
     3. Graded sensitivity across metabolism.yield_carbon with preserved
        molar closure and zero clipping across the sweep.
   ----------------------------------------------------------------------- */

#include "agent.h"
#include "chemical_field.h"
#include "domain.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

namespace {

Real total_species_moles(const Simulation& sim, Int species_idx) {
  const ChemicalField& chem = sim.chemical_field();
  const Real cell_vol = sim.domain().cell_volume();
  Real total = 0.0;
  for (Int c = 0; c < chem.global_ncells(); ++c) {
    total += chem.conc_global(species_idx, c) * cell_vol;
  }
  return total;
}

SimulationConfig delivery_growth_config(Real yield_carbon) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;
  cfg.seed = 20261121;
  cfg.time.bio_dt = 60.0;
  cfg.time.total_time = 20.0 * 60.0;  // 20 steps = 20 min
  cfg.time.output_interval = 600.0;

  // Match PatchTable carbon-ladder grid (4 um) and 10 um regularized delivery.
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {60.0e-6, 60.0e-6, 60.0e-6};
  cfg.domain.periodic = {true, true, false};
  cfg.domain.grid_dx = 4.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;

  // 2x crypt carbon supply matching PatchTable ladder arm.
  cfg.carbon_boundary_conc = 10.0e-3;
  cfg.chem_env.oxygen.enabled = false;
  cfg.chem_env.acetate.enabled = false;
  cfg.chem_env.siderophore.enabled = false;
  cfg.advection.crypts_enabled = false;
  cfg.advection.radial_turnover = 1.0e9;  // isolate growth/chemistry closure
  cfg.advection.distal_transit_time = 1.0e9;

  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = 1.0e-5;
  cfg.fixes.metabolism.yield_carbon = yield_carbon;

  cfg.closure.enforce_reaction_clip = true;
  cfg.closure.reaction_clip_tolerance_fraction = 0.0;
  cfg.closure.enforce_delivery_realization = true;
  cfg.dysbiosis_threshold = 0.0;

  for (auto& spec : cfg.chemicals) {
    spec.delivery_enabled = false;
    if (spec.name == species::CARBON) {
      spec.initial_conc = cfg.carbon_boundary_conc;
      spec.boundary_conc = cfg.carbon_boundary_conc;
      spec.delivery_enabled = true;
    }
  }

  cfg.enabled_fixes = {"metabolism", "mechanics"};
  cfg.initial_population.placement = "z_slab";
  cfg.initial_population.z_min = 4.0e-6;
  cfg.initial_population.z_max = 20.0e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 30;
  strain.mu_max = 5.0e-4;
  cfg.initial_strains.push_back(strain);
  return cfg;
}

void test_wet_mass_biomass_convention() {
  const Agent a = Agent::create_default(1, 1, Vec3{10e-6, 10e-6, 10e-6}, 5.0e-4);
  const Real expected_vol = (4.0 / 3.0) * PI * a.radius * a.radius * a.radius;
  const Real expected_wet_mass = expected_vol * CELL_DENSITY_DEFAULT;
  assert(std::abs(CELL_DENSITY_DEFAULT - 1100.0) < 1e-12);
  assert(std::abs(a.biomass - expected_wet_mass) <= 1e-15 * expected_wet_mass);
  assert(std::abs(a.mass - a.biomass) <= 1e-15 * expected_wet_mass);

  // Verify wet-mass identity persists after metabolic growth steps.
  SimulationConfig cfg = delivery_growth_config(0.5);
  cfg.time.total_time = 5.0 * cfg.time.bio_dt;
  Simulation sim;
  sim.init(cfg);
  for (int s = 0; s < 5; ++s) {
    sim.step(cfg.time.bio_dt);
  }
  assert(sim.agents().size() >= 30);
  for (const Agent& cell : sim.agents()) {
    if (cell.state == PhenoState::DEAD) continue;
    assert(cell.biomass > 0.0);
    assert(std::abs(cell.mass - cell.biomass) <= 1e-12 * cell.biomass);
    const Real vol = (4.0 / 3.0) * PI * cell.radius * cell.radius * cell.radius;
    const Real wet_from_radius = vol * CELL_DENSITY_DEFAULT;
    assert(std::abs(wet_from_radius - cell.biomass) <= 1e-11 * cell.biomass);
  }
  std::cout << "test_wet_mass_biomass_convention: PASSED\n";
}

void test_delivery_carbon_molar_balance_and_zero_clip() {
  // Part A: step-by-step exact carbon molar residual and zero reaction clip.
  {
    SimulationConfig cfg = delivery_growth_config(0.5);
    Simulation sim;
    sim.init(cfg);
    const Int carbon = sim.chemical_field().find(species::CARBON);
    assert(carbon >= 0);
    assert(sim.chemical_field().spec(carbon).delivery_enabled);

    Real max_abs_residual = 0.0;
    Real max_rel_residual = 0.0;
    Real max_clip = 0.0;
    Real total_uptake = 0.0;

    for (int step = 0; step < 60; ++step) {
      const Real before_mol = total_species_moles(sim, carbon);
      sim.step(cfg.time.bio_dt);
      const Real after_mol = total_species_moles(sim, carbon);

      const auto& flux = sim.chemical_field().flux_accounting();
      const auto idx = static_cast<size_t>(carbon);
      const Real clip = flux.reaction_clip_for_step(carbon);
      const Real residual =
          before_mol + flux.boundary_last_step[idx]
          + flux.gradient_source_last_step[idx]
          + flux.vbf_source_last_step[idx]
          - flux.agent_uptake_last_step[idx]
          - flux.maintenance_last_step[idx]
          - flux.vbf_sink_last_step[idx]
          + clip
          - after_mol;

      const Real scale = std::max({before_mol, after_mol, 1.0e-18});
      const Real rel = std::abs(residual) / scale;
      max_abs_residual = std::max(max_abs_residual, std::abs(residual));
      max_rel_residual = std::max(max_rel_residual, rel);
      max_clip = std::max(max_clip, clip);
      total_uptake += flux.agent_uptake_last_step[idx];

      // Report and verify clip separately from molar balance residual.
      assert(clip == 0.0);
      assert(std::isfinite(residual));
      assert(rel < 1.0e-12);
    }

    const auto& flux = sim.chemical_field().flux_accounting();
    const auto idx = static_cast<size_t>(carbon);
    assert(flux.reaction_clip_interval[idx] == 0.0);
    assert(total_uptake > 0.0);
    const Int total_divs =
        sim.step_events().divisions + sim.cumulative_events().divisions;
    assert(total_divs > 0);
    std::cout << "test_delivery_carbon_molar_balance_and_zero_clip(step): "
              << "divisions=" << total_divs
              << " uptake_mol=" << total_uptake
              << " max_rel_residual=" << max_rel_residual
              << " max_abs_residual=" << max_abs_residual
              << " max_clip=" << max_clip << '\n';
  }

  // Part B: run() loop with enforce_reaction_clip = true and zero tolerance
  // must reach HorizonReached without tripping ClosureViolation.
  {
    SimulationConfig cfg = delivery_growth_config(0.5);
    Simulation sim;
    sim.init(cfg);
    const int rc = sim.run();
    assert(rc == 0);
    assert(sim.termination_cause() == TerminationCause::HorizonReached);
    const Int carbon = sim.chemical_field().find(species::CARBON);
    assert(carbon >= 0);
    const auto& flux = sim.chemical_field().flux_accounting();
    const auto idx = static_cast<size_t>(carbon);
    assert(flux.reaction_clip_cumulative[idx] == 0.0);
    assert(flux.reaction_clip_interval[idx] == 0.0);
    assert(flux.agent_uptake_cumulative[idx] + flux.agent_uptake_interval[idx] >
           0.0);
  }
  std::cout << "test_delivery_carbon_molar_balance_and_zero_clip: PASSED\n";
}

void test_delivery_yield_carbon_graded_sensitivity() {
  const std::vector<Real> yields = {0.25, 0.50, 1.00};
  std::vector<Real> uptakes;
  uptakes.reserve(yields.size());

  for (const Real y : yields) {
    SimulationConfig cfg = delivery_growth_config(y);
    Simulation sim;
    sim.init(cfg);
    const Int carbon = sim.chemical_field().find(species::CARBON);
    assert(carbon >= 0);
    Real total_uptake = 0.0;
    for (int step = 0; step < 15; ++step) {
      const Real before_mol = total_species_moles(sim, carbon);
      sim.step(cfg.time.bio_dt);
      const Real after_mol = total_species_moles(sim, carbon);
      const auto& flux = sim.chemical_field().flux_accounting();
      const auto idx = static_cast<size_t>(carbon);
      const Real clip = flux.reaction_clip_for_step(carbon);
      const Real residual =
          before_mol + flux.boundary_last_step[idx]
          + flux.gradient_source_last_step[idx]
          + flux.vbf_source_last_step[idx]
          - flux.agent_uptake_last_step[idx]
          - flux.maintenance_last_step[idx]
          - flux.vbf_sink_last_step[idx]
          + clip
          - after_mol;
      const Real scale = std::max({before_mol, after_mol, 1.0e-18});
      assert(clip == 0.0);
      assert(std::abs(residual) / scale < 1.0e-12);
      total_uptake += flux.agent_uptake_last_step[idx];
    }
    assert(std::isfinite(total_uptake) && total_uptake > 0.0);
    uptakes.push_back(total_uptake);
  }

  // Higher yield_carbon (mol substrate / kg wet biomass) demands more carbon
  // per unit growth -> strictly increasing cumulative carbon uptake.
  assert(uptakes[0] < uptakes[1]);
  assert(uptakes[1] < uptakes[2]);
  assert(uptakes[2] > 1.5 * uptakes[0]);
  std::cout << "test_delivery_yield_carbon_graded_sensitivity: uptakes=["
            << uptakes[0] << ", " << uptakes[1] << ", " << uptakes[2]
            << "] PASSED\n";
}

}  // namespace

int main() {
  test_wet_mass_biomass_convention();
  test_delivery_carbon_molar_balance_and_zero_clip();
  test_delivery_yield_carbon_graded_sensitivity();
  std::cout << "test_metabolism: all passed\n";
  return 0;
}
