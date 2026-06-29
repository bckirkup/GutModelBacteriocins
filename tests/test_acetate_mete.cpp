/* -----------------------------------------------------------------------
   GutIBM – Tests for acetate inhibition of MetE pathway (issue #8)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace gutibm;

// Helper: build a minimal sim config with BtuB-downregulated agents
static SimulationConfig make_acetate_cfg(Real acetate_conc) {
  SimulationConfig cfg = InputParser::default_config();

  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;

  cfg.total_time      = 120.0;
  cfg.bio_dt          = 60.0;
  cfg.output_interval = 120.0;
  cfg.seed            = 42;
  cfg.hdf5.enabled    = false;

  cfg.advection.mucus_thickness    = 25e-6;
  cfg.advection.distal_length      = 50e-6;
  cfg.advection.radial_turnover    = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;

  cfg.qssa.toxin_cutoff    = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  // Override acetate initial/boundary concentration
  for (auto& spec : cfg.chemicals) {
    if (spec.name == "acetate") {
      spec.initial_conc  = acetate_conc;
      spec.boundary_conc = acetate_conc;
    }
  }

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 5;
  s.mu_max = 5e-4;
  s.plasmids = {};
  s.conjugative = false;
  cfg.initial_strains.push_back(s);

  return cfg;
}

void test_acetate_species_present() {
  SimulationConfig cfg = InputParser::default_config();
  bool found = false;
  for (const auto& spec : cfg.chemicals) {
    if (spec.name == "acetate") {
      found = true;
      assert(std::abs(spec.diff_coeff - 1.2e-9) < 1e-15);
      assert(std::abs(spec.retardation - 1.0) < 1e-15);
      assert(std::abs(spec.initial_conc - 80.0) < 1e-6);
      assert(std::abs(spec.boundary_conc - 80.0) < 1e-6);
      assert(std::abs(spec.decay_rate) < 1e-15);
    }
  }
  assert(found);
  std::cout << "  test_acetate_species_present: PASSED\n";
}

void test_metE_config_defaults() {
  MetabolismConfig mcfg;
  assert(std::abs(mcfg.metE_penalty - 0.05) < 1e-12);
  assert(std::abs(mcfg.metE_acetate_km - 40.0) < 1e-12);
  assert(std::abs(mcfg.metE_acetate_max_factor - 2.5) < 1e-12);
  assert(std::abs(mcfg.eut_km - 0.1e-3) < 1e-15);
  assert(std::abs(mcfg.eut_max_penalty - 0.10) < 1e-12);
  std::cout << "  test_metE_config_defaults: PASSED\n";
}

void test_acetate_increases_penalty() {
  // Run two simulations: zero acetate vs 80 mM acetate.
  // With BtuB-downregulated agents, higher acetate should yield
  // lower realized growth rates (larger MetE penalty).
  // Use high mu_max so agents survive the penalty + maintenance.

  SimulationConfig cfg_low  = make_acetate_cfg(0.0);
  SimulationConfig cfg_high = make_acetate_cfg(80.0);

  // Use high mu_max so net growth stays positive
  for (auto& s : cfg_low.initial_strains)  s.mu_max = 5e-3;
  for (auto& s : cfg_high.initial_strains) s.mu_max = 5e-3;

  Simulation sim_low;
  sim_low.init(cfg_low);
  for (Agent& a : sim_low.agents()) {
    a.receptor_expr[to_underlying(ReceptorType::BtuB)] = 0.1;
  }

  Simulation sim_high;
  sim_high.init(cfg_high);
  for (Agent& a : sim_high.agents()) {
    a.receptor_expr[to_underlying(ReceptorType::BtuB)] = 0.1;
  }

  // Record initial biomass
  Real bm_low_before = 0.0;
  Real bm_high_before = 0.0;
  for (const Agent& a : sim_low.agents())
    bm_low_before += a.biomass;
  for (const Agent& a : sim_high.agents())
    bm_high_before += a.biomass;

  sim_low.step(60.0);
  sim_high.step(60.0);

  // Collect total biomass after step (surviving agents)
  Real bm_low_after = 0.0;
  Real bm_high_after = 0.0;
  for (const Agent& a : sim_low.agents())
    bm_low_after += a.biomass;
  for (const Agent& a : sim_high.agents())
    bm_high_after += a.biomass;

  Real growth_low  = bm_low_after - bm_low_before;
  Real growth_high = bm_high_after - bm_high_before;

  // Higher acetate → larger MetE penalty → less biomass growth
  assert(growth_high < growth_low);

  std::cout << "  test_acetate_increases_penalty: PASSED"
            << " (growth_no_acetate=" << growth_low
            << " growth_80mM=" << growth_high << ")\n";
}

void test_acetate_penalty_scaling() {
  // Verify that the Michaelis-Menten scaling produces expected values.
  // At 80 mM with Km=40, factor = 1 + 1.5 * 80/(40+80) = 1 + 1.0 = 2.0
  // So effective penalty = 0.05 * 2.0 = 0.10 (10%)
  MetabolismConfig mcfg;
  Real acetate_conc = 80.0;
  Real factor = 1.0
      + (mcfg.metE_acetate_max_factor - 1.0)
        * acetate_conc / (mcfg.metE_acetate_km + acetate_conc);
  Real eff = mcfg.metE_penalty * factor;

  // Should be 10% at 80 mM (within the 8-12% spec from issue)
  assert(eff > 0.08 && eff < 0.12);
  assert(std::abs(eff - 0.10) < 1e-12);

  // At zero acetate, factor = 1.0, penalty stays at 5%
  Real factor_zero = 1.0
      + (mcfg.metE_acetate_max_factor - 1.0)
        * 0.0 / (mcfg.metE_acetate_km + 0.0);
  assert(std::abs(factor_zero - 1.0) < 1e-12);
  assert(std::abs(mcfg.metE_penalty * factor_zero - 0.05) < 1e-12);

  // At saturating acetate (>> Km), factor → 2.5, penalty → 12.5%
  Real sat_conc = 1e6;
  Real factor_sat = 1.0
      + (mcfg.metE_acetate_max_factor - 1.0)
        * sat_conc / (mcfg.metE_acetate_km + sat_conc);
  Real eff_sat = mcfg.metE_penalty * factor_sat;
  assert(std::abs(eff_sat - 0.125) < 0.001);

  std::cout << "  test_acetate_penalty_scaling: PASSED"
            << " (eff_80mM=" << eff
            << " eff_sat=" << eff_sat << ")\n";
}

void test_smoke_with_acetate() {
  // Full mini simulation with acetate species present
  SimulationConfig cfg = make_acetate_cfg(80.0);
  cfg.total_time = 600.0;
  cfg.output_interval = 300.0;

  // Add a producer strain too
  SimulationConfig::InitialStrain producer;
  producer.type = 2;
  producer.count = 5;
  producer.mu_max = 5e-4;
  producer.plasmids = {"ColE1"};
  producer.conjugative = false;
  cfg.initial_strains.push_back(producer);

  Simulation sim;
  sim.init(cfg);
  sim.run();

  assert(sim.time() > 0.0);
  assert(sim.step_count() > 0);

  Int alive = 0;
  for (const Agent& a : sim.agents()) {
    if (a.state != PhenoState::DEAD) alive++;
  }
  assert(alive > 0);

  // Verify acetate species exists in chemical field
  Int i_acetate = sim.chemical_field().find("acetate");
  assert(i_acetate >= 0);

  std::cout << "  test_smoke_with_acetate: PASSED"
            << " (alive=" << alive
            << " steps=" << sim.step_count() << ")\n";
}

int main() {
  std::cout << "=== Acetate MetE Inhibition Tests ===\n";
  test_acetate_species_present();
  test_metE_config_defaults();
  test_acetate_penalty_scaling();
  test_acetate_increases_penalty();
  test_smoke_with_acetate();
  std::cout << "All acetate MetE tests passed.\n";
  return 0;
}
