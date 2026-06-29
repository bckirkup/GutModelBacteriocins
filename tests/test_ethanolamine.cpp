/* -----------------------------------------------------------------------
   GutIBM -- Tests for nutrient-conditional ethanolamine utilization penalty
   Verifies that higher ambient [ethanolamine] produces a larger growth
   penalty when BtuB is downregulated (issue #7).
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace gutibm;

// Helper: build a tiny 1-agent simulation, set ethanolamine concentration,
// step once, and return the realized growth rate.
static Real measure_mu(Real eut_initial_conc) {
  SimulationConfig cfg = InputParser::default_config();

  // Tiny domain
  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {20e-6, 20e-6, 10e-6};
  cfg.domain.grid_dx = 10e-6;
  cfg.domain.hash_cell_size = 10e-6;

  cfg.total_time      = 60.0;
  cfg.bio_dt          = 60.0;
  cfg.output_interval = 60.0;
  cfg.seed            = 42;
  cfg.hdf5.enabled    = false;

  cfg.advection.mucus_thickness     = 10e-6;
  cfg.advection.distal_length       = 20e-6;
  cfg.advection.radial_turnover     = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff    = 10e-6;
  cfg.qssa.nutrient_cutoff = 10e-6;

  // Override ethanolamine concentration in the chemical spec
  for (auto& cs : cfg.chemicals) {
    if (cs.name == "ethanolamine") {
      cs.initial_conc  = eut_initial_conc;
      cs.boundary_conc = eut_initial_conc;
    }
  }

  // Single agent with BtuB downregulated
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type       = 1;
  s.count      = 1;
  s.mu_max     = 5.0e-4;
  s.plasmids   = {};
  s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);

  // Downregulate BtuB so the eut penalty path fires (expr < 0.5)
  Agent& a = sim.agents()[0];
  a.receptor_expr[to_underlying(ReceptorType::BtuB)] = 0.1;

  // Run one biological timestep
  sim.step(cfg.bio_dt);

  return sim.agents()[0].mu_realized;
}

void test_higher_eut_means_larger_penalty() {
  Real mu_low  = measure_mu(0.01e-3);   // 0.01 mM -- low ethanolamine
  Real mu_high = measure_mu(5.0e-3);    // 5.0  mM -- inflamed-gut level

  // Higher ethanolamine => bigger penalty => lower mu
  assert(mu_high < mu_low);

  std::cout << "  test_higher_eut_means_larger_penalty: PASSED"
            << " (mu_low=" << mu_low << " mu_high=" << mu_high << ")\n";
}

void test_zero_eut_means_no_eut_penalty() {
  Real mu_zero = measure_mu(0.0);       // no ethanolamine
  Real mu_some = measure_mu(1.0e-3);    // 1 mM

  // With zero ethanolamine, penalty term vanishes; mu should be higher
  assert(mu_zero > mu_some);

  std::cout << "  test_zero_eut_means_no_eut_penalty: PASSED"
            << " (mu_zero=" << mu_zero << " mu_some=" << mu_some << ")\n";
}

void test_ethanolamine_species_present() {
  SimulationConfig cfg = InputParser::default_config();
  bool found = false;
  for (const auto& cs : cfg.chemicals) {
    if (cs.name == "ethanolamine") {
      found = true;
      // Check default parameters from issue spec
      assert(std::abs(cs.diff_coeff - 1.0e-9) < 1e-15);
      assert(std::abs(cs.retardation - 1.0) < 1e-15);
      assert(std::abs(cs.initial_conc - 0.5e-3) < 1e-15);
      assert(std::abs(cs.boundary_conc - 0.5e-3) < 1e-15);
      assert(std::abs(cs.decay_rate) < 1e-15);
    }
  }
  assert(found);
  std::cout << "  test_ethanolamine_species_present: PASSED\n";
}

void test_metabolism_config_defaults() {
  MetabolismConfig mcfg;
  assert(std::abs(mcfg.eut_km - 0.1e-3) < 1e-15);
  assert(std::abs(mcfg.eut_max_penalty - 0.10) < 1e-15);
  std::cout << "  test_metabolism_config_defaults: PASSED\n";
}

void test_smoke_with_ethanolamine() {
  // Full mini-simulation with ethanolamine species present
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.total_time      = 300.0;
  cfg.bio_dt          = 60.0;
  cfg.output_interval = 300.0;
  cfg.seed            = 321;
  cfg.hdf5.enabled    = false;
  cfg.advection.mucus_thickness     = 50e-6;
  cfg.advection.distal_length       = 100e-6;
  cfg.advection.radial_turnover     = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff    = 50e-6;
  cfg.qssa.nutrient_cutoff = 25e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 20; s.mu_max = 5e-4;
  s.plasmids = {"ColE1"}; s.conjugative = false;
  cfg.initial_strains.push_back(s);

  SimulationConfig::InitialStrain s2;
  s2.type = 2; s2.count = 10; s2.mu_max = 5e-4;
  s2.plasmids = {}; s2.conjugative = false;
  cfg.initial_strains.push_back(s2);

  Simulation sim;
  sim.init(cfg);

  // Verify ethanolamine species was registered
  Int i_eut = sim.chemical_field().find("ethanolamine");
  assert(i_eut >= 0);

  sim.run();

  Int alive = 0;
  for (const Agent& a : sim.agents())
    if (a.state != PhenoState::DEAD) alive++;

  assert(alive > 0);
  std::cout << "  test_smoke_with_ethanolamine: PASSED (alive=" << alive << ")\n";
}

int main() {
  std::cout << "=== Ethanolamine Tests ===\n";
  test_ethanolamine_species_present();
  test_metabolism_config_defaults();
  test_higher_eut_means_larger_penalty();
  test_zero_eut_means_no_eut_penalty();
  test_smoke_with_ethanolamine();
  std::cout << "All ethanolamine tests passed.\n";
  return 0;
}
