#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"
#include "fix_metabolism.h"
#include "uptake_limit.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

using namespace gutibm;

namespace {

constexpr Real kDt = 60.0;

SimulationConfig base_config() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {100.0e-6, 100.0e-6, 40.0e-6};
  cfg.hdf5.enabled = false;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 1;
  strain.mu_max = 0.0;
  cfg.initial_strains.push_back(strain);
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0;
      chemical.boundary_conc = 0.0;
    }
  }
  return cfg;
}

struct MaintenanceResult {
  Real draw = 0.0;
  Real shortfall = 0.0;
  Real limited_agents = 0.0;
  Real biomass = 0.0;
  Real carbon_reaction = 0.0;
};

struct TrajectoryResult {
  Real biomass = 0.0;
  std::vector<Real> carbon;
};

struct LedgerResult {
  Real biomass = 0.0;
  Real maintenance = 0.0;
  Real shortfall = 0.0;
  Real limited_agents = 0.0;
  Real expected = 0.0;
};

MaintenanceResult carbon_maintenance_step(SimulationConfig cfg, Real rate,
                                          Real carbon_concentration,
                                          UptakeLimitMode mode =
                                              UptakeLimitMode::None) {
  cfg.fixes.metabolism.carbon_maintenance_rate = rate;
  cfg.fixes.metabolism.uptake_limit_mode = mode;
  cfg.fixes.metabolism.uptake_limit = mode == UptakeLimitMode::Voxel
      ? "voxel"
      : mode == UptakeLimitMode::Sherwood ? "sherwood" : "none";
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = carbon_concentration;
      chemical.boundary_conc = 0.0;
    }
  }
  Simulation sim;
  sim.init(cfg);
  const Real biomass = sim.agents()[0].biomass;
  FixMetabolism metabolism(sim, sim.config().fixes.metabolism);
  metabolism.compute(kDt);
  const Int carbon = sim.chemical_field().find(species::CARBON);
  assert(carbon >= 0);
  const auto& flux = sim.chemical_field().flux_accounting();
  const Real reaction = sim.chemical_field().reac_global(
      carbon, sim.agents()[0].grid_cell);
  return {
      flux.maintenance_step[static_cast<size_t>(carbon)],
      flux.maintenance_shortfall_step[static_cast<size_t>(carbon)],
      flux.maintenance_limited_agents_step[static_cast<size_t>(carbon)],
      biomass,
      reaction};
}

void test_non_growing_agent_is_charged() {
  const Real rate = 2.0e-5;
  SimulationConfig cfg = base_config();
  cfg.fixes.metabolism.maintenance_rate = 1.0e-3;
  const MaintenanceResult result =
      carbon_maintenance_step(cfg, rate, 1.0);
  const Real expected = rate * result.biomass * kDt;
  assert(result.draw > 0.0);
  assert(std::abs(result.draw - expected) <= expected * 1.0e-12);
  std::cout << "  test_non_growing_agent_is_charged: PASSED\n";
}

void test_maintenance_rate_is_independent() {
  SimulationConfig low_tax = base_config();
  SimulationConfig high_tax = low_tax;
  high_tax.fixes.metabolism.maintenance_rate = 1.0e-3;
  const Real rate = 2.0e-5;
  const Real low = carbon_maintenance_step(low_tax, rate, 1.0).draw;
  const Real high = carbon_maintenance_step(high_tax, rate, 1.0).draw;
  assert(std::abs(low - high) <= std::max(low, high) * 1.0e-12);
  std::cout << "  test_maintenance_rate_is_independent: PASSED\n";
}

void test_maintenance_sensitivity_and_shortfall() {
  const std::vector<Real> rates = {0.0, 1.0e-6, 1.0e-5, 1.0e-4};
  std::vector<Real> draws;
  draws.reserve(rates.size());
  for (const Real rate : rates) {
    const MaintenanceResult result =
        carbon_maintenance_step(base_config(), rate, 1.0);
    assert(std::isfinite(result.draw));
    assert(result.draw >= 0.0);
    draws.push_back(result.draw);
  }
  assert(std::is_sorted(draws.begin(), draws.end()));
  assert(draws.back() > draws.front());

  const MaintenanceResult result =
      carbon_maintenance_step(base_config(), 1.0e10, 1.0e-30,
                              UptakeLimitMode::Voxel);
  assert(result.draw >= 0.0);
  assert(1.0e-30 + result.carbon_reaction * kDt >= -1.0e-45);
  assert(result.shortfall > 0.0);
  assert(std::isfinite(result.shortfall));
  assert(std::isfinite(result.limited_agents));
  assert(result.limited_agents > 0.0);
  std::cout << "  test_maintenance_sensitivity_and_shortfall: PASSED\n";
}

TrajectoryResult run_trajectory(Real rate, bool set_rate = true) {
  SimulationConfig cfg = base_config();
  cfg.time.total_time = 600.0;
  if (set_rate) cfg.fixes.metabolism.carbon_maintenance_rate = rate;
  cfg.initial_strains[0].count = 2;
  cfg.initial_strains[0].mu_max = 5.0e-4;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) chemical.initial_conc = 1.0e-6;
  }
  Simulation sim;
  sim.init(cfg);
  sim.run();
  const Int carbon = sim.chemical_field().find(species::CARBON);
  assert(carbon >= 0);
  std::vector<Real> field(static_cast<size_t>(sim.chemical_field().ncells()));
  for (Int cell = 0; cell < sim.chemical_field().ncells(); ++cell) {
    field[static_cast<size_t>(cell)] = sim.chemical_field().conc(carbon, cell);
  }
  return {sim.agents()[0].biomass, std::move(field)};
}

LedgerResult run_ledger_probe() {
  SimulationConfig cfg = base_config();
  cfg.time.total_time = 600.0;
  cfg.fixes.metabolism.carbon_maintenance_rate = 1.0e-6;
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Voxel;
  cfg.fixes.metabolism.uptake_limit = "voxel";
  cfg.initial_strains[0].mu_max = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e-30;
      chemical.boundary_conc = 0.0;
    }
  }
  Simulation sim;
  sim.init(cfg);
  const Real biomass = sim.agents()[0].biomass;
  FixMetabolism metabolism(sim, sim.config().fixes.metabolism);
  Real expected = 0.0;
  for (int step = 0; step < 10; ++step) {
    expected += 1.0e-6 * sim.agents()[0].biomass * kDt;
    metabolism.compute(kDt);
    sim.chemical_field().flux_accounting().commit_agent_uptake_step();
    sim.chemical_field().flux_accounting().close_interval();
  }
  const Int carbon = sim.chemical_field().find(species::CARBON);
  assert(carbon >= 0);
  const auto& flux = sim.chemical_field().flux_accounting();
  return {biomass,
          flux.maintenance_cumulative[static_cast<size_t>(carbon)],
          flux.maintenance_shortfall_cumulative[static_cast<size_t>(carbon)],
          flux.maintenance_limited_agents_cumulative[
              static_cast<size_t>(carbon)],
          expected};
}

void test_trajectory_sensitivity_and_zero_compatibility() {
  const std::vector<Real> rates = {0.0, 1.0e-4, 1.0e-3, 1.0e-2};
  std::vector<TrajectoryResult> results;
  results.reserve(rates.size());
  for (const Real rate : rates) results.push_back(run_trajectory(rate));
  for (size_t i = 1; i < results.size(); ++i) {
    assert(results[i].biomass <= results[i - 1].biomass);
  }
  assert(results.front().biomass > results.back().biomass);
  const TrajectoryResult absent = run_trajectory(0.0, false);
  assert(results[0].carbon == absent.carbon);
  assert(results[0].biomass == absent.biomass);
  std::cout << "  test_trajectory_sensitivity_and_zero_compatibility: PASSED\n";
}

void test_maintenance_ledger_closure() {
  const LedgerResult result = run_ledger_probe();
  assert(result.maintenance > 0.0);
  assert(result.shortfall > 0.0);
  assert(result.limited_agents > 0.0);
  assert(std::isfinite(result.limited_agents));
  assert(std::abs(result.maintenance + result.shortfall - result.expected)
         <= result.expected * 1.0e-12);
  std::cout << "  test_maintenance_ledger_closure: PASSED\n";
}

void test_uptake_limit_mode_contrast() {
  constexpr Real rate = 1.0e10;
  constexpr Real concentration = 1.0e-30;
  const MaintenanceResult none = carbon_maintenance_step(
      base_config(), rate, concentration, UptakeLimitMode::None);
  const MaintenanceResult voxel = carbon_maintenance_step(
      base_config(), rate, concentration, UptakeLimitMode::Voxel);
  assert(none.draw > voxel.draw);
  assert(voxel.shortfall > 0.0);
  assert(std::isfinite(none.draw));
  assert(std::isfinite(voxel.draw));
  std::cout << "  test_uptake_limit_mode_contrast: PASSED\n";
}

}  // namespace

int main() {
  std::cout << "=== Carbon Maintenance Tests ===\n";
  test_non_growing_agent_is_charged();
  test_maintenance_rate_is_independent();
  test_maintenance_sensitivity_and_shortfall();
  test_trajectory_sensitivity_and_zero_compatibility();
  test_maintenance_ledger_closure();
  test_uptake_limit_mode_contrast();
  std::cout << "All carbon maintenance tests passed.\n";
  return 0;
}
