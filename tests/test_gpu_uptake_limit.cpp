/* -----------------------------------------------------------------------
   GutIBM – CPU/GPU parity for agent-side uptake limitation
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"
#include "gpu_test_support.h"
#include "gpu_diagnostic_format.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;
using gutibm::gpu_diagnostic::format_real;

namespace {

#ifdef GUTIBM_CUDA
SimulationConfig make_config(UptakeLimitMode mode) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = 9311;
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;
  cfg.time.total_time = 300.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 300.0;
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {40e-6, 40e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.fixes.metabolism.carbon_maintenance_rate = 1.0e-3;
  if (mode == UptakeLimitMode::Sherwood) {
    cfg.fixes.metabolism.uptake_limit = "sherwood";
  } else if (mode == UptakeLimitMode::Voxel) {
    cfg.fixes.metabolism.uptake_limit = "voxel";
  } else if (mode == UptakeLimitMode::Delivery) {
    cfg.fixes.metabolism.uptake_limit = "delivery";
    cfg.fixes.metabolism.delivery_far_field_radius = 0.0;
  } else {
    cfg.fixes.metabolism.uptake_limit = "none";
  }
  cfg.fixes.metabolism.uptake_limit_mode = mode;
  cfg.enabled_fixes = {"metabolism"};
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.z_gradient_enabled = false;
      chemical.retardation = 1.0e3;
      chemical.initial_conc = 1.0e-7;
      chemical.boundary_conc = 1.0e-7;
    }
  }
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 8;
  strain.mu_max = 5.0e-4;
  cfg.initial_strains.push_back(strain);
  return cfg;
}

Simulation run(const SimulationConfig& config, bool gpu) {
  SimulationConfig run_config = config;
  run_config.gpu.enabled = gpu;
  Simulation simulation;
  simulation.init(run_config);
  if (gpu) assert(simulation.gpu_active());
  simulation.run();
  return simulation;
}

Real report_residual(const char* field, Real cpu_value, Real gpu_value) {
  const Real absolute_difference = std::abs(cpu_value - gpu_value);
  const Real scale = std::max({1.0, std::abs(cpu_value), std::abs(gpu_value)});
  std::cerr << "[gpu_diag][gpu_uptake_limit][parity]"
            << " field=" << field
            << " cpu=" << format_real(cpu_value)
            << " gpu=" << format_real(gpu_value)
            << " abs_diff=" << format_real(absolute_difference)
            << " rel_diff=" << format_real(absolute_difference / scale) << "\n";
  return absolute_difference / scale;
}

Real total_mu(const Simulation& simulation) {
  Real value = 0.0;
  for (const Agent& agent : simulation.agents()) {
    if (agent.state == PhenoState::DEAD || agent.flags.is_ghost) continue;
    value += agent.mu_realized;
  }
  return value;
}

Real flux_value(const Simulation& simulation,
                const std::vector<Real>& cumulative,
                const std::vector<Real>& interval) {
  const Int carbon = simulation.chemical_field().find(species::CARBON);
  assert(carbon >= 0);
  return cumulative[static_cast<size_t>(carbon)]
      + interval[static_cast<size_t>(carbon)];
}

Real carbon_realized_sink(const Simulation& simulation) {
  const auto& chem = simulation.chemical_field();
  const Int carbon = chem.find(species::CARBON);
  assert(carbon >= 0);
  Real realized = 0.0;
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    realized += chem.sink_realized_global(carbon, cell);
  }
  return realized;
}

struct MaintenanceLedger {
  Real realized = 0.0;
  Real shortfall = 0.0;
  Real limited_agents = 0.0;
};

MaintenanceLedger maintenance_ledger(const Simulation& simulation) {
  const auto& flux = simulation.chemical_field().flux_accounting();
  return {
      flux_value(simulation, flux.maintenance_cumulative,
                 flux.maintenance_interval),
      flux_value(simulation, flux.maintenance_shortfall_cumulative,
                 flux.maintenance_shortfall_interval),
      flux_value(simulation, flux.maintenance_limited_agents_cumulative,
                 flux.maintenance_limited_agents_interval)};
}

void test_sherwood_parity() {
  const SimulationConfig config = make_config(UptakeLimitMode::Sherwood);
  Simulation cpu = run(config, false);
  Simulation gpu = run(config, true);
  assert(gpu.agents_gpu().metabolism_gpu_steps() > 0);

  const auto& cpu_flux = cpu.chemical_field().flux_accounting();
  const auto& gpu_flux = gpu.chemical_field().flux_accounting();
  const Real cpu_demand = flux_value(cpu, cpu_flux.uptake_demand_cumulative,
                                     cpu_flux.uptake_demand_interval);
  const Real gpu_demand = flux_value(gpu, gpu_flux.uptake_demand_cumulative,
                                     gpu_flux.uptake_demand_interval);
  const Real cpu_realized = flux_value(cpu, cpu_flux.agent_uptake_cumulative,
                                       cpu_flux.agent_uptake_interval);
  const Real gpu_realized = flux_value(gpu, gpu_flux.agent_uptake_cumulative,
                                       gpu_flux.agent_uptake_interval);
  const Real cpu_limited = flux_value(cpu, cpu_flux.uptake_limited_cumulative,
                                      cpu_flux.uptake_limited_interval);
  const Real gpu_limited = flux_value(gpu, gpu_flux.uptake_limited_cumulative,
                                      gpu_flux.uptake_limited_interval);
  const MaintenanceLedger cpu_maintenance = maintenance_ledger(cpu);
  const MaintenanceLedger gpu_maintenance = maintenance_ledger(gpu);

  constexpr Real tolerance = 5.0e-8;
  const Real mu_residual = report_residual("total_mu_realized", total_mu(cpu),
                                           total_mu(gpu));
  const Real demand_residual =
      report_residual("carbon_uptake_demand", cpu_demand, gpu_demand);
  const Real realized_residual =
      report_residual("carbon_uptake_realized", cpu_realized, gpu_realized);
  const Real limited_residual =
      report_residual("carbon_uptake_limited_agents", cpu_limited, gpu_limited);
  const Real maintenance_residual =
      report_residual("maintenance_realized", cpu_maintenance.realized,
                      gpu_maintenance.realized);
  const Real shortfall_residual =
      report_residual("maintenance_shortfall_mol", cpu_maintenance.shortfall,
                      gpu_maintenance.shortfall);
  const Real maintenance_limited_residual = report_residual(
      "maintenance_limited_agents", cpu_maintenance.limited_agents,
      gpu_maintenance.limited_agents);

  assert(mu_residual <= tolerance);
  assert(demand_residual <= tolerance);
  assert(realized_residual <= tolerance);
  assert(limited_residual <= tolerance);
  assert(maintenance_residual <= tolerance);
  assert(shortfall_residual <= tolerance);
  assert(maintenance_limited_residual <= tolerance);
  assert(cpu_limited > 0.0);
  assert(cpu_maintenance.shortfall > 0.0);
  assert(cpu_maintenance.limited_agents > 0.0);
  assert(cpu_realized <= cpu_demand);
  assert(gpu_realized <= gpu_demand);
  std::cout << "  test_sherwood_parity: PASSED\n";
}

void test_none_mode_parity_and_absence_of_limitation() {
  const SimulationConfig config = make_config(UptakeLimitMode::None);
  Simulation cpu = run(config, false);
  Simulation gpu = run(config, true);
  assert(gpu.agents_gpu().metabolism_gpu_steps() > 0);

  const auto& cpu_flux = cpu.chemical_field().flux_accounting();
  const auto& gpu_flux = gpu.chemical_field().flux_accounting();
  const Real cpu_limited = flux_value(cpu, cpu_flux.uptake_limited_cumulative,
                                      cpu_flux.uptake_limited_interval);
  const Real gpu_limited = flux_value(gpu, gpu_flux.uptake_limited_cumulative,
                                      gpu_flux.uptake_limited_interval);
  const Real cpu_demand = flux_value(cpu, cpu_flux.uptake_demand_cumulative,
                                     cpu_flux.uptake_demand_interval);
  const Real cpu_realized = flux_value(cpu, cpu_flux.agent_uptake_cumulative,
                                       cpu_flux.agent_uptake_interval);
  const MaintenanceLedger cpu_maintenance = maintenance_ledger(cpu);
  const MaintenanceLedger gpu_maintenance = maintenance_ledger(gpu);
  const Real residual = report_residual("total_mu_realized", total_mu(cpu),
                                        total_mu(gpu));
  assert(report_residual("maintenance_realized", cpu_maintenance.realized,
                         gpu_maintenance.realized) <= 5.0e-8);
  assert(report_residual("maintenance_shortfall_mol", cpu_maintenance.shortfall,
                         gpu_maintenance.shortfall) <= 5.0e-8);
  assert(report_residual("maintenance_limited_agents",
                         cpu_maintenance.limited_agents,
                         gpu_maintenance.limited_agents) <= 5.0e-8);
  std::cerr << "[gpu_diag][gpu_uptake_limit][none]"
            << " cpu_limited_agents=" << format_real(cpu_limited)
            << " gpu_limited_agents=" << format_real(gpu_limited)
            << " cpu_demand=" << format_real(cpu_demand)
            << " cpu_realized=" << format_real(cpu_realized) << "\n";
  assert(residual <= 5.0e-8);
  assert(cpu_limited == 0.0);
  assert(gpu_limited == 0.0);
  assert(cpu_maintenance.shortfall == 0.0);
  assert(gpu_maintenance.shortfall == 0.0);
  assert(cpu_maintenance.limited_agents == 0.0);
  assert(gpu_maintenance.limited_agents == 0.0);
  assert(std::abs(cpu_realized - cpu_demand) <= 1.0e-12 * cpu_demand);
  std::cout << "  test_none_mode_parity_and_absence_of_limitation: PASSED\n";
}

void test_delivery_forces_host_chemistry() {
  const SimulationConfig delivery_config =
      make_config(UptakeLimitMode::Delivery);
  const SimulationConfig none_config = make_config(UptakeLimitMode::None);
  Simulation delivery = run(delivery_config, true);
  Simulation none = run(none_config, true);

  assert(std::string(delivery.chemistry_placement())
         == "host_forced_delivery");
  const auto& delivery_flux = delivery.chemical_field().flux_accounting();
  const Real delivery_realized = flux_value(
      delivery, delivery_flux.agent_uptake_cumulative,
      delivery_flux.agent_uptake_interval);
  assert(delivery_realized > 0.0);
  const Real delivery_field_realized = carbon_realized_sink(delivery);
  assert(std::isfinite(delivery_field_realized));
  assert(delivery_field_realized > 0.0);

  const auto& none_flux = none.chemical_field().flux_accounting();
  const Real none_realized = flux_value(
      none, none_flux.agent_uptake_cumulative,
      none_flux.agent_uptake_interval);
  const Real none_field_realized = carbon_realized_sink(none);
  assert(std::isfinite(none_field_realized));
  assert(none_field_realized == 0.0);
  const Real contrast_scale =
      std::max(std::abs(delivery_realized), std::abs(none_realized));
  assert(std::abs(delivery_realized - none_realized)
         > 1.0e-6 * contrast_scale);
  std::cout << "    delivery carbon sink realized="
            << format_real(delivery_field_realized)
            << " none carbon sink realized="
            << format_real(none_field_realized) << "\n";
  std::cout << "  test_delivery_forces_host_chemistry: PASSED\n";
}
#endif

}  // namespace

int main() {
  std::cout << "=== GPU Uptake Limitation Parity ===\n";
  if (const int gpu_status = gutibm::test::require_gpu("gpu_uptake_limit");
      gpu_status != 0) {
    return gpu_status;
  }
#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  test_sherwood_parity();
  test_none_mode_parity_and_absence_of_limitation();
  test_delivery_forces_host_chemistry();
  std::cout << "GPU uptake limitation parity tests passed.\n";
  return 0;
#endif
}
