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
  const Real cpu_maintenance = flux_value(
      cpu, cpu_flux.maintenance_cumulative, cpu_flux.maintenance_interval);
  const Real gpu_maintenance = flux_value(
      gpu, gpu_flux.maintenance_cumulative, gpu_flux.maintenance_interval);
  const Real cpu_shortfall = flux_value(
      cpu, cpu_flux.maintenance_shortfall_cumulative,
      cpu_flux.maintenance_shortfall_interval);
  const Real gpu_shortfall = flux_value(
      gpu, gpu_flux.maintenance_shortfall_cumulative,
      gpu_flux.maintenance_shortfall_interval);
  const Real cpu_maintenance_limited = flux_value(
      cpu, cpu_flux.maintenance_limited_agents_cumulative,
      cpu_flux.maintenance_limited_agents_interval);
  const Real gpu_maintenance_limited = flux_value(
      gpu, gpu_flux.maintenance_limited_agents_cumulative,
      gpu_flux.maintenance_limited_agents_interval);

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
      report_residual("maintenance_realized", cpu_maintenance, gpu_maintenance);
  const Real shortfall_residual =
      report_residual("maintenance_shortfall_mol", cpu_shortfall, gpu_shortfall);
  const Real maintenance_limited_residual = report_residual(
      "maintenance_limited_agents", cpu_maintenance_limited,
      gpu_maintenance_limited);

  assert(mu_residual <= tolerance);
  assert(demand_residual <= tolerance);
  assert(realized_residual <= tolerance);
  assert(limited_residual <= tolerance);
  assert(maintenance_residual <= tolerance);
  assert(shortfall_residual <= tolerance);
  assert(maintenance_limited_residual <= tolerance);
  assert(cpu_limited > 0.0);
  assert(cpu_shortfall > 0.0);
  assert(cpu_maintenance_limited > 0.0);
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
  const Real cpu_maintenance = flux_value(
      cpu, cpu_flux.maintenance_cumulative, cpu_flux.maintenance_interval);
  const Real gpu_maintenance = flux_value(
      gpu, gpu_flux.maintenance_cumulative, gpu_flux.maintenance_interval);
  const Real cpu_shortfall = flux_value(
      cpu, cpu_flux.maintenance_shortfall_cumulative,
      cpu_flux.maintenance_shortfall_interval);
  const Real gpu_shortfall = flux_value(
      gpu, gpu_flux.maintenance_shortfall_cumulative,
      gpu_flux.maintenance_shortfall_interval);
  const Real cpu_maintenance_limited = flux_value(
      cpu, cpu_flux.maintenance_limited_agents_cumulative,
      cpu_flux.maintenance_limited_agents_interval);
  const Real gpu_maintenance_limited = flux_value(
      gpu, gpu_flux.maintenance_limited_agents_cumulative,
      gpu_flux.maintenance_limited_agents_interval);
  const Real residual = report_residual("total_mu_realized", total_mu(cpu),
                                        total_mu(gpu));
  assert(report_residual("maintenance_realized", cpu_maintenance,
                         gpu_maintenance) <= 5.0e-8);
  assert(report_residual("maintenance_shortfall_mol", cpu_shortfall,
                         gpu_shortfall) <= 5.0e-8);
  assert(report_residual("maintenance_limited_agents",
                         cpu_maintenance_limited,
                         gpu_maintenance_limited) <= 5.0e-8);
  std::cerr << "[gpu_diag][gpu_uptake_limit][none]"
            << " cpu_limited_agents=" << format_real(cpu_limited)
            << " gpu_limited_agents=" << format_real(gpu_limited)
            << " cpu_demand=" << format_real(cpu_demand)
            << " cpu_realized=" << format_real(cpu_realized) << "\n";
  assert(residual <= 5.0e-8);
  assert(cpu_limited == 0.0);
  assert(gpu_limited == 0.0);
  assert(cpu_shortfall == 0.0);
  assert(gpu_shortfall == 0.0);
  assert(cpu_maintenance_limited == 0.0);
  assert(gpu_maintenance_limited == 0.0);
  assert(std::abs(cpu_realized - cpu_demand) <= 1.0e-12 * cpu_demand);
  std::cout << "  test_none_mode_parity_and_absence_of_limitation: PASSED\n";
}
#endif

}  // namespace

int main() {
  std::cout << "=== GPU Uptake Limitation Parity ===\n";
  const int gpu_status = gutibm::test::require_gpu("gpu_uptake_limit");
  if (gpu_status != 0) return gpu_status;
#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  test_sherwood_parity();
  test_none_mode_parity_and_absence_of_limitation();
  std::cout << "GPU uptake limitation parity tests passed.\n";
  return 0;
#endif
}
