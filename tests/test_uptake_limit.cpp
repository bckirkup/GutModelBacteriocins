/* -----------------------------------------------------------------------
   GutIBM – Agent-side nutrient uptake limitation tests
   ----------------------------------------------------------------------- */

#include "fix_metabolism.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <format>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

using namespace gutibm;

namespace {

constexpr Real kDt = 60.0;

struct Measurement {
  Real mu_realized = 0.0;
  Real demand = 0.0;
  Real realized = 0.0;
  Real limited_agents = 0.0;
  Real field_removal = 0.0;
};

struct DeliveryStepMeasurement {
  Real demand = 0.0;
  Real funded = 0.0;
  Real ceiling = 0.0;
  Real field_removal = 0.0;
  Real min_concentration = 0.0;
  Real funded_fraction = 0.0;
  Real delivery_reduction = 0.0;
  Real delivery_retry_events = 0.0;
};

SimulationConfig base_config() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;
  cfg.initial_strains.clear();
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {20e-6, 20e-6, 20e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.chem_env.acetate.enabled = false;
  cfg.chem_env.oxygen.enabled = false;
  cfg.chem_env.siderophore.enabled = false;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 1;
  strain.mu_max = 5.0e-4;
  cfg.initial_strains.push_back(strain);
  return cfg;
}

const char* mode_name(UptakeLimitMode mode) {
  switch (mode) {
    case UptakeLimitMode::Sherwood: return "sherwood";
    case UptakeLimitMode::Voxel:    return "voxel";
    case UptakeLimitMode::Delivery: return "delivery";
    case UptakeLimitMode::None:     break;
  }
  return "none";
}

Measurement measure_single_agent(UptakeLimitMode mode, Real carbon_conc,
                                 Real radius, Real carbon_retardation) {
  SimulationConfig cfg = base_config();
  cfg.fixes.metabolism.uptake_limit = mode_name(mode);
  cfg.fixes.metabolism.uptake_limit_mode = mode;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.retardation = carbon_retardation;
      chemical.z_gradient_enabled = false;
      chemical.initial_conc = carbon_conc;
      chemical.boundary_conc = carbon_conc;
    }
  }

  Simulation sim;
  sim.init(cfg);
  auto& chem = sim.chemical_field();
  const Int carbon = chem.find(species::CARBON);
  const Int iron = chem.find(species::IRON);
  assert(carbon >= 0);
  assert(iron >= 0);
  assert(sim.agents().size() == 1);
  Agent& agent = sim.agents()[0];
  agent.radius = radius;
  agent.outer_radius = radius * 1.05;
  // Saturating carbon Monod term: demand becomes concentration-independent so
  // the sweeps below isolate the transport cap rather than the kinetics.
  agent.km.km_carbon = 1.0e-9;
  const Int cell = agent.grid_cell;
  for (Int c = 0; c < chem.global_ncells(); ++c) {
    chem.conc_global(carbon, c) = carbon_conc;
    chem.conc_global(iron, c) = 1.0e-3;
  }

  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  fix.compute(kDt);

  const auto& flux = chem.flux_accounting();
  Measurement result;
  result.mu_realized = sim.agents()[0].mu_realized;
  result.demand = flux.uptake_demand_step[static_cast<size_t>(carbon)];
  result.realized = flux.agent_uptake_step[static_cast<size_t>(carbon)];
  result.limited_agents =
      flux.uptake_limited_step[static_cast<size_t>(carbon)];
  result.field_removal =
      -chem.reac_global(carbon, cell) * sim.domain().cell_volume() * kDt;
  return result;
}

Real funded_fraction(UptakeLimitMode mode, Real carbon_conc, Real radius,
                     Real retardation) {
  const Measurement limited =
      measure_single_agent(mode, carbon_conc, radius, retardation);
  const Measurement unlimited = measure_single_agent(
      UptakeLimitMode::None, carbon_conc, radius, retardation);
  assert(unlimited.mu_realized > 0.0);
  return limited.mu_realized / unlimited.mu_realized;
}

DeliveryStepMeasurement measure_delivery_step(
    Real grid_dx, Real dt, Real concentration, Real mu_max = 1.0e-2,
    Real radius = 5.0e-7, Real boundary_concentration = -1.0,
    Real far_field_radius = 0.0) {
  SimulationConfig cfg = base_config();
  cfg.domain.grid_dx = grid_dx;
  cfg.domain.hi = {4.0 * grid_dx, 4.0 * grid_dx, 4.0 * grid_dx};
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = far_field_radius;
  cfg.fixes.metabolism.carbon_maintenance_rate = 0.0;
  cfg.initial_strains.front().mu_max = mu_max;
  cfg.time.total_time = dt;
  cfg.time.bio_dt = dt;
  cfg.time.output_interval = dt;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = concentration;
      chemical.boundary_conc = boundary_concentration >= 0.0
          ? boundary_concentration : concentration;
      chemical.diff_coeff = 5.0e-10;
      chemical.z_gradient_enabled = false;
    }
  }

  Simulation sim;
  sim.init(cfg);
  auto& agent = sim.agents()[0];
  agent.radius = radius;
  agent.outer_radius = radius * 1.05;
  agent.km.km_carbon = 1.0e-9;
  const Real initial_radius = agent.radius;
  const Int carbon = sim.chemical_field().find(species::CARBON);
  const Int cell = agent.grid_cell;
  const Real pre_step_concentration =
      sim.chemical_field().conc_global(carbon, cell);
  const Real effective_diffusivity =
      sim.chemical_field().spec(carbon).diff_coeff
      / sim.chemical_field().spec(carbon).retardation;
  sim.step(dt);

  const auto& chem = sim.chemical_field();
  const auto& flux = chem.flux_accounting();
  const auto index = static_cast<size_t>(carbon);
  DeliveryStepMeasurement result;
  result.demand = flux.uptake_demand_interval[index];
  result.funded = flux.agent_uptake_interval[index];
  result.ceiling = uptake::allowed_uptake_mol(
      to_underlying(UptakeLimitMode::Sherwood), pre_step_concentration,
      effective_diffusivity, initial_radius, sim.domain().cell_volume(), dt);
  for (Int global_cell = 0;
       global_cell < chem.global_ncells(); ++global_cell) {
    result.field_removal += chem.sink_realized_global(carbon, global_cell);
    result.min_concentration = global_cell == 0
        ? chem.conc_global(carbon, global_cell)
        : std::min(result.min_concentration,
                   chem.conc_global(carbon, global_cell));
  }
  result.funded_fraction = result.demand > 0.0
      ? result.funded / result.demand : 0.0;
  result.delivery_reduction =
      flux.delivery_reduction_interval[index]
      + flux.delivery_reduction_cumulative[index];
  result.delivery_retry_events =
      flux.delivery_retry_events_interval[index]
      + flux.delivery_retry_events_cumulative[index];
  return result;
}

Real run_resolution_probe(Real grid_dx, Real far_field_radius) {
  SimulationConfig cfg = base_config();
  cfg.domain.grid_dx = grid_dx;
  cfg.domain.hi = {120.0e-6, 120.0e-6, 120.0e-6};
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = far_field_radius;
  cfg.initial_strains.front().mu_max = 1.0e3;
  cfg.time.total_time = 120.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 120.0;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e-6;
      chemical.boundary_conc = 1.0e-6;
      chemical.diff_coeff = 5.0e-10;
      chemical.z_gradient_enabled = false;
    }
  }

  Simulation sim;
  sim.init(cfg);
  Agent& agent = sim.agents()[0];
  agent.x = {61.0e-6, 61.0e-6, 61.0e-6};
  agent.biomass = agent.mass;
  agent.radius = 5.0e-7;
  agent.outer_radius = agent.radius * 1.05;
  agent.km.km_carbon = 0.0;
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  sim.domain().pos_to_grid(agent.x, ix, iy, iz);
  agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  const Int carbon = sim.chemical_field().find(species::CARBON);
  const Real profile_width = 0.8e-6;
  for (Int cell = 0; cell < sim.chemical_field().global_ncells(); ++cell) {
    const Int cell_ix = cell % sim.domain().nx();
    const Int cell_iy = (cell / sim.domain().nx()) % sim.domain().ny();
    const Int cell_iz = cell / (sim.domain().nx() * sim.domain().ny());
    const Vec3 center = sim.domain().cell_center(cell_ix, cell_iy, cell_iz);
    const Real distance_sq = sim.domain().min_image_dist_sq(
        agent.x, center);
    sim.chemical_field().conc_global(carbon, cell) =
        1.0e-8 + (1.0e-6 - 1.0e-8)
        * (1.0 - std::exp(-distance_sq
                          / (2.0 * profile_width * profile_width)));
  }
  sim.step(60.0);
  const auto& flux = sim.chemical_field().flux_accounting();
  const auto index = static_cast<size_t>(carbon);
  return flux.agent_uptake_interval[index]
      / std::max(flux.uptake_demand_interval[index], 1.0e-30);
}

void test_far_field_resolution_invariance_and_anti_vacuity() {
  const std::vector<Real> resolutions = {2.0e-6, 4.0e-6, 6.0e-6};
  std::vector<Real> far_field;
  std::vector<Real> own_voxel;
  for (const Real resolution : resolutions) {
    far_field.push_back(run_resolution_probe(resolution, 10.0e-6));
    own_voxel.push_back(run_resolution_probe(resolution, 0.0));
  }
  const Real far_span = *std::max_element(
      far_field.begin(), far_field.end())
      - *std::min_element(far_field.begin(), far_field.end());
  const Real own_span = *std::max_element(
      own_voxel.begin(), own_voxel.end())
      - *std::min_element(own_voxel.begin(), own_voxel.end());
  const Real far_scale = std::max(
      *std::max_element(far_field.begin(), far_field.end()), 1.0e-30);
  const Real own_scale = std::max(
      *std::max_element(own_voxel.begin(), own_voxel.end()), 1.0e-30);
  std::cout << "    far-field fractions 2um/4um/6um="
            << far_field[0] << "/" << far_field[1] << "/" << far_field[2]
            << "\n";
  std::cout << "    own-voxel fractions 2um/4um/6um="
            << own_voxel[0] << "/" << own_voxel[1] << "/" << own_voxel[2]
            << "\n";
  assert(far_span / far_scale <= 0.10);
  assert(own_span / own_scale > 0.20);
  assert(own_span / own_scale > 2.0 * far_span / far_scale);
  std::cout << "  test_far_field_resolution_invariance_and_anti_vacuity:"
            << " PASSED\n";
}

void test_far_field_uniform_field_invariant() {
  const DeliveryStepMeasurement own = measure_delivery_step(
      5.0e-6, 60.0, 1.0e-4, 1.0e1, 5.0e-7, -1.0, 0.0);
  const DeliveryStepMeasurement near = measure_delivery_step(
      5.0e-6, 60.0, 1.0e-4, 1.0e1, 5.0e-7, -1.0, 5.0e-6);
  const DeliveryStepMeasurement far = measure_delivery_step(
      5.0e-6, 60.0, 1.0e-4, 1.0e1, 5.0e-7, -1.0, 10.0e-6);
  assert(std::abs(own.funded_fraction - near.funded_fraction) <= 1.0e-12);
  assert(std::abs(own.funded_fraction - far.funded_fraction) <= 1.0e-12);
  std::cout << "  test_far_field_uniform_field_invariant: PASSED\n";
}

void test_far_field_concentration_monotonicity() {
  const std::vector<Real> concentrations = {
      1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4};
  std::vector<Real> fractions;
  fractions.reserve(concentrations.size());
  for (const Real concentration : concentrations) {
    const DeliveryStepMeasurement result = measure_delivery_step(
        5.0e-6, 60.0, concentration, 1.0e1, 5.0e-7, -1.0, 10.0e-6);
    fractions.push_back(result.funded_fraction);
  }
  for (size_t i = 0; i + 1 < fractions.size(); ++i) {
    assert(fractions[i] <= fractions[i + 1]);
  }
  assert(fractions.back() > fractions.front());
  std::cout << "  test_far_field_concentration_monotonicity: PASSED\n";
}

void test_delivery_resolution_and_timestep_invariance() {
  const DeliveryStepMeasurement fine =
      measure_delivery_step(2.0e-6, 60.0, 1.0e2, 1.0e1);
  const DeliveryStepMeasurement coarse =
      measure_delivery_step(6.0e-6, 60.0, 1.0e2, 1.0e1);
  const DeliveryStepMeasurement short_step =
      measure_delivery_step(2.0e-6, 10.0, 1.0e2, 1.0e1);
  assert(fine.demand > 0.0);
  assert(coarse.demand > 0.0);
  assert(short_step.demand > 0.0);
  assert(fine.delivery_reduction == 0.0);
  assert(fine.delivery_retry_events == 0.0);
  assert(coarse.delivery_reduction == 0.0);
  assert(coarse.delivery_retry_events == 0.0);
  assert(short_step.delivery_reduction == 0.0);
  assert(short_step.delivery_retry_events == 0.0);
  assert(std::abs(fine.funded_fraction - coarse.funded_fraction) <= 0.05);
  assert(std::abs(fine.funded_fraction - short_step.funded_fraction) <= 0.05);
  std::cout << "    delivery fractions 2um/6um/10s="
            << fine.funded_fraction << "/" << coarse.funded_fraction << "/"
            << short_step.funded_fraction << "\n";
  std::cout << "  test_delivery_resolution_and_timestep_invariance: PASSED\n";
}

void test_delivery_funding_matches_demand_or_ceiling() {
  const DeliveryStepMeasurement demand_limited =
      measure_delivery_step(5.0e-6, 60.0, 1.0e2, 1.0e-6);
  // C* = demand / (4*pi*D*r*dt) ~= 0.63 mol/m^3 for the high-demand arm;
  // 0.1 mol/m^3 is below that transition and is therefore ceiling-limited.
  const DeliveryStepMeasurement ceiling_limited =
      measure_delivery_step(5.0e-6, 60.0, 1.0e-1, 1.0e1);
  assert(demand_limited.demand < demand_limited.ceiling);
  assert(std::abs(demand_limited.funded - demand_limited.demand)
         <= 1.0e-12 * std::max(demand_limited.demand, 1.0e-30));
  assert(ceiling_limited.demand > ceiling_limited.ceiling);
  assert(std::abs(ceiling_limited.funded - ceiling_limited.ceiling)
         <= 1.0e-12 * std::max(ceiling_limited.ceiling, 1.0e-30));
  std::cout << "  test_delivery_funding_matches_demand_or_ceiling: PASSED\n";
}

void test_delivery_concentration_response_and_positivity() {
  const std::vector concentrations = {
      1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4};
  std::vector<Real> fractions;
  fractions.reserve(concentrations.size());
  for (const Real concentration : concentrations) {
    const DeliveryStepMeasurement result =
        measure_delivery_step(5.0e-6, 60.0, concentration);
    assert(std::isfinite(result.funded));
    assert(std::isfinite(result.min_concentration));
    assert(result.min_concentration >= 0.0);
    fractions.push_back(result.funded_fraction);
  }
  for (size_t i = 0; i + 1 < fractions.size(); ++i) {
    assert(fractions[i] <= fractions[i + 1]);
  }
  assert(fractions.back() > fractions.front());
  std::cout << "  test_delivery_concentration_response_and_positivity: PASSED\n";
}

void test_delivery_exact_removal_and_depletion_cap() {
  const DeliveryStepMeasurement supplied =
      measure_delivery_step(5.0e-6, 60.0, 1.0e2);
  const DeliveryStepMeasurement starved =
      measure_delivery_step(5.0e-6, 60.0, 1.0e-12, 1.0e1,
                            5.0e-7, 0.0);
  assert(std::abs(supplied.field_removal - supplied.funded)
         <= 1.0e-12 * std::max(supplied.funded, 1.0e-30));
  assert(supplied.delivery_reduction == 0.0);
  assert(supplied.delivery_retry_events == 0.0);
  assert(std::abs(starved.field_removal - starved.funded)
         <= 1.0e-12 * std::max(starved.funded, 1.0e-30));
  assert(starved.funded < std::min(starved.demand, starved.ceiling));
  assert(starved.min_concentration >= 0.0);
  assert(starved.delivery_reduction > 0.0);
  assert(starved.delivery_retry_events > 0.0);
  assert(std::abs(starved.delivery_reduction
                  - (starved.ceiling - starved.funded))
         <= 1.0e-12 * std::max(starved.ceiling, 1.0e-30));
  std::cout << "  test_delivery_exact_removal_and_depletion_cap: PASSED\n";
}

std::pair<Real, Real> run_delivery_adjacent_agents(
    bool reverse_order, bool starved, Real& reduction) {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains[0].count = 2;
  cfg.initial_strains[0].mu_max = 1.0e-2;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = starved ? 1.0e-12 : 1.0e2;
      chemical.boundary_conc = starved ? 0.0 : 1.0e2;
      chemical.diff_coeff = 5.0e-10;
      chemical.z_gradient_enabled = false;
    }
  }

  Simulation sim;
  sim.init(cfg);
  const std::array<Real, 3> cell_one = {7.5e-6, 7.5e-6, 12.5e-6};
  const std::array<Real, 3> cell_two = {12.5e-6, 7.5e-6, 12.5e-6};
  for (size_t i = 0; i < sim.agents().size(); ++i) {
    Agent& agent = sim.agents()[i];
    const bool first_is_cell_two = reverse_order && i == 0;
    agent.x = first_is_cell_two || (!reverse_order && i == 1)
        ? cell_two : cell_one;
    agent.mu_max = agent.x[0] == cell_one[0] ? 1.0e-2 : 5.0e-3;
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid(agent.x, ix, iy, iz);
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  }
  sim.step(kDt);

  const auto& flux = sim.chemical_field().flux_accounting();
  const Int carbon = sim.chemical_field().find(species::CARBON);
  reduction = flux.delivery_reduction_interval[static_cast<size_t>(carbon)]
      + flux.delivery_reduction_cumulative[static_cast<size_t>(carbon)];
  std::pair<Real, Real> funding{};
  for (const Agent& agent : sim.agents()) {
    if (agent.x[0] < 10.0e-6) {
      funding.first = agent.pending_carbon_funding;
    } else {
      funding.second = agent.pending_carbon_funding;
    }
  }
  return funding;
}

void test_delivery_adjacent_agents_are_order_independent() {
  Real forward_reduction = 0.0;
  Real reverse_reduction = 0.0;
  const auto [forward_left, forward_right] = run_delivery_adjacent_agents(
      false, false, forward_reduction);
  const auto [reverse_left, reverse_right] = run_delivery_adjacent_agents(
      true, false, reverse_reduction);
  assert(forward_left > 0.0);
  assert(forward_right > 0.0);
  assert(std::abs(forward_left - reverse_left)
         <= 1.0e-12 * std::max(forward_left, 1.0e-30));
  assert(std::abs(forward_right - reverse_right)
         <= 1.0e-12 * std::max(forward_right, 1.0e-30));
  assert(forward_reduction == 0.0);
  assert(reverse_reduction == 0.0);
  const auto [starved_forward_left, starved_forward_right] =
      run_delivery_adjacent_agents(
      false, true, forward_reduction);
  const auto [starved_reverse_left, starved_reverse_right] =
      run_delivery_adjacent_agents(
      true, true, reverse_reduction);
  assert(forward_reduction > 0.0);
  assert(reverse_reduction > 0.0);
  assert(std::abs(starved_forward_left - starved_reverse_left)
         <= 1.0e-12 * std::max(starved_forward_left, 1.0e-30));
  assert(std::abs(starved_forward_right - starved_reverse_right)
         <= 1.0e-12 * std::max(starved_forward_right, 1.0e-30));
  std::cout << "  test_delivery_adjacent_agents_are_order_independent: PASSED\n";
}

bool finite_and_nonnegative(const Measurement& m) {
  return std::isfinite(m.mu_realized) && std::isfinite(m.demand)
      && std::isfinite(m.realized) && std::isfinite(m.field_removal)
      && m.demand >= 0.0 && m.realized >= 0.0 && m.limited_agents >= 0.0;
}

void test_funded_fraction_falls_with_local_concentration() {
  const std::vector<Real> concentrations = {1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7};
  std::vector<Real> fractions;
  fractions.reserve(concentrations.size());
  for (Real concentration : concentrations) {
    fractions.push_back(
        funded_fraction(UptakeLimitMode::Sherwood, concentration, 5.0e-7, 10.0));
  }
  for (size_t i = 0; i + 1 < fractions.size(); ++i) {
    std::cout << "    carbon=" << concentrations[i]
              << " funded_fraction=" << fractions[i] << "\n";
    assert(fractions[i] >= fractions[i + 1]);
  }
  std::cout << "    carbon=" << concentrations.back()
            << " funded_fraction=" << fractions.back() << "\n";
  assert(fractions.front() > fractions.back());
  assert(fractions.back() > 0.0);
  std::cout << "  test_funded_fraction_falls_with_local_concentration: PASSED\n";
}

void test_funded_fraction_falls_with_radius_and_diffusivity() {
  const Real carbon = 1.0e-7;
  const std::vector<Real> radii = {2.0e-6, 5.0e-7, 1.0e-7};
  std::vector<Real> by_radius;
  by_radius.reserve(radii.size());
  for (Real radius : radii) {
    by_radius.push_back(
        funded_fraction(UptakeLimitMode::Sherwood, carbon, radius, 10.0));
  }
  for (size_t i = 0; i + 1 < by_radius.size(); ++i) {
    std::cout << "    radius=" << radii[i]
              << " funded_fraction=" << by_radius[i] << "\n";
    assert(by_radius[i] >= by_radius[i + 1]);
  }
  assert(by_radius.front() > by_radius.back());

  const std::vector<Real> retardations = {1.0, 4.0, 16.0};
  std::vector<Real> by_diffusivity;
  by_diffusivity.reserve(retardations.size());
  for (Real retardation : retardations) {
    by_diffusivity.push_back(funded_fraction(
        UptakeLimitMode::Sherwood, carbon, 5.0e-7, retardation));
  }
  for (size_t i = 0; i + 1 < by_diffusivity.size(); ++i) {
    std::cout << "    retardation=" << retardations[i]
              << " funded_fraction=" << by_diffusivity[i] << "\n";
    assert(by_diffusivity[i] >= by_diffusivity[i + 1]);
  }
  assert(by_diffusivity.front() > by_diffusivity.back());
  std::cout << "  test_funded_fraction_falls_with_radius_and_diffusivity: PASSED\n";
}

void test_realized_uptake_respects_the_cap_and_the_field() {
  const Real carbon = 1.0e-7;
  const Real radius = 5.0e-7;
  SimulationConfig probe = base_config();
  Real effective_diffusivity = 0.0;
  for (const auto& chemical : probe.chemicals) {
    if (chemical.name == species::CARBON) {
      effective_diffusivity = chemical.diff_coeff / chemical.retardation;
    }
  }
  assert(effective_diffusivity > 0.0);
  const Real sherwood_ceiling = 4.0 * uptake::kPi * effective_diffusivity
      * radius * carbon * kDt;

  const Measurement capped =
      measure_single_agent(UptakeLimitMode::Sherwood, carbon, radius, 1.0);
  assert(finite_and_nonnegative(capped));
  assert(capped.realized <= capped.demand);
  assert(capped.realized <= sherwood_ceiling * (1.0 + 1.0e-12));
  assert(std::abs(capped.realized - capped.field_removal)
         <= 1.0e-12 * std::max(capped.realized, 1.0e-30));
  assert(capped.limited_agents == 1.0);
  assert(capped.mu_realized <= 5.0e-4);

  SimulationConfig voxel_probe = base_config();
  Real cell_volume = 0.0;
  {
    Simulation sim;
    sim.init(voxel_probe);
    cell_volume = sim.domain().cell_volume();
  }
  const Measurement voxel =
      measure_single_agent(UptakeLimitMode::Voxel, carbon, radius, 1.0);
  assert(finite_and_nonnegative(voxel));
  assert(voxel.realized <= carbon * cell_volume * (1.0 + 1.0e-12));
  assert(voxel.realized <= voxel.demand);
  std::cout << "    sherwood_realized=" << capped.realized
            << " voxel_realized=" << voxel.realized
            << " demand=" << capped.demand << "\n";
  std::cout << "  test_realized_uptake_respects_the_cap_and_the_field: PASSED\n";
}

void test_none_mode_leaves_growth_unfunded() {
  const Measurement none =
      measure_single_agent(UptakeLimitMode::None, 1.0e-9, 1.0e-7, 1.0);
  assert(none.limited_agents == 0.0);
  assert(none.realized == none.demand);
  assert(none.demand > 0.0);
  const Measurement capped =
      measure_single_agent(UptakeLimitMode::Sherwood, 1.0e-9, 1.0e-7, 1.0);
  assert(capped.realized < none.realized);
  assert(capped.mu_realized < none.mu_realized);
  std::cout << "    none_realized=" << none.realized
            << " sherwood_realized=" << capped.realized << "\n";
  std::cout << "  test_none_mode_leaves_growth_unfunded: PASSED\n";
}

struct DensityProbe {
  Real bound_fraction = 0.0;
  Real funded_ratio = 0.0;
};

// Transport-limited regime: carbon retardation is raised so the Sherwood
// ceiling sits near the demand, and the bulk concentration starts above the
// carbon Km so that depletion by the population feeds back on the cap.
DensityProbe density_probe(Int founders, Real retardation) {
  SimulationConfig cfg = base_config();
  cfg.seed = 4242;
  cfg.fixes.metabolism.uptake_limit = "sherwood";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Sherwood;
  cfg.time.total_time = 600.0;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = 600.0;
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = founders;
  strain.mu_max = 5.0e-4;
  cfg.initial_strains.push_back(strain);
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.z_gradient_enabled = false;
      chemical.retardation = retardation;
      chemical.initial_conc = 5.0e-2;
      chemical.boundary_conc = 5.0e-2;
    }
  }

  Simulation sim;
  sim.init(cfg);
  sim.run();

  const auto& flux = sim.chemical_field().flux_accounting();
  const Int carbon = sim.chemical_field().find(species::CARBON);
  const auto index = static_cast<size_t>(carbon);
  const Real limited = flux.uptake_limited_cumulative[index]
      + flux.uptake_limited_interval[index];
  const Real demanded = flux.uptake_demand_cumulative[index]
      + flux.uptake_demand_interval[index];
  const Real realized = flux.agent_uptake_cumulative[index]
      + flux.agent_uptake_interval[index];
  assert(demanded > 0.0);
  Int live = 0;
  for (const Agent& agent : sim.agents()) {
    if (agent.state != PhenoState::DEAD && !agent.flags.is_ghost) ++live;
  }
  assert(live > 0);
  const Real steps = cfg.time.total_time / cfg.time.bio_dt;
  DensityProbe probe;
  probe.bound_fraction = limited / (static_cast<Real>(live) * steps);
  probe.funded_ratio = realized / demanded;
  return probe;
}

void test_limitation_severity_rises_with_agent_density() {
  // Marginal transport regime: the ceiling sits close to the demand, so the
  // extra depletion caused by a denser population is what pushes agents over
  // the cap rather than the parameters alone.
  const Real retardation = 1.3e3;
  const DensityProbe sparse = density_probe(2, retardation);
  const DensityProbe dense = density_probe(48, retardation);
  std::cout << std::setprecision(6);
  std::cout << "    sparse_bound_fraction=" << sparse.bound_fraction
            << " funded_ratio=" << sparse.funded_ratio << "\n";
  std::cout << "    dense_bound_fraction=" << dense.bound_fraction
            << " funded_ratio=" << dense.funded_ratio << "\n";
  assert(std::isfinite(sparse.funded_ratio) && std::isfinite(dense.funded_ratio));
  assert(dense.bound_fraction > sparse.bound_fraction);
  assert(dense.funded_ratio < sparse.funded_ratio);
  assert(dense.funded_ratio > 0.0);
  assert(sparse.funded_ratio <= 1.0 + 1.0e-12);
  std::cout << "  test_limitation_severity_rises_with_agent_density: PASSED\n";
}

void test_delivery_is_positive_and_funds_only_removed_carbon() {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains[0].count = 2;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.initial_strains[0].mu_max = 5.0e-2;
  cfg.time.total_time = 120.0;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = 120.0;
  cfg.dysbiosis_threshold = 0.0;
  cfg.carbon_boundary_conc = 1.0e-7;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e-7;
      chemical.boundary_conc = 1.0e-7;
      chemical.z_gradient_enabled = false;
    }
  }
  Simulation sim;
  sim.init(cfg);
  std::vector<Real> initial_biomass_by_agent;
  initial_biomass_by_agent.reserve(sim.agents().size());
  for (Agent& agent : sim.agents()) {
    initial_biomass_by_agent.push_back(agent.biomass);
    agent.x = {7.5e-6, 7.5e-6, 12.5e-6};
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid(agent.x, ix, iy, iz);
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  }
  const Real initial_biomass_agent = sim.agents()[0].biomass;
  sim.run();
  const auto& chem = sim.chemical_field();
  const Int carbon = chem.find(species::CARBON);
  const auto& flux = chem.flux_accounting();
  const auto index = static_cast<size_t>(carbon);
  const Real removed = flux.agent_uptake_cumulative[index]
      + flux.agent_uptake_interval[index]
      + flux.maintenance_cumulative[index]
      + flux.maintenance_interval[index];
  const Real demanded = flux.uptake_demand_cumulative[index]
      + flux.uptake_demand_interval[index];
  assert(removed >= 0.0 && removed <= demanded + 1.0e-18);
  assert(flux.reaction_clip_cumulative[index]
         + flux.reaction_clip_interval[index] == 0.0);
  assert(flux.uptake_shortfall_cumulative[index]
         + flux.uptake_shortfall_interval[index] >= 0.0);
  assert(static_cast<size_t>(sim.agents().size())
         == initial_biomass_by_agent.size());
  Real committed_growth_carbon = 0.0;
  for (size_t i = 0; i < initial_biomass_by_agent.size(); ++i) {
    committed_growth_carbon +=
        (sim.agents()[i].biomass - initial_biomass_by_agent[i])
        * cfg.fixes.metabolism.yield_carbon;
  }
  assert(std::abs(committed_growth_carbon
                  + flux.maintenance_cumulative[index]
                  + flux.maintenance_interval[index]
                  - removed) <= 1.0e-18);
  assert(sim.agents()[0].biomass >= initial_biomass_agent);
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    assert(chem.conc_global(carbon, cell) >= 0.0);
  }
  std::cout << "  test_delivery_is_positive_and_funds_only_removed_carbon: PASSED\n";
}

void test_delivery_shared_voxel_funding_does_not_exceed_removal() {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains[0].count = 2;
  cfg.initial_strains[0].mu_max = 5.0e-2;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = 5.0e-6;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  cfg.dysbiosis_threshold = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e2;
      chemical.boundary_conc = 1.0e2;
      chemical.z_gradient_enabled = false;
    }
  }

  Simulation sim;
  sim.init(cfg);
  const std::array<Real, 3> shared_cell = {
      7.5e-6, 7.5e-6, 12.5e-6};
  for (Agent& agent : sim.agents()) {
    agent.x = shared_cell;
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid(agent.x, ix, iy, iz);
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  }
  sim.step(kDt);

  const auto& chem = sim.chemical_field();
  const Int carbon = chem.find(species::CARBON);
  Real funded = 0.0;
  for (const Agent& agent : sim.agents()) {
    if (agent.state != PhenoState::DEAD && !agent.flags.is_ghost) {
      funded += agent.pending_carbon_funding;
    }
  }
  Real realized = 0.0;
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    realized += chem.sink_realized_global(carbon, cell);
  }
  assert(funded > 0.0);
  assert(realized > 0.0);
  assert(funded <= realized
         + 1.0e-12 * std::max(realized, 1.0e-30));
  std::cout << "  test_delivery_shared_voxel_funding_does_not_exceed_removal:"
            << " PASSED\n";
}

void test_invalid_cell_agent_cannot_claim_delivery_funding() {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains[0].count = 2;
  cfg.initial_strains[0].mu_max = 5.0e-2;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = 5.0e-6;
  cfg.fixes.metabolism.division_threshold = 1.0e9;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e2;
      chemical.boundary_conc = 1.0e2;
      chemical.z_gradient_enabled = false;
    }
  }

  Simulation sim;
  sim.init(cfg);
  const Vec3 shared_position = {7.5e-6, 7.5e-6, 12.5e-6};
  Agent& valid = sim.agents()[0];
  Agent& invalid = sim.agents()[1];
  for (Agent& agent : sim.agents()) {
    agent.x = shared_position;
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid(agent.x, ix, iy, iz);
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  }
  invalid.grid_cell = -1;

  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  fix.compute(kDt);
  invalid.pending_growth_carbon = valid.pending_growth_carbon;
  invalid.pending_carbon_funding = valid.pending_carbon_funding;
  sim.chemical_field().apply_diffusion(sim.domain(), kDt);
  fix.post_chemistry(kDt);

  const Int carbon = sim.chemical_field().find(species::CARBON);
  const auto& flux = sim.chemical_field().flux_accounting();
  const auto index = static_cast<size_t>(carbon);
  Real realized = 0.0;
  for (Int cell = 0;
       cell < sim.chemical_field().global_ncells(); ++cell) {
    realized += sim.chemical_field().sink_realized_global(carbon, cell);
  }
  assert(invalid.pending_carbon_funding == 0.0);
  assert(flux.agent_uptake_interval[index] <= realized
         + 1.0e-12 * std::max(realized, 1.0e-30));
  assert(valid.pending_carbon_funding > 0.0);
  std::cout << "  test_invalid_cell_agent_cannot_claim_delivery_funding:"
            << " PASSED\n";
}

struct DeliveryRationingMeasurement {
  Real demand = 0.0;
  Real funded = 0.0;
  Real field_removal = 0.0;
  Real min_concentration = 0.0;
  Real rationing_factor = 1.0;
  Real delivery_reduction = 0.0;
};

DeliveryRationingMeasurement measure_delivery_rationing(Int agent_count) {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains.front().count = agent_count;
  cfg.initial_strains.front().mu_max = 1.0e-2;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = 10.0e-6;
  cfg.fixes.metabolism.division_threshold = 1.0e9;
  cfg.domain.grid_dx = 4.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e-7;
      chemical.boundary_conc = 1.0e-7;
      chemical.diff_coeff = 1.0e-12;
      chemical.z_gradient_enabled = false;
    }
  }

  Simulation sim;
  sim.init(cfg);
  const Vec3 center = {10.0e-6, 10.0e-6, 10.0e-6};
  for (Agent& agent : sim.agents()) {
    agent.x = center;
    agent.radius = 5.0e-6;
    agent.outer_radius = agent.radius * 1.05;
    agent.km.km_carbon = 1.0e-9;
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid(agent.x, ix, iy, iz);
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  }
  const Int carbon = sim.chemical_field().find(species::CARBON);
  sim.step(kDt);

  const auto& chem = sim.chemical_field();
  const auto& flux = chem.flux_accounting();
  const auto index = static_cast<size_t>(carbon);
  DeliveryRationingMeasurement result;
  result.demand = flux.uptake_demand_interval[index];
  result.funded = flux.agent_uptake_interval[index];
  result.rationing_factor =
      flux.delivery_rationing_factor_interval[index];
  result.delivery_reduction = flux.delivery_reduction_interval[index];
  result.min_concentration = std::numeric_limits<Real>::infinity();
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    result.min_concentration = std::min(
        result.min_concentration, chem.conc_global(carbon, cell));
    result.field_removal += chem.sink_realized_global(carbon, cell);
  }
  return result;
}

void test_delivery_positivity_rationing_is_sensitive_and_conservative() {
  const DeliveryRationingMeasurement sparse =
      measure_delivery_rationing(1);
  const DeliveryRationingMeasurement dense =
      measure_delivery_rationing(50);
  assert(sparse.min_concentration >= 0.0);
  assert(dense.min_concentration >= 0.0);
  assert(sparse.rationing_factor == 1.0);
  assert(sparse.delivery_reduction == 0.0);
  assert(dense.rationing_factor < 1.0);
  assert(dense.delivery_reduction > 0.0);
  assert(sparse.funded <= sparse.demand
         + 1.0e-12 * std::max(sparse.demand, 1.0e-30));
  assert(dense.funded <= dense.demand
         + 1.0e-12 * std::max(dense.demand, 1.0e-30));
  assert(std::abs(sparse.funded - sparse.field_removal)
         <= 1.0e-10 * std::max(sparse.field_removal, 1.0e-30));
  assert(std::abs(dense.funded - dense.field_removal)
         <= 1.0e-10 * std::max(dense.field_removal, 1.0e-30));
  std::cout << "    positivity sparse/dense factor="
            << sparse.rationing_factor << "/" << dense.rationing_factor
            << ", min concentration=" << sparse.min_concentration << "/"
            << dense.min_concentration << "\n";
  std::cout << "  test_delivery_positivity_rationing_is_sensitive_and_"
               "conservative: PASSED\n";
}

void test_delivery_rationing_is_local() {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains.front().count = 51;
  cfg.initial_strains.front().mu_max = 1.0e3;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = 10.0e-6;
  cfg.fixes.metabolism.division_threshold = 1.0e9;
  cfg.domain.hi = {100.0e-6, 20.0e-6, 20.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 0.0;
      chemical.boundary_conc = 0.0;
      chemical.diff_coeff = 1.0e-12;
      chemical.z_gradient_enabled = false;
      chemical.delivery_enabled = true;
      chemical.epithelial_boundary_mode = EpithelialBoundaryMode::Flux;
    }
  }

  Simulation sim;
  sim.init(cfg);
  const Vec3 starved_position = {20.0e-6, 10.0e-6, 15.0e-6};
  const Vec3 supplied_position = {80.0e-6, 10.0e-6, 15.0e-6};
  for (size_t i = 0; i < sim.agents().size(); ++i) {
    Agent& agent = sim.agents()[i];
    agent.x = i < 50 ? starved_position : supplied_position;
    agent.radius = 5.0e-6;
    agent.outer_radius = agent.radius * 1.05;
    agent.km.km_carbon = 1.0e-9;
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid(agent.x, ix, iy, iz);
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  }
  const Int carbon = sim.chemical_field().find(species::CARBON);
  for (Int cell = 0; cell < sim.domain().ncells(); ++cell) {
    sim.chemical_field().conc_global(carbon, cell) = 1.0e-7;
    const Int ix = cell % sim.domain().nx();
    if (ix >= 14) {
      sim.chemical_field().conc_global(carbon, cell) = 1.0e-4;
    }
  }

  sim.step(kDt);

  const auto& chem = sim.chemical_field();
  const auto& flux = chem.flux_accounting();
  const size_t index = static_cast<size_t>(carbon);
  Real minimum = std::numeric_limits<Real>::infinity();
  Real field_removal = 0.0;
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    minimum = std::min(minimum, chem.conc_global(carbon, cell));
    field_removal += chem.sink_realized_global(carbon, cell);
  }
  const Real starved_funding = sim.agents()[0].pending_carbon_funding;
  const Real supplied_funding = sim.agents()[50].pending_carbon_funding;
  const Real total_funding = flux.agent_uptake_interval[index];
  assert(flux.delivery_retry_events_interval[index] > 0.0);
  assert(flux.delivery_rationing_factor_interval[index] < 1.0);
  assert(minimum >= 0.0);
  assert(supplied_funding > starved_funding * 1.5);
  assert(total_funding <= flux.uptake_demand_interval[index]
         + 1.0e-12 * std::max(
             flux.uptake_demand_interval[index], 1.0e-30));
  assert(std::abs(field_removal - total_funding)
         <= 1.0e-10 * std::max(field_removal, 1.0e-30));
  std::cout << "    local funding starved/supplied="
            << starved_funding << "/" << supplied_funding
            << ", minimum concentration=" << minimum << "\n";
  std::cout << "  test_delivery_rationing_is_local: PASSED\n";
}

struct RegularizedDeliveryMeasurement {
  Real prescribed = 0.0;
  Real realized = 0.0;
  Int nonzero_cells = 0;
  bool wraps_x = false;
};

RegularizedDeliveryMeasurement measure_regularized_delivery(
    const Vec3& position) {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = 6.0e-6;
  cfg.initial_strains.front().mu_max = 1.0e1;
  cfg.fixes.metabolism.division_threshold = 1.0e9;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e2;
      chemical.boundary_conc = 1.0e2;
      chemical.z_gradient_enabled = false;
    }
  }

  Simulation sim;
  sim.init(cfg);
  Agent& agent = sim.agents()[0];
  agent.x = position;
  agent.radius = 5.0e-7;
  agent.outer_radius = agent.radius * 1.05;
  agent.km.km_carbon = 1.0e-9;
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  sim.domain().pos_to_grid(agent.x, ix, iy, iz);
  agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  const Int carbon = sim.chemical_field().find(species::CARBON);
  sim.step(kDt);

  RegularizedDeliveryMeasurement result;
  const auto& chem = sim.chemical_field();
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    const Real prescribed = chem.prescribed_sink_global(carbon, cell);
    const Real realized = chem.sink_realized_global(carbon, cell);
    result.prescribed += prescribed;
    result.realized += realized;
    result.nonzero_cells += prescribed > 0.0 ? 1 : 0;
    const Int cell_ix = cell % sim.domain().nx();
    if (prescribed > 0.0
        && (cell_ix == 0 || cell_ix == sim.domain().nx() - 1)) {
      result.wraps_x = true;
    }
  }
  return result;
}

void test_regularized_delivery_mass_is_conservative() {
  const RegularizedDeliveryMeasurement lower =
      measure_regularized_delivery({0.1e-6, 10.0e-6, 0.1e-6});
  const RegularizedDeliveryMeasurement upper =
      measure_regularized_delivery({0.1e-6, 10.0e-6, 19.9e-6});
  assert(lower.prescribed > 0.0);
  assert(upper.prescribed > 0.0);
  assert(lower.nonzero_cells > 1);
  assert(upper.nonzero_cells > 1);
  assert(lower.realized > 0.0);
  assert(upper.realized > 0.0);
  assert(lower.prescribed <= upper.prescribed
         + 1.0e-12 * std::max(upper.prescribed, 1.0e-30));
  assert(lower.realized <= lower.prescribed
         + 1.0e-12 * std::max(lower.prescribed, 1.0e-30));
  assert(upper.realized <= upper.prescribed
         + 1.0e-12 * std::max(upper.prescribed, 1.0e-30));
  assert(std::abs(lower.realized - lower.prescribed)
         <= 1.0e-12 * std::max(lower.prescribed, 1.0e-30));
  assert(std::abs(upper.realized - upper.prescribed)
         <= 1.0e-12 * std::max(upper.prescribed, 1.0e-30));
  assert(lower.wraps_x);
  std::cout << "  test_regularized_delivery_mass_is_conservative: PASSED\n";
}

std::pair<Real, Real> run_regularized_retry_case(Real support_radius) {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = support_radius;
  cfg.initial_strains.front().count = 1;
  cfg.initial_strains.front().mu_max = 1.0e1;
  cfg.domain.grid_dx = 2.0e-6;
  cfg.domain.hash_cell_size = 4.0e-6;
  cfg.fixes.metabolism.division_threshold = 1.0e9;
  constexpr Int kRetrySteps = 10;
  cfg.time.total_time = kRetrySteps * kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = cfg.time.total_time;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e-9;
      chemical.boundary_conc = 1.0e-9;
      chemical.diff_coeff = 1.0e-12;
      chemical.z_gradient_enabled = false;
    }
  }
  Simulation sim;
  sim.init(cfg);
  for (Agent& agent : sim.agents()) {
    agent.radius = 2.0e-6;
    agent.km.km_carbon = 1.0e-9;
    agent.x = {7.5e-6, 7.5e-6, 7.5e-6};
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid(agent.x, ix, iy, iz);
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  }
  for (Int step = 0; step < kRetrySteps; ++step) {
    sim.agents()[0].radius = 2.0e-6;
    sim.step(kDt);
    if (step + 1 < kRetrySteps) {
      const Int carbon = sim.chemical_field().find(species::CARBON);
      for (Int cell = 0;
           cell < sim.chemical_field().global_ncells(); ++cell) {
        sim.chemical_field().conc_global(carbon, cell) = 1.0e-9;
      }
      for (Agent& agent : sim.agents()) {
        agent.mu_realized = agent.mu_max;
      }
    }
  }
  const Int carbon = sim.chemical_field().find(species::CARBON);
  const auto index = static_cast<size_t>(carbon);
  const auto& flux = sim.chemical_field().flux_accounting();
  return {
      flux.delivery_reduction_interval[index]
          + flux.delivery_reduction_cumulative[index],
      flux.delivery_retry_events_interval[index]
          + flux.delivery_retry_events_cumulative[index]};
}

void test_regularized_delivery_reduces_depletion_retries() {
  const auto radius_zero = run_regularized_retry_case(0.0);
  const auto radius_ten = run_regularized_retry_case(10.0e-6);
  std::cout << "    retry reductions radius 0/10um="
            << radius_zero.first << "/" << radius_ten.first
            << ", events=" << radius_zero.second << "/" << radius_ten.second
            << "\n";
  assert(radius_zero.first > 0.0);
  assert(radius_zero.second > 0.0);
  assert(radius_ten.first >= 0.0);
  assert(radius_ten.second >= 0.0);
  assert(radius_zero.first >= 10.0 * std::max(radius_ten.first, 1.0e-30));
  assert(radius_zero.second >= 10.0 * std::max(radius_ten.second, 1.0));
  std::cout << "  test_regularized_delivery_reduces_depletion_retries: PASSED\n";
}

struct ResolutionFundingMeasurement {
  Real funded = 0.0;
  Real demanded = 0.0;
  Real delivery_reduction = 0.0;
};

ResolutionFundingMeasurement measure_regularized_resolution_funding(
    Real grid_dx, Real far_field_radius) {
  SimulationConfig cfg = base_config();
  cfg.domain.grid_dx = grid_dx;
  cfg.domain.hi = {120.0e-6, 120.0e-6, 120.0e-6};
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = far_field_radius;
  cfg.initial_strains.front().mu_max = 1.0e-4;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e-4;
      chemical.boundary_conc = 1.0e-4;
      chemical.diff_coeff = 1.0e-12;
      chemical.z_gradient_enabled = false;
    }
  }

  Simulation sim;
  sim.init(cfg);
  Agent& agent = sim.agents()[0];
  agent.x = {61.0e-6, 61.0e-6, 61.0e-6};
  agent.radius = 5.0e-6;
  agent.outer_radius = agent.radius * 1.05;
  agent.km.km_carbon = 1.0e-9;
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  sim.domain().pos_to_grid(agent.x, ix, iy, iz);
  agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  const Int carbon = sim.chemical_field().find(species::CARBON);
  for (Int cell = 0;
       cell < sim.chemical_field().global_ncells(); ++cell) {
    sim.chemical_field().conc_global(carbon, cell) = 1.0e-4;
  }
  sim.step(kDt);
  const auto& flux = sim.chemical_field().flux_accounting();
  const size_t index = static_cast<size_t>(carbon);
  return {flux.agent_uptake_interval[index],
          flux.uptake_demand_interval[index],
          flux.delivery_reduction_interval[index]
              + flux.delivery_reduction_cumulative[index]};
}

void test_regularized_delivery_funding_resolution_invariance() {
  const std::vector<Real> resolutions = {
      2.0e-6, 4.0e-6, 6.0e-6};
  std::vector<Real> regularized;
  std::vector<Real> own_voxel;
  std::cout << std::setprecision(12);
  for (const Real resolution : resolutions) {
    const ResolutionFundingMeasurement far_field =
        measure_regularized_resolution_funding(resolution, 10.0e-6);
    const ResolutionFundingMeasurement own =
        measure_regularized_resolution_funding(resolution, 0.0);
    assert(far_field.demanded > 0.0);
    assert(own.demanded > 0.0);
    const Real own_fraction = own.funded / own.demanded;
    assert(own_fraction > 0.0);
    assert(own_fraction < 1.0);
    assert(own.delivery_reduction > 0.0);
    std::cout << "    resolution " << resolution
              << " radius 10um funded/demand="
              << far_field.funded / far_field.demanded
              << " reduction=" << far_field.delivery_reduction
              << ", radius 0 funded/demand="
              << own_fraction
              << " reduction=" << own.delivery_reduction << "\n";
    regularized.push_back(far_field.funded / far_field.demanded);
    own_voxel.push_back(own_fraction);
  }
  const Real regularized_span = *std::max_element(
      regularized.begin(), regularized.end())
      - *std::min_element(regularized.begin(), regularized.end());
  const Real own_span = *std::max_element(
      own_voxel.begin(), own_voxel.end())
      - *std::min_element(own_voxel.begin(), own_voxel.end());
  const Real regularized_scale = std::max(
      *std::max_element(regularized.begin(), regularized.end()),
      1.0e-30);
  const Real own_scale = std::max(
      *std::max_element(own_voxel.begin(), own_voxel.end()), 1.0e-30);
  std::cout << "    realized funded fractions radius 10um 2/4/6um="
            << regularized[0] << "/" << regularized[1] << "/"
            << regularized[2] << "\n";
  std::cout << "    realized funded fractions radius 0 2/4/6um="
            << own_voxel[0] << "/" << own_voxel[1] << "/"
            << own_voxel[2] << "\n";
  assert(regularized_span / regularized_scale <= 0.05);
  assert(own_span / own_scale > 0.20);
  std::cout << "  test_regularized_delivery_funding_resolution_invariance:"
            << " PASSED\n";
}

struct PopulationResolutionMeasurement {
  Real funded = 0.0;
  Real demanded = 0.0;
  Real initial_biomass = 0.0;
  Real biomass = 0.0;
  Int divisions = 0;
  Int final_agents = 0;
};

PopulationResolutionMeasurement measure_population_resolution(
    Real grid_dx, Real far_field_radius) {
  SimulationConfig cfg = base_config();
  cfg.seed = 1001;
  cfg.enabled_fixes = {"metabolism"};
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {180.0e-6, 180.0e-6, 108.0e-6};
  cfg.domain.grid_dx = grid_dx;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.initial_population.placement = "anatomic";
  cfg.initial_population.anatomic_exclusion_floor = 20.0e-6;
  cfg.initial_population.anatomic_outer_extent = 100.0e-6;
  cfg.initial_strains.front().count = 20;
  cfg.initial_strains.front().mu_max = 5.0e-4;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = far_field_radius;
  cfg.fixes.metabolism.division_threshold = 2.0;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.dysbiosis_threshold = 0.0;
  cfg.advection.radial_turnover = 1.0e100;
  cfg.advection.distal_transit_time = 1.0e100;
  constexpr Int kPopulationSteps = 360;
  cfg.time.total_time = kPopulationSteps * kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = cfg.time.total_time;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 5.0e-3;
      chemical.boundary_conc = 5.0e-3;
      chemical.diff_coeff = 1.0e-12;
      chemical.z_gradient_enabled = false;
    }
  }

  PopulationResolutionMeasurement result;
  Simulation sim;
  sim.init(cfg);
  for (const Agent& agent : sim.agents()) {
    if (agent.state != PhenoState::DEAD && !agent.flags.is_ghost) {
      result.initial_biomass += agent.biomass;
    }
  }
  sim.run();
  const Int carbon = sim.chemical_field().find(species::CARBON);
  const auto index = static_cast<size_t>(carbon);
  const auto& flux = sim.chemical_field().flux_accounting();
  result.funded = flux.agent_uptake_cumulative[index]
      + flux.agent_uptake_interval[index];
  result.demanded = flux.uptake_demand_cumulative[index]
      + flux.uptake_demand_interval[index];
  for (const Agent& agent : sim.agents()) {
    if (agent.state != PhenoState::DEAD && !agent.flags.is_ghost) {
      result.biomass += agent.biomass;
    }
  }
  result.divisions = sim.cumulative_events().divisions;
  result.final_agents = sim.agents().size();
  return result;
}

Real relative_resolution_difference(Real lower, Real upper) {
  return std::abs(upper - lower)
      / std::max(0.5 * (std::abs(lower) + std::abs(upper)), 1.0e-30);
}

[[maybe_unused]] void test_population_scale_delivery_resolution_invariance() {
  const PopulationResolutionMeasurement regularized_2um =
      measure_population_resolution(2.0e-6, 10.0e-6);
  const PopulationResolutionMeasurement regularized_6um =
      measure_population_resolution(6.0e-6, 10.0e-6);
  const PopulationResolutionMeasurement voxel_2um =
      measure_population_resolution(2.0e-6, 0.0);
  const PopulationResolutionMeasurement voxel_6um =
      measure_population_resolution(6.0e-6, 0.0);

  const Real regularized_funded_difference =
      relative_resolution_difference(
          regularized_2um.funded, regularized_6um.funded);
  const Real regularized_demand_difference =
      relative_resolution_difference(
          regularized_2um.demanded, regularized_6um.demanded);
  const Real regularized_biomass_difference =
      relative_resolution_difference(
          regularized_2um.biomass, regularized_6um.biomass);
  const Real voxel_funded_difference =
      relative_resolution_difference(voxel_2um.funded, voxel_6um.funded);
  const Real voxel_demand_difference =
      relative_resolution_difference(voxel_2um.demanded, voxel_6um.demanded);
  const Real voxel_biomass_difference =
      relative_resolution_difference(voxel_2um.biomass, voxel_6um.biomass);

  assert(regularized_2um.funded > 0.0);
  assert(regularized_2um.demanded > regularized_2um.funded);
  assert(regularized_2um.biomass > regularized_2um.initial_biomass);
  assert(regularized_6um.biomass > regularized_6um.initial_biomass);
  assert(regularized_funded_difference <= 0.015);
  assert(regularized_demand_difference <= 0.015);
  assert(regularized_biomass_difference <= 0.012);

  assert(voxel_2um.funded > 0.0);
  assert(voxel_2um.demanded > voxel_2um.funded);
  assert(voxel_2um.biomass > voxel_2um.initial_biomass);
  assert(voxel_6um.biomass > voxel_6um.initial_biomass);
  assert(voxel_funded_difference >= 0.10);
  assert(voxel_demand_difference >= 0.10);
  assert(voxel_biomass_difference >= 0.10);
  assert(voxel_funded_difference > 3.0 * regularized_funded_difference);
  assert(voxel_demand_difference > 3.0 * regularized_demand_difference);
  assert(voxel_biomass_difference > 3.0 * regularized_biomass_difference);

  std::cout << std::format(
      "    population resolution measurements:\n"
      "      radius 10um, 2um grid: funded={:.12g} demanded={:.12g} "
      "initial_biomass={:.12g} biomass={:.12g} divisions={}\n"
      "      radius 10um, 6um grid: funded={:.12g} demanded={:.12g} "
      "initial_biomass={:.12g} biomass={:.12g} divisions={}\n"
      "      radius 0, 2um grid: funded={:.12g} demanded={:.12g} "
      "initial_biomass={:.12g} biomass={:.12g} divisions={}\n"
      "      radius 0, 6um grid: funded={:.12g} demanded={:.12g} "
      "initial_biomass={:.12g} biomass={:.12g} divisions={}\n"
      "      relative spreads funded/demanded/biomass: regularized={:.12g}/"
      "{:.12g}/{:.12g}, radius 0={:.12g}/{:.12g}/{:.12g}\n",
      regularized_2um.funded, regularized_2um.demanded,
      regularized_2um.initial_biomass, regularized_2um.biomass,
      regularized_2um.divisions, regularized_6um.funded,
      regularized_6um.demanded, regularized_6um.initial_biomass,
      regularized_6um.biomass, regularized_6um.divisions, voxel_2um.funded,
      voxel_2um.demanded, voxel_2um.initial_biomass, voxel_2um.biomass,
      voxel_2um.divisions, voxel_6um.funded, voxel_6um.demanded,
      voxel_6um.initial_biomass, voxel_6um.biomass, voxel_6um.divisions,
      regularized_funded_difference, regularized_demand_difference,
      regularized_biomass_difference, voxel_funded_difference,
      voxel_demand_difference, voxel_biomass_difference);
  std::cout << "  test_population_scale_delivery_resolution_invariance:"
            << " PASSED\n";
}

std::pair<Real, Real> run_delivery_maintenance_case(Real maintenance_rate) {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains[0].count = 2;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  // Keep growth demand large enough that maintenance competes for the
  // deliberately supply-limited carbon fixture.
  cfg.initial_strains[0].mu_max = 5.0e-3;
  cfg.fixes.metabolism.carbon_maintenance_rate = maintenance_rate;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.dysbiosis_threshold = 0.0;
  cfg.carbon_boundary_conc = 1.0e-6;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e-6;
      chemical.boundary_conc = 1.0e-6;
      chemical.diff_coeff = 1.0e-12;
      chemical.z_gradient_enabled = false;
    }
  }
  Simulation sim;
  sim.init(cfg);
  for (Agent& agent : sim.agents()) {
    agent.x = {7.5e-6, 7.5e-6, 12.5e-6};
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid(agent.x, ix, iy, iz);
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  }
  const Real initial_biomass = sim.agents()[0].biomass;
  sim.run();
  const auto& flux = sim.chemical_field().flux_accounting();
  const Int carbon = sim.chemical_field().find(species::CARBON);
  const auto index = static_cast<size_t>(carbon);
  const Real growth = sim.agents()[0].biomass - initial_biomass;
  return {flux.maintenance_interval[index], growth};
}

void test_delivery_maintenance_reduces_growth() {
  const auto [low_maintenance, low_growth] =
      run_delivery_maintenance_case(1.0e-9);
  const auto [high_maintenance, high_growth] =
      run_delivery_maintenance_case(2.0e-9);
  assert(high_maintenance > low_maintenance);
  assert(high_growth < low_growth);
  std::cout << "  test_delivery_maintenance_reduces_growth: PASSED\n";
}

Real delivery_density_funded_fraction(Int founders) {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains[0].count = founders;
  cfg.initial_strains[0].mu_max = 2.0e-2;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.carbon_maintenance_rate = 0.0;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  cfg.dysbiosis_threshold = 0.0;
  cfg.carbon_boundary_conc = 1.0e-6;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e-6;
      chemical.boundary_conc = 1.0e-6;
      chemical.diff_coeff = 1.0e-12;
      chemical.z_gradient_enabled = false;
    }
  }
  Simulation sim;
  sim.init(cfg);
  for (Agent& agent : sim.agents()) {
    agent.x = {7.5e-6, 7.5e-6, 12.5e-6};
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid(agent.x, ix, iy, iz);
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  }
  sim.run();
  const auto& flux = sim.chemical_field().flux_accounting();
  const Int carbon = sim.chemical_field().find(species::CARBON);
  const auto index = static_cast<size_t>(carbon);
  const Real demand = flux.uptake_demand_interval[index];
  const Real realized = flux.agent_uptake_interval[index];
  assert(demand > 0.0);
  return realized / demand;
}

void test_delivery_density_brake() {
  const Real sparse = delivery_density_funded_fraction(2);
  const Real medium = delivery_density_funded_fraction(8);
  const Real dense = delivery_density_funded_fraction(32);
  std::cout << "    delivery_density_funded_fraction=" << sparse << ", "
            << medium << ", " << dense << "\n";
  assert(std::isfinite(sparse) && std::isfinite(medium)
         && std::isfinite(dense));
  assert(sparse > medium);
  assert(medium > dense);
  assert(dense >= 0.0);
  std::cout << "  test_delivery_density_brake: PASSED\n";
}

void test_delivery_queues_noncarbon_chemistry_once() {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.division_threshold = 1.0e9;
  cfg.chem_env.acetate.enabled = true;
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.oxygen.metabolic_switch_enabled = true;
  cfg.chem_env.oxygen.tau_metabolic_switch = 1.0;
  cfg.chem_env.oxygen.ferm_acid_yield = 1.0;
  cfg.chem_env.oxygen.aerobic_carbon_cost_factor = 1.0;
  cfg.chem_env.oxygen.anaerobic_carbon_cost_factor = 1.0;
  cfg.chem_env.acetate.vbf_production = 0.0;
  cfg.chem_env.acetate.vbf_consumption = 0.0;
  cfg.chem_env.acetate.epithelial_uptake = 0.0;
  cfg.chem_env.acetate.overflow_threshold = 1.0e9;
  cfg.chem_env.acetate.scavenge_rate = 0.0;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.vbf.nutrient_sink = 0.0;
  cfg.initial_strains[0].mu_max = 5.0e-4;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON
        || chemical.name == species::IRON
        || chemical.name == species::B12) {
      chemical.z_gradient_enabled = false;
      chemical.initial_conc = 1.0e-6;
      chemical.boundary_conc = 1.0e-6;
    } else if (chemical.name == species::OXYGEN
               || chemical.name == species::ACETATE) {
      chemical.initial_conc = 0.0;
      chemical.boundary_conc = 0.0;
    }
  }
  Simulation sim;
  sim.init(cfg);
  const Int cell = sim.agents()[0].grid_cell;
  const Real initial_biomass = sim.agents()[0].biomass;
  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  fix.compute(kDt);
  sim.chemical_field().apply_diffusion(sim.domain(), kDt);
  fix.post_chemistry(kDt);
  const Real funded_biomass = sim.agents()[0].biomass - initial_biomass;
  assert(funded_biomass > 0.0);
  sim.chemical_field().zero_reactions();
  fix.compute(kDt);
  assert(funded_biomass > 0.0);
  const auto& chem = sim.chemical_field();
  const Int iron = chem.find(species::IRON);
  const Int acetate = chem.find(species::ACETATE);
  const Real volume = sim.domain().cell_volume();
  const Real iron_draw = -chem.reac_global(iron, cell) * volume * kDt;
  const Real acetate_production =
      chem.reac_global(acetate, cell) * volume * kDt;
  const Real expected_acid = funded_biomass
      * sim.agents()[0].realized_fermentation_fraction
      * cfg.fixes.metabolism.yield_carbon
      * cfg.chem_env.oxygen.ferm_acid_yield;
  assert(std::abs(iron_draw
                  - funded_biomass * cfg.fixes.metabolism.yield_iron)
         <= 1.0e-18);
  assert(std::abs(acetate_production - expected_acid) <= 1.0e-18);
  std::cout << "  test_delivery_queues_noncarbon_chemistry_once: PASSED\n";
}

void test_delivery_preserves_negative_growth() {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.maintenance_rate = 2.0e-5;
  cfg.initial_strains[0].mu_max = 1.0e-5;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON
        || chemical.name == species::IRON) {
      chemical.initial_conc = 1.0;
      chemical.boundary_conc = 1.0;
    }
  }
  Simulation sim;
  sim.init(cfg);
  const Real initial_biomass = sim.agents()[0].biomass;
  sim.step(kDt);
  assert(sim.agents()[0].biomass < initial_biomass);
  assert(sim.agents()[0].mu_realized < 0.0);
  std::cout << "  test_delivery_preserves_negative_growth: PASSED\n";
}

void test_delivery_negative_growth_books_maintenance() {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.maintenance_rate = 2.0e-5;
  cfg.fixes.metabolism.carbon_maintenance_rate = 1.0e-3;
  cfg.initial_strains[0].mu_max = 1.0e-5;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON
        || chemical.name == species::IRON) {
      chemical.z_gradient_enabled = false;
      chemical.initial_conc = 1.0e-7;
      chemical.boundary_conc = 1.0e-7;
    }
  }
  Simulation sim;
  sim.init(cfg);
  const Int carbon = sim.chemical_field().find(species::CARBON);
  const Real initial_biomass = sim.agents()[0].biomass;
  sim.step(kDt);

  const auto& chem = sim.chemical_field();
  const auto& flux = chem.flux_accounting();
  const auto index = static_cast<size_t>(carbon);
  Real realized = 0.0;
  for (Int global_cell = 0;
       global_cell < chem.global_ncells(); ++global_cell) {
    realized += chem.sink_realized_global(carbon, global_cell);
  }
  const Real maintenance = flux.maintenance_interval[index];
  const Real growth = flux.agent_uptake_interval[index];
  const Real maintenance_shortfall =
      flux.maintenance_shortfall_interval[index];
  assert(sim.agents()[0].biomass < initial_biomass);
  assert(sim.agents()[0].mu_realized < 0.0);
  assert(maintenance > 0.0);
  assert(maintenance_shortfall > 0.0);
  assert(std::abs(maintenance + growth - realized)
         <= 1.0e-12 * std::max(realized, 1.0e-30));
  std::cout << "  test_delivery_negative_growth_books_maintenance: PASSED\n";
}

void test_delivery_maintenance_limited_runs_to_horizon() {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains[0].count = 2;
  cfg.initial_strains[0].mu_max = 5.0e-4;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.carbon_maintenance_rate = 1.0e-3;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.time.total_time = 600.0;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = 600.0;
  cfg.closure.zero_realization_grace_steps = 1;
  cfg.dysbiosis_threshold = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON
        || chemical.name == species::IRON) {
      chemical.z_gradient_enabled = false;
      chemical.initial_conc = 1.0e-7;
      chemical.boundary_conc = 1.0e-7;
      chemical.diff_coeff = 1.0e-12;
    }
  }

  Simulation sim;
  sim.init(cfg);
  sim.run();

  const auto& flux = sim.chemical_field().flux_accounting();
  const Int carbon = sim.chemical_field().find(species::CARBON);
  const auto index = static_cast<size_t>(carbon);
  const Real maintenance = flux.maintenance_interval[index]
      + flux.maintenance_cumulative[index];
  const Real growth = flux.agent_uptake_interval[index]
      + flux.agent_uptake_cumulative[index];
  const Real shortfall = flux.maintenance_shortfall_interval[index]
      + flux.maintenance_shortfall_cumulative[index];
  assert(sim.termination_cause() == TerminationCause::HorizonReached);
  assert(sim.step_count() == 10);
  assert(maintenance > 0.0);
  assert(growth == 0.0);
  assert(shortfall > 0.0);
  std::cout << "    maintenance=" << maintenance
            << " growth=" << growth
            << " shortfall=" << shortfall << "\n";
  std::cout << "  test_delivery_maintenance_limited_runs_to_horizon: PASSED\n";
}

struct GradientDeliveryMeasurement {
  Real demand = 0.0;
  Real realized = 0.0;
  Real funded_fraction = 0.0;
  Real total_concentration = 0.0;
  Real background_concentration = 0.0;
};

struct GradientRepresentationMeasurement {
  Real demand = 0.0;
  Real funded = 0.0;
  Real pending_ceiling = 0.0;
  Real concentration = 0.0;
  Real analytic_ceiling = 0.0;
  Real effective_diffusivity = 0.0;
  Real radius = 0.0;
};

GradientRepresentationMeasurement run_gradient_representation(
    bool gradient_enabled) {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = 0.0;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.initial_strains.front().mu_max = 1.0e-2;
  cfg.fixes.metabolism.division_threshold = 1.0e9;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  constexpr Real initial_concentration = 1.0e-3;
  constexpr Real gradient_lambda = 25.0e-6;
  constexpr Real agent_z = 12.5e-6;
  const Real physical_concentration = initial_concentration
      * std::exp(-agent_z / gradient_lambda);
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = gradient_enabled
          ? initial_concentration : physical_concentration;
      chemical.boundary_conc = chemical.initial_conc;
      chemical.z_gradient_enabled = gradient_enabled;
      chemical.z_gradient_lambda = gradient_lambda;
    } else if (chemical.name == species::IRON) {
      chemical.initial_conc = initial_concentration;
      chemical.boundary_conc = initial_concentration;
      chemical.z_gradient_enabled = false;
    }
  }

  Simulation sim;
  sim.init(cfg);
  Agent& agent = sim.agents()[0];
  agent.x = {7.5e-6, 7.5e-6, agent_z};
  agent.radius = 5.0e-7;
  agent.outer_radius = agent.radius * 1.05;
  agent.km.km_carbon = 1.0e-9;
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  sim.domain().pos_to_grid(agent.x, ix, iy, iz);
  agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  const Int carbon = sim.chemical_field().find(species::CARBON);
  const Real concentration =
      sim.chemical_field().conc_global(carbon, agent.grid_cell);
  const Real effective_diffusivity = uptake::effective_diffusivity(
      sim.chemical_field().spec(carbon).diff_coeff,
      sim.chemical_field().spec(carbon).retardation);
  const Real analytic_ceiling = uptake::allowed_uptake_mol(
      to_underlying(UptakeLimitMode::Sherwood), concentration,
      effective_diffusivity, agent.radius, sim.domain().cell_volume(), kDt);
  sim.step(kDt);

  const auto& flux = sim.chemical_field().flux_accounting();
  const auto index = static_cast<size_t>(carbon);
  GradientRepresentationMeasurement result;
  result.demand = flux.uptake_demand_interval[index];
  result.funded = flux.agent_uptake_interval[index];
  result.pending_ceiling = sim.agents()[0].pending_carbon_ceiling;
  result.concentration = concentration;
  result.analytic_ceiling = analytic_ceiling;
  result.effective_diffusivity = effective_diffusivity;
  result.radius = agent.radius;
  return result;
}

void test_delivery_gradient_representation_invariance() {
  const GradientRepresentationMeasurement gradient =
      run_gradient_representation(true);
  const GradientRepresentationMeasurement flat =
      run_gradient_representation(false);
  const Real relative_tolerance = 1.0e-10;
  std::cerr << std::setprecision(17)
            << "    gradient detail funded/demand/pending/concentration="
            << gradient.funded << "/" << gradient.demand << "/"
            << gradient.pending_ceiling << "/" << gradient.concentration
            << ", flat=" << flat.funded << "/" << flat.demand << "/"
            << flat.pending_ceiling << "/" << flat.concentration << "\n";
  assert(gradient.demand > gradient.funded);
  assert(flat.demand > flat.funded);
  assert(std::abs(gradient.concentration - flat.concentration)
         <= relative_tolerance * flat.concentration);
  assert(std::abs(gradient.funded - flat.funded)
         <= relative_tolerance * flat.funded);
  assert(std::abs(gradient.pending_ceiling - flat.pending_ceiling)
         <= relative_tolerance * flat.pending_ceiling);

  assert(std::abs(gradient.pending_ceiling - gradient.analytic_ceiling)
         <= relative_tolerance * gradient.analytic_ceiling);
  const Real double_concentration_ceiling = uptake::allowed_uptake_mol(
      to_underlying(UptakeLimitMode::Sherwood),
      2.0 * gradient.concentration,
      gradient.effective_diffusivity, gradient.radius,
      0.0, kDt);
  // CHANGE DETECTOR: before this fix, the gradient representation priced
  // delivery at approximately twice the flat representation's concentration.
  assert(std::abs(gradient.pending_ceiling - double_concentration_ceiling)
         > relative_tolerance * double_concentration_ceiling);
  std::cout << "    gradient/flat funded=" << gradient.funded << "/"
            << flat.funded << " ratio=" << gradient.funded / flat.funded
            << ", pending ceiling=" << gradient.pending_ceiling << "\n";
  std::cout << "  test_delivery_gradient_representation_invariance: PASSED\n";
}

struct PopulationGradientDeliveryMeasurement {
  Real demand = 0.0;
  Real funded = 0.0;
  Real field_removal = 0.0;
  Real minimum_concentration = 0.0;
};

PopulationGradientDeliveryMeasurement measure_population_gradient_delivery() {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains.front().count = 128;
  cfg.initial_strains.front().mu_max = 1.0e-2;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.division_threshold = 1.0e9;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.domain.hi = {100.0e-6, 100.0e-6, 100.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.time.total_time = 3.0 * kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = cfg.time.total_time;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.advection.radial_turnover = 1.0e100;
  cfg.advection.distal_transit_time = 1.0e100;
  cfg.dysbiosis_threshold = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::IRON) {
      chemical.initial_conc = 1.0;
      chemical.boundary_conc = 1.0;
      chemical.z_gradient_enabled = false;
    }
  }

  Simulation sim;
  sim.init(cfg);
  const Int agent_count = static_cast<Int>(sim.agents().size());
  for (Int i = 0; i < agent_count; ++i) {
    Agent& agent = sim.agents()[static_cast<size_t>(i)];
    const Int ix = static_cast<Int>(i % 8) * 2 + 1;
    const Int iy = static_cast<Int>((i / 8) % 8) * 2 + 1;
    const Int iz = static_cast<Int>(i / 64) + 1;
    agent.x = sim.domain().cell_center(ix, iy, iz);
    agent.radius = 5.0e-7;
    agent.outer_radius = agent.radius * 1.05;
    agent.km.km_carbon = 1.0e-9;
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  }
  const auto& chem = sim.chemical_field();
  const Int carbon = chem.find(species::CARBON);
  const auto index = static_cast<size_t>(carbon);
  PopulationGradientDeliveryMeasurement result;
  constexpr Int population_steps = 3;
  for (Int step = 0; step < population_steps; ++step) {
    sim.step(kDt);
    for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
      result.field_removal += chem.sink_realized_global(carbon, cell);
    }
  }
  const auto& flux = chem.flux_accounting();
  result.demand = flux.uptake_demand_cumulative[index]
      + flux.uptake_demand_interval[index];
  result.funded = flux.agent_uptake_cumulative[index]
      + flux.agent_uptake_interval[index];
  result.minimum_concentration = std::numeric_limits<Real>::infinity();
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    result.minimum_concentration = std::min(
        result.minimum_concentration, chem.conc_global(carbon, cell));
  }
  return result;
}

void test_population_delivery_gradient_budget() {
  const PopulationGradientDeliveryMeasurement result =
      measure_population_gradient_delivery();
  std::cerr << std::setprecision(17)
            << "    population detail funded/demand/removal/min="
            << result.funded << "/" << result.demand << "/"
            << result.field_removal << "/" << result.minimum_concentration
            << "\n";
  assert(result.demand > 0.0);
  assert(result.funded > 0.0);
  assert(result.funded < result.demand);
  assert(result.minimum_concentration >= 0.0);
  assert(result.funded <= result.demand
         + 1.0e-12 * std::max(result.demand, 1.0e-30));
  assert(std::abs(result.field_removal - result.funded)
         <= 1.0e-10 * std::max(result.funded, 1.0e-30));
  std::cout << "    population gradient funded/demand="
            << result.funded << "/" << result.demand
            << ", minimum concentration=" << result.minimum_concentration
            << "\n";
  std::cout << "  test_population_delivery_gradient_budget: PASSED\n";
}

GradientDeliveryMeasurement run_gradient_delivery(
    bool gradient_enabled, Real gradient_lambda, Real radius) {
  SimulationConfig cfg = base_config();
  cfg.enabled_fixes = {"metabolism"};
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.initial_strains[0].mu_max = 5.0e-4;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.dysbiosis_threshold = 0.0;
  constexpr Real initial_concentration = 1.0e-3;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = initial_concentration;
      chemical.boundary_conc = initial_concentration;
      chemical.z_gradient_enabled = gradient_enabled;
      chemical.z_gradient_lambda = gradient_lambda;
    } else if (chemical.name == species::IRON) {
      chemical.z_gradient_enabled = false;
      chemical.initial_conc = initial_concentration;
      chemical.boundary_conc = initial_concentration;
    }
  }

  Simulation sim;
  sim.init(cfg);
  Agent& agent = sim.agents()[0];
  agent.radius = radius;
  agent.outer_radius = radius * 1.05;
  agent.km.km_carbon = 1.0e-9;
  const Int cell = agent.grid_cell;
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  sim.domain().pos_to_grid(agent.x, ix, iy, iz);
  const Real z_rel = (iz + 0.5) * sim.domain().dx_z();
  sim.step(kDt);

  const auto& flux = sim.chemical_field().flux_accounting();
  const Int carbon = sim.chemical_field().find(species::CARBON);
  const auto index = static_cast<size_t>(carbon);
  const Real demand = flux.uptake_demand_interval[index];
  const Real realized = flux.agent_uptake_interval[index];
  GradientDeliveryMeasurement result;
  result.demand = demand;
  result.realized = realized;
  result.funded_fraction = demand > 0.0 ? realized / demand : 0.0;
  result.total_concentration =
      sim.chemical_field().conc_global(carbon, cell);
  result.background_concentration = initial_concentration
      * std::exp(-z_rel / gradient_lambda);
  (void)ix;
  (void)iy;
  return result;
}

void test_delivery_gradient_realizes_and_funds() {
  const GradientDeliveryMeasurement result =
      run_gradient_delivery(true, 10.0e-6, 5.0e-7);
  assert(result.demand > 0.0);
  assert(result.realized > 0.0);
  assert(result.funded_fraction > 0.0);
  std::cout << "    gradient_realized=" << result.realized
            << " demand=" << result.demand << "\n";
  std::cout << "  test_delivery_gradient_realizes_and_funds: PASSED\n";
}

void test_delivery_gradient_large_lambda_matches_flat_profile() {
  const GradientDeliveryMeasurement gradient =
      run_gradient_delivery(true, 1.0e9, 5.0e-7);
  const GradientDeliveryMeasurement flat =
      run_gradient_delivery(false, 1.0e9, 5.0e-7);
  const Real tolerance = 1.0e-10;
  assert(std::abs(gradient.realized - flat.realized)
         <= tolerance * std::max(flat.realized, 1.0e-30));
  assert(std::abs(gradient.funded_fraction - flat.funded_fraction)
         <= tolerance * std::max(flat.funded_fraction, 1.0e-30));
  std::cout << "  test_delivery_gradient_large_lambda_matches_flat_profile: PASSED\n";
}

void test_delivery_gradient_depletes_below_background() {
  const GradientDeliveryMeasurement result =
      run_gradient_delivery(true, 10.0e-6, 2.0e-6);
  assert(result.realized > 0.0);
  assert(result.total_concentration < result.background_concentration);
  std::cout << "  test_delivery_gradient_depletes_below_background: PASSED\n";
}

void test_delivery_gradient_sensitivity() {
  const std::vector radii = {1.0e-9, 5.0e-9, 2.0e-8};
  std::vector<Real> realized;
  for (Real radius : radii) {
    realized.push_back(
        run_gradient_delivery(true, 10.0e-6, radius).realized);
  }
  assert(realized[0] > 0.0);
  assert(realized[0] < realized[1]);
  assert(realized[1] < realized[2]);
  std::cout << "    gradient_realized_by_radius=" << realized[0] << ", "
            << realized[1] << ", " << realized[2] << "\n";
  std::cout << "  test_delivery_gradient_sensitivity: PASSED\n";
}

void test_delivery_gradient_inertness_change_detectors() {
  const Measurement flat =
      measure_single_agent(UptakeLimitMode::Delivery, 1.0e-6, 5.0e-7, 1.0);
  assert(flat.realized == flat.field_removal);
  assert(flat.realized >= 0.0);

  SimulationConfig cfg = base_config();
  cfg.fixes.metabolism.uptake_limit = "none";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::None;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.z_gradient_enabled = true;
      chemical.z_gradient_lambda = 10.0e-6;
    }
  }
  Simulation sim;
  sim.init(cfg);
  const Real initial_biomass = sim.agents()[0].biomass;
  sim.step(kDt);
  const Int carbon = sim.chemical_field().find(species::CARBON);
  assert(sim.agents()[0].biomass > initial_biomass);
  assert(sim.chemical_field().flux_accounting().agent_uptake_step[
             static_cast<size_t>(carbon)] == 0.0);
  std::cout << "  test_delivery_gradient_inertness_change_detectors: PASSED\n";
}

void test_delivery_zero_realization_closure() {
  SimulationConfig cfg = base_config();
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.time.total_time = 600.0;
  cfg.time.output_interval = 600.0;
  cfg.closure.zero_realization_grace_steps = 1;
  cfg.initial_strains.front().count = 2;
  cfg.fixes.metabolism.carbon_maintenance_rate = 1.0e-3;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.vbf.carbon_sink_vmax = 0.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.z_gradient_enabled = false;
      chemical.initial_conc = 0.0;
      chemical.boundary_conc = 0.0;
    }
  }
  Simulation sim;
  sim.init(cfg);
  const int status = sim.run();
  assert(status != 0);
  assert(sim.termination_cause() == TerminationCause::ClosureViolation);
  assert(sim.step_count() == 2);
  std::cout << "  test_delivery_zero_realization_closure: PASSED\n";
}

void test_none_mode_clip_does_not_close_by_default() {
  SimulationConfig cfg = base_config();
  cfg.initial_strains.front().count = 2;
  cfg.fixes.metabolism.uptake_limit = "none";
  cfg.time.total_time = 60.0;
  cfg.time.output_interval = 60.0;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.z_gradient_enabled = false;
      chemical.initial_conc = 1.0e-12;
      chemical.boundary_conc = 1.0e-12;
    }
  }
  Simulation sim;
  sim.init(cfg);
  assert(sim.run() == 0);
  assert(sim.termination_cause() == TerminationCause::HorizonReached);
  std::cout << "  test_none_mode_clip_does_not_close_by_default: PASSED\n";
}

}  // namespace

int main() {
#ifdef GUTIBM_POPULATION_RESOLUTION_ONLY
  test_population_scale_delivery_resolution_invariance();
  test_population_delivery_gradient_budget();
#else
  std::cout << "=== Uptake Limitation Tests ===\n";
  test_delivery_resolution_and_timestep_invariance();
  test_far_field_resolution_invariance_and_anti_vacuity();
  test_far_field_uniform_field_invariant();
  test_far_field_concentration_monotonicity();
  test_delivery_funding_matches_demand_or_ceiling();
  test_delivery_concentration_response_and_positivity();
  test_delivery_exact_removal_and_depletion_cap();
  test_delivery_adjacent_agents_are_order_independent();
  test_funded_fraction_falls_with_local_concentration();
  test_funded_fraction_falls_with_radius_and_diffusivity();
  test_realized_uptake_respects_the_cap_and_the_field();
  test_none_mode_leaves_growth_unfunded();
  test_limitation_severity_rises_with_agent_density();
  test_delivery_is_positive_and_funds_only_removed_carbon();
  test_delivery_shared_voxel_funding_does_not_exceed_removal();
  test_invalid_cell_agent_cannot_claim_delivery_funding();
  test_delivery_positivity_rationing_is_sensitive_and_conservative();
  test_delivery_rationing_is_local();
  test_regularized_delivery_mass_is_conservative();
  test_regularized_delivery_reduces_depletion_retries();
  test_regularized_delivery_funding_resolution_invariance();
  test_delivery_maintenance_reduces_growth();
  test_delivery_density_brake();
  test_delivery_queues_noncarbon_chemistry_once();
  test_delivery_preserves_negative_growth();
  test_delivery_negative_growth_books_maintenance();
  test_delivery_maintenance_limited_runs_to_horizon();
  test_delivery_gradient_realizes_and_funds();
  test_delivery_gradient_large_lambda_matches_flat_profile();
  test_delivery_gradient_depletes_below_background();
  test_delivery_gradient_sensitivity();
  test_delivery_gradient_inertness_change_detectors();
  test_delivery_gradient_representation_invariance();
  test_population_delivery_gradient_budget();
  test_delivery_zero_realization_closure();
  test_none_mode_clip_does_not_close_by_default();
  std::cout << "All uptake limitation tests passed.\n";
#endif
  return 0;
}
