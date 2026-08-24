/* -----------------------------------------------------------------------
   GutIBM – Agent-side nutrient uptake limitation tests
   ----------------------------------------------------------------------- */

#include "fix_metabolism.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
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
    Real radius = 5.0e-7, Real boundary_concentration = -1.0) {
  SimulationConfig cfg = base_config();
  cfg.domain.grid_dx = grid_dx;
  cfg.domain.hi = {4.0 * grid_dx, 4.0 * grid_dx, 4.0 * grid_dx};
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
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
      sim.chemical_field().total_conc_global(
          carbon, cell, sim.domain());
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
  return result;
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
  assert(coarse.delivery_reduction == 0.0);
  assert(short_step.delivery_reduction == 0.0);
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
  assert(std::abs(starved.field_removal - starved.funded)
         <= 1.0e-12 * std::max(starved.funded, 1.0e-30));
  assert(starved.funded < std::min(starved.demand, starved.ceiling));
  assert(starved.min_concentration >= 0.0);
  assert(starved.delivery_reduction > 0.0);
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
  std::cout << "=== Uptake Limitation Tests ===\n";
  test_delivery_resolution_and_timestep_invariance();
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
  test_delivery_zero_realization_closure();
  test_none_mode_clip_does_not_close_by_default();
  std::cout << "All uptake limitation tests passed.\n";
  return 0;
}
