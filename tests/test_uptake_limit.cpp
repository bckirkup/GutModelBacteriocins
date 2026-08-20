/* -----------------------------------------------------------------------
   GutIBM – Agent-side nutrient uptake limitation tests
   ----------------------------------------------------------------------- */

#include "fix_metabolism.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
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

}  // namespace

int main() {
  std::cout << "=== Uptake Limitation Tests ===\n";
  test_funded_fraction_falls_with_local_concentration();
  test_funded_fraction_falls_with_radius_and_diffusivity();
  test_realized_uptake_respects_the_cap_and_the_field();
  test_none_mode_leaves_growth_unfunded();
  test_limitation_severity_rises_with_agent_density();
  std::cout << "All uptake limitation tests passed.\n";
  return 0;
}
