/* -----------------------------------------------------------------------
   GutIBM – Delivery support cache OpenMP regression test
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#ifdef GUTIBM_OPENMP
#include <omp.h>
#endif

using namespace gutibm;

namespace {

constexpr Real kDt = 60.0;
constexpr Int kAgentCount = 96;

#ifdef GUTIBM_OPENMP
struct DeliveryResult {
  std::vector<Real> deposits;
  Real funded = 0.0;
  Real demanded = 0.0;
};

SimulationConfig delivery_config(Real far_field_radius) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = 1729;
  cfg.hdf5.enabled = false;
  cfg.enabled_fixes = {"metabolism"};
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {50.0e-6, 50.0e-6, 50.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = kAgentCount;
  strain.mu_max = 1.0e-2;
  cfg.initial_strains.push_back(strain);
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.fixes.metabolism.carbon_maintenance_rate = 0.0;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  cfg.fixes.metabolism.delivery_far_field_radius = far_field_radius;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.chem_env.oxygen.enabled = false;
  cfg.chem_env.acetate.enabled = false;
  cfg.chem_env.siderophore.enabled = false;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e2;
      chemical.boundary_conc = 1.0e2;
      chemical.diff_coeff = 5.0e-10;
      chemical.z_gradient_enabled = false;
    }
  }
  return cfg;
}

void place_agents(Simulation& sim) {
  const Real dx = sim.config().domain.grid_dx;
  const Int side = 10;
  for (Int i = 0; i < sim.agents().size(); ++i) {
    const Int x = i % side;
    const Int y = (i / side) % side;
    const Int z = i / (side * side);
    Agent& agent = sim.agents()[i];
    agent.x = {(static_cast<Real>(x) + 0.5) * dx,
               (static_cast<Real>(y) + 0.5) * dx,
               (static_cast<Real>(z) + 0.5) * dx};
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid(agent.x, ix, iy, iz);
    agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
    agent.radius = 5.0e-7;
    agent.outer_radius = agent.radius * 1.05;
    agent.km.km_carbon = 1.0e-9;
  }
}

DeliveryResult run_delivery(Real far_field_radius, int thread_count) {
#ifdef GUTIBM_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(thread_count);
#else
  (void)thread_count;
#endif
  Simulation sim;
  sim.init(delivery_config(far_field_radius));
  place_agents(sim);
  sim.step(kDt);

  const ChemicalField& chem = sim.chemical_field();
  const Int carbon = chem.find(species::CARBON);
  assert(carbon >= 0);
  DeliveryResult result;
  result.deposits.reserve(static_cast<size_t>(chem.global_ncells()));
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    const Real deposit = chem.sink_realized_global(carbon, cell);
    assert(std::isfinite(deposit));
    assert(deposit >= 0.0);
    result.deposits.push_back(deposit);
  }
  const auto& flux = chem.flux_accounting();
  const size_t index = static_cast<size_t>(carbon);
  result.funded = flux.agent_uptake_interval[index]
      + flux.agent_uptake_cumulative[index];
  result.demanded = flux.uptake_demand_interval[index]
      + flux.uptake_demand_cumulative[index];
  assert(std::isfinite(result.funded));
  assert(std::isfinite(result.demanded));
  assert(result.funded >= 0.0);
  assert(result.demanded >= 0.0);
  return result;
}

void assert_same_delivery(const DeliveryResult& serial,
                          const DeliveryResult& parallel) {
  assert(serial.deposits.size() == parallel.deposits.size());
  for (size_t cell = 0; cell < serial.deposits.size(); ++cell) {
    const Real scale = std::max({1.0, std::abs(serial.deposits[cell]),
                                 std::abs(parallel.deposits[cell])});
    assert(std::abs(serial.deposits[cell] - parallel.deposits[cell])
           <= 1.0e-12 * scale);
  }
  assert(std::abs(serial.funded - parallel.funded)
         <= 1.0e-12 * std::max({1.0, std::abs(serial.funded),
                                std::abs(parallel.funded)}));
  assert(std::abs(serial.demanded - parallel.demanded)
         <= 1.0e-12 * std::max({1.0, std::abs(serial.demanded),
                                std::abs(parallel.demanded)}));
}
#endif

void test_delivery_support_serial_parallel_parity() {
#ifdef GUTIBM_OPENMP
  const int parallel_threads = 4;
  for (const Real radius : {0.0, 1.0e-5}) {
    const DeliveryResult serial = run_delivery(radius, 1);
    const DeliveryResult parallel = run_delivery(radius, parallel_threads);
    assert(serial.funded > 0.0);
    assert(serial.demanded > 0.0);
    assert(std::any_of(
        serial.deposits.begin(), serial.deposits.end(),
        [](const Real deposit) { return deposit > 0.0; }));
    assert_same_delivery(serial, parallel);
    std::cout << "  radius=" << radius
              << " serial_funded=" << serial.funded
              << " parallel_funded=" << parallel.funded
              << " serial_demanded=" << serial.demanded
              << " parallel_demanded=" << parallel.demanded << "\n";
  }
  std::cout << "  test_delivery_support_serial_parallel_parity: PASSED\n";
#else
  std::cout << "  test_delivery_support_serial_parallel_parity: SKIPPED"
            << " (OpenMP unavailable)\n";
#endif
}

}  // namespace

int main() {
  test_delivery_support_serial_parallel_parity();
  return 0;
}
