/* -----------------------------------------------------------------------
   GutIBM – Fur-regulated receptor expression tests (Spec 3)
   ----------------------------------------------------------------------- */

#include "fix_metabolism.h"
#include "fix_receptor.h"
#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

static Simulation make_fur_sim(uint64_t seed = 42) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.cell_bio.fur.enabled = true;
  cfg.cell_bio.fur.Km = 1.0e-5;
  cfg.cell_bio.fur.upregulation_max = 4.0;
  cfg.cell_bio.fur.receptor_max = 5.0;
  cfg.seed = seed;

  Simulation sim;
  sim.init(cfg);
  return sim;
}

static Agent make_agent_at_center(Simulation& sim) {
  Vec3 center = {
    0.5 * (sim.domain().lo()[0] + sim.domain().hi()[0]),
    0.5 * (sim.domain().lo()[1] + sim.domain().hi()[1]),
    0.5 * (sim.domain().lo()[2] + sim.domain().hi()[2]),
  };
  Agent a = Agent::create_default(sim.agents().next_tag(), 1, center, 5e-4);
  Int ix;
  Int iy;
  Int iz;
  sim.domain().pos_to_grid(a.x, ix, iy, iz);
  a.grid_cell = sim.domain().cell_index(ix, iy, iz);
  return a;
}

static void set_iron(Simulation& sim, Int cell, Real conc) {
  Int i_iron = sim.chemical_field().find("iron");
  assert(i_iron >= 0);
  sim.chemical_field().conc(i_iron, cell) = conc;
}

void test_fur_upregulation() {
  auto sim = make_fur_sim();
  Agent a = make_agent_at_center(sim);
  set_iron(sim, a.grid_cell, 0.0);
  sim.agents().push_back(std::move(a));

  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  fix.compute(1.0);

  const int ri = to_underlying(ReceptorType::FepA);
  assert(sim.agents()[0].receptor_expr[ri] > sim.agents()[0].receptor_expr_base[ri]);
  assert(sim.agents()[0].receptor_expr[ri] > 4.0);

  std::cout << "  test_fur_upregulation: PASSED\n";
}

void test_fur_repression() {
  auto sim = make_fur_sim();
  Agent a = make_agent_at_center(sim);
  set_iron(sim, a.grid_cell, 1.0);
  sim.agents().push_back(std::move(a));

  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  fix.compute(1.0);

  const int ri = to_underlying(ReceptorType::FepA);
  assert(std::abs(sim.agents()[0].receptor_expr[ri]
                  - sim.agents()[0].receptor_expr_base[ri]) < 0.05);

  std::cout << "  test_fur_repression: PASSED\n";
}

void test_fur_increases_susceptibility() {
  ReceptorConfig rcfg;
  rcfg.kill_rate_colicin = 0.02;

  int deaths_low = 0;
  int deaths_high = 0;
  const int trials = 40;

  for (int trial = 0; trial < trials; ++trial) {
    auto sim_low = make_fur_sim(1000 + static_cast<uint64_t>(trial));
    Agent low_iron = make_agent_at_center(sim_low);
    low_iron.receptor_expr_base[to_underlying(ReceptorType::FepA)] = 1.0;
    low_iron.receptor_expr[to_underlying(ReceptorType::CirA)] = 0.0;
    low_iron.receptor_expr[to_underlying(ReceptorType::BtuB)] = 0.0;
    Int cell_low = low_iron.grid_cell;
    set_iron(sim_low, cell_low, 0.0);
    Int i_tox = sim_low.chemical_field().find("bacteriocin");
    Int i_b12 = sim_low.chemical_field().find("b12");
    sim_low.chemical_field().conc(i_tox, cell_low) = 1.0e-4;
    sim_low.chemical_field().conc(i_b12, cell_low) = 0.0;
    sim_low.agents().push_back(std::move(low_iron));

    auto sim_high = make_fur_sim(2000 + static_cast<uint64_t>(trial));
    Agent high_iron = make_agent_at_center(sim_high);
    high_iron.receptor_expr_base[to_underlying(ReceptorType::FepA)] = 1.0;
    high_iron.receptor_expr[to_underlying(ReceptorType::CirA)] = 0.0;
    high_iron.receptor_expr[to_underlying(ReceptorType::BtuB)] = 0.0;
    Int cell_high = high_iron.grid_cell;
    set_iron(sim_high, cell_high, 1.0);
    sim_high.chemical_field().conc(i_tox, cell_high) = 1.0e-4;
    sim_high.chemical_field().conc(i_b12, cell_high) = 0.0;
    sim_high.agents().push_back(std::move(high_iron));

    FixMetabolism metab_low(sim_low, sim_low.config().fixes.metabolism);
    FixMetabolism metab_high(sim_high, sim_high.config().fixes.metabolism);
    metab_low.compute(1.0);
    metab_high.compute(1.0);

    assert(sim_low.agents()[0].receptor_expr[to_underlying(ReceptorType::FepA)]
           > sim_high.agents()[0].receptor_expr[to_underlying(ReceptorType::FepA)]);

    FixReceptor receptor_low(sim_low, rcfg);
    FixReceptor receptor_high(sim_high, rcfg);
    receptor_low.compute(60.0);
    receptor_high.compute(60.0);

    if (sim_low.agents()[0].state == PhenoState::DEAD) ++deaths_low;
    if (sim_high.agents()[0].state == PhenoState::DEAD) ++deaths_high;
  }

  assert(deaths_low > deaths_high);

  std::cout << "  test_fur_increases_susceptibility: PASSED\n";
}

int main() {
  std::cout << "=== Fur Regulation Tests ===\n";
  test_fur_upregulation();
  test_fur_repression();
  test_fur_increases_susceptibility();
  std::cout << "All fur regulation tests passed.\n";
  return 0;
}
