/* -----------------------------------------------------------------------
   GutIBM – O2 growth boost near epithelium (Spec 1)
   ----------------------------------------------------------------------- */

#include "fix_metabolism.h"
#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

static Simulation make_o2_sim() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {40e-6, 40e-6, 80e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.oxygen.enabled = true;
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  return sim;
}

static Agent make_agent_at_z(Simulation& sim, Real z, Int type = 1) {
  Vec3 pos = {
    0.5 * (sim.domain().lo()[0] + sim.domain().hi()[0]),
    0.5 * (sim.domain().lo()[1] + sim.domain().hi()[1]),
    z,
  };
  Agent a = Agent::create_default(sim.agents().next_tag(), type, pos, 5e-4);
  Int ix;
  Int iy;
  Int iz;
  sim.domain().pos_to_grid(a.x, ix, iy, iz);
  a.grid_cell = sim.domain().cell_index(ix, iy, iz);
  return a;
}

void test_o2_growth_boost_epithelial_vs_lumen() {
  auto sim = make_o2_sim();
  MetabolismConfig mcfg;
  FixMetabolism fix(sim, mcfg);

  Agent near_epithelium = make_agent_at_z(sim, 2.5e-6);
  Agent near_lumen = make_agent_at_z(sim, sim.domain().hi()[2] - 2.5e-6);
  sim.agents().push_back(near_epithelium);
  sim.agents().push_back(near_lumen);

  fix.compute(60.0);

  const Real mu_ep = sim.agents()[0].mu_realized;
  const Real mu_lumen = sim.agents()[1].mu_realized;

  assert(mu_ep > mu_lumen);
  assert(mu_ep > 0.0);

  std::cout << "  test_o2_growth_boost_epithelial_vs_lumen: PASSED"
            << " (mu_ep=" << mu_ep
            << " mu_lumen=" << mu_lumen << ")\n";
}

void test_local_o2_accessor() {
  auto sim = make_o2_sim();
  Agent a = make_agent_at_z(sim, 2.5e-6);
  sim.agents().push_back(a);
  const Real o2 = sim.local_O2(sim.agents()[0]);
  assert(o2 > 0.0);
  assert(o2 <= sim.config().oxygen.epithelial_conc * 1.01);

  std::cout << "  test_local_o2_accessor: PASSED (O2=" << o2 << ")\n";
}

int main() {
  std::cout << "=== O2 Growth Boost Tests ===\n";
  test_o2_growth_boost_epithelial_vs_lumen();
  test_local_o2_accessor();
  std::cout << "All O2 growth boost tests passed.\n";
  return 0;
}
