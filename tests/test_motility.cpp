/* -----------------------------------------------------------------------
   GutIBM – Active motility tests (Spec 3)
   ----------------------------------------------------------------------- */

#include "fix_motility.h"
#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

static Simulation make_motility_sim(bool chemotaxis = false) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.cell_bio.motility.enabled = true;
  cfg.cell_bio.motility.swim_speed = 7.76e-6;
  cfg.cell_bio.motility.run_mean_duration = 5.0;
  cfg.cell_bio.motility.stop_probability = 0.0;
  cfg.cell_bio.motility.chemotaxis_enabled = chemotaxis;
  cfg.cell_bio.motility.chi_carbon = 1.0;
  cfg.cell_bio.motility.chi_oxygen = 0.0;

  Simulation sim;
  sim.init(cfg);
  return sim;
}

static Agent make_motile_agent(Simulation& sim) {
  Vec3 center = {
    0.5 * (sim.domain().lo()[0] + sim.domain().hi()[0]),
    0.5 * (sim.domain().lo()[1] + sim.domain().hi()[1]),
    0.5 * (sim.domain().lo()[2] + sim.domain().hi()[2]),
  };
  Agent a = Agent::create_default(sim.agents().next_tag(), 1, center, 5e-4);
  FixMotility::init_agent_motility(a, sim.config().cell_bio.motility, sim.rng());
  a.motility.swim_direction = {1.0, 0.0, 0.0};
  a.motility.is_stopped = false;
  a.motility.run_timer = 100.0;
  Int ix;
  Int iy;
  Int iz;
  sim.domain().pos_to_grid(a.x, ix, iy, iz);
  a.grid_cell = sim.domain().cell_index(ix, iy, iz);
  return a;
}

static Real displacement(const Vec3& start, const Vec3& end) {
  Real dx = end[0] - start[0];
  Real dy = end[1] - start[1];
  Real dz = end[2] - start[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void test_motility_displacement() {
  auto sim_motile = make_motility_sim();
  Agent motile = make_motile_agent(sim_motile);
  Vec3 start = motile.x;
  sim_motile.agents().push_back(std::move(motile));

  Simulation sim_static;
  SimulationConfig static_cfg = InputParser::default_config();
  static_cfg.initial_strains.clear();
  static_cfg.hdf5.enabled = false;
  static_cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  static_cfg.domain.grid_dx = 5e-6;
  static_cfg.cell_bio.motility.enabled = false;
  sim_static.init(static_cfg);
  Agent still = make_motile_agent(sim_static);
  Vec3 still_start = still.x;
  sim_static.agents().push_back(std::move(still));

  const Real dt = 1.0;
  const int steps = 20;
  FixMotility mot_fix(sim_motile, sim_motile.config().cell_bio.motility);
  for (int i = 0; i < steps; ++i) {
    mot_fix.pre_step(dt);
    Agent& a = sim_motile.agents()[0];
    a.x[0] += a.motility.swim_direction[0] * sim_motile.config().cell_bio.motility.swim_speed * dt;
    a.x[1] += a.motility.swim_direction[1] * sim_motile.config().cell_bio.motility.swim_speed * dt;
    a.x[2] += a.motility.swim_direction[2] * sim_motile.config().cell_bio.motility.swim_speed * dt;
  }

  const Real motile_disp = displacement(start, sim_motile.agents()[0].x);
  const Real static_disp = displacement(still_start, sim_static.agents()[0].x);
  assert(motile_disp > static_disp);
  assert(motile_disp > 1.0e-6);

  std::cout << "  test_motility_displacement: PASSED\n";
}

void test_chemotaxis_bias() {
  auto sim = make_motility_sim(true);
  Agent a = make_motile_agent(sim);
  a.motility.swim_direction = {1.0, 0.0, 0.0};
  a.motility.prev_carbon = 0.0;
  Int cell = a.grid_cell;
  Int i_carbon = sim.chemical_field().find("carbon");
  assert(i_carbon >= 0);
  sim.chemical_field().conc(i_carbon, cell) = 0.0;
  sim.agents().push_back(std::move(a));

  FixMotility fix(sim, sim.config().cell_bio.motility);
  const Real dt = 1.0;
  fix.pre_step(dt);
  const Real run_timer_after_low = sim.agents()[0].motility.run_timer;
  sim.chemical_field().conc(i_carbon, cell) = 1.0;
  fix.pre_step(dt);

  assert(sim.agents()[0].motility.run_timer > run_timer_after_low);

  std::cout << "  test_chemotaxis_bias: PASSED\n";
}

int main() {
  std::cout << "=== Motility Tests ===\n";
  test_motility_displacement();
  test_chemotaxis_bias();
  std::cout << "All motility tests passed.\n";
  return 0;
}
