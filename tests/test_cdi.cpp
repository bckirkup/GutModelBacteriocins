/* -----------------------------------------------------------------------
   GutIBM – Contact-dependent inhibition tests (Spec 3)
   ----------------------------------------------------------------------- */

#include "fix_cdi.h"
#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

static void rebuild_hash(Simulation& sim) {
  auto& hash = sim.domain().spatial_hash();
  hash.clear();
  for (Int i = 0; i < sim.agents().size(); ++i) {
    const Agent& a = sim.agents()[i];
    if (a.state != PhenoState::DEAD) {
      hash.insert(i, a.x);
      continue;
    }
    if (sim.config().cell_bio.cdi.enabled && a.timers.death_time >= 0.0
        && (sim.time() - a.timers.death_time) < sim.config().cell_bio.cdi.corpse_persistence) {
      hash.insert(i, a.x);
    }
  }
}

static Simulation make_cdi_sim() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.cell_bio.cdi.enabled = true;
  cfg.cell_bio.cdi.kill_rate = 1.0;
  cfg.cell_bio.cdi.contact_radius = 2.0e-6;
  cfg.cell_bio.cdi.corpse_persistence = 300.0;

  Simulation sim;
  sim.init(cfg);
  return sim;
}

static Vec3 domain_center(const Simulation& sim) {
  return {
    0.5 * (sim.domain().lo()[0] + sim.domain().hi()[0]),
    0.5 * (sim.domain().lo()[1] + sim.domain().hi()[1]),
    0.5 * (sim.domain().lo()[2] + sim.domain().hi()[2]),
  };
}

static Agent make_live_agent(Simulation& sim, Vec3 pos, uint16_t cdi_type,
                             uint16_t cdi_immunity) {
  Agent a = Agent::create_default(sim.agents().next_tag(), 1, pos, 5e-4);
  a.genome.cdi_type = cdi_type;
  a.genome.cdi_immunity = cdi_immunity;
  Int ix;
  Int iy;
  Int iz;
  sim.domain().pos_to_grid(a.x, ix, iy, iz);
  a.grid_cell = sim.domain().cell_index(ix, iy, iz);
  return a;
}

void test_cdi_kills_neighbor() {
  auto sim = make_cdi_sim();
  Vec3 center = domain_center(sim);
  Real offset = 1.0e-6;

  Agent attacker = make_live_agent(sim, center, 1, 1);
  Agent victim = make_live_agent(sim,
      {center[0] + offset, center[1], center[2]}, 0, 0);
  sim.agents().push_back(std::move(attacker));
  sim.agents().push_back(std::move(victim));

  rebuild_hash(sim);

  CdiConfig cfg = sim.config().cell_bio.cdi;
  FixCdi fix(sim, cfg);
  for (int step = 0; step < 50; ++step) {
    fix.compute(0.1);
  }

  assert(sim.agents()[1].state == PhenoState::DEAD);
  assert(sim.agents()[1].timers.death_time >= 0.0);

  std::cout << "  test_cdi_kills_neighbor: PASSED\n";
}

void test_cdi_immunity_protects() {
  auto sim = make_cdi_sim();
  Vec3 center = domain_center(sim);
  Real offset = 1.0e-6;

  Agent attacker = make_live_agent(sim, center, 1, 1);
  Agent immune = make_live_agent(sim,
      {center[0] + offset, center[1], center[2]}, 0, 1);
  sim.agents().push_back(std::move(attacker));
  sim.agents().push_back(std::move(immune));

  rebuild_hash(sim);

  FixCdi fix(sim, sim.config().cell_bio.cdi);
  for (int step = 0; step < 50; ++step) {
    fix.compute(0.1);
  }

  assert(sim.agents()[1].state != PhenoState::DEAD);

  std::cout << "  test_cdi_immunity_protects: PASSED\n";
}

static Int count_cdi_kills(bool with_barrier, Real offset) {
  Simulation trial = make_cdi_sim();
  const Vec3 c = domain_center(trial);

  Agent attacker = make_live_agent(trial, c, 1, 1);
  trial.agents().push_back(std::move(attacker));

  if (with_barrier) {
    Agent corpse = make_live_agent(trial, {c[0] + offset, c[1], c[2]}, 0, 0);
    corpse.state = PhenoState::DEAD;
    corpse.timers.death_time = 0.0;
    trial.agents().push_back(std::move(corpse));
  }

  Agent victim = make_live_agent(trial,
      {c[0] + (with_barrier ? 2.0 : 1.0) * offset, c[1], c[2]}, 0, 0);
  trial.agents().push_back(std::move(victim));

  const Int victim_idx = trial.agents().size() - 1;
  rebuild_hash(trial);

  FixCdi fix(trial, trial.config().cell_bio.cdi);
  Int kills = 0;
  for (int step = 0; step < 30; ++step) {
    if (trial.agents()[victim_idx].state == PhenoState::DEAD) ++kills;
    fix.compute(0.1);
    rebuild_hash(trial);
  }
  if (trial.agents()[victim_idx].state == PhenoState::DEAD) ++kills;
  return kills;
}

void test_cdi_corpse_barrier() {
  const Real offset = 1.0e-6;
  const Int kills_barrier = count_cdi_kills(true, offset);
  const Int kills_clear = count_cdi_kills(false, offset);
  assert(kills_clear > kills_barrier);

  std::cout << "  test_cdi_corpse_barrier: PASSED\n";
}

int main() {
  std::cout << "=== CDI Tests ===\n";
  test_cdi_kills_neighbor();
  test_cdi_immunity_protects();
  test_cdi_corpse_barrier();
  std::cout << "All CDI tests passed.\n";
  return 0;
}
