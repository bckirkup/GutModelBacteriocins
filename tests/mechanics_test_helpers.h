/* -----------------------------------------------------------------------
   GutIBM – Shared mechanics test helpers
   ----------------------------------------------------------------------- */

#pragma once

#include "agent.h"
#include "fix_mechanics.h"
#include "input_parser.h"
#include "simulation.h"

namespace gutibm {
namespace test {

inline Simulation make_two_agent_sim(Vec3 pos_a, Vec3 pos_b,
                                     const MechanicsConfig& mcfg = {},
                                     bool gpu_enabled = false) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.domain.hi = {100e-6, 100e-6, 100e-6};
  cfg.domain.grid_dx = 10e-6;
  cfg.domain.hash_cell_size = 20e-6;
  cfg.fixes.mechanics = mcfg;
  cfg.hdf5.enabled = false;
  cfg.gpu.enabled = gpu_enabled;

  Simulation sim;
  sim.init(cfg);

  Agent a = Agent::create_default(sim.agents().next_tag(), 1, pos_a, 5e-4);
  Agent b = Agent::create_default(sim.agents().next_tag(), 1, pos_b, 5e-4);
  sim.agents().push_back(std::move(a));
  sim.agents().push_back(std::move(b));

  sim.domain().spatial_hash().clear();
  sim.domain().spatial_hash().insert(0, sim.agents()[0].x);
  sim.domain().spatial_hash().insert(1, sim.agents()[1].x);

  return sim;
}

}  // namespace test
}  // namespace gutibm
