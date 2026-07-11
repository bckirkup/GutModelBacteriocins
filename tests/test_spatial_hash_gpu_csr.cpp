/* -----------------------------------------------------------------------
   Spatial hash GPU CSR build smoke test
   ----------------------------------------------------------------------- */

#include "spatial_hash_gpu.h"
#include "agent_pool_gpu.h"
#include "dispatch.h"
#include "simulation.h"
#include "input_parser.h"

#include <cassert>
#include <iostream>

using namespace gutibm;

int main() {
  std::cout << "=== Spatial Hash GPU CSR ===\n";

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  GpuConfig gcfg;
  gcfg.enabled = true;
  gpu_set_config(gcfg);
  if (!gpu_init_for_rank(0, 1)) {
    std::cout << "  SKIPPED (no CUDA device)\n";
    return 0;
  }

  SimulationConfig cfg = InputParser::default_config();
  cfg.hdf5.enabled = false;
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 40;
  strain.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(strain);

  Simulation sim;
  sim.init(cfg);

  AgentPoolGpu ag;
  ag.sync_from_host(sim.agents());

  SpatialHashGpu hash;
  const Real cell_size = 10.0e-6;
  const Vec3 lo = sim.domain().lo();
  const Vec3 hi = sim.domain().hi();
  assert(gpu_build_spatial_hash(ag, ag.size(), lo, hi, cell_size, hash));
  assert(hash.active());
  assert(hash.cell_offsets.size() > 0);
  assert(hash.cell_counts.size() > 0);
  assert(hash.sorted_agent_indices.size() == static_cast<size_t>(ag.size()));

  std::cout << "  test_spatial_hash_gpu_csr: PASSED\n";
  return 0;
#endif
}
