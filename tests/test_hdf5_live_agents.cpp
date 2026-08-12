/* -----------------------------------------------------------------------
   GutIBM – analysis output excludes retained dead agents
   ----------------------------------------------------------------------- */

#include "hdf5_writer.h"
#include "input_parser.h"
#include "path_utils.h"
#include "simulation.h"

#include <cassert>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

using namespace gutibm;

int main() {
#ifndef GUTIBM_HDF5
  std::cout << "HDF5 disabled — skipping live-agent output test.\n";
  return 0;
#else
  const std::string filename =
      resolve_test_h5_path("GUTIBM_LIVE_AGENTS_H5", "live_agents");

  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {40e-6, 40e-6, 20e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.hdf5.enabled = false;
  cfg.initial_strains[0].count = 3;
  cfg.initial_strains.resize(1);

  Simulation sim;
  sim.init(cfg);
  assert(sim.agents().size() == 3);
  const TagID dead_id = sim.agents()[0].identity.tag;
  const Real dead_biomass = sim.agents()[0].biomass;
  sim.agents()[0].state = PhenoState::DEAD;

  Real live_biomass = 0.0;
  std::vector<TagID> live_ids;
  for (const Agent& agent : sim.agents()) {
    if (agent.state == PhenoState::DEAD) continue;
    live_biomass += agent.biomass;
    live_ids.push_back(agent.identity.tag);
  }
  assert(live_ids.size() == 2);

  HDF5Config hdf5;
  hdf5.filename = filename;
  hdf5.schedule.summary = 0;
  hdf5.schedule.agents = 1;
  hdf5.schedule.grid = 0;
  hdf5.schedule.lineage = 1;
  hdf5.schedule.genome = 1;

  HDF5Writer writer;
  writer.init(hdf5, sim.domain());
  assert(writer.is_enabled());
  writer.write_step(sim, 0, 0.0, 0.0);
  writer.finalize();

  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  hid_t dataset =
      H5Dopen2(file, "agents/step_000000/id", H5P_DEFAULT);
  assert(dataset >= 0);
  hid_t space = H5Dget_space(dataset);
  hsize_t count = 0;
  H5Sget_simple_extent_dims(space, &count, nullptr);
  assert(count == live_ids.size());
  std::vector<int64_t> ids(static_cast<size_t>(count));
  H5Dread(dataset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids.data());
  H5Sclose(space);
  H5Dclose(dataset);

  dataset = H5Dopen2(file, "agents/step_000000/biomass", H5P_DEFAULT);
  assert(dataset >= 0);
  std::vector<double> biomass(live_ids.size());
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          biomass.data());
  H5Dclose(dataset);

  double serialized_biomass = 0.0;
  for (const int64_t id : ids) {
    assert(id != dead_id);
    assert(std::find(live_ids.begin(), live_ids.end(), id) != live_ids.end());
  }
  for (const double value : biomass) serialized_biomass += value;
  assert(std::abs(serialized_biomass - live_biomass) < 1.0e-30);
  assert(std::abs(serialized_biomass - (live_biomass + dead_biomass)) > 1.0e-30);

  const std::array<const char*, 7> aligned_datasets = {
      "lineage/step_000000/btuB_expression",
      "lineage/step_000000/fepA_expression",
      "lineage/step_000000/num_bi_loci",
      "lineage/step_000000/generation",
      "genome/step_000000/parent_id",
      "genome/step_000000/mutations",
      "genome/step_000000/has_conjugative_plasmid",
  };
  for (const char* path : aligned_datasets) {
    dataset = H5Dopen2(file, path, H5P_DEFAULT);
    assert(dataset >= 0);
    space = H5Dget_space(dataset);
    hsize_t extent = 0;
    H5Sget_simple_extent_dims(space, &extent, nullptr);
    assert(extent == live_ids.size());
    H5Sclose(space);
    H5Dclose(dataset);
  }

  H5Fclose(file);
  std::cout << "All HDF5 live-agent output tests passed.\n";
  return 0;
#endif
}
