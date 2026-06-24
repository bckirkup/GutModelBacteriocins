/* -----------------------------------------------------------------------
   GutIBM – Entry point
   
   Usage:
     gut_ibm [config.json]
     mpirun -np 4 gut_ibm config.json
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"

#include <iostream>
#include <string>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

int main(int argc, char** argv) {
#ifdef GUTIBM_MPI
  MPI_Init(&argc, &argv);
#endif

  gutibm::SimulationConfig cfg;

  if (argc > 1) {
    std::string config_file = argv[1];
    cfg = gutibm::InputParser::parse(config_file);
  } else {
    cfg = gutibm::InputParser::default_config();
  }

  gutibm::Simulation sim;
  if (!cfg.checkpoint.file.empty()) {
    sim.init_from_checkpoint(cfg, cfg.checkpoint.file, cfg.checkpoint.step);
  } else {
    sim.init(cfg);
  }
  sim.run();

#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif

  return 0;
}
