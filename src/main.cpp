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

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

#ifdef GUTIBM_MPI
#include <mpi.h>
#include "error.h"
#endif

namespace {

// OpenMPI on some stacks (notably WSL2) probes CUDA after MPI_Init and can leave
// the runtime without a visible device unless the driver is touched first.
// gpu_smoke avoids this because it never calls MPI_Init.
void cuda_runtime_probe_before_mpi() {
#ifdef GUTIBM_CUDA
  int count = 0;
  if (cudaError_t err = cudaGetDeviceCount(&count); err != cudaSuccess) {
    std::cerr << "Warning: CUDA probe before MPI_Init failed: "
              << cudaGetErrorString(err) << "\n";
    return;
  }
  if (count > 0) {
    (void)cudaFree(nullptr);
  }
#endif
}

}  // namespace

int main(int argc, char** argv) {
  cuda_runtime_probe_before_mpi();

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
