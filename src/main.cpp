/* -----------------------------------------------------------------------
   GutIBM – Entry point
   
   Usage:
     gut_ibm [config.json]
     mpirun -np 4 gut_ibm config.json
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "segment.h"
#include "input_parser.h"
#include "stop_signal.h"
#include "error.h"

#include <exception>
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

  try {
    gutibm::SimulationConfig cfg;

    if (argc > 1) {
      std::string config_file = argv[1];
      cfg = gutibm::InputParser::parse(config_file);
    } else {
      cfg = gutibm::InputParser::default_config();
    }

    gutibm::install_stop_signal_handlers();

    int exit_code = 0;
    int n_ranks = 1;
#ifdef GUTIBM_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
#endif
    if (cfg.layer2.enabled && n_ranks > 1) {
      throw gutibm::ConfigError(
          "layer2.enabled requires a single rank in Phase 1 "
          "(concurrent-CPU patches; MPI segment decomposition is a "
          "later phase)");
    }
    if (cfg.layer2.enabled) {
      // Spec 13 Phase 1: the Layer-2 mucus segment drives the run;
      // Layer-1 Simulation instances appear only inside audit slots.
      gutibm::MucusSegment segment;
      if (!cfg.checkpoint.file.empty()) {
        segment.init_from_checkpoint(cfg, cfg.checkpoint.file);
      } else {
        segment.init(cfg);
      }
      exit_code = segment.run();
    } else {
      gutibm::Simulation sim;
      if (!cfg.checkpoint.file.empty()) {
        sim.init_from_checkpoint(cfg, cfg.checkpoint.file, cfg.checkpoint.step);
      } else {
        sim.init(cfg);
      }
      exit_code = sim.run();
    }

#ifdef GUTIBM_MPI
    MPI_Finalize();
#endif

    return exit_code;
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << "\n";
#ifdef GUTIBM_MPI
    MPI_Abort(MPI_COMM_WORLD, 2);
#endif
    return 2;
  }
}
