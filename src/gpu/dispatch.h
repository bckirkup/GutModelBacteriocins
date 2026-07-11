#ifndef GUTIBM_GPU_DISPATCH_H
#define GUTIBM_GPU_DISPATCH_H

#include "gpu_config.h"
#include "device.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

namespace gutibm {

// Global GPU runtime state (set once during Simulation::init).
void gpu_set_config(const GpuConfig& cfg);
const GpuConfig& gpu_config();
const DeviceContext& gpu_device();

// Initialize device for this MPI rank; returns true if GPU path is active.
bool gpu_init_for_rank(int mpi_rank, int mpi_nprocs);

// True when CUDA was compiled in, config enables GPU, and a device is active.
bool gpu_runtime_enabled();

#ifdef GUTIBM_CUDA
// Default async stream for chemistry/field kernels (nullptr before init).
cudaStream_t gpu_compute_stream();
#endif

// Block until the default compute stream completes.
void gpu_sync_compute();

}  // namespace gutibm

#endif  // GUTIBM_GPU_DISPATCH_H
