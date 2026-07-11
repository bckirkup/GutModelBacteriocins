/* -----------------------------------------------------------------------
   GutIBM – CUDA-aware MPI runtime detection (issue #156)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_CUDA_AWARE_MPI_H
#define GUTIBM_CUDA_AWARE_MPI_H

namespace gutibm {

// True when the linked MPI library reports CUDA buffer support (Open MPI
// MPIX_Query_cuda_support) or GUTIBM_CUDA_AWARE_MPI_FORCE=1 for manual HPC
// clusters where auto-detection is unavailable.
bool cuda_aware_mpi_runtime_available();

}  // namespace gutibm

#endif  // GUTIBM_CUDA_AWARE_MPI_H
