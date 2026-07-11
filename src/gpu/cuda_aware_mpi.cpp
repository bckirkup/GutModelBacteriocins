/* -----------------------------------------------------------------------
   GutIBM – CUDA-aware MPI runtime detection (issue #156)
   ----------------------------------------------------------------------- */

#include "cuda_aware_mpi.h"

#include <cstdlib>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

#if defined(GUTIBM_CUDA) && defined(GUTIBM_MPI) \
    && (defined(OPEN_MPI) || defined(OMPI_MAJOR_VERSION))
#include <mpi-ext.h>
#endif

namespace gutibm {

bool cuda_aware_mpi_runtime_available() {
#ifndef GUTIBM_CUDA
  return false;
#else
  const char* force = std::getenv("GUTIBM_CUDA_AWARE_MPI_FORCE");
  if (force != nullptr
      && (force[0] == '1' || force[0] == 't' || force[0] == 'T')) {
    return true;
  }

#ifdef GUTIBM_MPI
#if defined(OPEN_MPI) || defined(OMPI_MAJOR_VERSION)
  return MPIX_Query_cuda_support() > 0;
#else
  return false;
#endif
#else
  return false;
#endif
#endif
}

}  // namespace gutibm
