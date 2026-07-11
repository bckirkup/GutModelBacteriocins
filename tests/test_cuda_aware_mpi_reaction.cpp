/* -----------------------------------------------------------------------
   GutIBM – CUDA-aware MPI reaction reduce validation (issue #156)
   Run with: mpirun -np 2 test_cuda_aware_mpi_reaction
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "chemical_field_gpu.h"
#include "cuda_aware_mpi.h"
#include "dispatch.h"
#include "device.h"
#include "domain.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

using namespace gutibm;

namespace {

#ifdef GUTIBM_MPI

void require_two_ranks() {
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  assert(nprocs == 2);
}

void test_device_reduce_gated_without_env() {
  require_two_ranks();

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifndef GUTIBM_CUDA
  if (rank == 0) {
    std::cout << "  test_device_reduce_gated_without_env: SKIPPED (no CUDA)\n";
  }
  return;
#else
  if (DeviceContext::device_count() <= 0) {
    if (rank == 0) {
      std::cout << "  test_device_reduce_gated_without_env: SKIPPED (no device)\n";
    }
    return;
  }

  unsetenv("GUTIBM_CUDA_AWARE_MPI");

  DomainConfig domain_cfg;
  domain_cfg.lo = {0.0, 0.0, 0.0};
  domain_cfg.hi = {20e-6, 15e-6, 15e-6};
  domain_cfg.grid_dx = 5e-6;
  domain_cfg.periodic = {true, true, false};
  Domain domain;
  domain.init(domain_cfg);

  ChemicalSpec oxygen;
  oxygen.name = "oxygen";
  oxygen.diff_coeff = 2.1e-9;
  oxygen.retardation = 1.0;
  oxygen.initial_conc = 0.0;
  oxygen.boundary_conc = 1.0;
  oxygen.diffusion_enabled = true;

  ChemicalField chem;
  chem.init(domain, {oxygen});
  ChemicalFieldGpu chem_gpu;
  chem_gpu.init(chem);

  const Int reaction_cell = domain.cell_index(1, 1, 1);
  chem.reac(0, reaction_cell) = static_cast<Real>(rank + 1);
  chem_gpu.sync_reactions_to_device(chem);

  assert(!chem_gpu.try_sum_reactions_on_device(chem));

  if (rank == 0) {
    std::cout << "  test_device_reduce_gated_without_env: PASSED\n";
  }
#endif
}

void test_device_reduce_matches_host_when_available() {
  require_two_ranks();

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifndef GUTIBM_CUDA
  if (rank == 0) {
    std::cout << "  test_device_reduce_matches_host_when_available: SKIPPED (no CUDA)\n";
  }
  return;
#else
  if (DeviceContext::device_count() <= 0) {
    if (rank == 0) {
      std::cout << "  test_device_reduce_matches_host_when_available: SKIPPED (no device)\n";
    }
    return;
  }

  if (!cuda_aware_mpi_runtime_available()) {
    if (rank == 0) {
      std::cout << "  test_device_reduce_matches_host_when_available: SKIPPED"
                << " (CUDA-aware MPI not available; set GUTIBM_CUDA_AWARE_MPI_FORCE=1 on HPC)\n";
    }
    return;
  }

  setenv("GUTIBM_CUDA_AWARE_MPI", "1", 1);

  DomainConfig domain_cfg;
  domain_cfg.lo = {0.0, 0.0, 0.0};
  domain_cfg.hi = {20e-6, 15e-6, 15e-6};
  domain_cfg.grid_dx = 5e-6;
  domain_cfg.periodic = {true, true, false};
  Domain domain;
  domain.init(domain_cfg);

  ChemicalSpec oxygen;
  oxygen.name = "oxygen";
  oxygen.diff_coeff = 2.1e-9;
  oxygen.retardation = 1.0;
  oxygen.initial_conc = 0.0;
  oxygen.boundary_conc = 1.0;
  oxygen.diffusion_enabled = true;

  ChemicalField host_ref;
  host_ref.init(domain, {oxygen});
  ChemicalField device_field;
  device_field.init(domain, {oxygen});
  ChemicalFieldGpu chem_gpu;
  chem_gpu.init(device_field);

  const Int reaction_cell = domain.cell_index(1, 1, 1);
  const Real local_value = static_cast<Real>(rank + 1);
  host_ref.reac(0, reaction_cell) = local_value;
  device_field.reac(0, reaction_cell) = local_value;
  chem_gpu.sync_reactions_to_device(device_field);

  host_ref.sum_reactions_across_ranks();
  assert(std::abs(host_ref.reac(0, reaction_cell) - 3.0) < 1e-15);

  assert(chem_gpu.try_sum_reactions_on_device(device_field));
  assert(std::abs(device_field.reac(0, reaction_cell) - 3.0) < 1e-12);

  if (rank == 0) {
    std::cout << "  test_device_reduce_matches_host_when_available: PASSED\n";
  }

  unsetenv("GUTIBM_CUDA_AWARE_MPI");
#endif
}

void test_runtime_detection_reports_status() {
  require_two_ranks();

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const bool available = cuda_aware_mpi_runtime_available();
  if (rank == 0) {
    std::cout << "  test_runtime_detection_reports_status: PASSED"
              << " (cuda_aware_mpi_runtime_available=" << (available ? "true" : "false")
              << ")\n";
  }
}

#endif  // GUTIBM_MPI

}  // namespace

int main(int argc, char** argv) {
#ifdef GUTIBM_MPI
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "=== CUDA-Aware MPI Reaction Tests ===\n";
  }

  test_device_reduce_gated_without_env();
  test_device_reduce_matches_host_when_available();
  test_runtime_detection_reports_status();

  if (rank == 0) {
    std::cout << "All CUDA-aware MPI reaction tests passed.\n";
  }

  MPI_Finalize();
#else
  (void)argc;
  (void)argv;
  std::cout << "MPI disabled at build time — skipping CUDA-aware MPI tests.\n";
#endif
  return 0;
}
