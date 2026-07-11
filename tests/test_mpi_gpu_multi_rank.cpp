/* -----------------------------------------------------------------------
   GutIBM – MPI + GPU integration (Spec 9 PR5 / issue #33)
   Run with: mpirun -np 2 test_mpi_gpu_multi_rank
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "dispatch.h"
#include "device.h"

#include <cassert>
#include <cmath>
#include <iostream>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

using namespace gutibm;

namespace {

SimulationConfig make_mpi_gpu_config() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.domain.ghost_width = 10e-6;
  cfg.domain.periodic = {false, true, false};
  cfg.time.total_time = 180.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 180.0;
  cfg.seed = 5150;
  cfg.hdf5.enabled = false;
  cfg.gpu.enabled = true;
  cfg.gpu.device_id = -1;
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.acetate.enabled = true;
  cfg.chem_env.mucin.enabled = true;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 30;
  s.mu_max = 5e-4;
  s.plasmids = {};
  s.conjugative = false;
  cfg.initial_strains.push_back(s);
  return cfg;
}

#ifdef GUTIBM_MPI

void require_two_ranks() {
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  assert(nprocs == 2);
}

Real chemistry_checksum(const Simulation& sim) {
  Real sum = 0.0;
  const auto& chem = sim.chemical_field();
  for (Int s = 0; s < chem.num_species(); ++s) {
    for (Int c = 0; c < chem.ncells(); ++c) {
      sum += chem.conc(s, c) * static_cast<Real>(c + 1 + s * chem.ncells());
    }
  }
  return sum;
}

void test_mpi_gpu_chemistry_identical_across_ranks() {
  require_two_ranks();

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifndef GUTIBM_CUDA
  if (rank == 0) {
    std::cout << "  test_mpi_gpu_chemistry_identical_across_ranks: SKIPPED (no CUDA)\n";
  }
  return;
#else
  if (DeviceContext::device_count() <= 0) {
    if (rank == 0) {
      std::cout << "  test_mpi_gpu_chemistry_identical_across_ranks: SKIPPED (no device)\n";
    }
    return;
  }

  SimulationConfig cfg = make_mpi_gpu_config();
  Simulation sim;
  sim.init(cfg);
  assert(sim.gpu_active());
  sim.run();

  Real local = chemistry_checksum(sim);
  Real minimum = 0.0;
  Real maximum = 0.0;
  MPI_Allreduce(&local, &minimum, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local, &maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  assert(std::abs(maximum - minimum) < 1e-9);

  if (rank == 0) {
    std::cout << "  test_mpi_gpu_chemistry_identical_across_ranks: PASSED\n";
  }
#endif
}

#endif  // GUTIBM_MPI

}  // namespace

int main() {
#ifdef GUTIBM_MPI
  MPI_Init(nullptr, nullptr);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout << "=== MPI GPU Multi-Rank Tests ===\n";
  test_mpi_gpu_chemistry_identical_across_ranks();
  if (rank == 0) {
    std::cout << "All MPI GPU multi-rank tests passed.\n";
  }
  MPI_Finalize();
  return 0;
#else
  std::cout << "=== MPI GPU Multi-Rank Tests ===\n";
  std::cout << "  SKIPPED (MPI not enabled)\n";
  std::cout << "All MPI GPU multi-rank tests passed.\n";
  return 0;
#endif
}
