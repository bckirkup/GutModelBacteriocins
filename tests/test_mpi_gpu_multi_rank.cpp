/* -----------------------------------------------------------------------
   GutIBM – MPI + GPU integration (Spec 9 PR5 / issue #33)
   Run with: mpirun -np 2 test_mpi_gpu_multi_rank
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "dispatch.h"
#include "device.h"
#include "gpu_test_support.h"
#include "fix_metabolism.h"

#include <cassert>
#include <cmath>
#include <iostream>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

using namespace gutibm;

namespace {

[[maybe_unused]] SimulationConfig make_mpi_gpu_config() {
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
  cfg.cell_bio.fur.enabled = true;
  cfg.chem_env.siderophore.enabled = true;
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
  cfg.enabled_fixes = {"metabolism"};
  return cfg;
}

#ifdef GUTIBM_MPI

void require_two_ranks() {
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  assert(nprocs == 2);
}

[[maybe_unused]] Real chemistry_checksum(const Simulation& sim) {
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
  const int gpu_status = test::require_gpu("mpi_gpu_multi_rank");
  if (gpu_status != 0) return;

#ifdef GUTIBM_CUDA
  SimulationConfig cfg = make_mpi_gpu_config();
  Simulation sim;
  sim.init(cfg);
  assert(sim.gpu_active());
  bool positioned_boundary_agent = false;
  for (Agent& agent : sim.agents()) {
    if (agent.state == PhenoState::DEAD) continue;
    agent.x[0] = rank == 0
        ? sim.domain().local_hi_x() - 0.5 * sim.domain().ghost_width()
        : sim.domain().local_lo_x() + 0.5 * sim.domain().ghost_width();
    positioned_boundary_agent = true;
    break;
  }
  assert(positioned_boundary_agent);
  sim.run();
  if (const auto metabolism_gpu_steps =
          sim.agents_gpu().metabolism_gpu_steps();
      metabolism_gpu_steps <= 0) {
    const Real metabolism_abs_diff =
        std::abs(static_cast<Real>(metabolism_gpu_steps));
    std::cerr << "[gpu_diag][mpi_gpu_multi_rank] rank=" << rank
              << " measured_metabolism_gpu_steps=" << metabolism_gpu_steps
              << " reference=0 abs_diff=" << metabolism_abs_diff
              << " rel_diff=" << metabolism_abs_diff
              << " tolerance=0 (strictly positive required)"
              << " fur_enabled="
              << sim.config().cell_bio.fur.enabled
              << " siderophore_enabled="
              << sim.config().chem_env.siderophore.enabled << "\n";
  }
  assert(sim.agents_gpu().metabolism_gpu_steps() > 0);

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

void test_mpi_gpu_ghost_receptor_parity() {
  require_two_ranks();
  const SimulationConfig cfg = make_mpi_gpu_config();
  Simulation cpu;
  Simulation gpu;
  SimulationConfig cpu_cfg = cfg;
  cpu_cfg.gpu.enabled = false;
  cpu.init(cpu_cfg);
  gpu.init(cfg);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for (Simulation* simulation : {&cpu, &gpu}) {
    for (Agent& agent : simulation->agents()) {
      if (agent.state == PhenoState::DEAD) continue;
      agent.x[0] = rank == 0
          ? simulation->domain().local_hi_x()
              - 0.5 * simulation->domain().ghost_width()
          : simulation->domain().local_lo_x()
              + 0.5 * simulation->domain().ghost_width();
      break;
    }
    simulation->exchange_ghost_agents();
    FixMetabolism metabolism(*simulation, cfg.fixes.metabolism);
    metabolism.compute(60.0);
  }

  assert(cpu.agents().size() == gpu.agents().size());
  bool found_ghost = false;
  for (Int i = 0; i < cpu.agents().size(); ++i) {
    if (!cpu.agents()[i].flags.is_ghost) continue;
    found_ghost = true;
    for (int receptor = 0; receptor < NUM_RECEPTORS; ++receptor) {
      assert(std::abs(cpu.agents()[i].receptor_expr[receptor]
                      - gpu.agents()[i].receptor_expr[receptor]) < 1.0e-12);
    }
  }
  assert(found_ghost);
}

#endif  // GUTIBM_MPI

}  // namespace

int main() {
  const int gpu_status = test::require_gpu("mpi_gpu_multi_rank");
  if (gpu_status != 0) return gpu_status;
#ifdef GUTIBM_MPI
  MPI_Init(nullptr, nullptr);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout << "=== MPI GPU Multi-Rank Tests ===\n";
  test_mpi_gpu_chemistry_identical_across_ranks();
  test_mpi_gpu_ghost_receptor_parity();
  if (rank == 0) {
    std::cout << "All MPI GPU multi-rank tests passed.\n";
  }
  MPI_Finalize();
  return 0;
#else
  std::cout << "=== MPI GPU Multi-Rank Tests ===\n";
  return 77;
#endif
}
