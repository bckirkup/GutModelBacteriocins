#include "input_parser.h"
#include "mpi_test_helpers.h"
#include "simulation.h"
#include "species_names.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

using namespace gutibm;
using gutibm::test::make_mpi_config;
using gutibm::test::require_mpi_ranks;

namespace {

constexpr Real kCoupling = 1.0e-20;
constexpr Real kFieldTolerance = 1.0e-12;
constexpr Real kContrastMargin = 1.0e-8;
constexpr Int kAgentCount = 512;

struct ProbeResult {
  Int global_agents = 0;
  Int occupied_cells = 0;
  std::vector<Int> occupancy;
  std::vector<Real> carbon;
};

SimulationConfig make_density_config(Real coupling) {
  SimulationConfig cfg = make_mpi_config(73017, kAgentCount);
  cfg.initial_strains[0].count = 0;
  cfg.time.total_time = cfg.time.bio_dt;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.carbon_sink_km = 1.0e-4;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.vbf.agent_carbon_coupling = coupling;
  cfg.advection.radial_turnover = 1.0e30;
  cfg.advection.distal_transit_time = 1.0e30;
  cfg.disabled_fixes = {
      "metabolism", "quorum_sensing", "bacteriocin", "receptor",
      "motility", "conjugation", "cdi", "mutation", "mechanics"};
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.initial_conc = 1.0e-3;
      chemical.boundary_conc = 1.0e-3;
      chemical.diffusion_enabled = false;
    }
  }
  return cfg;
}

void populate_deterministic(Simulation& sim) {
  sim.agents().reserve(kAgentCount);
  const Vec3 lo = sim.domain().lo();
  for (Int index = 0; index < kAgentCount; ++index) {
    const Int ix = index % sim.domain().nx();
    const Int iy = (index / sim.domain().nx()) % sim.domain().ny();
    const Int iz = (index / (sim.domain().nx() * sim.domain().ny()))
        % sim.domain().nz();
    const Vec3 position = {
        lo[0] + (static_cast<Real>(ix) + 0.5) * sim.domain().dx_x(),
        lo[1] + (static_cast<Real>(iy) + 0.5) * sim.domain().dx_y(),
        lo[2] + (static_cast<Real>(iz) + 0.5) * sim.domain().dx_z()};
    if (!sim.domain().is_local(position)) continue;
    Agent agent = Agent::create_default(
        sim.agents().next_tag(), 1, position, 0.0);
    agent.identity.owner_rank = sim.domain().rank();
    sim.agents().push_back(std::move(agent));
  }
}

ProbeResult run_probe(Real coupling) {
  Simulation sim;
  sim.init(make_density_config(coupling));
  populate_deterministic(sim);
  const Int carbon = sim.chemical_field().find(species::CARBON);
  assert(carbon >= 0);
  sim.step(sim.config().time.bio_dt);

  ProbeResult result;
  result.global_agents = sim.global_agent_count();
  result.occupancy.assign(
      static_cast<size_t>(sim.domain().ncells()), 0);
  for (const Agent& agent : sim.agents()) {
    if (agent.state != PhenoState::DEAD && !agent.flags.is_ghost
        && agent.grid_cell >= 0
        && agent.grid_cell < sim.domain().ncells()) {
      ++result.occupancy[static_cast<size_t>(agent.grid_cell)];
    }
  }
#ifdef GUTIBM_MPI
  MPI_Allreduce(MPI_IN_PLACE, result.occupancy.data(),
                sim.domain().ncells(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  for (const Int count : result.occupancy) {
    if (count > 0) ++result.occupied_cells;
  }

  result.carbon.resize(static_cast<size_t>(sim.domain().ncells()));
  for (Int cell = 0; cell < sim.domain().ncells(); ++cell) {
    result.carbon[static_cast<size_t>(cell)] =
        sim.chemical_field().conc_global(carbon, cell);
  }
  return result;
}

void write_mpi_result(const std::string& path, const ProbeResult& result) {
  std::ofstream output(path);
  assert(output);
  output << std::setprecision(std::numeric_limits<Real>::max_digits10);
  output << result.global_agents << " " << result.occupied_cells << "\n";
  output << result.carbon.size() << "\n";
  for (const Real value : result.carbon) {
    output << value << "\n";
  }
}

ProbeResult read_mpi_result(const std::string& path) {
  std::ifstream input(path);
  assert(input);
  ProbeResult result;
  size_t cell_count = 0;
  input >> result.global_agents >> result.occupied_cells;
  input >> cell_count;
  result.carbon.resize(cell_count);
  for (Real& value : result.carbon) {
    input >> value;
  }
  assert(input);
  return result;
}

#ifdef GUTIBM_MPI

void run_replicated_mpi(const std::string& path) {
  require_mpi_ranks(2);
  const ProbeResult result = run_probe(kCoupling);

  std::vector<Real> minimum = result.carbon;
  std::vector<Real> maximum = result.carbon;
  MPI_Allreduce(MPI_IN_PLACE, minimum.data(),
                static_cast<int>(minimum.size()), MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, maximum.data(),
                static_cast<int>(maximum.size()), MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  Real rank_field_difference = 0.0;
  for (size_t cell = 0; cell < result.carbon.size(); ++cell) {
    rank_field_difference = std::max(
        rank_field_difference, maximum[cell] - minimum[cell]);
  }
  assert(rank_field_difference <= kFieldTolerance);
  assert(result.occupied_cells >= 8);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    write_mpi_result(path, result);
    std::cout << "  replicated MPI field consistency: max_rank_difference="
              << rank_field_difference
              << " occupied_cells=" << result.occupied_cells
              << " global_agents=" << result.global_agents << "\n";
  }
}

#endif

void run_serial(const std::string& path) {
  const ProbeResult expected = read_mpi_result(path);
  const ProbeResult coupled = run_probe(kCoupling);
  assert(coupled.global_agents == expected.global_agents);
  assert(coupled.occupied_cells == expected.occupied_cells);
  assert(coupled.carbon.size() == expected.carbon.size());

  Real max_serial_difference = 0.0;
  for (size_t cell = 0; cell < coupled.carbon.size(); ++cell) {
    max_serial_difference = std::max(
        max_serial_difference,
        std::abs(coupled.carbon[cell] - expected.carbon[cell]));
  }
  assert(max_serial_difference <= kFieldTolerance);

  const ProbeResult uncoupled = run_probe(0.0);
  Real max_occupied_contrast = 0.0;
  for (size_t cell = 0; cell < coupled.occupancy.size(); ++cell) {
    if (coupled.occupancy[cell] > 0) {
      max_occupied_contrast = std::max(
          max_occupied_contrast,
          std::abs(coupled.carbon[cell] - uncoupled.carbon[cell]));
    }
  }
  assert(max_occupied_contrast > kContrastMargin);
  std::cout << "  serial agreement: max_field_difference="
            << max_serial_difference
            << " max_occupied_contrast=" << max_occupied_contrast << "\n";
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GUTIBM_MPI
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
  return 77;
#endif

  assert(argc == 3);
  const std::string mode = argv[1];
  const std::string path = argv[2];
#ifdef GUTIBM_MPI
  if (mode == "mpi") {
    run_replicated_mpi(path);
  } else if (mode == "serial") {
    int ranks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    assert(ranks == 1);
    run_serial(path);
  } else {
    assert(false);
  }
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "VBF replicated density regression passed.\n";
  }
  MPI_Finalize();
#endif
  return 0;
}
