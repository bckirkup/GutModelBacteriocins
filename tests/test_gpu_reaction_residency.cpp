#include "gpu_test_support.h"
#include "input_parser.h"
#include "species_names.h"
#include "simulation.h"
#include "gpu_profile.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <unistd.h>
#include <vector>

using namespace gutibm;

namespace {

#ifdef GUTIBM_CUDA
GpuTransferSiteProfile find_site(
    const std::vector<GpuTransferSiteProfile>& profiles,
    const std::string& label) {
  for (const auto& profile : profiles) {
    if (profile.label == label) return profile;
  }
  // A site that issued no transfer is never registered, which is a legitimate
  // outcome for the resident mode rather than a test failure.
  GpuTransferSiteProfile empty;
  empty.label = label;
  return empty;
}

SimulationConfig reaction_config(uint64_t seed) {
  SimulationConfig config = InputParser::default_config();
  config.seed = seed;
  config.time.total_time = 180.0;
  config.time.bio_dt = 60.0;
  config.time.output_interval = 180.0;
  config.hdf5.enabled = false;
  config.gpu.enabled = true;
  config.quorum_sensing.enabled = true;
  config.chem_env.siderophore.enabled = true;
  config.chem_env.oxygen.enabled = true;
  config.domain.hi = {80.0e-6, 20.0e-6, 20.0e-6};
  config.domain.grid_dx = 5.0e-6;
  config.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 12;
  strain.mu_max = 5.0e-4;
  strain.plasmids = {"ColE1"};
  config.initial_strains.push_back(strain);
  return config;
}

void run_child(const char* executable, const char* mode, const char* output) {
  const std::string command = "GUTIBM_GPU_REACTION_RESIDENCY=" + std::string(mode)
      + " " + executable + " --child " + output;
  assert(std::system(command.c_str()) == 0);
}

void write_child_result(const char* output) {
  Simulation simulation;
  gpu_transfer_profile_set_enabled(true);
  simulation.init(reaction_config(9137));
  gpu_transfer_profile_reset();
  std::vector<Real> oxygen_ledger;
  const Int oxygen = simulation.chemical_field().find(species::OXYGEN);
  assert(oxygen >= 0);
  for (int step = 0; step < 3; ++step) {
    simulation.step(simulation.config().time.bio_dt);
    const auto& flux = simulation.chemical_field().flux_accounting();
    oxygen_ledger.push_back(
        flux.vbf_sink_for_step(oxygen)
        + flux.agent_uptake_for_step(oxygen)
        + flux.reaction_clip_for_step(oxygen));
  }
  const auto chem_reactions = find_site(
      gpu_transfer_site_profiles(), "chem_reactions");
  const auto full_field_bytes = static_cast<unsigned long long>(
      simulation.chemical_field().num_species())
      * static_cast<unsigned long long>(simulation.chemical_field().ncells())
      * sizeof(Real);
  std::ofstream stream(output);
  assert(stream.good());
  stream << chem_reactions.h2d_bytes << '\n'
         << chem_reactions.d2h_bytes << '\n'
         << full_field_bytes << '\n'
         << std::scientific << std::setprecision(17);
  const auto& field = simulation.chemical_field();
  for (const auto& row : field.conc_data()) {
    for (const Real value : row) stream << value << '\n';
  }
  for (const Real value : oxygen_ledger) stream << value << '\n';
}

void compare_results(const char* lhs_path, const char* rhs_path) {
  std::ifstream lhs(lhs_path);
  std::ifstream rhs(rhs_path);
  assert(lhs.good() && rhs.good());
  unsigned long long lhs_h2d = 0;
  unsigned long long lhs_d2h = 0;
  unsigned long long rhs_h2d = 0;
  unsigned long long rhs_d2h = 0;
  unsigned long long lhs_full_field = 0;
  unsigned long long rhs_full_field = 0;
  assert(lhs >> lhs_h2d >> lhs_d2h >> lhs_full_field);
  assert(rhs >> rhs_h2d >> rhs_d2h >> rhs_full_field);
  const unsigned long long lhs_bytes = lhs_h2d + lhs_d2h;
  const unsigned long long rhs_bytes = rhs_h2d + rhs_d2h;
  Real lhs_value = 0.0;
  Real rhs_value = 0.0;
  while (lhs >> lhs_value && rhs >> rhs_value) {
    const Real scale = std::max({1.0, std::abs(lhs_value), std::abs(rhs_value)});
    assert(std::isfinite(lhs_value) && std::isfinite(rhs_value));
    assert(std::abs(lhs_value - rhs_value) <= 1.0e-12 * scale);
  }
  assert(lhs.eof() && rhs.eof());
  assert(lhs_full_field == rhs_full_field);
  assert(rhs_bytes >= rhs_full_field);
  assert(lhs_bytes * 2ULL < rhs_bytes);
}
#endif

}  // namespace

int main(int argc, char** argv) {
#ifndef GUTIBM_CUDA
  (void)argc;
  (void)argv;
#endif
#ifdef GUTIBM_CUDA
  if (argc == 3 && std::string(argv[1]) == "--child") {
    write_child_result(argv[2]);
    return 0;
  }
#endif
  if (const int gpu_status = test::require_gpu("gpu_reaction_residency");
      gpu_status != 0) {
    return gpu_status;
  }
#ifndef GUTIBM_CUDA
  std::cout << "SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  const std::string base = "gpu_reaction_residency_" + std::to_string(
      static_cast<long long>(::getpid()));
  const std::string resident = base + "_resident.txt";
  const std::string legacy = base + "_legacy.txt";
  run_child(argv[0], "1", resident.c_str());
  run_child(argv[0], "0", legacy.c_str());
  compare_results(resident.c_str(), legacy.c_str());
  std::remove(resident.c_str());
  std::remove(legacy.c_str());
  std::cout << "GPU reaction residency equivalence passed.\n";
  return 0;
#endif
}
