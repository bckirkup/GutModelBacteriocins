#include "device.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"
#include "gpu_test_support.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace gutibm;

namespace {

struct FeedbackResult {
  Real initial_mass = 0.0;
  Real final_mass = 0.0;
  Real uptake = 0.0;
  Real boundary = 0.0;
  Real vbf_source = 0.0;
  Real vbf_sink = 0.0;
  Real clipped = 0.0;
  std::vector<Real> first_field;
  std::vector<Real> second_field;
};

SimulationConfig feedback_config(Int count, bool gpu) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.gpu.enabled = gpu;
  cfg.gpu.device_id = 0;
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;
  cfg.time.bio_dt = 60.0;
  cfg.time.total_time = 120.0;
  cfg.time.output_interval = 120.0;
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {40e-6, 40e-6, 40e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.chemistry_decomposition = "replicated";
  cfg.chem_env.siderophore.enabled = true;
  cfg.cell_bio.fur.enabled = true;
  cfg.vbf.nutrient_sink = 0.0;
  cfg.vbf.mucin_liberation = 0.0;
  cfg.vbf.use_dynamic_mucin = false;
  cfg.initial_strains.clear();

  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = count;
  strain.mu_max = 5.5e-4;
  strain.plasmids = {"ColE1", "ColB"};
  strain.conjugative = true;
  cfg.initial_strains.push_back(strain);
  return cfg;
}

Real carbon_mass(const Simulation& sim, Int carbon) {
  const auto& chem = sim.chemical_field();
  Real total = 0.0;
  for (const Real value : chem.conc_data()[static_cast<size_t>(carbon)]) {
    total += value;
  }
  return total * sim.domain().cell_volume();
}

std::vector<Real> carbon_field(const Simulation& sim, Int carbon) {
  const auto& values =
      sim.chemical_field().conc_data()[static_cast<size_t>(carbon)];
  return values;
}

FeedbackResult run_feedback(Int count, bool gpu) {
  // Fur and siderophore chemistry force CPU metabolism while oxygen depletion
  // remains eligible for the GPU, exercising the host-reaction/device-kernel
  // handoff that previously discarded carbon and iron reactions.
  Simulation sim;
  sim.init(feedback_config(count, gpu));
  assert(sim.gpu_active() == gpu);
  assert(sim.config().cell_bio.fur.enabled);
  assert(sim.config().chem_env.siderophore.enabled);
  const Int carbon = sim.chemical_field().find(species::CARBON);
  assert(carbon >= 0);

  FeedbackResult result;
  result.initial_mass = carbon_mass(sim, carbon);
  sim.step(60.0);
  result.first_field = carbon_field(sim, carbon);
  sim.step(60.0);
  result.second_field = carbon_field(sim, carbon);
  result.final_mass = carbon_mass(sim, carbon);

  auto& accounting = sim.chemical_field().flux_accounting();
  accounting.close_interval();
  const auto index = static_cast<size_t>(carbon);
  result.uptake = accounting.agent_uptake_cumulative[index];
  result.boundary = accounting.boundary_cumulative[index];
  result.vbf_source = accounting.vbf_source_cumulative[index];
  result.vbf_sink = accounting.vbf_sink_cumulative[index];
  result.clipped = accounting.reaction_clip_cumulative[index];
  return result;
}

void assert_closure(const FeedbackResult& result, const std::string& path) {
  const Real expected = result.initial_mass + result.boundary
      + result.vbf_source - result.uptake - result.vbf_sink + result.clipped;
  const Real scale = std::max({std::abs(result.initial_mass),
                               std::abs(result.final_mass),
                               std::abs(result.uptake), 1.0e-18});
  const Real tolerance = std::max(2.0e-4 * scale, 2.0e-19);
  assert(std::abs(expected - result.final_mass) <= tolerance);
  const Real realized_uptake = result.initial_mass + result.boundary
      + result.vbf_source - result.vbf_sink + result.clipped
      - result.final_mass;
  assert(std::abs(realized_uptake - result.uptake) <= tolerance);
  if (result.uptake > 0.0) {
    assert(result.first_field != result.second_field);
  }
  std::cout << "  " << path << " closure residual="
            << expected - result.final_mass
            << " uptake=" << result.uptake << "\n";
}

void run_path_checks(bool gpu) {
  const std::vector counts = {0, 4, 8};
  std::vector<Real> depletion;
  for (const Int count : counts) {
    const FeedbackResult result = run_feedback(count, gpu);
    assert_closure(result, gpu ? "GPU" : "CPU");
    if (count > 0) {
      assert(result.uptake > 0.0);
    }
    const FeedbackResult control = run_feedback(0, gpu);
    depletion.push_back(control.final_mass - result.final_mass);
  }
  assert(depletion[1] > depletion[0]);
  assert(depletion[2] > depletion[1]);
}
}  // namespace

int main() {
  run_path_checks(false);
  const int gpu_status = test::require_gpu("gpu_nutrient_feedback");
  if (gpu_status != 0) return gpu_status;
#ifndef GUTIBM_CUDA
  std::cout << "GPU nutrient feedback checks passed (CPU only; CUDA not compiled in).\n";
  return 0;
#else
  run_path_checks(true);
  std::cout << "GPU nutrient feedback tests passed.\n";
  return 0;
#endif
}
