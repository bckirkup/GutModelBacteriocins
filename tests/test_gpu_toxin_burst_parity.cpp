/* -----------------------------------------------------------------------
   GutIBM – deterministic GPU parity for toxin bursts

   The burst source is injected through Simulation::add_toxin_burst(), so
   these checks do not consume stochastic SOS or receptor RNG draws.
   ----------------------------------------------------------------------- */

#include "device.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

namespace {

struct BurstRun {
  std::vector<Real> field;
  Real maximum = 0.0;
  Real mean = 0.0;
  Int colicin_kills = 0;
};

SimulationConfig make_config(bool gpu_enabled) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {100.0e-6, 100.0e-6, 50.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.time.total_time = 60.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 60.0;
  cfg.seed = 4701;
  cfg.hdf5.enabled = false;
  cfg.gpu.enabled = gpu_enabled;
  cfg.gpu.device_id = 0;
  cfg.qssa.toxin_cutoff = 200.0e-6;
  cfg.qssa.use_fmm = false;
  cfg.chem_env.siderophore.enabled = false;
  for (auto& chemical : cfg.chemicals) {
    if (chemical.name == species::B12) {
      chemical.initial_conc = 0.0;
    }
  }
  cfg.fixes.bacteriocin.sos_basal_rate = 0.0;
  cfg.fixes.bacteriocin.sos_lysis_prob = 0.0;
  cfg.fixes.receptor.kill_rate_colicin = 1.0e6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain target;
  target.type = 2;
  target.count = 32;
  target.mu_max = 0.0;
  target.plasmids = {};
  target.conjugative = false;
  cfg.initial_strains.push_back(target);
  return cfg;
}

void add_burst(Simulation& sim, Real strength) {
  ToxinBurstSource burst;
  burst.pos = {
      0.5 * (sim.domain().lo()[0] + sim.domain().hi()[0]),
      0.5 * (sim.domain().lo()[1] + sim.domain().hi()[1]),
      0.5 * (sim.domain().lo()[2] + sim.domain().hi()[2])};
  burst.params.diff_coeff = 4.0e-11;
  burst.params.retardation = 1.0;
  burst.params.pI = 9.0;
  burst.params.source_rate = strength / burst.release_tau;
  burst.params.decay_rate = 0.0;
  burst.creation_time = 0.0;
  burst.target = ReceptorType::BtuB;
  sim.add_toxin_burst(burst);
}

BurstRun run_case(bool gpu_enabled, Real strength) {
  Simulation sim;
  SimulationConfig cfg = make_config(gpu_enabled);
  sim.init(cfg);
  add_burst(sim, strength);
  sim.step(cfg.time.bio_dt);

  const Int toxin = sim.chemical_field().find(species::BACTERIOCIN_BTUB);
  assert(toxin >= 0);
  BurstRun result;
  result.field.resize(static_cast<size_t>(sim.chemical_field().ncells()));
  for (Int cell = 0; cell < sim.chemical_field().ncells(); ++cell) {
    const Real value = sim.chemical_field().conc(toxin, cell);
    assert(std::isfinite(value));
    assert(value >= 0.0);
    result.field[static_cast<size_t>(cell)] = value;
    result.maximum = std::max(result.maximum, value);
    result.mean += value;
  }
  result.mean /= static_cast<Real>(result.field.size());
  result.colicin_kills = sim.step_events().mortality_colicin;
  return result;
}

#ifdef GUTIBM_CUDA
void compare_fields(const BurstRun& cpu, const BurstRun& gpu) {
  assert(cpu.field.size() == gpu.field.size());
  const Real max_scale = std::max(cpu.maximum, 1.0e-30);
  const Real mean_scale = std::max(cpu.mean, 1.0e-30);
  assert(std::abs(cpu.maximum - gpu.maximum) / max_scale < 1.0e-10);
  assert(std::abs(cpu.mean - gpu.mean) / mean_scale < 1.0e-10);

  Real max_relative = 0.0;
  for (size_t cell = 0; cell < cpu.field.size(); ++cell) {
    const Real scale = std::max(std::abs(cpu.field[cell]), 1.0e-30);
    max_relative = std::max(max_relative,
                            std::abs(cpu.field[cell] - gpu.field[cell]) / scale);
  }
  assert(max_relative < 1.0e-10);
  std::cout << "  field parity: max=" << cpu.maximum
            << " mean=" << cpu.mean
            << " max_relative=" << max_relative << "\n";
}
#endif

}  // namespace

int main() {
  std::cout << "=== GPU Toxin Burst Parity ===\n";

#ifndef GUTIBM_CUDA
  const BurstRun cpu = run_case(false, 1.0e-18);
  assert(cpu.maximum > 0.0);
  std::cout << "  CPU reference: max=" << cpu.maximum
            << " kills=" << cpu.colicin_kills << "\n";
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  if (DeviceContext::device_count() <= 0) {
    const BurstRun cpu = run_case(false, 1.0e-18);
    assert(cpu.maximum > 0.0);
    std::cout << "  CPU reference: max=" << cpu.maximum
              << " kills=" << cpu.colicin_kills << "\n";
    std::cout << "  SKIPPED (no CUDA device)\n";
    return 0;
  }

  const Real strength = 1.0e-18;
  const BurstRun cpu = run_case(false, strength);
  const BurstRun gpu = run_case(true, strength);
  compare_fields(cpu, gpu);

  const std::vector<Real> strengths{0.0, 0.1 * strength, strength, 10.0 * strength};
  std::vector<Int> cpu_kills;
  std::vector<Int> gpu_kills;
  for (const Real dose : strengths) {
    cpu_kills.push_back(run_case(false, dose).colicin_kills);
    gpu_kills.push_back(run_case(true, dose).colicin_kills);
  }
  assert(std::is_sorted(cpu_kills.begin(), cpu_kills.end()));
  assert(std::is_sorted(gpu_kills.begin(), gpu_kills.end()));
  constexpr Int target_count = 32;
  assert(cpu_kills.front() == 0);
  assert(gpu_kills.front() == 0);
  assert(cpu_kills.back() == target_count);
  assert(gpu_kills.back() == target_count);
  for (size_t i = 1; i + 1 < cpu_kills.size(); ++i) {
    // ULP-level reduction-order differences can flip one borderline stochastic
    // kill, so intermediate CPU/GPU counts need not match exactly.
    assert(std::abs(cpu_kills[i] - gpu_kills[i]) <= 1);
  }
  std::cout << "  deterministic kill parity: " << cpu_kills[0] << ", "
            << cpu_kills[1] << ", " << cpu_kills[2] << ", " << cpu_kills[3]
            << "\n";
  std::cout << "All GPU toxin burst parity tests passed.\n";
  return 0;
#endif
}
