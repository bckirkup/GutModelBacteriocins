/* -----------------------------------------------------------------------
   GutIBM – GPU smoke integration test
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "dispatch.h"
#include "device.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

#ifdef GUTIBM_CUDA
static Real fingerprint(const Simulation& sim) {
  Real fp = 0.0;
  const auto& agents = sim.agents();
  for (const Agent& a : agents) {
    if (a.state == PhenoState::DEAD) continue;
    fp += a.biomass + a.mu_realized * 1e3 + a.x[0] * 1e6;
  }
  fp += static_cast<Real>(sim.global_agent_count());
  return fp;
}

static std::vector<Real> chemical_fingerprints(const Simulation& sim) {
  const auto& chem = sim.chemical_field();
  std::vector<Real> fingerprints;
  fingerprints.reserve(static_cast<size_t>(chem.num_species()));
  const auto cells = static_cast<Real>(chem.ncells());
  for (const auto& species : chem.conc_data()) {
    Real fp = 0.0;
    for (Int cell = 0; cell < chem.ncells(); ++cell) {
      fp += species[static_cast<size_t>(cell)]
          * static_cast<Real>(cell + 1);
    }
    fingerprints.push_back(fp / (cells * cells));
  }
  return fingerprints;
}
#endif

int main() {
  std::cout << "=== GPU Smoke Test ===\n";

  SimulationConfig cfg = InputParser::default_config();
  cfg.time.total_time = 300.0;
  cfg.time.bio_dt = 60.0;
  cfg.hdf5.enabled = false;
  cfg.initial_strains.clear();

  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 20;
  s.mu_max = 5e-4;
  s.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(s);

#ifndef GUTIBM_CUDA
  std::cout << "  test_gpu_smoke: SKIPPED (CUDA not compiled in)\n";
  std::cout << "All GPU smoke tests passed.\n";
  return 0;
#else
  // CPU baseline
  cfg.gpu.enabled = false;
  Simulation sim_cpu;
  sim_cpu.init(cfg);
  sim_cpu.run();
  Real fp_cpu = fingerprint(sim_cpu);
  const auto chem_fp_cpu = chemical_fingerprints(sim_cpu);

  if (DeviceContext::device_count() <= 0) {
    std::cout << "  test_gpu_smoke: SKIPPED (no CUDA device)\n";
    std::cout << "All GPU smoke tests passed.\n";
    return 0;
  }

  cfg.gpu.enabled = true;
  cfg.gpu.device_id = 0;
  Simulation sim_gpu;
  sim_gpu.init(cfg);
  if (!sim_gpu.gpu_active()) {
    std::cerr << "  test_gpu_smoke: FAILED (GPU init failed";
#ifdef GUTIBM_CUDA
    std::cerr << ": " << DeviceContext::last_error();
#endif
    std::cerr << ")\n";
    return 1;
  }
  sim_gpu.run();
  Real fp_gpu = fingerprint(sim_gpu);
  const auto chem_fp_gpu = chemical_fingerprints(sim_gpu);

  Real rel = std::abs(fp_cpu - fp_gpu) / std::max(std::abs(fp_cpu), 1e-30);
  Real chem_rel = 0.0;
  for (size_t species = 0; species < chem_fp_cpu.size(); ++species) {
    const Real species_rel = std::abs(chem_fp_cpu[species] - chem_fp_gpu[species])
        / std::max(std::abs(chem_fp_cpu[species]), 1e-30);
    chem_rel = std::max(chem_rel, species_rel);
  }
  if (rel > 0.05 || chem_rel > 0.05) {
    std::cerr << "  test_gpu_smoke: FAILED (agent_rel=" << rel
              << " chem_rel=" << chem_rel << ")\n";
    return 1;
  }

  std::cout << "  test_gpu_smoke: PASSED (agent_rel=" << rel
            << " chem_rel=" << chem_rel << ")\n";
  std::cout << "All GPU smoke tests passed.\n";
  return 0;
#endif
}
