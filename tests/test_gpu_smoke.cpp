/* -----------------------------------------------------------------------
   GutIBM – GPU smoke integration test
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "dispatch.h"
#include "device.h"
#include "gpu_diagnostic_format.h"
#include "gpu_test_support.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;
using gutibm::gpu_diagnostic::format_real;

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

static void print_species_diagnostics(const Simulation& gpu,
                                      const std::vector<Real>& cpu_fp,
                                      const std::vector<Real>& gpu_fp) {
  const auto& chem = gpu.chemical_field();
  std::cerr << "[gpu_diag][gpu_smoke][species]\n";
  for (Int species = 0; species < chem.num_species(); ++species) {
    const Real absolute_difference =
        std::abs(cpu_fp[static_cast<size_t>(species)]
                 - gpu_fp[static_cast<size_t>(species)]);
    const Real scale = std::max(
        std::abs(cpu_fp[static_cast<size_t>(species)]), 1.0e-30);
    const Real relative_difference = absolute_difference / scale;
    std::cerr << "  name=" << chem.spec(species).name
              << " cpu_fingerprint="
              << format_real(cpu_fp[static_cast<size_t>(species)])
              << " gpu_fingerprint="
              << format_real(gpu_fp[static_cast<size_t>(species)])
              << " abs_diff=" << format_real(absolute_difference)
              << " rel_diff=" << format_real(relative_difference) << "\n";
  }
  std::cerr << "[gpu_diag][gpu_smoke][dispatch]"
            << " metabolism_gpu_steps="
            << gpu.agents_gpu().metabolism_gpu_steps()
            << " gpu_metabolism_active="
            << gpu.gpu_metabolism_active()
            << " chemistry_steps=unavailable"
            << " reaction_application_steps=unavailable"
            << " diffusion_steps=unavailable"
            << " qssa_superposition_steps=unavailable"
            << " mechanics_steps=unavailable\n";
}
#endif

int main() {
  std::cout << "=== GPU Smoke Test ===\n";
  const int gpu_status = test::require_gpu("gpu_smoke");
  if (gpu_status != 0) return gpu_status;

  SimulationConfig cfg = InputParser::default_config();
  cfg.time.total_time = 300.0;
  cfg.time.bio_dt = 60.0;
  cfg.hdf5.enabled = false;
  // Keep this parity fixture on the GPU metabolism path; siderophore chemistry
  // is CPU-authoritative and intentionally tested by a separate fallback case.
  cfg.chem_env.siderophore.enabled = false;
  cfg.initial_strains.clear();
  cfg.domain.hi = {200.0e-6, 200.0e-6, 100.0e-6};
  cfg.domain.grid_dx = 4.0e-6;

  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 20;
  s.mu_max = 5e-4;
  s.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(s);

  // Default Fur upregulation makes iron receptors differ from BtuB. GPU
  // metabolism must read receptor_expr as SoA (r * n + i) or agent_rel blows up
  // while chem_rel stays ~1e-9 (see test_gpu_receptor_layout).

#ifndef GUTIBM_CUDA
  std::cout << "  test_gpu_smoke: SKIPPED (CUDA not compiled in)\n";
  std::cout << "All GPU smoke tests passed.\n";
  return 0;
#else
  // CPU baseline
  cfg.gpu.enabled = false;
  Real fp_cpu = 0.0;
  std::vector<Real> chem_fp_cpu;
  {
    Simulation sim_cpu;
    sim_cpu.init(cfg);
    if (sim_cpu.domain().ncells() != 62500) {
      std::cerr << "  test_gpu_smoke: FAILED (expected 62500 grid cells, got "
                << sim_cpu.domain().ncells() << ")\n";
      return 1;
    }
    sim_cpu.run();
    fp_cpu = fingerprint(sim_cpu);
    chem_fp_cpu = chemical_fingerprints(sim_cpu);
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
  Real chem_abs_diff = 0.0;
  Real chem_scale = 1.0e-30;
  for (size_t species = 0; species < chem_fp_cpu.size(); ++species) {
    const Real species_abs_diff =
        std::abs(chem_fp_cpu[species] - chem_fp_gpu[species]);
    const Real species_scale = std::max(std::abs(chem_fp_cpu[species]),
                                        1.0e-30);
    const Real species_rel = species_abs_diff / species_scale;
    chem_abs_diff = std::max(chem_abs_diff, species_abs_diff);
    chem_scale = std::max(chem_scale, species_scale);
    chem_rel = std::max(chem_rel, species_rel);
  }
  if (rel > 0.05 || chem_rel > 0.05) {
    std::cerr << "[gpu_diag][gpu_smoke] fp_cpu=" << format_real(fp_cpu)
              << " fp_gpu=" << format_real(fp_gpu)
              << " agent_abs_diff=" << format_real(std::abs(fp_cpu - fp_gpu))
              << " agent_rel=" << format_real(rel)
              << " agent_tolerance=0.05"
              << " chem_abs_diff=" << format_real(chem_abs_diff)
              << " chem_scale=" << format_real(chem_scale)
              << " chem_rel=" << format_real(chem_rel)
              << " chem_tolerance=0.05\n";
    print_species_diagnostics(sim_gpu, chem_fp_cpu, chem_fp_gpu);
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
