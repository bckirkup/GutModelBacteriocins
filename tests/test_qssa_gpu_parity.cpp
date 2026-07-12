/* -----------------------------------------------------------------------
   QSSA near-field CPU vs GPU parity test
   ----------------------------------------------------------------------- */

#include "qssa_solver.h"
#include "greens_function_gpu.h"
#include "dispatch.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

#ifdef GUTIBM_CUDA
static Simulation make_qssa_sim(bool gpu_enabled, Real microcin_secretion,
                                uint64_t seed = 77) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.hdf5.enabled = false;
  cfg.gpu.enabled = gpu_enabled;
  cfg.domain.hi = {120.0e-6, 120.0e-6, 60.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.qssa.toxin_cutoff = 60.0e-6;
  cfg.qssa.use_fmm = false;
  cfg.qssa.microcin_secretion = microcin_secretion;
  cfg.fixes.bacteriocin.sos_basal_rate = 0.0;
  cfg.fixes.bacteriocin.sos_lysis_prob = 0.0;
  cfg.seed = seed;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 12;
  strain.mu_max = 5.5e-4;
  strain.plasmids = {"ColE1", "MccV"};
  strain.conjugative = true;
  cfg.initial_strains.push_back(strain);

  Simulation sim;
  sim.init(cfg);
  return sim;
}

static bool has_continuous_cira_source(const Simulation& sim) {
  if (sim.agents().size() == 0) return false;
  for (const auto& locus : sim.agents()[0].genome.bi_loci) {
    if (locus.target == ReceptorType::CirA
        && locus.release_mode == ReleaseMode::CONTINUOUS) {
      return true;
    }
  }
  return false;
}

static bool capture_cpu_field(Real microcin_secretion,
                              std::vector<Real>& values,
                              Real& total) {
  Simulation sim = make_qssa_sim(false, microcin_secretion, 77);
  if (!has_continuous_cira_source(sim)) return false;
  const Int idx = sim.chemical_field().find(species::BACTERIOCIN_CIRA);
  if (idx < 0) return false;

  sim.step(sim.config().time.bio_dt);
  const auto& chem = sim.chemical_field();
  values.resize(static_cast<size_t>(chem.ncells()));
  total = 0.0;
  for (Int c = 0; c < chem.ncells(); ++c) {
    const Real concentration = chem.conc(idx, c);
    values[static_cast<size_t>(c)] = concentration;
    total += concentration;
  }
  return std::isfinite(total) && total > 0.0;
}
#endif

int main() {
  std::cout << "=== QSSA GPU Near-Field Parity ===\n";

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  GpuConfig gcfg;
  gcfg.enabled = true;
  gcfg.device_id = 0;
  gpu_set_config(gcfg);

  if (!gpu_init_for_rank(0, 1)) {
    std::cout << "  SKIPPED (no CUDA device)\n";
    return 0;
  }

  constexpr Real kMicrocinSecretion = 1.0e-20;
  std::vector<Real> cpu_values;
  Real cpu_total = 0.0;
  if (!capture_cpu_field(kMicrocinSecretion, cpu_values, cpu_total)) {
    std::cout << "  FAIL: CPU continuous CirA reference missing\n";
    return 1;
  }

  std::vector<Real> stronger_cpu_values;
  Real stronger_cpu_total = 0.0;
  if (!capture_cpu_field(2.0 * kMicrocinSecretion, stronger_cpu_values,
                         stronger_cpu_total)
      || stronger_cpu_total <= cpu_total) {
    std::cout << "  FAIL: microcin secretion sensitivity missing\n";
    return 1;
  }

  Simulation sim_gpu = make_qssa_sim(true, kMicrocinSecretion, 77);
  if (!has_continuous_cira_source(sim_gpu)) {
    std::cout << "  FAIL: GPU continuous CirA source missing\n";
    return 1;
  }
  if (!sim_gpu.gpu_active() || !gpu_runtime_enabled()) {
    std::cout << "  FAIL: GPU dispatch inactive before GPU step\n";
    return 1;
  }
  const Int idx_gpu =
      sim_gpu.chemical_field().find(species::BACTERIOCIN_CIRA);
  if (idx_gpu < 0) {
    std::cout << "  FAIL: GPU bacteriocin species missing\n";
    return 1;
  }
  sim_gpu.step(sim_gpu.config().time.bio_dt);

  const auto& chem_gpu = sim_gpu.chemical_field();
  const Int ncells = chem_gpu.ncells();
  if (cpu_values.size() != static_cast<size_t>(ncells)) {
    std::cout << "  FAIL: CPU/GPU field sizes differ\n";
    return 1;
  }

  double max_abs = 0.0;
  double max_rel = 0.0;
  int compared = 0;
  for (Int c = 0; c < ncells; ++c) {
    const double cpu_v = cpu_values[static_cast<size_t>(c)];
    const double gpu_v = chem_gpu.conc(idx_gpu, c);
    if (cpu_v <= 0.0 && gpu_v <= 0.0) continue;
    ++compared;
    const double abs_err = std::abs(cpu_v - gpu_v);
    const double denom = std::max(std::abs(cpu_v), 1.0e-30);
    max_abs = std::max(max_abs, abs_err);
    max_rel = std::max(max_rel, abs_err / denom);
  }

  if (compared < 10) {
    std::cout << "  FAIL: insufficient nonzero cells (" << compared << ")\n";
    return 1;
  }

  constexpr double kAbsTol = 1.0e-12;
  constexpr double kRelTol = 1.0e-4;
  if (max_abs > kAbsTol && max_rel > kRelTol) {
    std::cout << "  FAIL: max_abs=" << max_abs << " max_rel=" << max_rel << "\n";
    return 1;
  }

  std::cout << "  test_qssa_gpu_near_field_parity: PASSED"
            << " (cells=" << compared << " max_rel=" << max_rel << ")\n";
  return 0;
#endif
}
