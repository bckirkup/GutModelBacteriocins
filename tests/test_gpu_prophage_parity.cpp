/* -----------------------------------------------------------------------
   GutIBM – CPU vs GPU parity test for Spec 14 temperate prophage module
   ----------------------------------------------------------------------- */

#include "fix_prophage.h"
#include "gpu_test_support.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

namespace {

SimulationConfig make_parity_config(bool gpu_enabled) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {60.0e-6, 60.0e-6, 40.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.time.bio_dt = 60.0;
  cfg.time.total_time = 600.0;
  cfg.time.output_interval = 600.0;
  cfg.seed = 8401;
  cfg.hdf5.enabled = false;
  cfg.dysbiosis_threshold = 0.0;
  cfg.gpu.enabled = gpu_enabled;
  cfg.gpu.device_id = 0;

  cfg.advection.radial_turnover = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.advection.peristaltic_enabled = false;
  cfg.advection.crypts_enabled = false;
  cfg.fixes.bacteriocin.sos_basal_rate = 0.0;
  cfg.fixes.bacteriocin.sos_lysis_prob = 0.0;
  cfg.fixes.bacteriocin.sos_cross_induction_rate = 0.0;
  cfg.enabled_fixes = {"prophage"};

  cfg.prophage.enabled = true;
  cfg.prophage.k_spont = 1.0e-4;
  cfg.prophage.k_max_damage = 5.0e-3;
  cfg.prophage.K_damage = 1.0e-6;
  cfg.prophage.hill_n = 2.0;
  cfg.prophage.damage_stimulus = 1.0e-6;
  cfg.prophage.virion_burst_size = 15.0;
  cfg.prophage.D_virion = 2.0e-11;
  cfg.prophage.virion_retardation = 1.0;
  cfg.prophage.k_virion_decay = 5.0e-4;
  cfg.prophage.virion_washout_scale = 1.0;
  cfg.prophage.k_adsorb = 2.0e-14;
  cfg.prophage.p_lysogeny = 0.4;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain active;
  active.type = 1;
  active.count = 60;
  active.mu_max = 0.0;
  active.lysogen_prevalence = 1.0;
  active.active_inducible_fraction = 1.0;

  SimulationConfig::InitialStrain cryptic;
  cryptic.type = 2;
  cryptic.count = 40;
  cryptic.mu_max = 0.0;
  cryptic.lysogen_prevalence = 1.0;
  cryptic.active_inducible_fraction = 0.0;

  SimulationConfig::InitialStrain susceptible;
  susceptible.type = 3;
  susceptible.count = 100;
  susceptible.mu_max = 0.0;
  susceptible.lysogen_prevalence = 0.0;
  susceptible.active_inducible_fraction = 0.0;

  cfg.initial_strains = {active, cryptic, susceptible};
  return cfg;
}

StepEvents cumulative_total(const Simulation& sim) {
  StepEvents total = sim.cumulative_events();
  total.add(sim.step_events());
  return total;
}

}  // namespace

int main() {
  if (const int skip = test::require_gpu("gpu_prophage_parity"); skip != 0) {
    return skip;
  }

  std::cout << "=== GPU Prophage Parity Tests ===\n";

  Simulation cpu_sim;
  cpu_sim.init(make_parity_config(false));
  assert(!cpu_sim.gpu_active());

  Simulation gpu_sim;
  gpu_sim.init(make_parity_config(true));
  assert(gpu_sim.gpu_active());

  for (int step = 0; step < 10; ++step) {
    cpu_sim.step(60.0);
    gpu_sim.step(60.0);
  }

  const StepEvents cpu_ev = cumulative_total(cpu_sim);
  const StepEvents gpu_ev = cumulative_total(gpu_sim);

  // Exact event count parity (same RNG stream, identical induction/adsorption math)
  assert(cpu_ev.prophage_inductions > 0);
  assert(cpu_ev.prophage_inductions == gpu_ev.prophage_inductions);
  assert(cpu_ev.prophage_lytic_kills == gpu_ev.prophage_lytic_kills);
  assert(cpu_ev.prophage_lysogen_conversions == gpu_ev.prophage_lysogen_conversions);
  assert(cpu_sim.agents().size() == gpu_sim.agents().size());

  // Virion ledger parity & conservation closure on GPU
  const Real cpu_v = FixProphage::total_virions_in_field(cpu_sim);
  const Real gpu_v = FixProphage::total_virions_in_field(gpu_sim);
  const Real v_scale = std::max(1.0, cpu_ev.virions_released);
  assert(std::abs(cpu_v - gpu_v) / v_scale < 1.0e-9);
  assert(std::abs(cpu_ev.virions_released - gpu_ev.virions_released) / v_scale < 1.0e-9);
  assert(std::abs(cpu_ev.virions_adsorbed - gpu_ev.virions_adsorbed) / v_scale < 1.0e-9);
  assert(std::abs(cpu_ev.virions_decayed - gpu_ev.virions_decayed) / v_scale < 1.0e-9);

  const Real gpu_lhs = gpu_v + gpu_ev.virions_adsorbed + gpu_ev.virions_decayed;
  const Real gpu_rhs = gpu_ev.virions_released;
  assert(std::abs(gpu_lhs - gpu_rhs) / std::max(1.0, gpu_rhs) < 1.0e-9);

  // Cryptic carriers must remain 100% immune and unlysed on GPU
  Int gpu_cryptic = 0;
  for (const Agent& a : gpu_sim.agents()) {
    if (a.identity.type == 2) {
      assert(a.flags.has_prophage_immunity);
      assert(!a.flags.has_inducible_prophage);
      ++gpu_cryptic;
    }
  }
  assert(gpu_cryptic == 40);

  std::cout << "  CPU vs GPU prophage parity: PASSED (inductions="
            << gpu_ev.prophage_inductions << ", conversions="
            << gpu_ev.prophage_lysogen_conversions << ", lytic_kills="
            << gpu_ev.prophage_lytic_kills << ")\n";
  return 0;
}
