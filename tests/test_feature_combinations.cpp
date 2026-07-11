/* -----------------------------------------------------------------------
   GutIBM – Feature-combination integration tests (Spec 8)

   Short full Simulation::init() → Simulation::run() exercises that turn on
   multiple optional modules together and assert biologically directional or
   stability outcomes. Complements per-mechanism unit tests.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "species_names.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

namespace {

SimulationConfig make_combo_config(uint64_t seed) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = seed;
  cfg.time.total_time = 300.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 300.0;
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;

  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {80e-6, 80e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;

  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length = 80e-6;
  cfg.qssa.toxin_cutoff = 40e-6;
  cfg.qssa.nutrient_cutoff = 20e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = 20;
  resident.mu_max = 5.5e-4;
  resident.plasmids = {"ColE1", "ColB"};
  resident.conjugative = true;
  cfg.initial_strains.push_back(resident);

  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = 10;
  immigrant.mu_max = 5.0e-4;
  immigrant.plasmids = {};
  immigrant.conjugative = false;
  cfg.initial_strains.push_back(immigrant);

  return cfg;
}

Simulation run_combo(const SimulationConfig& cfg) {
  Simulation sim;
  sim.init(cfg);
  sim.run();
  return sim;
}

void assert_chemistry_sane(const Simulation& sim) {
  const auto& chem = sim.chemical_field();
  for (Int s = 0; s < chem.num_species(); ++s) {
    for (Int c = 0; c < chem.ncells(); ++c) {
      const Real val = chem.conc(s, c);
      assert(std::isfinite(val));
      assert(val >= 0.0);
    }
  }
}

void assert_population_sane(const Simulation& sim, Int min_expected = 1) {
  Int live = 0;
  for (const Agent& a : sim.agents()) {
    if (a.state == PhenoState::DEAD) continue;
    assert(std::isfinite(a.biomass));
    assert(a.biomass > 0.0);
    assert(std::isfinite(a.mu_realized));
    ++live;
  }
  assert(live >= min_expected);
}

void zero_acetate_pool(SimulationConfig& cfg) {
  for (auto& spec : cfg.chemicals) {
    if (spec.name == species::ACETATE) {
      spec.initial_conc = 0.0;
      spec.boundary_conc = 0.0;
      return;
    }
  }
}

void fill_acetate_field(Simulation& sim, Real conc) {
  const Int i_acetate = sim.chemical_field().find(species::ACETATE);
  if (i_acetate < 0) return;
  for (Int c = 0; c < sim.chemical_field().ncells(); ++c) {
    sim.chemical_field().conc(i_acetate, c) = conc;
  }
}

void downregulate_btuB(Agent& agent) {
  const Int idx = to_underlying(ReceptorType::BtuB);
  agent.receptor_expr[idx] = 0.1;
  agent.receptor_expr_base[idx] = 0.1;
  agent.genome.receptor_expression[idx] = 0.1;
}

}  // namespace

void test_aerobic_growth_advantage() {
  SimulationConfig cfg = make_combo_config(2001);
  cfg.chem_env.oxygen.enabled = true;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain near_epi;
  near_epi.type = 1;
  near_epi.count = 10;
  near_epi.mu_max = 5e-4;
  cfg.initial_strains.push_back(near_epi);

  Simulation sim;
  sim.init(cfg);

  Int i = 0;
  std::vector<bool> near_group;
  near_group.reserve(sim.agents().size());
  for (Agent& a : sim.agents()) {
    if (i < 5) {
      a.x[2] = 5e-6;
      near_group.push_back(true);
    } else {
      a.x[2] = 45e-6;
      near_group.push_back(false);
    }
    ++i;
  }

  sim.step(cfg.time.bio_dt);

  Real mu_near = 0.0;
  Real mu_far = 0.0;
  Int n_near = 0;
  Int n_far = 0;
  for (Int idx = 0; idx < sim.agents().size(); ++idx) {
    const Agent& a = sim.agents()[idx];
    if (near_group[idx]) {
      mu_near += a.mu_realized;
      ++n_near;
    } else {
      mu_far += a.mu_realized;
      ++n_far;
    }
  }
  assert(n_near > 0 && n_far > 0);
  mu_near /= n_near;
  mu_far /= n_far;

  assert(mu_near > mu_far);

  assert_chemistry_sane(sim);
  std::cout << "  test_aerobic_growth_advantage: PASSED"
            << " (mu_near=" << mu_near << " mu_far=" << mu_far << ")\n";
}

void test_metabolic_trap_with_acetate() {
  SimulationConfig cfg_on = make_combo_config(2002);
  cfg_on.chem_env.oxygen.enabled = true;
  cfg_on.chem_env.acetate.enabled = true;

  Simulation sim_on;
  sim_on.init(cfg_on);

  Agent& target_on = sim_on.agents()[0];
  downregulate_btuB(target_on);
  fill_acetate_field(sim_on, 80.0);

  sim_on.step(cfg_on.time.bio_dt);
  const Real mu_after_on = sim_on.agents()[0].mu_realized;

  SimulationConfig cfg_off = make_combo_config(2002);
  cfg_off.chem_env.oxygen.enabled = true;
  cfg_off.chem_env.acetate.enabled = false;
  zero_acetate_pool(cfg_off);

  Simulation sim_off;
  sim_off.init(cfg_off);

  Agent& target_off = sim_off.agents()[0];
  downregulate_btuB(target_off);

  sim_off.step(cfg_off.time.bio_dt);
  const Real mu_after_off = sim_off.agents()[0].mu_realized;

  assert(mu_after_on < mu_after_off);

  assert_chemistry_sane(sim_on);
  assert_chemistry_sane(sim_off);
  std::cout << "  test_metabolic_trap_with_acetate: PASSED"
            << " (mu_acetate_on=" << mu_after_on
            << " mu_acetate_off=" << mu_after_off << ")\n";
}

void test_full_chemical_environment() {
  SimulationConfig cfg = make_combo_config(2003);
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.acetate.enabled = true;
  cfg.chem_env.mucin.enabled = true;
  cfg.chem_env.siderophore.enabled = true;

  Simulation sim = run_combo(cfg);

  assert_chemistry_sane(sim);
  assert_population_sane(sim);

  const auto& chem = sim.chemical_field();
  for (Int s = 0; s < chem.num_species(); ++s) {
    Real max_conc = 0.0;
    for (Int c = 0; c < chem.ncells(); ++c) {
      max_conc = std::max(max_conc, chem.conc(s, c));
    }
    assert(max_conc < 1e6);
  }

  const Int i_carbon = chem.find(species::CARBON);
  if (i_carbon >= 0) {
    Real sum = 0.0;
    for (Int c = 0; c < chem.ncells(); ++c) {
      sum += chem.conc(i_carbon, c);
    }
    assert(sum > 0.0);
  }

  std::cout << "  test_full_chemical_environment: PASSED\n";
}

void test_adaptive_dt_with_crypts() {
  SimulationConfig cfg = make_combo_config(2004);
  cfg.chem_env.oxygen.enabled = true;
  cfg.adaptive_dt.enabled = true;
  cfg.adaptive_dt.min = 1.0;
  cfg.adaptive_dt.max = 120.0;
  cfg.adaptive_dt.safety = 0.8;
  cfg.adaptive_dt.growth_limit = 0.1;
  cfg.advection.crypts_enabled = true;
  cfg.advection.crypt_depth = 10e-6;
  cfg.advection.crypt_entry_rate = 1e-3;
  cfg.advection.crypt_exit_rate = 1e-4;
  cfg.advection.crypt_carrying_capacity = 10;

  Simulation sim = run_combo(cfg);

  assert_chemistry_sane(sim);
  assert_population_sane(sim);

  Int crypt_count = 0;
  for (const Agent& a : sim.agents()) {
    if (a.state != PhenoState::DEAD && a.flags.in_crypt) {
      ++crypt_count;
    }
  }
  assert(crypt_count > 0);

  const Int steps = sim.step_count();
  const auto min_steps = static_cast<Int>(cfg.time.total_time / cfg.adaptive_dt.max);
  const auto max_steps = static_cast<Int>(cfg.time.total_time / cfg.adaptive_dt.min);
  assert(steps >= min_steps);
  assert(steps <= max_steps);

  std::cout << "  test_adaptive_dt_with_crypts: PASSED"
            << " (crypt_count=" << crypt_count << " steps=" << steps << ")\n";
}

void test_full_biology_with_fmm() {
  SimulationConfig cfg = make_combo_config(2005);
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.acetate.enabled = true;
  cfg.qssa.use_fmm = true;
  cfg.qssa.fmm_theta = 0.5;
  cfg.qssa.fmm_expansion_order = 2;
  // One bio step is enough to exercise the FMM toxin path while keeping
  // enough residents alive for a directional retention assertion.
  cfg.time.total_time = 60.0;

  Simulation sim = run_combo(cfg);

  assert_chemistry_sane(sim);
  assert_population_sane(sim);

  Int type1_count = 0;
  for (const Agent& a : sim.agents()) {
    if (a.state != PhenoState::DEAD && a.identity.type == 1) {
      ++type1_count;
    }
  }
  assert(type1_count > 0);

  std::cout << "  test_full_biology_with_fmm: PASSED"
            << " (type1_survivors=" << type1_count << ")\n";
}

void test_kitchen_sink() {
  SimulationConfig cfg = make_combo_config(2006);

  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.acetate.enabled = true;
  cfg.chem_env.mucin.enabled = true;
  cfg.chem_env.siderophore.enabled = true;

  cfg.adaptive_dt.enabled = true;
  cfg.adaptive_dt.min = 1.0;
  cfg.adaptive_dt.max = 120.0;
  cfg.adaptive_dt.safety = 0.8;

  cfg.advection.crypts_enabled = true;
  cfg.advection.crypt_depth = 10e-6;
  cfg.advection.peristaltic_enabled = true;
  cfg.advection.peristaltic_period = 15.0;
  cfg.advection.peristaltic_amplitude = 0.3;
  cfg.advection.peristaltic_wavelength = 40e-6;

  cfg.qssa.use_fmm = true;
  cfg.qssa.fmm_theta = 0.5;

  cfg.dysbiosis_threshold = 1e12;

  cfg.initial_strains[0].cdi_type = 1;
  cfg.initial_strains[0].cdi_immunity = 1;

  Simulation sim = run_combo(cfg);

  assert_chemistry_sane(sim);

  const auto& chem = sim.chemical_field();
  for (Int s = 0; s < chem.num_species(); ++s) {
    for (Int c = 0; c < chem.ncells(); ++c) {
      const Real val = chem.conc(s, c);
      assert(std::isfinite(val));
      assert(val >= 0.0);
      assert(val < 1e8);
    }
  }

  std::cout << "  test_kitchen_sink: PASSED\n";
}

int main() {
  std::cout << "=== Feature Combination Integration Tests (Spec 8) ===\n";
  test_aerobic_growth_advantage();
  test_metabolic_trap_with_acetate();
  test_full_chemical_environment();
  test_adaptive_dt_with_crypts();
  test_full_biology_with_fmm();
  test_kitchen_sink();
  std::cout << "All feature combination tests passed.\n";
  return 0;
}
