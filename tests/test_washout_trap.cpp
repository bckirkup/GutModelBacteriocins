/* -----------------------------------------------------------------------
   GutIBM – Metabolic washout trap long-horizon regression (issue #160)
   VADI combinatorial trap: sustained mu_realized < gamma_flow expels cells.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "sim_fingerprint.h"
#include "advection.h"
#include "domain.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

namespace {

constexpr Real kTol = 1e-12;

SimulationConfig make_washout_horizon_config(unsigned seed, Real total_time) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {80e-6, 80e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.time.total_time = total_time;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = total_time;
  cfg.seed = seed;
  cfg.hdf5.enabled = false;
  cfg.advection.crypts_enabled = false;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.radial_turnover = 5400.0;
  cfg.advection.washout_trap = WashoutTrapMode::IMPOSED;
  cfg.advection.distal_length = 80e-6;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;
  cfg.cell_bio.fur.enabled = false;
  cfg.cell_bio.cdi.enabled = false;
  cfg.cell_bio.motility.enabled = false;
  return cfg;
}

void apply_trap_immigrant_profile(Agent& agent, Real z_lumen) {
  agent.x[2] = z_lumen;
  agent.flags.in_crypt = false;
  agent.mu_max = 1e-6;
  agent.km.km_carbon = 500.0;
  for (Real& expr : agent.receptor_expr) {
    expr = 0.01;
  }
  for (Real& expr : agent.genome.receptor_expression) {
    expr = 0.01;
  }
}

void apply_resident_profile(Agent& agent, Real z_epithelium) {
  agent.x[2] = z_epithelium;
  agent.flags.in_crypt = false;
}

Int count_live_agents(const Simulation& sim) {
  Int alive = 0;
  for (const Agent& agent : sim.agents()) {
    if (agent.state != PhenoState::DEAD) {
      ++alive;
    }
  }
  return alive;
}

Simulation run_trap_cohort(unsigned seed, int count, Real total_time) {
  SimulationConfig cfg = make_washout_horizon_config(seed, total_time);
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = count;
  immigrant.mu_max = 5e-4;
  immigrant.plasmids = {};
  immigrant.conjugative = false;
  cfg.initial_strains.push_back(immigrant);

  Simulation sim;
  sim.init(cfg);
  assert(sim.global_agent_count() == count);

  Real idx = 0.0;
  for (Agent& agent : sim.agents()) {
    agent.x[0] = 20e-6 + std::fmod(idx * 7e-6, 40e-6);
    agent.x[1] = 20e-6 + std::fmod(idx * 5e-6, 40e-6);
    apply_trap_immigrant_profile(agent, 44e-6);
    idx += 1.0;
  }

  sim.run();
  return sim;
}

Simulation run_resident_cohort(unsigned seed, int count, Real total_time) {
  SimulationConfig cfg = make_washout_horizon_config(seed, total_time);
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain resident;
  resident.type = 1;
  resident.count = count;
  resident.mu_max = 5e-4;
  resident.plasmids = {};
  resident.conjugative = false;
  cfg.initial_strains.push_back(resident);

  Simulation sim;
  sim.init(cfg);
  assert(sim.global_agent_count() == count);

  Real idx = 0.0;
  for (Agent& agent : sim.agents()) {
    agent.x[0] = 20e-6 + std::fmod(idx * 7e-6, 40e-6);
    agent.x[1] = 20e-6 + std::fmod(idx * 5e-6, 40e-6);
    apply_resident_profile(agent, 6e-6);
    idx += 1.0;
  }

  sim.run();
  return sim;
}

void test_long_horizon_trap_expels_cohort() {
  constexpr int kInitialCount = 24;
  constexpr Real kHorizon = 3600.0;
  Simulation sim = run_trap_cohort(16001, kInitialCount, kHorizon);

  assert(sim.global_agent_count() <= 1);
  assert(sim.step_count() == 1);
  assert(std::abs(sim.time() - 60.0) < kTol);

  Int cumulative_washout = sim.step_events().outflow_washout;
  assert(cumulative_washout == kInitialCount);

  std::cout << "  test_long_horizon_trap_expels_cohort: PASSED"
            << " (final_global=" << sim.global_agent_count()
            << " steps=" << sim.step_count() << ")\n";
}

void test_long_horizon_resident_persists() {
  constexpr int kInitialCount = 24;
  constexpr Real kHorizon = 3600.0;
  Simulation sim = run_resident_cohort(16011, kInitialCount, kHorizon);

  assert(sim.global_agent_count() == 24);
  assert(sim.step_count() == 60);
  assert(std::abs(sim.time() - kHorizon) < kTol);

  for (const Agent& agent : sim.agents()) {
    if (agent.state == PhenoState::DEAD) continue;
    assert(agent.mu_realized
           > sim.advection().washout_rate(agent.x[2]));
  }

  std::cout << "  test_long_horizon_resident_persists: PASSED"
            << " (final_global=" << sim.global_agent_count()
            << " steps=" << sim.step_count() << ")\n";
}

void test_trap_monotonic_live_decline() {
  constexpr Int kInitialCount = 12;
  SimulationConfig cfg = make_washout_horizon_config(16021, 1800.0);
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain immigrant;
  immigrant.type = 2;
  immigrant.count = kInitialCount;
  immigrant.mu_max = 5e-4;
  cfg.initial_strains.push_back(immigrant);

  Simulation sim;
  sim.init(cfg);
  Real idx = 0.0;
  for (Agent& agent : sim.agents()) {
    agent.x[0] = 20e-6 + std::fmod(idx * 7e-6, 40e-6);
    agent.x[1] = 20e-6 + std::fmod(idx * 5e-6, 40e-6);
    apply_trap_immigrant_profile(agent, 44e-6);
    idx += 1.0;
  }

  Int prev_live = sim.global_agent_count();
  Int cumulative_washout = 0;

  while (sim.time() < cfg.time.total_time && prev_live > 1) {
    sim.step(cfg.time.bio_dt);
    cumulative_washout += sim.step_events().outflow_washout;
    const Int live = count_live_agents(sim);
    assert(live <= prev_live);
    prev_live = live;
  }

  assert(prev_live <= 1);
  assert(cumulative_washout == kInitialCount - prev_live);
  assert(cumulative_washout > 0);

  std::cout << "  test_trap_monotonic_live_decline: PASSED"
            << " (outflow_washout=" << cumulative_washout << ")\n";
}

void test_trap_vs_resident_fingerprints_differ() {
  constexpr Real kHorizon = 1200.0;
  const uint64_t fp_trap =
      test_util::simulation_fingerprint(run_trap_cohort(16031, 16, kHorizon));
  const uint64_t fp_resident =
      test_util::simulation_fingerprint(run_resident_cohort(16031, 16, kHorizon));
  assert(fp_trap != fp_resident);

  std::cout << "  test_trap_vs_resident_fingerprints_differ: PASSED\n";
}

void test_washout_trap_reproducible() {
  constexpr Real kHorizon = 1200.0;
  const uint64_t fp_a =
      test_util::simulation_fingerprint(run_trap_cohort(16041, 20, kHorizon));
  const uint64_t fp_b =
      test_util::simulation_fingerprint(run_trap_cohort(16041, 20, kHorizon));
  assert(fp_a == fp_b);

  std::cout << "  test_washout_trap_reproducible: PASSED\n";
}

void test_radial_turnover_controls_washout_threshold() {
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {80e-6, 80e-6, 50e-6};
  Domain domain;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.mucus_thickness = 50e-6;
  acfg.crypts_enabled = false;

  AdvectionField fast_adv;
  acfg.radial_turnover = 3600.0;
  fast_adv.init(acfg, domain);

  AdvectionField slow_adv;
  acfg.radial_turnover = 120000.0;
  slow_adv.init(acfg, domain);

  const Real z_lumen = 42e-6;
  const Real gamma_fast = fast_adv.washout_rate(z_lumen);
  const Real gamma_slow = slow_adv.washout_rate(z_lumen);
  assert(gamma_fast > gamma_slow);
  assert(gamma_fast > 1.0e-4);
  assert(gamma_slow < 1.0e-4);

  std::cout << "  test_radial_turnover_controls_washout_threshold: PASSED"
            << " (gamma_fast=" << gamma_fast
            << " gamma_slow=" << gamma_slow << ")\n";
}

Int live_count_for_profile(bool extreme_trap, int count, Real dt) {
  SimulationConfig cfg = make_washout_horizon_config(16052, dt);
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 2;
  strain.count = count;
  strain.mu_max = 5e-4;
  cfg.initial_strains.push_back(strain);

  Simulation sim;
  sim.init(cfg);
  Real idx = 0.0;
  for (Agent& agent : sim.agents()) {
    agent.x[0] = 20e-6 + std::fmod(idx * 7e-6, 40e-6);
    agent.x[1] = 20e-6 + std::fmod(idx * 5e-6, 40e-6);
    if (extreme_trap) {
      apply_trap_immigrant_profile(agent, 44e-6);
    } else {
      agent.x[2] = 6e-6;
      agent.flags.in_crypt = false;
      agent.mu_max = 2e-4;
      for (Real& expr : agent.receptor_expr) {
        expr = 0.55;
      }
      for (Real& expr : agent.genome.receptor_expression) {
        expr = 0.55;
      }
    }
    idx += 1.0;
  }

  sim.step(dt);
  return count_live_agents(sim);
}

void test_trap_profile_more_lethal_than_mild_downregulation() {
  constexpr int kCount = 6;
  constexpr Real kDt = 60.0;

  const Int live_extreme = live_count_for_profile(true, kCount, kDt);
  const Int live_mild = live_count_for_profile(false, kCount, kDt);
  assert(live_extreme == 0);
  assert(live_mild == kCount);

  std::cout << "  test_trap_profile_more_lethal_than_mild_downregulation: PASSED\n";
}

void test_population_stocks_are_instantaneous_and_ordered() {
  constexpr Int kCount = 12;
  const std::vector turnovers = {3600.0, 120000.0, 1.0e6};
  std::vector<Int> trapped;

  for (const Real turnover : turnovers) {
    SimulationConfig cfg = make_washout_horizon_config(16053, 60.0);
    cfg.initial_strains.clear();
    cfg.fixes.metabolism.bacteriostasis_threshold = 2.0e-4;
    cfg.advection.radial_turnover = turnover;
    SimulationConfig::InitialStrain strain;
    strain.type = 2;
    strain.count = kCount;
    strain.mu_max = 5e-4;
    cfg.initial_strains.push_back(strain);

    Simulation sim;
    sim.init(cfg);
    for (Agent& agent : sim.agents()) {
      agent.x[2] = 44e-6;
      agent.mu_realized = 1.0e-4;
    }
    sim.prepare_population_stocks_for_summary();
    assert(sim.population_stocks().bacteriostatic_live == kCount);
    trapped.push_back(sim.population_stocks().washout_trapped_live);
  }

  assert(std::is_sorted(trapped.rbegin(), trapped.rend()));
  assert(trapped.front() > trapped.back());

  SimulationConfig cfg = make_washout_horizon_config(16054, 60.0);
  cfg.initial_strains.clear();
  cfg.fixes.metabolism.bacteriostasis_threshold = 2.0e-4;
  SimulationConfig::InitialStrain strain;
  strain.type = 2;
  strain.count = kCount;
  strain.mu_max = 5e-4;
  cfg.initial_strains.push_back(strain);
  Simulation sim;
  sim.init(cfg);
  for (Agent& agent : sim.agents()) {
    agent.x[2] = 6e-6;
    agent.mu_realized = 1.0e-4;
  }
  sim.prepare_population_stocks_for_summary();
  const Int initial_starving = sim.population_stocks().bacteriostatic_live;
  for (Agent& agent : sim.agents()) {
    agent.mu_realized = 5.0e-4;
  }
  sim.prepare_population_stocks_for_summary();
  assert(initial_starving > sim.population_stocks().bacteriostatic_live);
  assert(sim.population_stocks().bacteriostatic_live <= sim.global_agent_count());

  std::cout << "  test_population_stocks_are_instantaneous_and_ordered: PASSED\n";
}

Int run_threshold_removal(WashoutTrapMode mode, Real turnover) {
  SimulationConfig cfg = make_washout_horizon_config(16055, 60.0);
  cfg.advection.washout_trap = mode;
  cfg.advection.radial_turnover = turnover;
  cfg.initial_strains.clear();
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains.push_back({2, 8, 5e-4, {}, false, 0, 0, {}});

  Simulation sim;
  sim.init(cfg);
  for (Agent& agent : sim.agents()) {
    agent.x[2] = 44e-6;
    apply_trap_immigrant_profile(agent, 44e-6);
  }
  sim.step(60.0);
  return sim.step_events().outflow_washout;
}

std::pair<Int, Int> run_transport_loss(WashoutTrapMode mode, Real domain_z,
                                       Real turnover, Real total_time) {
  SimulationConfig cfg = make_washout_horizon_config(16057, total_time);
  cfg.domain.hi[2] = domain_z;
  cfg.advection.washout_trap = mode;
  cfg.advection.radial_turnover = turnover;
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains.clear();
  cfg.initial_strains.push_back({2, 8, 1e-9, {}, false, 0, 0, {}});

  Simulation sim;
  sim.init(cfg);
  for (Agent& agent : sim.agents()) {
    apply_trap_immigrant_profile(agent, 44e-6);
  }
  Int total_outflow = 0;
  const auto steps = static_cast<Int>(std::ceil(total_time / 60.0));
  for (Int step = 0; step < steps && sim.global_agent_count() > 0; ++step) {
    sim.step(60.0);
    const StepEvents& events = sim.step_events();
    total_outflow += events.outflow_boundary + events.outflow_washout;
  }
  return {total_outflow,
          sim.global_agent_count()};
}

void test_washout_modes_have_distinct_low_flow_behavior() {
  const Int imposed = run_threshold_removal(WashoutTrapMode::IMPOSED, 1.0e7);
  const Int emergent = run_threshold_removal(WashoutTrapMode::EMERGENT, 1.0e7);
  assert(imposed > 0);
  assert(emergent == 0);

  std::cout << "  test_washout_modes_have_distinct_low_flow_behavior: PASSED\n";
}

void test_washout_modes_converge_at_high_flow() {
  const auto [imposed_survivors, imposed_washout] =
      run_transport_loss(WashoutTrapMode::IMPOSED, 50e-6, 1.0, 60.0);
  const auto [emergent_survivors, emergent_washout] =
      run_transport_loss(WashoutTrapMode::EMERGENT, 50e-6, 1.0, 60.0);
  constexpr Int initial_agents = 8;
  assert(imposed_survivors == emergent_survivors);
  assert(imposed_survivors == initial_agents);
  assert(emergent_survivors == initial_agents);
  assert(imposed_washout == 0);
  assert(emergent_washout == 0);
  std::cout << "  test_washout_modes_converge_at_high_flow: PASSED\n";
}

void test_emergent_washout_is_geometry_sensitive() {
  const auto [short_box_survivors, short_box_washout] = run_transport_loss(
      WashoutTrapMode::EMERGENT, 50e-6, 5400.0, 1200.0);
  const auto [tall_box_survivors, tall_box_washout] = run_transport_loss(
      WashoutTrapMode::EMERGENT, 100e-6, 5400.0, 1200.0);
  assert(short_box_survivors > 0);
  assert(short_box_survivors > tall_box_survivors);
  assert(short_box_washout < tall_box_washout);
  std::cout << "  test_emergent_washout_is_geometry_sensitive: PASSED\n";
}

void test_bacteriostasis_threshold_default_is_reachable() {
  SimulationConfig cfg = InputParser::default_config();
  const Real threshold = cfg.fixes.metabolism.bacteriostasis_threshold;
  assert(threshold > -cfg.fixes.metabolism.maintenance_rate);
  assert(threshold < cfg.initial_strains.front().mu_max);

  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {80e-6, 80e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.initial_strains.clear();
  cfg.initial_strains.push_back({2, 6, 5e-4, {}, false, 0, 0, {}});
  cfg.hdf5.enabled = false;
  Simulation sim;
  sim.init(cfg);
  for (Agent& agent : sim.agents()) {
    agent.mu_realized = 0.0;
  }
  sim.prepare_population_stocks_for_summary();
  assert(sim.population_stocks().bacteriostatic_live == 6);

  std::cout << "  test_bacteriostasis_threshold_default_is_reachable: PASSED\n";
}

void test_bacteriostasis_threshold_is_sensitive() {
  const std::vector thresholds = {-4.0e-6, 2.0e-6, 1.0e-5};
  const std::vector realized_mu = {
      -5.0e-6, 0.0, 2.0e-6, 8.0e-6, 2.0e-5, 1.0e-4};
  std::vector<Int> stock_counts;

  for (const Real threshold : thresholds) {
    SimulationConfig cfg = make_washout_horizon_config(16058, 60.0);
    cfg.enabled_fixes = {"metabolism"};
    cfg.fixes.metabolism.bacteriostasis_threshold = threshold;
    cfg.initial_strains.clear();
    cfg.initial_strains.push_back({2, 6, 5e-4, {}, false, 0, 0, {}});

    Simulation sim;
    sim.init(cfg);
    Int index = 0;
    for (Agent& agent : sim.agents()) {
      agent.mu_realized = realized_mu[static_cast<size_t>(index)];
      ++index;
    }
    sim.prepare_population_stocks_for_summary();
    stock_counts.push_back(sim.population_stocks().bacteriostatic_live);
  }

  assert(stock_counts[0] < stock_counts[1]);
  assert(stock_counts[1] < stock_counts[2]);
  assert(stock_counts[0] == 1);
  assert(stock_counts[1] == 2);
  assert(stock_counts[2] == 4);

  std::cout << "  test_bacteriostasis_threshold_is_sensitive: PASSED\n";
}

void test_starvation_is_viable_bacteriostasis() {
  SimulationConfig cfg = make_washout_horizon_config(16056, 60.0);
  cfg.enabled_fixes = {"metabolism"};
  cfg.fixes.metabolism.bacteriostasis_threshold = 1.0;
  cfg.initial_strains.clear();
  cfg.initial_strains.push_back({2, 8, 5e-4, {}, false, 0, 0, {}});

  Simulation sim;
  sim.init(cfg);
  for (Agent& agent : sim.agents()) {
    agent.timers.age = 7200.0;
    agent.mu_realized = 0.0;
  }
  sim.step(60.0);
  sim.prepare_population_stocks_for_summary();
  assert(count_live_agents(sim) == 8);
  assert(sim.population_stocks().bacteriostatic_live == 8);
  assert(sim.step_events().outflow_boundary == 0);
  assert(sim.step_events().outflow_washout == 0);
  assert(sim.step_events().mortality_lysis == 0);
  assert(sim.step_events().mortality_colicin == 0);
  assert(sim.step_events().mortality_cdi == 0);

  std::cout << "  test_starvation_is_viable_bacteriostasis: PASSED\n";
}

}  // namespace

int main() {
  std::cout << "=== Washout Trap Long-Horizon Tests (issue #160) ===\n";
  test_long_horizon_trap_expels_cohort();
  test_long_horizon_resident_persists();
  test_trap_monotonic_live_decline();
  test_trap_vs_resident_fingerprints_differ();
  test_washout_trap_reproducible();
  test_radial_turnover_controls_washout_threshold();
  test_trap_profile_more_lethal_than_mild_downregulation();
  test_population_stocks_are_instantaneous_and_ordered();
  test_washout_modes_have_distinct_low_flow_behavior();
  test_washout_modes_converge_at_high_flow();
  test_emergent_washout_is_geometry_sensitive();
  test_bacteriostasis_threshold_default_is_reachable();
  test_bacteriostasis_threshold_is_sensitive();
  test_starvation_is_viable_bacteriostasis();
  std::cout << "All washout trap regression tests passed.\n";
  return 0;
}
