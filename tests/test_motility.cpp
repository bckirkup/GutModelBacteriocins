/* -----------------------------------------------------------------------
   GutIBM – Active motility tests (Spec 3 / Spec 10v2)
   Golden anchors + config-sensitivity probes. Kept short for CI.
   ----------------------------------------------------------------------- */

#include "fix_motility.h"
#include "simulation.h"
#include "input_parser.h"
#include "species_names.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

namespace {

SimulationConfig base_motility_cfg(uint64_t seed = 42) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = seed;
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.enabled_fixes = {"motility"};
  cfg.domain.hi = {50e-6, 50e-6, 100e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.periodic = {false, false, false};
  cfg.advection.radial_turnover = 1.0e30;
  cfg.advection.distal_transit_time = 1.0e30;
  cfg.cell_bio.motility.enabled = true;
  cfg.cell_bio.motility.swim_speed = 7.76e-6;
  cfg.cell_bio.motility.run_mean_duration = 5.0;
  cfg.cell_bio.motility.stop_probability = 0.0;
  cfg.cell_bio.motility.aerotaxis_enabled = false;
  cfg.cell_bio.motility.chemotaxis_enabled = false;
  cfg.cell_bio.motility.energy_taxis_enabled = false;
  cfg.cell_bio.motility.surface_sensing_enabled = false;
  cfg.cell_bio.motility.mucin_drag_enabled = false;
  return cfg;
}

void sync_grid_cell(Simulation& sim, Agent& agent) {
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  sim.domain().pos_to_grid(agent.x, ix, iy, iz);
  agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
}

Agent make_motile_agent(Simulation& sim, Vec3 pos, Vec3 direction) {
  Agent a = Agent::create_default(sim.agents().next_tag(), 1, pos, 5e-4);
  FixMotility::init_agent_motility(a, sim.config().cell_bio.motility, sim.rng());
  a.motility.swim_direction = direction;
  a.motility.is_stopped = false;
  a.motility.run_timer = 100.0;
  a.motility.step_displacement = {0.0, 0.0, 0.0};
  sync_grid_cell(sim, a);
  return a;
}

void apply_displacement(Simulation& sim, Agent& agent) {
  for (int d = 0; d < 3; ++d) {
    agent.x[d] += agent.motility.step_displacement[d];
  }
  // Keep agent inside domain for grid lookups
  const auto& lo = sim.domain().lo();
  const auto& hi = sim.domain().hi();
  for (int d = 0; d < 3; ++d) {
    agent.x[d] = std::clamp(agent.x[d], lo[d] + 1e-9, hi[d] - 1e-9);
  }
  sync_grid_cell(sim, agent);
}

Real displacement_norm(const Vec3& d) {
  return std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
}

Real partitioned_displacement(Real swim_speed, int steps) {
  SimulationConfig cfg = base_motility_cfg();
  cfg.domain.hi = {2.0e-3, 100e-6, 100e-6};
  cfg.domain.grid_dx = 10e-6;
  cfg.cell_bio.motility.swim_speed = swim_speed;
  cfg.cell_bio.motility.run_mean_duration = 1.0e9;
  cfg.cell_bio.motility.stop_probability = 0.0;

  Simulation sim;
  sim.init(cfg);
  Agent agent = make_motile_agent(sim, {1.0e-3, 50e-6, 50e-6}, {1.0, 0.0, 0.0});
  agent.motility.swim_speed = swim_speed;
  agent.motility.is_stopped = true;
  agent.motility.stop_timer = 59.0;
  agent.motility.run_timer = 0.0;
  agent.mu_realized = agent.mu_max;
  const Real start_x = agent.x[0];
  sim.agents().push_back(std::move(agent));

  const Real dt = 60.0 / static_cast<Real>(steps);
  for (int step = 0; step < steps; ++step) {
    sim.step(dt);
  }
  return sim.agents()[0].x[0] - start_x;
}

}  // namespace

void test_biological_timestep_subcycling() {
  constexpr Real swim_speed = 7.76e-6;
  const Real displacement_coarse = partitioned_displacement(swim_speed, 1);
  const Real displacement_fine = partitioned_displacement(swim_speed, 60);
  const Real displacement_fast = partitioned_displacement(2.0 * swim_speed, 1);

  assert(std::abs(displacement_coarse - swim_speed) < 1.0e-12);
  assert(std::abs(displacement_fine - displacement_coarse) < 1.0e-12);
  assert(std::abs(displacement_fast - 2.0 * displacement_coarse) < 1.0e-12);

  std::cout << "  test_biological_timestep_subcycling: PASSED\n";
}

void test_motility_displacement() {
  auto cfg = base_motility_cfg();
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  Simulation sim_motile;
  sim_motile.init(cfg);
  Agent motile = make_motile_agent(
      sim_motile, {25e-6, 25e-6, 12e-6}, {1.0, 0.0, 0.0});
  const Vec3 start = motile.x;
  sim_motile.agents().push_back(std::move(motile));

  SimulationConfig static_cfg = cfg;
  static_cfg.cell_bio.motility.enabled = false;
  Simulation sim_static;
  sim_static.init(static_cfg);
  Agent still = make_motile_agent(
      sim_static, {25e-6, 25e-6, 12e-6}, {1.0, 0.0, 0.0});
  const Vec3 still_start = still.x;
  sim_static.agents().push_back(std::move(still));

  FixMotility mot_fix(sim_motile, sim_motile.config().cell_bio.motility);
  constexpr Real dt = 1.0;
  constexpr int steps = 20;
  for (int i = 0; i < steps; ++i) {
    mot_fix.pre_step(dt);
    apply_displacement(sim_motile, sim_motile.agents()[0]);
  }

  const Real motile_disp = displacement_norm({
      sim_motile.agents()[0].x[0] - start[0],
      sim_motile.agents()[0].x[1] - start[1],
      sim_motile.agents()[0].x[2] - start[2]});
  const Real static_disp = displacement_norm({
      sim_static.agents()[0].x[0] - still_start[0],
      sim_static.agents()[0].x[1] - still_start[1],
      sim_static.agents()[0].x[2] - still_start[2]});
  assert(motile_disp > static_disp);
  assert(motile_disp > 1.0e-6);

  std::cout << "  test_motility_displacement: PASSED\n";
}

// Golden: Weber–Fechner carbon chemotaxis lengthens runs up-gradient.
// Sensitivity: chi_carbon=0 leaves timer unchanged for the same ΔC/C.
void test_weber_fechner_carbon_chemotaxis() {
  auto cfg = base_motility_cfg();
  cfg.cell_bio.motility.chemotaxis_enabled = true;
  cfg.cell_bio.motility.chi_carbon = 2.0;
  cfg.cell_bio.motility.chemotaxis_threshold = 1.0e-6;

  Simulation sim;
  sim.init(cfg);
  Agent a = make_motile_agent(sim, {25e-6, 25e-6, 50e-6}, {1.0, 0.0, 0.0});
  const Int i_carbon = sim.chemical_field().find(species::CARBON);
  assert(i_carbon >= 0);
  sim.chemical_field().conc(i_carbon, a.grid_cell) = 2.0e-3;
  a.motility.prev_carbon = 2.0e-3;
  a.motility.run_timer = 10.0;
  sim.agents().push_back(std::move(a));

  // Up-gradient: +8% fractional rate → modifier = 1 + 2*0.08 = 1.16
  sim.chemical_field().conc(i_carbon, sim.agents()[0].grid_cell) = 2.16e-3;
  FixMotility fix(sim, sim.config().cell_bio.motility);
  fix.pre_step(1.0);
  const Real timer_up = sim.agents()[0].motility.run_timer;
  // After pre_step, ~1s of run consumed from the modified timer
  // Initial 10 * 1.16 = 11.6, then minus ~1s of swimming → ~10.6
  assert(timer_up > 10.0);

  // Sensitivity: chi_carbon = 0 → no lengthening for same step
  auto cfg0 = cfg;
  cfg0.cell_bio.motility.chi_carbon = 0.0;
  Simulation sim0;
  sim0.init(cfg0);
  const Int i_carbon0 = sim0.chemical_field().find(species::CARBON);
  assert(i_carbon0 >= 0);
  Agent b = make_motile_agent(sim0, {25e-6, 25e-6, 50e-6}, {1.0, 0.0, 0.0});
  sim0.chemical_field().conc(i_carbon0, b.grid_cell) = 2.0e-3;
  b.motility.prev_carbon = 2.0e-3;
  b.motility.run_timer = 10.0;
  sim0.agents().push_back(std::move(b));
  sim0.chemical_field().conc(i_carbon0, sim0.agents()[0].grid_cell) = 2.16e-3;
  FixMotility fix0(sim0, sim0.config().cell_bio.motility);
  fix0.pre_step(1.0);
  const Real timer_chi0 = sim0.agents()[0].motility.run_timer;
  assert(timer_up > timer_chi0 + 0.5);

  std::cout << "  test_weber_fechner_carbon_chemotaxis: PASSED\n";
}

// Golden: aerotaxis lengthens runs toward rising O₂.
// Sensitivity: aerotaxis_sensitivity scales the timer response.
void test_aerotaxis_bias_toward_oxygen() {
  auto cfg = base_motility_cfg();
  cfg.chem_env.oxygen.enabled = true;
  cfg.cell_bio.motility.aerotaxis_enabled = true;
  cfg.cell_bio.motility.aerotaxis_sensitivity = 4.0;
  cfg.cell_bio.motility.run_mean_duration = 1.0;
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  const Int i_o2 = sim.chemical_field().find(species::OXYGEN);
  assert(i_o2 >= 0);

  Agent a = make_motile_agent(sim, {25e-6, 25e-6, 20e-6}, {0.0, 0.0, -1.0});
  const Real o2_mid = sim.chemical_field().conc(i_o2, a.grid_cell);
  a.motility.prev_oxygen = o2_mid;
  a.motility.run_timer = 10.0;
  // Synthetic up-gradient step matching Spec 10v2 (~+7.7%/s fractional)
  sim.chemical_field().conc(i_o2, a.grid_cell) = o2_mid * 1.077;
  sim.agents().push_back(std::move(a));

  FixMotility fix(sim, sim.config().cell_bio.motility);
  fix.pre_step(1.0);
  const Real timer_sensitive = sim.agents()[0].motility.run_timer;
  assert(timer_sensitive > 10.0);

  auto cfg0 = cfg;
  cfg0.cell_bio.motility.aerotaxis_sensitivity = 0.0;
  Simulation sim0;
  sim0.init(cfg0);
  Agent b = make_motile_agent(sim0, {25e-6, 25e-6, 20e-6}, {0.0, 0.0, -1.0});
  const Real o2_mid0 = sim0.chemical_field().conc(i_o2, b.grid_cell);
  b.motility.prev_oxygen = o2_mid0;
  b.motility.run_timer = 10.0;
  sim0.chemical_field().conc(i_o2, b.grid_cell) = o2_mid0 * 1.077;
  sim0.agents().push_back(std::move(b));
  FixMotility fix0(sim0, sim0.config().cell_bio.motility);
  fix0.pre_step(1.0);
  const Real timer_zero = sim0.agents()[0].motility.run_timer;
  assert(timer_sensitive > timer_zero + 1.0);

  // Spatial bias: run-and-reverse under O₂ z-gradient accumulates -z drift
  auto cfg_space = cfg;
  cfg_space.cell_bio.motility.aerotaxis_sensitivity = 4.0;
  cfg_space.seed = 99;
  Simulation sim_space;
  sim_space.init(cfg_space);
  Agent c = make_motile_agent(
      sim_space, {25e-6, 25e-6, 40e-6}, {0.0, 0.0, -1.0});
  c.motility.prev_oxygen =
      sim_space.chemical_field().conc(i_o2, c.grid_cell);
  c.motility.run_timer = 1.0;
  const Real z0 = c.x[2];
  sim_space.agents().push_back(std::move(c));
  FixMotility fix_space(sim_space, sim_space.config().cell_bio.motility);
  for (int i = 0; i < 40; ++i) {
    fix_space.pre_step(1.0);
    apply_displacement(sim_space, sim_space.agents()[0]);
  }
  const Real dz_aero = sim_space.agents()[0].x[2] - z0;

  auto cfg_blind = cfg_space;
  cfg_blind.cell_bio.motility.aerotaxis_sensitivity = 0.0;
  Simulation sim_blind;
  sim_blind.init(cfg_blind);
  Agent d = make_motile_agent(
      sim_blind, {25e-6, 25e-6, 40e-6}, {0.0, 0.0, -1.0});
  d.motility.prev_oxygen =
      sim_blind.chemical_field().conc(i_o2, d.grid_cell);
  d.motility.run_timer = 1.0;
  sim_blind.agents().push_back(std::move(d));
  FixMotility fix_blind(sim_blind, sim_blind.config().cell_bio.motility);
  for (int i = 0; i < 40; ++i) {
    fix_blind.pre_step(1.0);
    apply_displacement(sim_blind, sim_blind.agents()[0]);
  }
  const Real dz_blind = sim_blind.agents()[0].x[2] - z0;
  // Aerotaxis should produce more negative (epithelium-ward) net z than blind
  assert(dz_aero < dz_blind);

  std::cout << "  test_aerotaxis_bias_toward_oxygen: PASSED\n";
}

// Combined aerotaxis + carbon modifier > either alone (parameter wiring).
void test_aerotaxis_dominates_carbon() {
  auto make_timer = [](bool aero, bool carbon, Real aero_s, Real chi) {
    auto cfg = base_motility_cfg();
    cfg.chem_env.oxygen.enabled = true;
    cfg.cell_bio.motility.aerotaxis_enabled = aero;
    cfg.cell_bio.motility.aerotaxis_sensitivity = aero_s;
    cfg.cell_bio.motility.chemotaxis_enabled = carbon;
    cfg.cell_bio.motility.chi_carbon = chi;
    InputParser::finalize_config(cfg);
    Simulation sim;
    sim.init(cfg);
    const Int i_o2 = sim.chemical_field().find(species::OXYGEN);
    const Int i_c = sim.chemical_field().find(species::CARBON);
    assert(i_o2 >= 0 && i_c >= 0);
    Agent a = make_motile_agent(sim, {25e-6, 25e-6, 20e-6}, {0.0, 0.0, -1.0});
    const Real o2 = sim.chemical_field().conc(i_o2, a.grid_cell);
    a.motility.prev_oxygen = o2;
    a.motility.prev_carbon = 2.0e-3;
    a.motility.run_timer = 10.0;
    sim.chemical_field().conc(i_o2, a.grid_cell) = o2 * 1.077;
    sim.chemical_field().conc(i_c, a.grid_cell) = 2.16e-3;
    sim.agents().push_back(std::move(a));
    FixMotility fix(sim, sim.config().cell_bio.motility);
    fix.pre_step(1.0);
    return sim.agents()[0].motility.run_timer;
  };

  const Real both = make_timer(true, true, 4.0, 2.0);
  const Real aero_only = make_timer(true, false, 4.0, 2.0);
  const Real carbon_only = make_timer(false, true, 4.0, 2.0);
  assert(both > aero_only);
  assert(aero_only > carbon_only);

  std::cout << "  test_aerotaxis_dominates_carbon: PASSED\n";
}

// Golden + sensitivity: starved agents slow to energy_taxis_floor.
void test_energy_taxis_slows_starved_agents() {
  auto cfg = base_motility_cfg();
  cfg.cell_bio.motility.energy_taxis_enabled = true;
  cfg.cell_bio.motility.energy_taxis_floor = 0.1;
  cfg.cell_bio.motility.run_mean_duration = 1.0e9;

  Simulation sim_full;
  sim_full.init(cfg);
  Agent full = make_motile_agent(
      sim_full, {25e-6, 25e-6, 50e-6}, {1.0, 0.0, 0.0});
  full.mu_realized = full.mu_max;
  full.motility.run_timer = 100.0;
  sim_full.agents().push_back(std::move(full));
  FixMotility fix_full(sim_full, sim_full.config().cell_bio.motility);
  fix_full.pre_step(1.0);
  const Real disp_full =
      displacement_norm(sim_full.agents()[0].motility.step_displacement);

  Simulation sim_starved;
  sim_starved.init(cfg);
  Agent starved = make_motile_agent(
      sim_starved, {25e-6, 25e-6, 50e-6}, {1.0, 0.0, 0.0});
  starved.mu_realized = 0.0;
  starved.motility.run_timer = 100.0;
  sim_starved.agents().push_back(std::move(starved));
  FixMotility fix_starved(sim_starved, sim_starved.config().cell_bio.motility);
  fix_starved.pre_step(1.0);
  const Real disp_starved =
      displacement_norm(sim_starved.agents()[0].motility.step_displacement);

  assert(std::abs(disp_full - cfg.cell_bio.motility.swim_speed) < 1e-12);
  assert(std::abs(disp_starved - 0.1 * disp_full) < 1e-12);

  // Sensitivity: raising floor increases starved speed
  auto cfg_hi = cfg;
  cfg_hi.cell_bio.motility.energy_taxis_floor = 0.5;
  Simulation sim_hi;
  sim_hi.init(cfg_hi);
  Agent starved_hi = make_motile_agent(
      sim_hi, {25e-6, 25e-6, 50e-6}, {1.0, 0.0, 0.0});
  starved_hi.mu_realized = 0.0;
  starved_hi.motility.run_timer = 100.0;
  sim_hi.agents().push_back(std::move(starved_hi));
  FixMotility fix_hi(sim_hi, sim_hi.config().cell_bio.motility);
  fix_hi.pre_step(1.0);
  const Real disp_hi =
      displacement_norm(sim_hi.agents()[0].motility.step_displacement);
  assert(disp_hi > disp_starved * 2.0);

  std::cout << "  test_energy_taxis_slows_starved_agents: PASSED\n";
}

// Near-epithelium agents swim slower; toggling the flag changes outcome.
void test_surface_sensing_near_epithelium() {
  auto cfg = base_motility_cfg();
  cfg.cell_bio.motility.surface_sensing_enabled = true;
  cfg.cell_bio.motility.surface_sensing_depth = 10e-6;
  cfg.cell_bio.motility.surface_sensing_floor = 0.3;
  cfg.cell_bio.motility.run_mean_duration = 1.0e9;

  auto run_disp = [&](Real z) {
    Simulation sim;
    sim.init(cfg);
    Agent a = make_motile_agent(sim, {25e-6, 25e-6, z}, {1.0, 0.0, 0.0});
    a.motility.run_timer = 100.0;
    sim.agents().push_back(std::move(a));
    FixMotility fix(sim, sim.config().cell_bio.motility);
    fix.pre_step(1.0);
    return displacement_norm(sim.agents()[0].motility.step_displacement);
  };

  const Real near = run_disp(3e-6);   // z_rel = 0.3 → factor = 0.3 + 0.7*0.3 = 0.51
  const Real far = run_disp(50e-6);   // outside depth → full speed
  assert(far > near);
  assert(std::abs(far - cfg.cell_bio.motility.swim_speed) < 1e-12);
  assert(std::abs(near / far - 0.51) < 1e-6);

  // Sensitivity: disabled → near == far
  cfg.cell_bio.motility.surface_sensing_enabled = false;
  Simulation sim_off;
  sim_off.init(cfg);
  Agent a_off = make_motile_agent(
      sim_off, {25e-6, 25e-6, 3e-6}, {1.0, 0.0, 0.0});
  a_off.motility.run_timer = 100.0;
  sim_off.agents().push_back(std::move(a_off));
  FixMotility fix_off(sim_off, sim_off.config().cell_bio.motility);
  fix_off.pre_step(1.0);
  const Real near_off =
      displacement_norm(sim_off.agents()[0].motility.step_displacement);
  assert(std::abs(near_off - far) < 1e-12);

  std::cout << "  test_surface_sensing_near_epithelium: PASSED\n";
}

// Multiplicative composition of energy × surface factors.
void test_speed_modifiers_compose() {
  auto cfg = base_motility_cfg();
  cfg.cell_bio.motility.energy_taxis_enabled = true;
  cfg.cell_bio.motility.energy_taxis_floor = 0.0;  // pure mu_frac scaling
  cfg.cell_bio.motility.surface_sensing_enabled = true;
  cfg.cell_bio.motility.surface_sensing_depth = 10e-6;
  cfg.cell_bio.motility.surface_sensing_floor = 0.0;  // pure z_rel scaling
  cfg.cell_bio.motility.run_mean_duration = 1.0e9;

  Simulation sim;
  sim.init(cfg);
  // mu_frac = 0.5, z_rel = 0.6 → speed factor = 0.5 * 0.6 = 0.3
  Agent a = make_motile_agent(sim, {25e-6, 25e-6, 6e-6}, {1.0, 0.0, 0.0});
  a.mu_realized = 0.5 * a.mu_max;
  a.motility.run_timer = 100.0;
  sim.agents().push_back(std::move(a));
  FixMotility fix(sim, sim.config().cell_bio.motility);
  fix.pre_step(1.0);
  const Real disp =
      displacement_norm(sim.agents()[0].motility.step_displacement);
  const Real expected = 0.3 * cfg.cell_bio.motility.swim_speed;
  assert(std::abs(disp - expected) < 1e-12);

  std::cout << "  test_speed_modifiers_compose: PASSED\n";
}

// Mucin drag: 2× reference → speed = 1/3; disabling restores full speed.
void test_mucin_drag() {
  auto cfg = base_motility_cfg();
  cfg.chem_env.mucin.enabled = true;
  cfg.cell_bio.motility.mucin_drag_enabled = true;
  cfg.cell_bio.motility.mucin_drag_reference = 1.0e-2;
  cfg.cell_bio.motility.run_mean_duration = 1.0e9;
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  const Int i_mucin = sim.chemical_field().find(species::MUCIN);
  assert(i_mucin >= 0);
  Agent a = make_motile_agent(sim, {25e-6, 25e-6, 50e-6}, {1.0, 0.0, 0.0});
  sim.chemical_field().conc(i_mucin, a.grid_cell) = 2.0e-2;  // 2× reference
  a.motility.run_timer = 100.0;
  sim.agents().push_back(std::move(a));
  FixMotility fix(sim, sim.config().cell_bio.motility);
  fix.pre_step(1.0);
  const Real disp_drag =
      displacement_norm(sim.agents()[0].motility.step_displacement);
  const Real expected = cfg.cell_bio.motility.swim_speed / 3.0;
  assert(std::abs(disp_drag - expected) < 1e-12);

  // Sensitivity: mucin_drag_enabled = false → full speed at same mucin
  auto cfg_off = cfg;
  cfg_off.cell_bio.motility.mucin_drag_enabled = false;
  Simulation sim_off;
  sim_off.init(cfg_off);
  Agent b = make_motile_agent(sim_off, {25e-6, 25e-6, 50e-6}, {1.0, 0.0, 0.0});
  sim_off.chemical_field().conc(i_mucin, b.grid_cell) = 2.0e-2;
  b.motility.run_timer = 100.0;
  sim_off.agents().push_back(std::move(b));
  FixMotility fix_off(sim_off, sim_off.config().cell_bio.motility);
  fix_off.pre_step(1.0);
  const Real disp_off =
      displacement_norm(sim_off.agents()[0].motility.step_displacement);
  assert(std::abs(disp_off - cfg.cell_bio.motility.swim_speed) < 1e-12);
  assert(disp_off > disp_drag * 2.0);

  std::cout << "  test_mucin_drag: PASSED\n";
}

int main() {
  std::cout << "=== Motility Tests (Spec 10v2) ===\n";
  test_biological_timestep_subcycling();
  test_motility_displacement();
  test_weber_fechner_carbon_chemotaxis();
  test_aerotaxis_bias_toward_oxygen();
  test_aerotaxis_dominates_carbon();
  test_energy_taxis_slows_starved_agents();
  test_surface_sensing_near_epithelium();
  test_speed_modifiers_compose();
  test_mucin_drag();
  std::cout << "All motility tests passed.\n";
  return 0;
}
