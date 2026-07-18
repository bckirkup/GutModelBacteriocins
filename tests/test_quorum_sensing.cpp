/* -----------------------------------------------------------------------
   GutIBM – AI-2 quorum sensing tests (Spec 11)
   Golden anchors + config-sensitivity probes. Kept short for CI.
   ----------------------------------------------------------------------- */

#include "fix_quorum_sensing.h"
#include "fix_motility.h"
#include "simulation.h"
#include "input_parser.h"
#include "species_names.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

namespace {

bool has_ai2_species(const SimulationConfig& cfg) {
  return std::ranges::any_of(cfg.chemicals, [](const ChemicalSpec& s) {
    return s.name == species::AI2;
  });
}

SimulationConfig qs_cfg(uint64_t seed = 42) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = seed;
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 50e-6};
  cfg.domain.grid_dx = 10e-6;
  cfg.domain.periodic = {false, false, false};
  cfg.advection.radial_turnover = 1.0e30;
  cfg.advection.distal_transit_time = 1.0e30;
  cfg.quorum_sensing.enabled = true;
  cfg.quorum_sensing.ai2_basal_rate = 1.0e-20;
  cfg.quorum_sensing.ai2_growth_coupled = 1.0e-16;
  cfg.quorum_sensing.lsr_vmax = 1.0e-18;
  cfg.quorum_sensing.lsr_km = 1.0e-7;
  cfg.quorum_sensing.ai2_decay_rate = 0.0;  // isolate production/import
  cfg.quorum_sensing.ai2_chemotaxis_enabled = false;
  cfg.cell_bio.motility.enabled = false;
  cfg.enabled_fixes = {"quorum_sensing"};
  return cfg;
}

void sync_grid_cell(Simulation& sim, Agent& agent) {
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  sim.domain().pos_to_grid(agent.x, ix, iy, iz);
  agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
}

Agent place_agent(Simulation& sim, Vec3 pos, Real mu_realized) {
  Agent a = Agent::create_default(sim.agents().next_tag(), 1, pos, 5e-4);
  a.mu_realized = mu_realized;
  sync_grid_cell(sim, a);
  return a;
}

Real peak_ai2(const Simulation& sim, Int i_ai2) {
  Real peak = 0.0;
  for (Int c = 0; c < sim.chemical_field().ncells(); ++c) {
    peak = std::max(peak, sim.chemical_field().conc(i_ai2, c));
  }
  return peak;
}

}  // namespace

void test_qs_disabled_has_no_effect() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.quorum_sensing.enabled = false;
  InputParser::finalize_config(cfg);
  assert(!has_ai2_species(cfg));

  // default_config already finalized with QS off → no AI-2 species
  SimulationConfig cfg_default = InputParser::default_config();
  assert(!has_ai2_species(cfg_default));

  std::cout << "  test_qs_disabled_has_no_effect: PASSED\n";
}

void test_ai2_species_registered_when_enabled() {
  auto cfg = qs_cfg();
  InputParser::finalize_config(cfg);
  bool found = false;
  for (const auto& s : cfg.chemicals) {
    if (s.name == species::AI2) {
      found = true;
      assert(s.diffusion_enabled);
      assert(!s.z_gradient_enabled);
      assert(std::abs(s.diff_coeff - cfg.quorum_sensing.ai2_D_free) < 1e-20);
    }
  }
  assert(found);
  std::cout << "  test_ai2_species_registered_when_enabled: PASSED\n";
}

void test_ai2_production_deposits_reactions() {
  auto cfg = qs_cfg();
  Simulation sim;
  sim.init(cfg);
  const Int i_ai2 = sim.chemical_field().find(species::AI2);
  assert(i_ai2 >= 0);

  Agent growing = place_agent(sim, {25e-6, 25e-6, 25e-6}, 5e-4);
  const Int cell = growing.grid_cell;
  sim.agents().push_back(std::move(growing));

  sim.chemical_field().zero_reactions();
  FixQuorumSensing fix(sim, sim.config().quorum_sensing);
  fix.compute(60.0);
  const Real reac_growing = sim.chemical_field().reac(i_ai2, cell);
  assert(reac_growing > 0.0);

  // Sensitivity: starving agent (mu=0) produces only basal → lower rate
  Simulation sim0;
  sim0.init(cfg);
  Agent starved = place_agent(sim0, {25e-6, 25e-6, 25e-6}, 0.0);
  const Int cell0 = starved.grid_cell;
  sim0.agents().push_back(std::move(starved));
  sim0.chemical_field().zero_reactions();
  FixQuorumSensing fix0(sim0, sim0.config().quorum_sensing);
  fix0.compute(60.0);
  const Real reac_starved = sim0.chemical_field().reac(i_ai2, cell0);
  assert(reac_starved > 0.0);
  // basal 1e-20 vs basal + 1e-16*5e-4 (=6e-20) → ~6× ratio
  assert(reac_growing > reac_starved * 5.0);

  // Sensitivity: doubling basal rate increases starved production
  auto cfg_hi = cfg;
  cfg_hi.quorum_sensing.ai2_basal_rate = 2.0e-20;
  Simulation sim_hi;
  sim_hi.init(cfg_hi);
  Agent starved_hi = place_agent(sim_hi, {25e-6, 25e-6, 25e-6}, 0.0);
  const Int cell_hi = starved_hi.grid_cell;
  sim_hi.agents().push_back(std::move(starved_hi));
  sim_hi.chemical_field().zero_reactions();
  FixQuorumSensing fix_hi(sim_hi, sim_hi.config().quorum_sensing);
  fix_hi.compute(60.0);
  assert(sim_hi.chemical_field().reac(i_ai2, cell_hi) > reac_starved * 1.5);

  std::cout << "  test_ai2_production_deposits_reactions: PASSED\n";
}

void test_lsr_import_creates_local_sink() {
  auto cfg = qs_cfg();
  cfg.quorum_sensing.ai2_basal_rate = 0.0;
  cfg.quorum_sensing.ai2_growth_coupled = 0.0;
  cfg.quorum_sensing.lsr_vmax = 1.0e-18;

  Simulation sim;
  sim.init(cfg);
  const Int i_ai2 = sim.chemical_field().find(species::AI2);
  assert(i_ai2 >= 0);

  for (Int c = 0; c < sim.chemical_field().ncells(); ++c) {
    sim.chemical_field().conc(i_ai2, c) = 1.0e-6;
  }

  Agent a = place_agent(sim, {25e-6, 25e-6, 25e-6}, 0.0);
  const Int agent_cell = a.grid_cell;
  const Int empty_cell = (agent_cell + 1) % sim.chemical_field().ncells();
  sim.agents().push_back(std::move(a));

  sim.chemical_field().zero_reactions();
  FixQuorumSensing fix(sim, sim.config().quorum_sensing);
  fix.compute(60.0);

  const Real reac_agent = sim.chemical_field().reac(i_ai2, agent_cell);
  const Real reac_empty = sim.chemical_field().reac(i_ai2, empty_cell);
  assert(reac_agent < 0.0);
  assert(reac_empty == 0.0);
  assert(reac_agent < reac_empty);

  // Sensitivity: higher lsr_vmax → deeper sink
  auto cfg_hi = cfg;
  cfg_hi.quorum_sensing.lsr_vmax = 5.0e-18;
  Simulation sim_hi;
  sim_hi.init(cfg_hi);
  for (Int c = 0; c < sim_hi.chemical_field().ncells(); ++c) {
    sim_hi.chemical_field().conc(i_ai2, c) = 1.0e-6;
  }
  Agent b = place_agent(sim_hi, {25e-6, 25e-6, 25e-6}, 0.0);
  const Int cell_hi = b.grid_cell;
  sim_hi.agents().push_back(std::move(b));
  sim_hi.chemical_field().zero_reactions();
  FixQuorumSensing fix_hi(sim_hi, sim_hi.config().quorum_sensing);
  fix_hi.compute(60.0);
  assert(sim_hi.chemical_field().reac(i_ai2, cell_hi) < reac_agent);

  std::cout << "  test_lsr_import_creates_local_sink: PASSED\n";
}

void test_ai2_diffusion_spreads_from_source() {
  auto cfg = qs_cfg();
  cfg.quorum_sensing.lsr_vmax = 0.0;
  cfg.quorum_sensing.ai2_decay_rate = 0.0;
  cfg.quorum_sensing.ai2_basal_rate = 1.0e-18;
  cfg.quorum_sensing.ai2_growth_coupled = 0.0;
  cfg.enabled_fixes = {"quorum_sensing"};

  Simulation sim;
  sim.init(cfg);
  const Int i_ai2 = sim.chemical_field().find(species::AI2);
  assert(i_ai2 >= 0);

  // Place above the epithelial Dirichlet plane (z=0 is forced to boundary_conc)
  Agent a = place_agent(sim, {25e-6, 25e-6, 25e-6}, 0.0);
  const Int source_cell = a.grid_cell;
  sim.agents().push_back(std::move(a));

  for (int step = 0; step < 5; ++step) {
    sim.step(10.0);
  }

  const Real c_source = sim.chemical_field().conc(i_ai2, source_cell);
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  sim.domain().pos_to_grid({35e-6, 25e-6, 25e-6}, ix, iy, iz);
  const Int neighbor_cell = sim.domain().cell_index(ix, iy, iz);
  const Real c_neighbor = sim.chemical_field().conc(i_ai2, neighbor_cell);
  assert(c_source > 0.0);
  assert(c_neighbor > 0.0);
  assert(c_source > c_neighbor);

  // Sensitivity: decay lowers the peak concentration
  auto cfg_decay = cfg;
  cfg_decay.quorum_sensing.ai2_decay_rate = 1.0;  // strong decay
  Simulation sim_decay;
  sim_decay.init(cfg_decay);
  Agent b = place_agent(sim_decay, {25e-6, 25e-6, 25e-6}, 0.0);
  const Int src_decay = b.grid_cell;
  sim_decay.agents().push_back(std::move(b));
  for (int step = 0; step < 5; ++step) {
    sim_decay.step(10.0);
  }
  const Real c_decay = sim_decay.chemical_field().conc(i_ai2, src_decay);
  assert(c_source > c_decay);

  std::cout << "  test_ai2_diffusion_spreads_from_source: PASSED\n";
}

void test_ai2_chemotaxis_bias() {
  auto cfg = qs_cfg();
  cfg.quorum_sensing.ai2_chemotaxis_enabled = true;
  cfg.quorum_sensing.chi_ai2 = 3.0;
  cfg.cell_bio.motility.enabled = true;
  cfg.cell_bio.motility.aerotaxis_enabled = false;
  cfg.cell_bio.motility.chemotaxis_enabled = false;
  cfg.cell_bio.motility.energy_taxis_enabled = false;
  cfg.cell_bio.motility.stop_probability = 0.0;
  cfg.cell_bio.motility.run_mean_duration = 1.0e9;
  cfg.enabled_fixes = {"motility"};

  Simulation sim;
  sim.init(cfg);
  const Int i_ai2 = sim.chemical_field().find(species::AI2);
  assert(i_ai2 >= 0);

  // Concentrations above chemotaxis_threshold so Weber–Fechner uses ΔC/C
  constexpr Real ai2_prev = 2.0e-6;
  constexpr Real ai2_now = 2.2e-6;  // +10% → modifier = 1 + chi*0.1
  Agent a = place_agent(sim, {25e-6, 25e-6, 25e-6}, 5e-4);
  FixMotility::init_agent_motility(a, sim.config().cell_bio.motility, sim.rng());
  a.motility.swim_direction = {-1.0, 0.0, 0.0};
  a.motility.is_stopped = false;
  a.motility.run_timer = 10.0;
  a.motility.prev_ai2 = ai2_prev;
  sim.chemical_field().conc(i_ai2, a.grid_cell) = ai2_now;
  sim.agents().push_back(std::move(a));

  FixMotility fix(sim, sim.config().cell_bio.motility);
  fix.pre_step(1.0);
  const Real timer_hi = sim.agents()[0].motility.run_timer;
  // 10 * (1 + 3*0.1) - 1s swim ≈ 12.0
  assert(timer_hi > 11.0);

  // Sensitivity: chi_ai2 = 0 → no lengthening (timer ≈ 9 after 1s swim)
  auto cfg0 = cfg;
  cfg0.quorum_sensing.chi_ai2 = 0.0;
  Simulation sim0;
  sim0.init(cfg0);
  Agent b = place_agent(sim0, {25e-6, 25e-6, 25e-6}, 5e-4);
  FixMotility::init_agent_motility(b, sim0.config().cell_bio.motility, sim0.rng());
  b.motility.swim_direction = {-1.0, 0.0, 0.0};
  b.motility.is_stopped = false;
  b.motility.run_timer = 10.0;
  b.motility.prev_ai2 = ai2_prev;
  sim0.chemical_field().conc(i_ai2, b.grid_cell) = ai2_now;
  sim0.agents().push_back(std::move(b));
  FixMotility fix0(sim0, sim0.config().cell_bio.motility);
  fix0.pre_step(1.0);
  const Real timer_zero = sim0.agents()[0].motility.run_timer;
  assert(timer_hi > timer_zero + 2.0);

  std::cout << "  test_ai2_chemotaxis_bias: PASSED\n";
}

void test_clustering_increases_local_ai2() {
  auto make_peak = [](bool clustered) {
    auto cfg = qs_cfg(77);
    cfg.quorum_sensing.lsr_vmax = 0.0;
    cfg.quorum_sensing.ai2_decay_rate = 0.0;
    cfg.quorum_sensing.ai2_basal_rate = 1.0e-19;
    cfg.quorum_sensing.ai2_growth_coupled = 0.0;
    cfg.enabled_fixes = {"quorum_sensing"};
    Simulation sim;
    sim.init(cfg);
    const Int i_ai2 = sim.chemical_field().find(species::AI2);
    assert(i_ai2 >= 0);

    if (clustered) {
      for (int i = 0; i < 5; ++i) {
        sim.agents().push_back(
            place_agent(sim, {25e-6, 25e-6, 25e-6}, 0.0));
      }
    } else {
      const std::vector<Vec3> positions = {
          {5e-6, 5e-6, 5e-6},
          {45e-6, 5e-6, 5e-6},
          {5e-6, 45e-6, 5e-6},
          {45e-6, 45e-6, 5e-6},
          {25e-6, 25e-6, 45e-6},
      };
      for (const Vec3& p : positions) {
        sim.agents().push_back(place_agent(sim, p, 0.0));
      }
    }

    for (int step = 0; step < 4; ++step) {
      sim.step(10.0);
    }
    return peak_ai2(sim, i_ai2);
  };

  const Real peak_cluster = make_peak(true);
  const Real peak_spread = make_peak(false);
  assert(peak_cluster > 0.0);
  assert(peak_cluster > peak_spread);

  std::cout << "  test_clustering_increases_local_ai2: PASSED\n";
}

int main() {
  std::cout << "=== Quorum Sensing / AI-2 Tests (Spec 11) ===\n";
  test_qs_disabled_has_no_effect();
  test_ai2_species_registered_when_enabled();
  test_ai2_production_deposits_reactions();
  test_lsr_import_creates_local_sink();
  test_ai2_diffusion_spreads_from_source();
  test_ai2_chemotaxis_bias();
  test_clustering_increases_local_ai2();
  std::cout << "All quorum sensing tests passed.\n";
  return 0;
}
