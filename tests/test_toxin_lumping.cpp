/* -----------------------------------------------------------------------
   GutIBM – Lumped bacteriocin field modelling variant tests
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"
#include "fix_receptor.h"
#include "error.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

namespace {

Simulation make_sim(const std::string& lumping,
                    const std::string& evaluation = "grid",
                    bool use_fmm = false) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {120e-6, 120e-6, 120e-6};
  cfg.domain.grid_dx = 10e-6;
  cfg.qssa.toxin_lumping = lumping;
  cfg.qssa.toxin_evaluation = evaluation;
  cfg.qssa.use_fmm = use_fmm;
  Simulation sim;
  sim.init(cfg);
  return sim;
}

ToxinBurstSource source_at(Vec3 position, ReceptorType target) {
  ToxinBurstSource source;
  source.pos = position;
  source.creation_time = 0.0;
  source.release_tau = 300.0;
  source.target = target;
  source.params.diff_coeff = 4e-11;
  source.params.retardation = 1.0;
  source.params.source_rate = 1e-18;
  source.params.decay_rate = 0.0;
  return source;
}

Agent add_agent(Simulation& sim, Vec3 position) {
  Agent agent = Agent::create_default(sim.agents().next_tag(), 2, position,
                                      5e-4);
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  sim.domain().pos_to_grid(position, ix, iy, iz);
  agent.grid_cell = sim.domain().cell_index(ix, iy, iz);
  sim.agents().push_back(agent);
  return agent;
}

void solve(Simulation& sim, const std::vector<ToxinBurstSource>& sources,
           bool materialize) {
  sim.qssa().solve_all_bacteriocin_fields(
      sim.agents(), sources, 0.0, sim.config().chem_env.protease,
      sim.advection(), sim.chemical_field(), nullptr, materialize);
}

}  // namespace

void test_single_target_equivalence() {
  const auto source = source_at({60e-6, 60e-6, 60e-6}, ReceptorType::BtuB);
  auto per_receptor = make_sim("per_receptor");
  auto lumped = make_sim("lumped");
  solve(per_receptor, {source}, true);
  solve(lumped, {source}, true);
  const Int per_idx = per_receptor.chemical_field().find(
      species::BACTERIOCIN_BTUB);
  const Int lumped_idx = lumped.chemical_field().find(
      species::BACTERIOCIN_LUMPED);
  assert(per_idx >= 0 && lumped_idx >= 0);
  for (Int cell = 0; cell < per_receptor.chemical_field().ncells(); ++cell) {
    assert(std::abs(per_receptor.chemical_field().conc(per_idx, cell)
                    - lumped.chemical_field().conc(lumped_idx, cell))
           <= 1e-12);
  }
  std::cout << "  test_single_target_equivalence: PASSED\n";
}

void test_single_target_kill_equivalence() {
  auto per_receptor = make_sim("per_receptor");
  auto lumped = make_sim("lumped");
  add_agent(per_receptor, {60e-6, 60e-6, 60e-6});
  add_agent(lumped, {60e-6, 60e-6, 60e-6});
  for (Simulation* sim : {&per_receptor, &lumped}) {
    Agent& agent = (*sim).agents()[0];
    agent.receptor_expr[to_underlying(ReceptorType::BtuB)] = 1.0;
    agent.receptor_expr[to_underlying(ReceptorType::FepA)] = 0.0;
    agent.receptor_expr[to_underlying(ReceptorType::CirA)] = 0.0;
    agent.receptor_expr[to_underlying(ReceptorType::FhuA)] = 0.0;
  }
  const Int per_idx = per_receptor.chemical_field().find(
      species::BACTERIOCIN_BTUB);
  const Int lumped_idx = lumped.chemical_field().find(
      species::BACTERIOCIN_LUMPED);
  per_receptor.chemical_field().conc(
      per_idx, per_receptor.agents()[0].grid_cell) = 1.0e-4;
  lumped.chemical_field().conc(
      lumped_idx, lumped.agents()[0].grid_cell) = 1.0e-4;
  ReceptorConfig cfg;
  cfg.kill_rate_colicin = 1.0;
  FixReceptor per_fix(per_receptor, cfg);
  FixReceptor lumped_fix(lumped, cfg);
  per_fix.compute(60.0);
  lumped_fix.compute(60.0);
  assert(per_receptor.agents()[0].state == lumped.agents()[0].state);
  assert(per_receptor.agents()[0].state == PhenoState::DEAD);
  std::cout << "  test_single_target_kill_equivalence: PASSED\n";
}

void test_lumped_field_is_per_target_sum() {
  const std::vector<ToxinBurstSource> sources = {
      source_at({30e-6, 60e-6, 60e-6}, ReceptorType::BtuB),
      source_at({60e-6, 30e-6, 60e-6}, ReceptorType::FepA),
      source_at({90e-6, 60e-6, 60e-6}, ReceptorType::CirA),
      source_at({60e-6, 90e-6, 60e-6}, ReceptorType::FhuA)};
  auto per_receptor = make_sim("per_receptor");
  auto lumped = make_sim("lumped");
  solve(per_receptor, sources, true);
  solve(lumped, sources, true);
  const Int lumped_idx = lumped.chemical_field().find(
      species::BACTERIOCIN_LUMPED);
  assert(lumped_idx >= 0);
  for (Int cell = 0; cell < lumped.chemical_field().ncells(); ++cell) {
    Real sum = 0.0;
    for (const ReceptorType target : species::BACTERIOCIN_RECEPTOR_TARGETS) {
      const Int idx = per_receptor.chemical_field().find(
          species::bacteriocin_species_for(target));
      sum += per_receptor.chemical_field().conc(idx, cell);
    }
    assert(std::abs(sum - lumped.chemical_field().conc(lumped_idx, cell))
           <= 1e-12 * std::max(1.0, std::abs(sum)));
  }
  std::cout << "  test_lumped_field_is_per_target_sum: PASSED\n";
}

void test_lumped_agent_sampling_matches_analytic_sum() {
  auto sim = make_sim("lumped", "agents");
  const Agent target = add_agent(sim, {55e-6, 65e-6, 65e-6});
  const std::vector<ToxinBurstSource> sources = {
      source_at({65e-6, 65e-6, 65e-6}, ReceptorType::BtuB),
      source_at({45e-6, 65e-6, 65e-6}, ReceptorType::FepA)};
  solve(sim, sources, false);
  const Int toxin = sim.chemical_field().find(species::BACTERIOCIN_LUMPED);
  const Real sampled = sim.qssa().sampled_toxin_conc(0, toxin);
  const Real analytic = sim.qssa().point_concentration(
      target.x, {sources[0].pos, sources[1].pos},
      {sources[0].params, sources[1].params}, {1.0, 1.0});
  assert(sampled > 0.0);
  assert(std::isfinite(sampled));
  assert(std::abs(sampled - analytic)
         <= 1e-12 * std::max({1.0, std::abs(sampled), std::abs(analytic)}));
  std::cout << "  test_lumped_agent_sampling_matches_analytic_sum: PASSED\n";
}

void test_lumped_agent_sampling_fmm_is_finite() {
  auto sim = make_sim("lumped", "agents", true);
  add_agent(sim, {55e-6, 65e-6, 65e-6});
  const std::vector<ToxinBurstSource> sources = {
      source_at({65e-6, 65e-6, 65e-6}, ReceptorType::BtuB),
      source_at({45e-6, 65e-6, 65e-6}, ReceptorType::FepA)};
  solve(sim, sources, false);
  const Int toxin = sim.chemical_field().find(species::BACTERIOCIN_LUMPED);
  const Real sampled = sim.qssa().sampled_toxin_conc(0, toxin);
  assert(std::isfinite(sampled));
  assert(sampled > 0.0);
  std::cout << "  test_lumped_agent_sampling_fmm_is_finite: PASSED\n";
}

void test_lumped_rejects_mismatched_specs() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.qssa.toxin_lumping = "lumped";
  Int idx = -1;
  for (Int i = 0; i < static_cast<Int>(cfg.chemicals.size()); ++i) {
    if (cfg.chemicals[static_cast<size_t>(i)].name
        == species::BACTERIOCIN_FEPA) {
      idx = i;
      break;
    }
  }
  assert(idx >= 0);
  cfg.chemicals[static_cast<size_t>(idx)].retardation += 1.0;
  bool threw = false;
  try {
    Simulation sim;
    sim.init(cfg);
  } catch (const ConfigError& error) {
    threw = true;
    assert(std::string(error.what()).find(
               "matching diff_coeff, retardation, and decay_rate")
           != std::string::npos);
  }
  assert(threw);
  std::cout << "  test_lumped_rejects_mismatched_specs: PASSED\n";
}

int main() {
  test_single_target_equivalence();
  test_single_target_kill_equivalence();
  test_lumped_field_is_per_target_sum();
  test_lumped_agent_sampling_matches_analytic_sum();
  test_lumped_agent_sampling_fmm_is_finite();
  test_lumped_rejects_mismatched_specs();
  std::cout << "All toxin lumping tests passed.\n";
  return 0;
}
