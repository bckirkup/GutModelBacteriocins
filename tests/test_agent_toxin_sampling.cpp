/* -----------------------------------------------------------------------
   GutIBM – Agent-position toxin sampling invariants
   ----------------------------------------------------------------------- */

#include "agent.h"
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

Simulation make_sim(const std::string& mode, bool use_fmm = false,
                    Real domain_size = 120e-6, Real grid_dx = 10e-6,
                    Real toxin_cutoff = 80e-6) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {domain_size, domain_size, domain_size};
  cfg.domain.grid_dx = grid_dx;
  cfg.qssa.toxin_cutoff = toxin_cutoff;
  cfg.qssa.toxin_evaluation = mode;
  cfg.qssa.use_fmm = use_fmm;
  Simulation sim;
  sim.init(cfg);
  return sim;
}

Agent add_target(Simulation& sim, Vec3 position) {
  Agent target = Agent::create_default(sim.agents().next_tag(), 2, position,
                                        5e-4);
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  sim.domain().pos_to_grid(position, ix, iy, iz);
  target.grid_cell = sim.domain().cell_index(ix, iy, iz);
  sim.agents().push_back(target);
  return target;
}

ToxinBurstSource source_at(Vec3 position) {
  ToxinBurstSource source;
  source.pos = position;
  source.creation_time = 0.0;
  source.release_tau = 300.0;
  source.target = ReceptorType::BtuB;
  source.params.diff_coeff = 4e-11;
  source.params.retardation = 1.0;
  source.params.source_rate = 1e-18;
  source.params.decay_rate = 0.0;
  return source;
}

void solve(Simulation& sim, const std::vector<ToxinBurstSource>& sources,
           bool materialize_grid = false) {
  sim.qssa().solve_all_bacteriocin_fields(
      sim.agents(), sources, 0.0, sim.config().chem_env.protease,
      sim.advection(), sim.chemical_field(), nullptr, materialize_grid);
}

void solve(Simulation& sim, const ToxinBurstSource& source,
           bool materialize_grid = false) {
  solve(sim, std::vector<ToxinBurstSource>{source}, materialize_grid);
}

void solve_without_sources(Simulation& sim) {
  sim.qssa().solve_all_bacteriocin_fields(
      sim.agents(), {}, 0.0, sim.config().chem_env.protease,
      sim.advection(), sim.chemical_field(), nullptr, false);
}

}  // namespace

void test_cell_center_matches_grid() {
  auto grid_sim = make_sim("grid");
  auto agent = add_target(grid_sim, {55e-6, 65e-6, 65e-6});
  const ToxinBurstSource source = source_at({65e-6, 65e-6, 65e-6});
  solve(grid_sim, source, true);
  const Int toxin = grid_sim.chemical_field().find(
      species::BACTERIOCIN_BTUB);
  const Real grid_value =
      grid_sim.chemical_field().conc(toxin, agent.grid_cell);

  auto agent_sim = make_sim("agents");
  add_target(agent_sim, agent.x);
  solve(agent_sim, source, true);
  const Real sample = agent_sim.qssa().sampled_toxin_conc(
      0, agent_sim.chemical_field().find(species::BACTERIOCIN_BTUB));
  assert(std::isfinite(sample));
  assert(std::abs(sample - grid_value)
         <= 1e-12 * std::max({1.0, std::abs(sample), std::abs(grid_value)}));
  for (Int cell = 0; cell < grid_sim.chemical_field().ncells(); ++cell) {
    assert(std::abs(agent_sim.chemical_field().conc(toxin, cell)
                    - grid_sim.chemical_field().conc(toxin, cell))
           <= 1e-12 * std::max(
               {1.0, std::abs(agent_sim.chemical_field().conc(toxin, cell)),
                std::abs(grid_sim.chemical_field().conc(toxin, cell))}));
  }
  std::cout << "  test_cell_center_matches_grid: PASSED\n";
}

void test_sampling_exactness_and_bounds() {
  auto sim = make_sim("agents");
  add_target(sim, {55e-6, 65e-6, 65e-6});
  const ToxinBurstSource source = source_at({65e-6, 65e-6, 65e-6});
  solve(sim, source);
  const Int toxin = sim.chemical_field().find(species::BACTERIOCIN_BTUB);
  const Real sample = sim.qssa().sampled_toxin_conc(0, toxin);
  std::vector<Vec3> sources{source.pos};
  std::vector<GreensFunctionParams> params{source.params};
  std::vector<Real> strengths{1.0};
  const Real analytic =
      sim.qssa().point_concentration(sim.agents()[0].x, sources, params,
                                     strengths);
  assert(sample > 0.0);
  assert(std::isfinite(sample));
  assert(std::abs(sample - analytic)
         <= 1e-12 * std::max({1.0, std::abs(sample), std::abs(analytic)}));

  auto empty = make_sim("agents");
  add_target(empty, {55e-6, 65e-6, 65e-6});
  solve_without_sources(empty);
  const Real empty_sample = empty.qssa().sampled_toxin_conc(
      0, empty.chemical_field().find(species::BACTERIOCIN_BTUB));
  assert(empty_sample == 0.0);
  std::cout << "  test_sampling_exactness_and_bounds: PASSED\n";
}

void test_nonmaterialized_grid_is_zeroed() {
  auto sim = make_sim("agents");
  add_target(sim, {55e-6, 65e-6, 65e-6});
  const ToxinBurstSource source = source_at({65e-6, 65e-6, 65e-6});
  solve(sim, source, true);
  solve(sim, source, false);
  for (const ReceptorType target : species::BACTERIOCIN_RECEPTOR_TARGETS) {
    const char* name = species::bacteriocin_species_for(target);
    assert(name != nullptr);
    const Int index = sim.chemical_field().find(name);
    assert(index >= 0);
    for (Int cell = 0; cell < sim.chemical_field().ncells(); ++cell) {
      assert(sim.chemical_field().conc(index, cell) == 0.0);
    }
  }
  std::cout << "  test_nonmaterialized_grid_is_zeroed: PASSED\n";
}

void test_distance_ordering_and_fmm() {
  auto sim = make_sim("agents");
  add_target(sim, {55e-6, 65e-6, 65e-6});
  add_target(sim, {45e-6, 65e-6, 65e-6});
  add_target(sim, {35e-6, 65e-6, 65e-6});
  const ToxinBurstSource source = source_at({65e-6, 65e-6, 65e-6});
  solve(sim, source);
  const Int toxin = sim.chemical_field().find(species::BACTERIOCIN_BTUB);
  const Real near = sim.qssa().sampled_toxin_conc(0, toxin);
  const Real middle = sim.qssa().sampled_toxin_conc(1, toxin);
  const Real far = sim.qssa().sampled_toxin_conc(2, toxin);
  assert(near > middle);
  assert(middle > far);

  auto fmm_sim = make_sim("agents", true);
  add_target(fmm_sim, {55e-6, 65e-6, 65e-6});
  solve(fmm_sim, source);
  const Real fmm_sample = fmm_sim.qssa().sampled_toxin_conc(
      0, fmm_sim.chemical_field().find(species::BACTERIOCIN_BTUB));
  assert(std::isfinite(fmm_sample));
  assert(fmm_sample > 0.0);
  std::cout << "  test_distance_ordering_and_fmm: PASSED\n";
}

void test_sampling_caches_appended_agents() {
  auto sim = make_sim("agents");
  add_target(sim, {55e-6, 65e-6, 65e-6});
  const ToxinBurstSource source = source_at({65e-6, 65e-6, 65e-6});
  solve(sim, source);

  const Agent appended = add_target(sim, {45e-6, 65e-6, 65e-6});
  const Int toxin = sim.chemical_field().find(
      species::BACTERIOCIN_BTUB);
  const Real sampled = sim.qssa().sampled_toxin_conc(1, toxin);
  const Real analytic = sim.qssa().point_concentration(
      appended.x, {source.pos}, {source.params}, {1.0});
  assert(sampled > 0.0);
  assert(std::abs(sampled - analytic)
         <= 1e-12 * std::max({1.0, std::abs(sampled), std::abs(analytic)}));
  std::cout << "  test_sampling_caches_appended_agents: PASSED\n";
}

int main() {
  test_cell_center_matches_grid();
  test_sampling_exactness_and_bounds();
  test_nonmaterialized_grid_is_zeroed();
  test_distance_ordering_and_fmm();
  test_sampling_caches_appended_agents();
  std::cout << "All agent toxin sampling tests passed.\n";
  return 0;
}
