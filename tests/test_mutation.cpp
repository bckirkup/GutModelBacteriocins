/* -----------------------------------------------------------------------
   GutIBM – Unit tests for fix_mutation
   ----------------------------------------------------------------------- */

#include "fix_mutation.h"
#include "plasmid.h"
#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <array>
#include <cmath>
#include <iostream>

using namespace gutibm;

static Simulation make_empty_sim(uint64_t seed = 42) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.seed = seed;

  Simulation sim;
  sim.init(cfg);
  return sim;
}

static Agent make_dividing_agent(Simulation& sim) {
  Vec3 center = {
    0.5 * (sim.domain().lo()[0] + sim.domain().hi()[0]),
    0.5 * (sim.domain().lo()[1] + sim.domain().hi()[1]),
    0.5 * (sim.domain().lo()[2] + sim.domain().hi()[2]),
  };
  Agent a = Agent::create_default(sim.agents().next_tag(), 1, center, 5e-4);
  a.age = 0.0;
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_B());
  return a;
}

void test_bi_duplication_increases_loci() {
  MutationConfig cfg;
  cfg.bi_duplication_rate = 1.0;

  auto sim = make_empty_sim(7007);
  Agent a = make_dividing_agent(sim);
  auto n_before = static_cast<Int>(a.genome.bi_loci.size());
  sim.agents().push_back(std::move(a));

  FixMutation fix(sim, cfg);
  fix.compute(60.0);

  auto n_after = static_cast<Int>(sim.agents()[0].genome.bi_loci.size());
  assert(n_after == n_before + 1);
  assert(sim.agents()[0].genome.mutations > 0);

  std::cout << "  test_bi_duplication_increases_loci: PASSED\n";
}

void test_receptor_downregulation_reduces_expression() {
  MutationConfig cfg;
  cfg.receptor_mutation_rate = 1.0;
  cfg.receptor_reduction = 0.25;

  auto sim = make_empty_sim(8008);
  Agent a = make_dividing_agent(sim);
  std::array<Real, NUM_RECEPTORS> expr_before = a.receptor_expr_base;
  sim.agents().push_back(std::move(a));

  FixMutation fix(sim, cfg);
  fix.compute(60.0);

  Real total_before = 0.0;
  Real total_after = 0.0;
  for (int ri = 0; ri < NUM_RECEPTORS; ++ri) {
    total_before += expr_before[ri];
    total_after += sim.agents()[0].receptor_expr_base[ri];
  }
  assert(total_after < total_before);
  assert(std::abs(total_after - (total_before - cfg.receptor_reduction)) < 1e-12);

  std::cout << "  test_receptor_downregulation_reduces_expression: PASSED\n";
}

void test_partial_resistance_sets_affinity() {
  MutationConfig cfg;
  cfg.partial_resistance_rate = 1.0;

  auto sim = make_empty_sim(9009);
  Agent a = make_dividing_agent(sim);
  sim.agents().push_back(std::move(a));

  FixMutation fix(sim, cfg);
  fix.compute(60.0);

  bool reduced = false;
  for (int ri = 0; ri < NUM_RECEPTORS; ++ri) {
    if (sim.agents()[0].genome.toxin_affinity[ri] < 1.0) {
      reduced = true;
      assert(sim.agents()[0].genome.toxin_affinity[ri] >= 0.01);
      assert(sim.agents()[0].genome.toxin_affinity[ri] <= 0.1);
    }
  }
  assert(reduced);

  std::cout << "  test_partial_resistance_sets_affinity: PASSED\n";
}

void test_super_killer_immunity_escape() {
  MutationConfig cfg;
  cfg.super_killer_rate = 1.0;
  cfg.immunity_escape_prob = 1.0;
  cfg.escape_affinity_lo = 0.05;
  cfg.escape_affinity_hi = 0.2;

  auto sim = make_empty_sim(9010);
  Agent a = make_dividing_agent(sim);
  auto n_before = static_cast<Int>(a.genome.bi_loci.size());
  sim.agents().push_back(std::move(a));

  FixMutation fix(sim, cfg);
  fix.compute(60.0);

  auto n_after = static_cast<Int>(sim.agents()[0].genome.bi_loci.size());
  assert(n_after == n_before + 1);

  const BICluster& novel = sim.agents()[0].genome.bi_loci.back();
  assert(novel.immunity_binding_affinity < 1.0);
  assert(novel.immunity_binding_affinity >= cfg.escape_affinity_lo);
  assert(novel.immunity_binding_affinity <= cfg.escape_affinity_hi);

  std::cout << "  test_super_killer_immunity_escape: PASSED\n";
}

void test_mutation_skips_aged_agents() {
  MutationConfig cfg;
  cfg.bi_duplication_rate = 1.0;

  auto sim = make_empty_sim(9011);
  Agent a = make_dividing_agent(sim);
  a.age = 120.0;
  auto n_before = static_cast<Int>(a.genome.bi_loci.size());
  sim.agents().push_back(std::move(a));

  FixMutation fix(sim, cfg);
  fix.compute(60.0);

  assert(static_cast<Int>(sim.agents()[0].genome.bi_loci.size()) == n_before);

  std::cout << "  test_mutation_skips_aged_agents: PASSED\n";
}

int main() {
  std::cout << "=== Mutation Fix Tests ===\n";
  test_bi_duplication_increases_loci();
  test_receptor_downregulation_reduces_expression();
  test_partial_resistance_sets_affinity();
  test_super_killer_immunity_escape();
  test_mutation_skips_aged_agents();
  std::cout << "All mutation fix tests passed.\n";
  return 0;
}
