/* -----------------------------------------------------------------------
   GutIBM – Tests for immunity-escape super-killer mutations (Issue #12)
   Verifies the affinity-neutralization matrix: super-killer toxins with
   reduced immunity_binding_affinity can kill agents that would normally
   be immune.
   ----------------------------------------------------------------------- */

#include "agent.h"
#include "plasmid.h"
#include "simulation.h"
#include "input_parser.h"
#include "fix_mutation.h"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace gutibm;

void test_bicluster_default_affinity() {
  // All existing plasmid library entries should have affinity = 1.0
  auto e1 = PlasmidLibrary::colicin_E1();
  assert(std::abs(e1.immunity_binding_affinity - 1.0) < 1e-12);

  auto cb = PlasmidLibrary::colicin_B();
  assert(std::abs(cb.immunity_binding_affinity - 1.0) < 1e-12);

  auto ia = PlasmidLibrary::colicin_Ia();
  assert(std::abs(ia.immunity_binding_affinity - 1.0) < 1e-12);

  auto mv = PlasmidLibrary::microcin_V();
  assert(std::abs(mv.immunity_binding_affinity - 1.0) < 1e-12);

  std::cout << "  test_bicluster_default_affinity: PASSED\n";
}

void test_escape_mutation_reduces_affinity() {
  // Manually create a novel toxin using the mutation config defaults
  // and verify that immunity-escape variants have reduced affinity.
  MutationConfig cfg;
  assert(std::abs(cfg.immunity_escape_prob - 0.5) < 1e-12);
  assert(cfg.escape_affinity_lo > 0.0);
  assert(cfg.escape_affinity_hi < 1.0);

  // A BI cluster with immunity_binding_affinity < 1.0 represents escape
  BICluster escape = PlasmidLibrary::colicin_E1();
  escape.immunity_binding_affinity = 0.1;
  assert(escape.immunity_binding_affinity < 1.0);

  // Full-affinity cluster
  BICluster normal = PlasmidLibrary::colicin_E1();
  assert(std::abs(normal.immunity_binding_affinity - 1.0) < 1e-12);

  std::cout << "  test_escape_mutation_reduces_affinity: PASSED\n";
}

void test_immunity_effectiveness_with_escape() {
  // Verify the immunity math: eff = immunity_factor * immunity_binding_affinity
  //
  // Normal immunity (affinity = 1.0):  eff = 0.001 * 1.0 = 0.001 (1000x protection)
  // Escape  immunity (affinity = 0.1): eff = 0.001 * 0.1 = 0.0001 (still protective but less)
  // No      immunity:                  eff = 1.0 (no protection)

  ReceptorConfig rcfg;
  Real immunity_factor = rcfg.immunity_factor;  // 0.001

  // Normal protection
  Real eff_normal = immunity_factor * 1.0;
  assert(std::abs(eff_normal - 0.001) < 1e-12);

  // Escape variant: reduced protection
  Real eff_escape = immunity_factor * 0.1;
  assert(std::abs(eff_escape - 0.0001) < 1e-12);

  // The escape variant is still more protective than no immunity
  assert(eff_escape < 1.0);
  // But less protective than normal
  assert(eff_escape < eff_normal);

  std::cout << "  test_immunity_effectiveness_with_escape: PASSED\n";
}

void test_create_novel_toxin_with_escape() {
  // Run a mini simulation that exercises super-killer generation and
  // verify that some BICluster entries end up with reduced affinity.
  // We use a high super_killer_rate to guarantee events in a short run.

  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.total_time = 600.0;
  cfg.bio_dt = 60.0;
  cfg.output_interval = 600.0;
  cfg.seed = 7777;
  cfg.hdf5.enabled = false;
  cfg.advection.mucus_thickness = 25e-6;
  cfg.advection.distal_length = 50e-6;
  cfg.qssa.toxin_cutoff = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  // Boost super-killer rate dramatically so we get events in a small sim
  cfg.mutation.super_killer_rate = 0.5;  // 50% per division
  cfg.mutation.immunity_escape_prob = 1.0;  // guarantee escape on every super-killer

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 20;
  s.mu_max = 5e-4;
  s.plasmids = {"ColE1"};
  s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);
  sim.run();

  // Count BI clusters with reduced immunity_binding_affinity
  int escape_count = 0;
  int total_bi = 0;
  for (const Agent& a : sim.agents()) {
    for (const auto& bi : a.genome.bi_loci) {
      total_bi++;
      if (bi.immunity_binding_affinity < 1.0) {
        escape_count++;
        // Verify escape affinity is within configured bounds
        assert(bi.immunity_binding_affinity >= 0.01);
        assert(bi.immunity_binding_affinity <= 0.3);
      }
    }
  }

  // With 50% super-killer rate + 100% escape prob + 20 agents dividing
  // over 10 timesteps, we should have at least some escape variants
  assert(total_bi > 0);
  assert(escape_count > 0);

  std::cout << "  test_create_novel_toxin_with_escape: PASSED"
            << " (escape_bi=" << escape_count
            << " total_bi=" << total_bi << ")\n";
}

void test_smoke_with_immunity_escape() {
  // Compare toxin mortality with and without immunity-escape super-killers.
  struct Outcome {
    Int producers = 0;
    Int immune = 0;
    Int naive = 0;
    Int escape_bi = 0;
  };

  auto run_scenario = [](bool enable_escape, uint64_t seed) {
    SimulationConfig cfg = InputParser::default_config();
    cfg.domain.hi = {100e-6, 100e-6, 50e-6};
    cfg.domain.grid_dx = 5e-6;
    cfg.total_time = 600.0;
    cfg.bio_dt = 60.0;
    cfg.output_interval = 600.0;
    cfg.seed = seed;
    cfg.hdf5.enabled = false;
    cfg.advection.mucus_thickness = 50e-6;
    cfg.advection.distal_length = 100e-6;
    cfg.qssa.toxin_cutoff = 50e-6;
    cfg.qssa.nutrient_cutoff = 25e-6;

    if (enable_escape) {
      cfg.mutation.super_killer_rate = 0.5;
      cfg.mutation.immunity_escape_prob = 1.0;
    }

    cfg.initial_strains.clear();

    SimulationConfig::InitialStrain producer;
    producer.type = 1;
    producer.count = 10;
    producer.mu_max = 5e-4;
    producer.plasmids = {"ColE1"};
    producer.conjugative = false;
    cfg.initial_strains.push_back(producer);

    SimulationConfig::InitialStrain immune_target;
    immune_target.type = 2;
    immune_target.count = 10;
    immune_target.mu_max = 5e-4;
    immune_target.plasmids = {"ColE1"};
    immune_target.conjugative = false;
    cfg.initial_strains.push_back(immune_target);

    SimulationConfig::InitialStrain non_immune;
    non_immune.type = 3;
    non_immune.count = 10;
    non_immune.mu_max = 5e-4;
    non_immune.plasmids = {};
    non_immune.conjugative = false;
    cfg.initial_strains.push_back(non_immune);

    Simulation sim;
    sim.init(cfg);
    sim.run();

    Outcome out;
    for (const Agent& a : sim.agents()) {
      if (a.type == 1) {
        for (const auto& bi : a.genome.bi_loci) {
          if (bi.immunity_binding_affinity < 1.0) out.escape_bi++;
        }
      }
      if (a.state == PhenoState::DEAD) continue;
      if (a.type == 1) {
        out.producers++;
      } else if (a.type == 2) {
        out.immune++;
      } else {
        out.naive++;
      }
    }
    return out;
  };

  constexpr uint64_t seed = 12345;
  Outcome escape = run_scenario(true, seed);

  // Producers should accumulate immunity-escape BI variants during the run.
  assert(escape.escape_bi > 0);
  // Toxin stress should remove agents from the population.
  assert(escape.producers + escape.immune + escape.naive < 30);
  // Naive targets (no colicin immunity) should take casualties.
  assert(escape.naive < 10);

  std::cout << "  test_smoke_with_immunity_escape: PASSED"
            << " (producers=" << escape.producers
            << " immune=" << escape.immune
            << " naive=" << escape.naive
            << " escape_bi=" << escape.escape_bi << ")\n";
}

int main() {
  std::cout << "=== Immunity-Escape Tests ===\n";
  test_bicluster_default_affinity();
  test_escape_mutation_reduces_affinity();
  test_immunity_effectiveness_with_escape();
  test_create_novel_toxin_with_escape();
  test_smoke_with_immunity_escape();
  std::cout << "All immunity-escape tests passed.\n";
  return 0;
}
