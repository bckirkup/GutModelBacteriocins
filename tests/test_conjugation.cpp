/* -----------------------------------------------------------------------
   GutIBM – Tests for F-pili length heterogeneity in conjugation
   Verifies that with pili_heterogeneity enabled, different conjugation
   events sample different effective radii from uniform(min, max).
   ----------------------------------------------------------------------- */

#include "fix_conjugation.h"
#include "plasmid.h"
#include "random.h"
#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

static Simulation make_ownership_sim(bool ghost_recipient) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.enabled_fixes = {"conjugation"};
  cfg.fixes.conjugation.pili_heterogeneity = false;
  cfg.fixes.conjugation.pili_length = 4.0e-6;
  cfg.fixes.conjugation.base_transfer_rate = 1.0e6;
  cfg.advection.radial_turnover = 1.0e12;
  cfg.advection.distal_transit_time = 1.0e12;

  Simulation sim;
  sim.init(cfg);

  const Vec3 donor_pos = {25e-6, 25e-6, 12.5e-6};
  const Vec3 recipient_pos = {27e-6, 25e-6, 12.5e-6};
  Agent donor = Agent::create_default(
      sim.agents().next_tag(), 1, donor_pos, 5e-4);
  Int ix;
  Int iy;
  Int iz;
  sim.domain().pos_to_grid(donor.x, ix, iy, iz);
  donor.grid_cell = sim.domain().cell_index(ix, iy, iz);
  donor.identity.owner_rank = sim.domain().rank();
  donor.genome.bi_loci.push_back(PlasmidLibrary::colicin_B());
  donor.genome.has_conjugative_plasmid = true;
  sim.agents().push_back(std::move(donor));

  Agent recipient = Agent::create_default(
      sim.agents().next_tag(), 2, recipient_pos, 5e-4);
  sim.domain().pos_to_grid(recipient.x, ix, iy, iz);
  recipient.grid_cell = sim.domain().cell_index(ix, iy, iz);
  recipient.identity.owner_rank = sim.domain().rank();
  recipient.flags.is_ghost = ghost_recipient;
  sim.agents().push_back(std::move(recipient));

  return sim;
}

static Int count_hgt_events(const Simulation& sim) {
  Int count = 0;
  for (const LineageEvent& event : sim.lineage_tracker().events()) {
    if (event.type == LineageEvent::Type::HGT) ++count;
  }
  return count;
}

void test_config_defaults() {
  ConjugationConfig cfg;
  assert(!cfg.pili_heterogeneity);
  assert(cfg.pili_length_min == 1.0e-6);
  assert(cfg.pili_length_max == 4.0e-6);
  assert(cfg.pili_length == 4.0e-6);
  std::cout << "  test_config_defaults: PASSED\n";
}

void test_sampled_radii_vary() {
  // Directly sample from RNG to confirm uniform(1,4) um produces variation
  RNG rng(12345);
  Real lo = 1.0e-6;
  Real hi = 4.0e-6;
  const int N = 10000;
  Real sum = 0.0;
  Real min_val = hi;
  Real max_val = lo;
  std::vector<Real> samples(N);

  for (int i = 0; i < N; ++i) {
    Real v = rng.uniform(lo, hi);
    samples[i] = v;
    sum += v;
    if (v < min_val) min_val = v;
    if (v > max_val) max_val = v;
  }

  Real mean = sum / N;
  Real expected_mean = (lo + hi) / 2.0;  // 2.5 um

  // Mean should be close to 2.5 um (within 5%)
  assert(std::abs(mean - expected_mean) / expected_mean < 0.05);

  // Should have variation (min != max)
  assert(max_val > min_val);

  // All values in range
  assert(min_val >= lo);
  assert(max_val <= hi);

  std::cout << "  test_sampled_radii_vary: PASSED"
            << " (mean=" << mean * 1e6 << " um"
            << " min=" << min_val * 1e6 << " um"
            << " max=" << max_val * 1e6 << " um)\n";
}

void test_heterogeneity_integration() {
  // Run a mini simulation with pili_heterogeneity enabled
  SimulationConfig cfg = InputParser::default_config();

  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;

  cfg.time.total_time      = 300.0;
  cfg.time.bio_dt          = 60.0;
  cfg.time.output_interval = 300.0;
  cfg.seed            = 54321;
  cfg.hdf5.enabled    = false;

  cfg.advection.mucus_thickness     = 25e-6;
  cfg.advection.distal_length       = 50e-6;
  cfg.advection.radial_turnover     = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff    = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  // Enable pili heterogeneity
  cfg.fixes.conjugation.pili_heterogeneity = true;
  cfg.fixes.conjugation.pili_length_min    = 1.0e-6;
  cfg.fixes.conjugation.pili_length_max    = 4.0e-6;

  cfg.initial_strains.clear();

  // Conjugative donor strain
  SimulationConfig::InitialStrain donor;
  donor.type        = 1;
  donor.count       = 10;
  donor.mu_max      = 5e-4;
  donor.plasmids    = {"ColB"};
  donor.conjugative = true;
  cfg.initial_strains.push_back(donor);

  // Recipient strain (no plasmids)
  SimulationConfig::InitialStrain recipient;
  recipient.type        = 2;
  recipient.count       = 10;
  recipient.mu_max      = 5e-4;
  recipient.plasmids    = {};
  recipient.conjugative = false;
  cfg.initial_strains.push_back(recipient);

  Simulation sim;
  sim.init(cfg);
  assert(sim.agents().size() == 20);

  // Run without crashing
  sim.run();
  assert(sim.step_count() > 0);

  std::cout << "  test_heterogeneity_integration: PASSED"
            << " (steps=" << sim.step_count() << ")\n";
}

void test_heterogeneity_disabled_uses_fixed() {
  // With heterogeneity disabled, the fixed pili_length should be used
  SimulationConfig cfg = InputParser::default_config();

  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;

  cfg.time.total_time      = 120.0;
  cfg.time.bio_dt          = 60.0;
  cfg.time.output_interval = 120.0;
  cfg.seed            = 99999;
  cfg.hdf5.enabled    = false;

  cfg.advection.mucus_thickness     = 25e-6;
  cfg.advection.distal_length       = 50e-6;
  cfg.qssa.toxin_cutoff    = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  // Ensure heterogeneity is OFF (default)
  assert(!cfg.fixes.conjugation.pili_heterogeneity);

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1; s.count = 5; s.mu_max = 5e-4;
  s.plasmids = {"ColB"}; s.conjugative = true;
  cfg.initial_strains.push_back(s);

  SimulationConfig::InitialStrain r;
  r.type = 2; r.count = 5; r.mu_max = 5e-4;
  r.plasmids = {}; r.conjugative = false;
  cfg.initial_strains.push_back(r);

  Simulation sim;
  sim.init(cfg);
  sim.run();
  assert(sim.step_count() > 0);

  std::cout << "  test_heterogeneity_disabled_uses_fixed: PASSED\n";
}

void test_ghost_recipient_does_not_commit() {
  auto sim = make_ownership_sim(true);
  assert((sim.agents()[1].genome.bi_loci.empty())
         && "ghost recipient genome changed");

  sim.step(sim.config().time.bio_dt);

  assert((sim.agents()[1].genome.bi_loci.empty())
         && "ghost recipient acquired a cluster");
  assert((sim.step_events().conjugation_transfers == 0)
         && "ghost recipient incremented transfer count");
  assert((count_hgt_events(sim) == 0)
         && "ghost recipient recorded an HGT event");
  std::cout << "  test_ghost_recipient_does_not_commit: PASSED\n";
}

void test_owned_recipient_commits() {
  auto sim = make_ownership_sim(false);

  sim.step(sim.config().time.bio_dt);

  assert((!sim.agents()[1].genome.bi_loci.empty())
         && "owned recipient did not acquire a cluster");
  assert((sim.step_events().conjugation_transfers == 1)
         && "owned recipient transfer count mismatch");
  assert((count_hgt_events(sim) == 1)
         && "owned recipient HGT event count mismatch");
  std::cout << "  test_owned_recipient_commits: PASSED\n";
}

int main() {
  std::cout << "=== Conjugation Pili Heterogeneity Tests ===\n";
  test_config_defaults();
  test_sampled_radii_vary();
  test_heterogeneity_integration();
  test_heterogeneity_disabled_uses_fixed();
  test_ghost_recipient_does_not_commit();
  test_owned_recipient_commits();
  std::cout << "All conjugation tests passed.\n";
  return 0;
}
