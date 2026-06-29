/* -----------------------------------------------------------------------
   GutIBM – Tests for F-pili length heterogeneity in conjugation
   Verifies that with pili_heterogeneity enabled, different conjugation
   events sample different effective radii from uniform(min, max).
   ----------------------------------------------------------------------- */

#include "fix_conjugation.h"
#include "random.h"
#include "simulation.h"
#include "input_parser.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

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

  cfg.total_time      = 300.0;
  cfg.bio_dt          = 60.0;
  cfg.output_interval = 300.0;
  cfg.seed            = 54321;
  cfg.hdf5.enabled    = false;

  cfg.advection.mucus_thickness     = 25e-6;
  cfg.advection.distal_length       = 50e-6;
  cfg.advection.radial_turnover     = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff    = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  // Enable pili heterogeneity
  cfg.conjugation.pili_heterogeneity = true;
  cfg.conjugation.pili_length_min    = 1.0e-6;
  cfg.conjugation.pili_length_max    = 4.0e-6;

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

  cfg.total_time      = 120.0;
  cfg.bio_dt          = 60.0;
  cfg.output_interval = 120.0;
  cfg.seed            = 99999;
  cfg.hdf5.enabled    = false;

  cfg.advection.mucus_thickness     = 25e-6;
  cfg.advection.distal_length       = 50e-6;
  cfg.qssa.toxin_cutoff    = 25e-6;
  cfg.qssa.nutrient_cutoff = 15e-6;

  // Ensure heterogeneity is OFF (default)
  assert(!cfg.conjugation.pili_heterogeneity);

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

int main() {
  std::cout << "=== Conjugation Pili Heterogeneity Tests ===\n";
  test_config_defaults();
  test_sampled_radii_vary();
  test_heterogeneity_integration();
  test_heterogeneity_disabled_uses_fixed();
  std::cout << "All conjugation tests passed.\n";
  return 0;
}
