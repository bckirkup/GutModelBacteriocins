/* -----------------------------------------------------------------------
   GutIBM – Configuration diversity integration tests (issue #76)
   Ensures parsed/programmatic configs reach Simulation::init and produce
   distinct deterministic fingerprints — catches silent overrides to defaults.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "sim_fingerprint.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#ifndef GUTIBM_SOURCE_DIR
#define GUTIBM_SOURCE_DIR "."
#endif

using namespace gutibm;

namespace {

std::string fixture_path(const char* name) {
  return std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/" + name;
}

// Keep CI runs short while preserving config-specific fields.
void shrink_for_ci(SimulationConfig& cfg) {
  cfg.total_time = std::min(cfg.total_time, 180.0);
  cfg.output_interval = cfg.total_time;
  cfg.bio_dt = 60.0;
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;
  cfg.checkpoint.file.clear();

  constexpr Real kMaxSpan = 80e-6;
  for (int d = 0; d < 3; ++d) {
    if (cfg.domain.hi[d] > kMaxSpan) {
      cfg.domain.hi[d] = kMaxSpan;
    }
  }
  if (cfg.domain.grid_dx <= 0.0) {
    cfg.domain.grid_dx = 5e-6;
  }
  if (cfg.domain.hash_cell_size <= 0.0) {
    cfg.domain.hash_cell_size = 10e-6;
  }

  // Coarsen grid if the shrunk domain still implies a huge cell count.
  Int nx = static_cast<Int>(std::ceil(cfg.domain.hi[0] / cfg.domain.grid_dx));
  Int ny = static_cast<Int>(std::ceil(cfg.domain.hi[1] / cfg.domain.grid_dx));
  Int nz = static_cast<Int>(std::ceil(cfg.domain.hi[2] / cfg.domain.grid_dx));
  while (nx * ny * nz > 20000 && cfg.domain.grid_dx < 20e-6) {
    cfg.domain.grid_dx *= 2.0;
    nx = static_cast<Int>(std::ceil(cfg.domain.hi[0] / cfg.domain.grid_dx));
    ny = static_cast<Int>(std::ceil(cfg.domain.hi[1] / cfg.domain.grid_dx));
    nz = static_cast<Int>(std::ceil(cfg.domain.hi[2] / cfg.domain.grid_dx));
  }

  Int total_agents = 0;
  for (const auto& strain : cfg.initial_strains) {
    total_agents += strain.count;
  }
  constexpr Int kMaxAgents = 40;
  if (total_agents > kMaxAgents) {
    const Real scale = static_cast<Real>(kMaxAgents) / static_cast<Real>(total_agents);
    for (auto& strain : cfg.initial_strains) {
      strain.count = std::max<Int>(1, static_cast<Int>(std::llround(strain.count * scale)));
    }
  }
}

uint64_t run_fingerprint(const SimulationConfig& in_cfg) {
  SimulationConfig cfg = in_cfg;
  shrink_for_ci(cfg);

  Simulation sim;
  sim.init(cfg);
  sim.run();
  return test_util::simulation_fingerprint(sim);
}

SimulationConfig growth_baseline(uint64_t seed) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = seed;
  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {80e-6, 80e-6, 40e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.advection.mucus_thickness = 40e-6;
  cfg.advection.distal_length = 80e-6;
  cfg.qssa.toxin_cutoff = 40e-6;
  cfg.qssa.nutrient_cutoff = 20e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 12;
  s.mu_max = 4.5e-4;
  s.plasmids = {};
  s.conjugative = false;
  cfg.initial_strains.push_back(s);
  return cfg;
}

void report_duplicate_fingerprints(
    const std::string& name,
    uint64_t fp,
    const std::vector<std::pair<std::string, uint64_t>>& labeled) {
  std::cerr << "ERROR: duplicate fingerprint for scenario '" << name << "'\n";
  for (const auto& [other_name, other_fp] : labeled) {
    if (other_fp != fp) {
      continue;
    }
    std::cerr << "  matches '" << other_name << "'\n";
  }
}

void assert_all_distinct(const std::vector<std::pair<std::string, uint64_t>>& labeled) {
  std::set<uint64_t> seen;
  for (const auto& [name, fp] : labeled) {
    auto [it, inserted] = seen.insert(fp);
    if (inserted) {
      continue;
    }
    report_duplicate_fingerprints(name, fp, labeled);
    assert(false && "configurations should produce distinct fingerprints");
  }
}

}  // namespace

void test_fixture_configs_produce_distinct_fingerprints() {
  struct Scenario {
    const char* label;
    const char* fixture;
  };
  const Scenario scenarios[] = {
      {"strains", "parser_strains.json"},
      {"fix_subset", "parser_fixes.json"},
      {"fix_tunables", "parser_fix_tunables.json"},
      {"fmm_peristaltic", "parser_fmm_peristaltic.json"},
  };

  std::vector<std::pair<std::string, uint64_t>> fingerprints;
  for (const auto& sc : scenarios) {
    SimulationConfig cfg = InputParser::parse(fixture_path(sc.fixture));
    fingerprints.emplace_back(sc.label, run_fingerprint(cfg));
  }

  // Baseline without fixture-specific knobs.
  fingerprints.emplace_back("growth_baseline", run_fingerprint(growth_baseline(6001)));

  assert_all_distinct(fingerprints);

  std::cout << "  test_fixture_configs_produce_distinct_fingerprints: PASSED ("
            << fingerprints.size() << " scenarios)\n";
}

void test_example_configs_differ() {
  SimulationConfig single =
      InputParser::parse(std::string(GUTIBM_SOURCE_DIR) + "/examples/single_colony/input.json");
  SimulationConfig diversity =
      InputParser::parse(std::string(GUTIBM_SOURCE_DIR) + "/examples/diversity_paradox/input.json");

  uint64_t fp_single = run_fingerprint(single);
  uint64_t fp_diversity = run_fingerprint(diversity);
  assert(fp_single != fp_diversity);

  std::cout << "  test_example_configs_differ: PASSED\n";
}

void test_seed_and_fix_subset_change_outcomes() {
  uint64_t fp_a = run_fingerprint(growth_baseline(7001));
  uint64_t fp_b = run_fingerprint(growth_baseline(7002));
  assert(fp_a != fp_b);

  SimulationConfig subset = InputParser::parse(fixture_path("parser_fixes.json"));
  subset.seed = 7001;
  uint64_t fp_subset = run_fingerprint(subset);
  assert(fp_subset != fp_a);

  std::cout << "  test_seed_and_fix_subset_change_outcomes: PASSED\n";
}

void test_parsed_fix_list_is_honored() {
  SimulationConfig cfg = InputParser::parse(fixture_path("parser_fixes.json"));
  shrink_for_ci(cfg);

  Simulation sim_subset;
  sim_subset.init(cfg);
  auto names = sim_subset.fix_names();
  assert(names.size() == 2);
  assert(names[0] == "metabolism");
  assert(names[1] == "mechanics");

  for (int step = 0; step < 3; ++step) {
    sim_subset.step(cfg.bio_dt);
  }
  uint64_t fp_subset = test_util::simulation_fingerprint(sim_subset);

  SimulationConfig full = cfg;
  full.enabled_fixes.clear();

  Simulation sim_full;
  sim_full.init(full);
  for (int step = 0; step < 3; ++step) {
    sim_full.step(full.bio_dt);
  }
  uint64_t fp_full = test_util::simulation_fingerprint(sim_full);
  assert(fp_subset != fp_full);

  std::cout << "  test_parsed_fix_list_is_honored: PASSED\n";
}

void test_fix_tunables_reach_simulation() {
  SimulationConfig tuned = InputParser::parse(fixture_path("parser_fix_tunables.json"));
  assert(std::abs(tuned.receptor.kill_rate_colicin - 2e-3) < 1e-12);
  assert(tuned.conjugation.pili_heterogeneity == true);
  assert(tuned.mutation.max_bi_loci == 6);

  SimulationConfig baseline = tuned;
  baseline.receptor = InputParser::default_config().receptor;
  baseline.conjugation = InputParser::default_config().conjugation;
  baseline.mutation = InputParser::default_config().mutation;
  baseline.seed = tuned.seed;

  uint64_t fp_tuned = run_fingerprint(tuned);
  uint64_t fp_baseline = run_fingerprint(baseline);
  assert(fp_tuned != fp_baseline);

  std::cout << "  test_fix_tunables_reach_simulation: PASSED\n";
}

void test_peristaltic_toggle_changes_fingerprint() {
  SimulationConfig with_peri = growth_baseline(8001);
  with_peri.advection.peristaltic_enabled = true;
  with_peri.advection.peristaltic_period = 12.0;
  with_peri.advection.peristaltic_amplitude = 0.4;
  with_peri.advection.peristaltic_wavelength = 40e-6;

  SimulationConfig without_peri = with_peri;
  without_peri.advection.peristaltic_enabled = false;

  uint64_t fp_on = run_fingerprint(with_peri);
  uint64_t fp_off = run_fingerprint(without_peri);
  assert(fp_on != fp_off);

  std::cout << "  test_peristaltic_toggle_changes_fingerprint: PASSED\n";
}

void test_same_config_is_reproducible() {
  SimulationConfig cfg = growth_baseline(9001);
  const uint64_t fp1 = run_fingerprint(cfg);
  if (const uint64_t fp2 = run_fingerprint(cfg); fp1 != fp2) {
    std::cerr << "ERROR: reproducibility failure fp1=" << fp1 << " fp2=" << fp2 << "\n";
    assert(false);
  }

  std::cout << "  test_same_config_is_reproducible: PASSED\n";
}

int main() {
  std::cout << "=== Config Diversity Tests ===\n";
  test_fixture_configs_produce_distinct_fingerprints();
  test_example_configs_differ();
  test_seed_and_fix_subset_change_outcomes();
  test_parsed_fix_list_is_honored();
  test_fix_tunables_reach_simulation();
  test_peristaltic_toggle_changes_fingerprint();
  test_same_config_is_reproducible();
  std::cout << "All config diversity tests passed.\n";
  return 0;
}
