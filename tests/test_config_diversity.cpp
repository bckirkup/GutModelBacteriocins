/* -----------------------------------------------------------------------
   GutIBM – Configuration diversity integration tests (issue #76)
   Ensures parsed/programmatic configs reach Simulation::init and produce
   distinct deterministic fingerprints — catches silent overrides to defaults.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "sim_fingerprint.h"
#include "species_names.h"
#include "config_json.h"

#include <array>
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
  cfg.time.total_time = std::min(cfg.time.total_time, 180.0);
  cfg.time.output_interval = cfg.time.total_time;
  cfg.time.bio_dt = 60.0;
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
  auto nx = static_cast<Int>(std::ceil(cfg.domain.hi[0] / cfg.domain.grid_dx));
  auto ny = static_cast<Int>(std::ceil(cfg.domain.hi[1] / cfg.domain.grid_dx));
  auto nz = static_cast<Int>(std::ceil(cfg.domain.hi[2] / cfg.domain.grid_dx));
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
  constexpr std::array scenarios = {
      Scenario{"strains", "parser_strains.json"},
      Scenario{"fix_subset", "parser_fixes.json"},
      Scenario{"fix_tunables", "parser_fix_tunables.json"},
      Scenario{"fmm_peristaltic", "parser_fmm_peristaltic.json"},
      Scenario{"immigration", "parser_immigration.json"},
      Scenario{"initial_population", "parser_initial_population.json"},
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

void test_ros_driver_config_diversity() {
  SimulationConfig ambient = growth_baseline(7011);
  SimulationConfig funded = ambient;
  funded.chem_env.oxygen.enabled = true;
  funded.chem_env.oxygen.ros_driver = "funded";
  funded.chem_env.oxygen.delivery_uptake_enabled = true;
  funded.chem_env.oxygen.k_ROS_respiratory = 1.0;
  funded.fixes.metabolism.uptake_limit = "delivery";
  InputParser::finalize_config(funded);

  assert(ambient.chem_env.oxygen.ros_driver
         != funded.chem_env.oxygen.ros_driver);
  assert(ConfigJson::serialize_document(ambient)
         != ConfigJson::serialize_document(funded));
  assert(run_fingerprint(ambient) != run_fingerprint(funded));
  std::cout << "  test_ros_driver_config_diversity: PASSED\n";
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

  for (int step = 0; step < 6; ++step) {
    sim_subset.step(cfg.time.bio_dt);
  }
  uint64_t fp_subset = test_util::simulation_fingerprint(sim_subset);

  SimulationConfig full = cfg;
  full.enabled_fixes.clear();

  Simulation sim_full;
  sim_full.init(full);
  for (int step = 0; step < 6; ++step) {
    sim_full.step(full.time.bio_dt);
  }
  uint64_t fp_full = test_util::simulation_fingerprint(sim_full);
  assert(fp_subset != fp_full);

  std::cout << "  test_parsed_fix_list_is_honored: PASSED\n";
}

void test_fix_tunables_reach_simulation() {
  SimulationConfig tuned = InputParser::parse(fixture_path("parser_fix_tunables.json"));
  assert(std::abs(tuned.fixes.receptor.kill_rate_colicin - 2e-3) < 1e-12);
  assert(tuned.fixes.conjugation.pili_heterogeneity == true);
  assert(tuned.fixes.mutation.max_bi_loci == 6);
  assert(tuned.cell_bio.fur.enabled == true);

  SimulationConfig baseline = tuned;
  baseline.cell_bio.fur.Km = InputParser::default_config().cell_bio.fur.Km;
  baseline.fixes.receptor = InputParser::default_config().fixes.receptor;
  baseline.fixes.conjugation = InputParser::default_config().fixes.conjugation;
  baseline.fixes.mutation = InputParser::default_config().fixes.mutation;
  baseline.seed = tuned.seed;

  shrink_for_ci(tuned);
  shrink_for_ci(baseline);

  Simulation sim_tuned;
  Simulation sim_baseline;
  sim_tuned.init(tuned);
  sim_baseline.init(baseline);
  sim_tuned.step(tuned.time.bio_dt);
  sim_baseline.step(baseline.time.bio_dt);

  assert(sim_tuned.agents().size() > 0);
  assert(sim_baseline.agents().size() > 0);
  const Real mu_tuned = sim_tuned.agents()[0].mu_realized;
  const Real mu_baseline = sim_baseline.agents()[0].mu_realized;
  assert(std::abs(mu_tuned - mu_baseline) > 1e-10);

  auto count_bi_loci = [](SimulationConfig cfg) {
    shrink_for_ci(cfg);
    cfg.enabled_fixes = {"metabolism", "mutation"};
    cfg.initial_strains[0].count = 4;
    Simulation sim;
    sim.init(cfg);
    for (int step = 0; step < 10; ++step) {
      sim.step(cfg.time.bio_dt);
    }
    int total = 0;
    for (const Agent& a : sim.agents()) {
      if (a.state == PhenoState::DEAD) continue;
      total += static_cast<int>(a.genome.bi_loci.size());
    }
    return total;
  };

  SimulationConfig high_dup = tuned;
  high_dup.fixes.mutation.bi_duplication_rate = 0.5;
  high_dup.seed = 8801;

  SimulationConfig no_dup = tuned;
  no_dup.fixes.mutation.bi_duplication_rate = 0.0;
  no_dup.seed = 8801;

  const int bi_high = count_bi_loci(high_dup);
  const int bi_none = count_bi_loci(no_dup);
  assert(bi_high > bi_none);

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

void test_washout_trap_modes_change_outcome() {
  SimulationConfig imposed = growth_baseline(8011);
  imposed.advection.washout_trap = WashoutTrapMode::IMPOSED;
  imposed.advection.radial_turnover = 1.0e7;
  imposed.initial_strains.clear();
  imposed.initial_strains.push_back({2, 12, 1e-9, {}, false, 0, 0, {}});
  SimulationConfig emergent = imposed;
  emergent.advection.washout_trap = WashoutTrapMode::EMERGENT;
  assert(run_fingerprint(imposed) != run_fingerprint(emergent));
  std::cout << "  test_washout_trap_modes_change_outcome: PASSED\n";
}

void test_carbon_epithelial_boundary_changes_outcome() {
  SimulationConfig dirichlet = growth_baseline(8021);
  for (auto& chemical : dirichlet.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.z_gradient_enabled = false;
    }
  }

  SimulationConfig robin_low = dirichlet;
  robin_low.carbon_epithelial_boundary = "robin";
  robin_low.carbon_epithelial_transfer_coeff = 1.0e-6;

  SimulationConfig robin_high = robin_low;
  robin_high.carbon_epithelial_transfer_coeff = 1.0e-4;

  SimulationConfig flux = dirichlet;
  flux.carbon_epithelial_boundary = "flux";
  flux.carbon_epithelial_flux = 1.0e-10;

  const uint64_t fp_dirichlet = run_fingerprint(dirichlet);
  const uint64_t fp_robin_low = run_fingerprint(robin_low);
  const uint64_t fp_robin_high = run_fingerprint(robin_high);
  const uint64_t fp_flux = run_fingerprint(flux);
  assert(fp_dirichlet != fp_robin_low);
  assert(fp_robin_low != fp_robin_high);
  assert(fp_dirichlet != fp_flux);

  std::cout << "  test_carbon_epithelial_boundary_changes_outcome: PASSED\n";
}

void test_uptake_limit_changes_outcome() {
  SimulationConfig none = growth_baseline(8031);
  none.initial_strains[0].count = 24;
  none.fixes.metabolism.maintenance_rate = 0.0;
  for (auto& chemical : none.chemicals) {
    if (chemical.name == species::CARBON) {
      chemical.z_gradient_enabled = false;
      chemical.retardation = 3.0e3;
      chemical.initial_conc = 5.0e-2;
      chemical.boundary_conc = 5.0e-2;
    }
  }
  SimulationConfig sherwood = none;
  sherwood.fixes.metabolism.uptake_limit = "sherwood";
  sherwood.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Sherwood;
  SimulationConfig voxel = none;
  voxel.fixes.metabolism.uptake_limit = "voxel";
  voxel.fixes.metabolism.uptake_limit_mode = UptakeLimitMode::Voxel;

  const uint64_t fp_none = run_fingerprint(none);
  const uint64_t fp_sherwood = run_fingerprint(sherwood);
  const uint64_t fp_voxel = run_fingerprint(voxel);
  assert(fp_none != fp_sherwood);
  assert(fp_sherwood != fp_voxel);

  std::cout << "  test_uptake_limit_changes_outcome: PASSED\n";
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
  test_ros_driver_config_diversity();
  test_parsed_fix_list_is_honored();
  test_fix_tunables_reach_simulation();
  test_peristaltic_toggle_changes_fingerprint();
  test_washout_trap_modes_change_outcome();
  test_carbon_epithelial_boundary_changes_outcome();
  test_uptake_limit_changes_outcome();
  test_same_config_is_reproducible();
  std::cout << "All config diversity tests passed.\n";
  return 0;
}
