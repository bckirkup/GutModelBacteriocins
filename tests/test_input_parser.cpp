/* -----------------------------------------------------------------------
   GutIBM – Example input file parsing tests
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "simulation.h"
#include "error.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

#ifndef GUTIBM_SOURCE_DIR
#define GUTIBM_SOURCE_DIR "."
#endif

using namespace gutibm;

void test_single_colony_example() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/examples/single_colony/input.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(std::abs(cfg.time.total_time - 86400.0) < 1e-6);
  assert(std::abs(cfg.time.bio_dt - 60.0) < 1e-6);
  assert(cfg.seed == 12345);
  assert(std::abs(cfg.domain.hi[0] - 0.001) < 1e-12);
  assert(cfg.hdf5.filename == "single_colony_output.h5");
  assert(cfg.hdf5.schedule.agents == 60);
  assert(cfg.hdf5.schedule.provenance == 0);
  std::cout << "  test_single_colony_example: PASSED\n";
}

void test_diversity_paradox_example() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/examples/diversity_paradox/input.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(std::abs(cfg.time.total_time - 604800.0) < 1e-6);
  assert(cfg.seed == 42);
  assert(std::abs(cfg.domain.hi[0] - 0.002) < 1e-12);
  assert(cfg.hdf5.filename == "diversity_paradox_output.h5");
  assert(cfg.qssa.use_fmm == true);
  assert(std::abs(cfg.qssa.fmm_theta - 0.5) < 1e-12);
  assert(std::abs(cfg.qssa.toxin_cutoff - 200e-6) < 1e-15);
  assert(std::abs(cfg.qssa.nutrient_cutoff - 50e-6) < 1e-15);
  std::cout << "  test_diversity_paradox_example: PASSED\n";
}

void test_single_colony_peristaltic_keys() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/examples/single_colony/input.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(cfg.advection.peristaltic_enabled == true);
  assert(std::abs(cfg.advection.peristaltic_period - 20.0) < 1e-12);
  assert(std::abs(cfg.advection.peristaltic_amplitude - 0.5) < 1e-12);
  assert(std::abs(cfg.advection.peristaltic_wavelength - 0.001) < 1e-15);
  std::cout << "  test_single_colony_peristaltic_keys: PASSED\n";
}

void test_fmm_peristaltic_fixture() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/parser_fmm_peristaltic.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(std::abs(cfg.time.total_time - 3600.0) < 1e-6);
  assert(cfg.seed == 99);
  assert(cfg.advection.peristaltic_enabled == true);
  assert(std::abs(cfg.advection.peristaltic_period - 15.0) < 1e-12);
  assert(std::abs(cfg.advection.peristaltic_amplitude - 0.3) < 1e-12);
  assert(std::abs(cfg.advection.peristaltic_wavelength - 0.0005) < 1e-15);
  assert(cfg.qssa.use_fmm == true);
  assert(std::abs(cfg.qssa.fmm_theta - 0.3) < 1e-12);
  assert(cfg.qssa.fmm_expansion_order == 2);
  assert(std::abs(cfg.qssa.toxin_cutoff - 100e-6) < 1e-15);
  assert(std::abs(cfg.qssa.nutrient_cutoff - 25e-6) < 1e-15);
  std::cout << "  test_fmm_peristaltic_fixture: PASSED\n";
}

void test_strain_fixture() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/parser_strains.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(cfg.initial_strains.size() == 2);

  const auto& resident = cfg.initial_strains[0];
  assert(resident.type == 1);
  assert(resident.count == 12);
  assert(std::abs(resident.mu_max - 5.5e-4) < 1e-12);
  assert(resident.plasmids.size() == 1);
  assert(resident.plasmids[0] == "ColE1");
  assert(resident.conjugative == false);

  const auto& immigrant = cfg.initial_strains[1];
  assert(immigrant.type == 2);
  assert(immigrant.count == 4);
  assert(immigrant.plasmids.empty());
  std::cout << "  test_strain_fixture: PASSED\n";
}

void test_immigration_fixture() {
  const std::string path = std::string(GUTIBM_SOURCE_DIR) +
                           "/tests/fixtures/parser_immigration.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(cfg.immigration.enabled);
  assert(cfg.immigration.count == 2);
  assert(cfg.immigration.strain_index == 1);
  assert(cfg.immigration.placement == "z_slab");
  assert(cfg.immigration.schedule == "pulse");
  assert(std::abs(cfg.immigration.z_min - 1e-6) < 1e-15);
  assert(std::abs(cfg.immigration.z_max - 20e-6) < 1e-15);
  std::cout << "  test_immigration_fixture: PASSED\n";
}

void test_diversity_paradox_strains() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/examples/diversity_paradox/input.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(cfg.initial_strains.size() == 2);
  assert(cfg.initial_strains[0].plasmids.size() == 2);
  assert(cfg.initial_strains[0].plasmids[0] == "ColE1");
  assert(cfg.initial_strains[0].plasmids[1] == "ColB");
  assert(cfg.initial_strains[1].count == 100);
  std::cout << "  test_diversity_paradox_strains: PASSED\n";
}

void test_strain_spawn_integration() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/parser_strains.json";
  SimulationConfig cfg = InputParser::parse(path);
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.hdf5.enabled = false;
  cfg.time.total_time = 1.0;

  Simulation sim;
  sim.init(cfg);

  Int with_bi = 0;
  for (const Agent& a : sim.agents()) {
    if (!a.genome.bi_loci.empty()) ++with_bi;
  }
  assert(with_bi > 0);
  std::cout << "  test_strain_spawn_integration: PASSED\n";
}

void test_fixes_fixture() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/parser_fixes.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(cfg.enabled_fixes.size() == 2);
  assert(cfg.enabled_fixes[0] == "metabolism");
  assert(cfg.enabled_fixes[1] == "mechanics");
  assert(cfg.initial_strains.size() == 1);
  assert(cfg.initial_strains[0].count == 10);

  cfg.hdf5.enabled = false;
  Simulation sim;
  sim.init(cfg);
  auto names = sim.fix_names();
  assert(names.size() == 2);
  assert(names[0] == "metabolism");
  assert(names[1] == "mechanics");
  std::cout << "  test_fixes_fixture: PASSED\n";
}

void test_json_document_parser() {
  const std::string json = R"({
    "_comment": "inline JSON document test \u2014 batch rewrite \ud83d\ude80",
    "total_time": 1234,
    "seed": 99,
    "hdf5_file": "batch\u0020output.h5",
    "peristaltic_enabled": true,
    "initial_strains": [
      {"type": 1, "count": 3, "mu_max": 5e-4, "plasmids": ["ColE1"], "conjugative": false}
    ],
    "fixes": ["metabolism"]
  })";

  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/_inline_json_doc.json";
  {
    std::ofstream out(path);
    out << json;
  }

  SimulationConfig cfg = InputParser::parse(path);
  assert(std::abs(cfg.time.total_time - 1234.0) < 1e-6);
  assert(cfg.seed == 99);
  assert(cfg.hdf5.filename == "batch output.h5");
  assert(cfg.advection.peristaltic_enabled == true);
  assert(cfg.initial_strains.size() == 1);
  assert(cfg.initial_strains[0].plasmids[0] == "ColE1");
  assert(cfg.enabled_fixes.size() == 1);
  assert(cfg.enabled_fixes[0] == "metabolism");
  std::remove(path.c_str());
  std::cout << "  test_json_document_parser: PASSED\n";
}

void test_malformed_numeric_warnings_json() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/parser_bad_numeric.json";

  std::stringstream err;
  std::streambuf* old_err = std::cerr.rdbuf(err.rdbuf());
  SimulationConfig cfg = InputParser::parse(path);
  std::cerr.rdbuf(old_err);

  assert(std::abs(cfg.domain.hi[0]) < 1e-15);
  assert(cfg.seed == 0);

  const std::string warnings = err.str();
  assert(warnings.find("domain_x") != std::string::npos);
  assert(warnings.find("1mm") != std::string::npos);
  assert(warnings.find("seed") != std::string::npos);
  assert(warnings.find("not_a_number") != std::string::npos);
  std::cout << "  test_malformed_numeric_warnings_json: PASSED\n";
}

void test_malformed_numeric_warnings_legacy() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/parser_bad_numeric.legacy";

  std::stringstream err;
  std::streambuf* old_err = std::cerr.rdbuf(err.rdbuf());
  SimulationConfig cfg = InputParser::parse(path);
  std::cerr.rdbuf(old_err);

  assert(std::abs(cfg.domain.hi[0]) < 1e-15);
  assert(cfg.seed == 0);

  const std::string warnings = err.str();
  assert(warnings.find("domain_x") != std::string::npos);
  assert(warnings.find("1mm") != std::string::npos);
  std::cout << "  test_malformed_numeric_warnings_legacy: PASSED\n";
}

void test_chem_env_fixture() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/parser_chem_env.json";
  SimulationConfig cfg = InputParser::parse(path);

  assert(cfg.chem_env.oxygen.enabled == true);
  assert(std::abs(cfg.chem_env.oxygen.epithelial_conc - 60e-6) < 1e-15);
  assert(std::abs(cfg.chem_env.oxygen.boost_max - 1.8) < 1e-12);

  assert(cfg.chem_env.acetate.enabled == true);
  assert(std::abs(cfg.chem_env.acetate.vbf_production - 2e-3) < 1e-12);
  assert(std::abs(cfg.chem_env.acetate.overflow_rate - 2e-15) < 1e-18);

  assert(cfg.chem_env.mucin.enabled == true);
  assert(std::abs(cfg.chem_env.mucin.initial_conc - 2e-2) < 1e-12);
  assert(std::abs(cfg.chem_env.mucin.k_liberation - 2e-4) < 1e-12);

  assert(cfg.chem_env.protease.enabled == false);
  assert(std::abs(cfg.chem_env.protease.default_half_life - 1200.0) < 1e-6);
  assert(std::abs(cfg.chem_env.protease.dilution_rate - 2e-4) < 1e-12);

  assert(std::abs(cfg.chem_env.siderophore.secretion_rate - 1e-5) < 1e-12);
  assert(std::abs(cfg.chem_env.siderophore.Vmax_reimport - 1e-5) < 1e-12);

  assert(cfg.chem_env.ferrichrome.enabled == true);
  assert(std::abs(cfg.chem_env.ferrichrome.initial_conc - 2e-6) < 1e-15);
  assert(std::abs(cfg.chem_env.ferrichrome.boundary_conc - 4e-6) < 1e-15);
  assert(std::abs(cfg.carbon_boundary_conc - 7e-3) < 1e-15);
  assert(cfg.chemistry_decomposition == "slab");

  bool has_oxygen = false;
  bool has_acetate = false;
  bool has_mucin = false;
  bool has_ferrichrome = false;
  for (const auto& spec : cfg.chemicals) {
    if (spec.name == "oxygen") {
      has_oxygen = true;
      assert(std::abs(spec.boundary_conc - 60e-6) < 1e-15);
    } else if (spec.name == "acetate") {
      has_acetate = true;
    } else if (spec.name == "mucin") {
      has_mucin = true;
      assert(std::abs(spec.initial_conc - 2e-2) < 1e-12);
    } else if (spec.name == "ferrichrome") {
      has_ferrichrome = true;
      assert(std::abs(spec.initial_conc - 2e-6) < 1e-15);
      assert(std::abs(spec.boundary_conc - 4e-6) < 1e-15);
    } else if (spec.name == "carbon") {
      assert(std::abs(spec.boundary_conc - 7e-3) < 1e-15);
    }
  }
  assert(has_oxygen);
  assert(has_acetate);
  assert(has_mucin);
  assert(has_ferrichrome);

  std::cout << "  test_chem_env_fixture: PASSED\n";
}

void test_operating_envelope_fixture() {
  const std::string path = std::string(GUTIBM_SOURCE_DIR) +
                           "/tests/fixtures/parser_operating_envelope.json";
  const SimulationConfig cfg = InputParser::parse(path);
  assert(std::abs(cfg.dysbiosis_threshold - 1.0e8) < 1.0);
  assert(std::abs(cfg.dysbiosis_sampling_interval - 300.0) < 1e-12);
  assert(cfg.dysbiosis_sample_count == 7);
  std::cout << "  test_operating_envelope_fixture: PASSED\n";
}

void test_fix_tunables_fixture() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/parser_fix_tunables.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(std::abs(cfg.fixes.receptor.kd_colicinE_btuB - 1e-9) < 1e-15);
  assert(std::abs(cfg.fixes.receptor.kill_rate_colicin - 2e-3) < 1e-12);
  assert(std::abs(cfg.fixes.receptor.immunity_factor - 0.0005) < 1e-12);
  assert(std::abs(cfg.fixes.conjugation.base_transfer_rate - 2e-4) < 1e-12);
  assert(std::abs(cfg.fixes.conjugation.plasmid_copy_cost - 0.15) < 1e-12);
  assert(cfg.fixes.conjugation.pili_heterogeneity == true);
  assert(std::abs(cfg.fixes.conjugation.pili_length_min - 2e-6) < 1e-15);
  assert(std::abs(cfg.fixes.conjugation.pili_length_max - 3e-6) < 1e-15);
  assert(std::abs(cfg.fixes.mutation.bi_duplication_rate - 1e-4) < 1e-12);
  assert(cfg.fixes.mutation.max_bi_loci == 6);
  assert(std::abs(cfg.fixes.mutation.immunity_escape_prob - 0.75) < 1e-12);
  assert(cfg.hdf5.schedule.provenance == 7);
  std::cout << "  test_fix_tunables_fixture: PASSED\n";
}

void test_unknown_key_warning_json() {
  const std::string json = R"({
    "_comment": "unknown key warning test",
    "total_time": 100,
    "bogus_key_xyz": 5,
    "another.unknown_key": true,
    "siderophore.recapture_fraction": 0.5
  })";

  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/_unknown_key_doc.json";
  {
    std::ofstream out(path);
    out << json;
  }

  std::stringstream err;
  std::streambuf* old_err = std::cerr.rdbuf(err.rdbuf());
  SimulationConfig cfg = InputParser::parse(path);
  std::cerr.rdbuf(old_err);
  std::remove(path.c_str());

  // Known key still applied.
  assert(std::abs(cfg.time.total_time - 100.0) < 1e-6);

  const std::string warnings = err.str();
  // Unknown keys are surfaced.
  assert(warnings.find("bogus_key_xyz") != std::string::npos);
  assert(warnings.find("another.unknown_key") != std::string::npos);
  assert(warnings.find("siderophore.recapture_fraction") != std::string::npos);
  // Comment keys and recognized keys are not flagged.
  assert(warnings.find("_comment") == std::string::npos);
  assert(warnings.find("'total_time'") == std::string::npos);
  std::cout << "  test_unknown_key_warning_json: PASSED\n";
}

void test_unknown_key_warning_legacy() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/_unknown_key.legacy";
  {
    std::ofstream out(path);
    out << "total_time: 200\n";
    out << "made_up_key: 3\n";
    out << "_comment: ignore me\n";
  }

  std::stringstream err;
  std::streambuf* old_err = std::cerr.rdbuf(err.rdbuf());
  SimulationConfig cfg = InputParser::parse(path);
  std::cerr.rdbuf(old_err);
  std::remove(path.c_str());

  assert(std::abs(cfg.time.total_time - 200.0) < 1e-6);
  const std::string warnings = err.str();
  assert(warnings.find("made_up_key") != std::string::npos);
  assert(warnings.find("_comment") == std::string::npos);
  std::cout << "  test_unknown_key_warning_legacy: PASSED\n";
}

void test_gpu_enabled_fixture() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/parser_gpu.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(cfg.gpu.enabled == true);
  assert(cfg.gpu.device_id == 0);
  std::cout << "  test_gpu_enabled_fixture: PASSED\n";
}

void test_strict_config_aborts_on_bad_numeric() {
  std::string path = std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/parser_bad_numeric.json";

  const char* previous = std::getenv("GUTIBM_STRICT_CONFIG");
  std::string saved;
  if (previous) saved = previous;
  setenv("GUTIBM_STRICT_CONFIG", "1", 1);

  bool threw = false;
  try {
    (void)InputParser::parse(path);
  } catch (const ConfigError&) {
    threw = true;
  }

  if (saved.empty()) {
    unsetenv("GUTIBM_STRICT_CONFIG");
  } else {
    setenv("GUTIBM_STRICT_CONFIG", saved.c_str(), 1);
  }

  assert(threw);
  std::cout << "  test_strict_config_aborts_on_bad_numeric: PASSED\n";
}

void test_burst_release_tau_must_be_positive() {
  SimulationConfig cfg = InputParser::default_config();
  bool threw = false;
  try {
    (void)InputParser::apply_flat_key(cfg, "burst_release_tau", "0");
  } catch (const ConfigError&) {
    threw = true;
  }
  assert(threw);
  std::cout << "  test_burst_release_tau_must_be_positive: PASSED\n";
}

void test_chemistry_stride_requires_positive_integer() {
  for (const std::string value : {"0", "-2", "1.5", "abc"}) {
    SimulationConfig cfg = InputParser::default_config();
    bool threw = false;
    try {
      (void)InputParser::apply_flat_key(cfg, "chemistry_stride_x", value);
    } catch (const ConfigError&) {
      threw = true;
    }
    assert(threw);
  }
  std::cout << "  test_chemistry_stride_requires_positive_integer: PASSED\n";
}

void test_grid_halo_width_fixture() {
  const std::string path = std::string(GUTIBM_SOURCE_DIR)
      + "/tests/fixtures/parser_grid_halo_width.json";
  const SimulationConfig cfg = InputParser::parse(path);
  assert(cfg.domain.grid_halo_width == 2);
  std::cout << "  test_grid_halo_width_fixture: PASSED\n";
}

void test_grid_halo_width_requires_positive_integer() {
  for (const std::string value : {"0", "-2", "1.5", "abc"}) {
    SimulationConfig cfg = InputParser::default_config();
    bool threw = false;
    try {
      (void)InputParser::apply_flat_key(cfg, "grid_halo_width", value);
    } catch (const ConfigError&) {
      threw = true;
    }
    assert(threw);
  }
  std::cout << "  test_grid_halo_width_requires_positive_integer: PASSED\n";
}

void test_json_grid_halo_width_errors_are_not_swallowed() {
  const std::string path = std::string(GUTIBM_SOURCE_DIR)
      + "/tests/fixtures/parser_bad_grid_halo_width.json";
  bool threw = false;
  try {
    (void)InputParser::parse(path);
  } catch (const ConfigError&) {
    threw = true;
  }
  assert(threw);
  std::cout << "  test_json_grid_halo_width_errors_are_not_swallowed: PASSED\n";
}

void test_json_chemistry_stride_errors_are_not_swallowed() {
  const std::string path = std::string(GUTIBM_SOURCE_DIR)
      + "/tests/fixtures/parser_bad_chemistry_stride.json";
  bool threw = false;
  try {
    (void)InputParser::parse(path);
  } catch (const ConfigError&) {
    threw = true;
  }
  assert(threw);
  std::cout << "  test_json_chemistry_stride_errors_are_not_swallowed: PASSED\n";
}

void test_chemistry_decomposition_rejects_unimplemented_modes() {
  SimulationConfig cfg = InputParser::default_config();
  assert(cfg.chemistry_decomposition == "replicated");
  assert(InputParser::apply_flat_key(cfg, "chemistry_decomposition", "slab"));
  assert(cfg.chemistry_decomposition == "slab");
  for (const std::string mode : {"interface", "unknown"}) {
    bool threw = false;
    try {
      (void)InputParser::apply_flat_key(cfg, "chemistry_decomposition", mode);
    } catch (const ConfigError&) {
      threw = true;
    }
    assert(threw);
  }
  std::cout << "  test_chemistry_decomposition_rejects_unimplemented_modes: PASSED\n";
}

void test_toxin_evaluation_rejects_unknown_modes() {
  SimulationConfig cfg = InputParser::default_config();
  assert(cfg.qssa.toxin_evaluation == "grid");
  assert(InputParser::apply_flat_key(cfg, "toxin_evaluation", "agents"));
  assert(cfg.qssa.toxin_evaluation == "agents");
  bool threw = false;
  try {
    (void)InputParser::apply_flat_key(
        cfg, "chemistry.toxin_evaluation", "unknown");
  } catch (const ConfigError&) {
    threw = true;
  }
  assert(threw);

  const std::string path = std::string(GUTIBM_SOURCE_DIR)
      + "/tests/fixtures/parser_bad_toxin_evaluation.json";
  threw = false;
  try {
    (void)InputParser::parse(path);
  } catch (const ConfigError&) {
    threw = true;
  }
  assert(threw);
  std::cout << "  test_toxin_evaluation_rejects_unknown_modes: PASSED\n";
}

void test_toxin_lumping_rejects_unknown_modes() {
  SimulationConfig cfg = InputParser::default_config();
  assert(cfg.qssa.toxin_lumping == "per_receptor");
  assert(InputParser::apply_flat_key(cfg, "toxin_lumping", "lumped"));
  assert(cfg.qssa.toxin_lumping == "lumped");
  bool threw = false;
  try {
    (void)InputParser::apply_flat_key(
        cfg, "chemistry.toxin_lumping", "unknown");
  } catch (const ConfigError&) {
    threw = true;
  }
  assert(threw);

  const std::string path = std::string(GUTIBM_SOURCE_DIR)
      + "/tests/fixtures/parser_bad_toxin_lumping.json";
  threw = false;
  try {
    (void)InputParser::parse(path);
  } catch (const ConfigError&) {
    threw = true;
  }
  assert(threw);
  std::cout << "  test_toxin_lumping_rejects_unknown_modes: PASSED\n";
}

void test_species_subset_rejects_unknown_modes() {
  SimulationConfig cfg = InputParser::default_config();
  assert(cfg.species_subset == "full");
  assert(InputParser::apply_flat_key(
      cfg, "chemistry.species_subset", "nutrient_only"));
  assert(cfg.species_subset == "nutrient_only");
  bool threw = false;
  try {
    (void)InputParser::apply_flat_key(
        cfg, "chemistry.species_subset", "unknown");
  } catch (const ConfigError&) {
    threw = true;
  }
  assert(threw);

  const std::string path = std::string(GUTIBM_SOURCE_DIR)
      + "/tests/fixtures/parser_bad_species_subset.json";
  threw = false;
  try {
    (void)InputParser::parse(path);
  } catch (const ConfigError&) {
    threw = true;
  }
  assert(threw);
  std::cout << "  test_species_subset_rejects_unknown_modes: PASSED\n";
}

int main() {
  std::cout << "=== Input Parser Example Tests ===\n";
  test_single_colony_example();
  test_single_colony_peristaltic_keys();
  test_diversity_paradox_example();
  test_fmm_peristaltic_fixture();
  test_strain_fixture();
  test_immigration_fixture();
  test_diversity_paradox_strains();
  test_strain_spawn_integration();
  test_fixes_fixture();
  test_chem_env_fixture();
  test_operating_envelope_fixture();
  test_fix_tunables_fixture();
  test_json_document_parser();
  test_malformed_numeric_warnings_json();
  test_malformed_numeric_warnings_legacy();
  test_unknown_key_warning_json();
  test_unknown_key_warning_legacy();
  test_gpu_enabled_fixture();
  test_strict_config_aborts_on_bad_numeric();
  test_burst_release_tau_must_be_positive();
  test_chemistry_stride_requires_positive_integer();
  test_grid_halo_width_fixture();
  test_grid_halo_width_requires_positive_integer();
  test_json_grid_halo_width_errors_are_not_swallowed();
  test_json_chemistry_stride_errors_are_not_swallowed();
  test_chemistry_decomposition_rejects_unimplemented_modes();
  test_toxin_evaluation_rejects_unknown_modes();
  test_toxin_lumping_rejects_unknown_modes();
  test_species_subset_rejects_unknown_modes();
  std::cout << "All input parser example tests passed.\n";
  return 0;
}
