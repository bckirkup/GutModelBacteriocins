/* -----------------------------------------------------------------------
   GutIBM – JSON configuration parser
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_INPUT_PARSER_H
#define GUTIBM_INPUT_PARSER_H

#include "types.h"
#include "domain.h"
#include "advection.h"
#include "vbf.h"
#include "chemical_field.h"
#include "chem_environment_config.h"
#include "fur_config.h"
#include "cdi_config.h"
#include "motility_config.h"
#include "quorum_sensing_config.h"
#include "qssa_solver.h"
#include "fix_metabolism.h"
#include "fix_receptor.h"
#include "fix_bacteriocin.h"
#include "fix_conjugation.h"
#include "fix_mutation.h"
#include "fix_mechanics.h"
#include "hdf5_writer.h"

#include <functional>
#include "hdf5_reader.h"
#include "gpu_config.h"
#include "immigration_config.h"

#include <string>
#include <string_view>
#include <map>
#include <optional>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>

namespace gutibm {

struct TimeControlConfig {
  Real total_time = 86400.0;
  Real bio_dt = 60.0;
  Real output_interval = 3600.0;
};

struct AdaptiveTimestepConfig {
  bool enabled = false;
  Real min = 1.0;
  Real max = 300.0;
  Real safety = 0.8;
  Real growth_limit = 0.1;
};

struct ChemicalEnvironmentConfig {
  OxygenConfig oxygen;
  AcetateConfig acetate;
  MucinConfig mucin;
  ProteaseConfig protease;
  SiderophoreConfig siderophore;
  FerrichromeConfig ferrichrome;
};

struct CellBiologyConfig {
  FurConfig fur;
  CdiConfig cdi;
  MotilityConfig motility;
};

struct InitialPopulationConfig {
  std::string placement = "legacy";
  Real z_min = 0.0;
  Real z_max = 1.0e-6;
  Real anatomic_exclusion_floor = 20.0e-6;
  Real anatomic_exponential_scale = 40.0e-6;
  Real anatomic_outer_extent = 150.0e-6;
};

struct FixPluginsConfig {
  MetabolismConfig metabolism;
  ReceptorConfig receptor;
  BacteriocinConfig bacteriocin;
  ConjugationConfig conjugation;
  MutationConfig mutation;
  MechanicsConfig mechanics;
};

struct ClosureConfig {
  bool enforce_delivery_realization = true;
  Int zero_realization_grace_steps = 5;
  bool enforce_reaction_clip = false;
  Real reaction_clip_tolerance_fraction = 0.0;
};

struct SimulationConfig {
  TimeControlConfig time;
  AdaptiveTimestepConfig adaptive_dt;
  ClosureConfig closure;

  DomainConfig domain;
  // Chemistry layout selector; stage 2a retains global storage in both modes.
  std::string chemistry_decomposition = "replicated";
  // Scientific chemistry model variant; full retains the historical species set.
  std::string species_subset = "full";
  AdvectionConfig advection;
  VBFConfig vbf;
  Real carbon_boundary_conc = 5.0e-3;
  Real carbon_z_amplitude = 0.0;
  Real b12_initial_conc = 1.0e-3;
  std::string carbon_epithelial_boundary = "dirichlet";
  Real carbon_epithelial_transfer_coeff = 0.0;
  Real carbon_epithelial_flux = 0.0;
  std::string oxygen_epithelial_boundary = "dirichlet";
  Real oxygen_epithelial_transfer_coeff = 0.0;
  Real oxygen_epithelial_flux = 0.0;
  bool oxygen_z_gradient_enabled = true;
  std::vector<ChemicalSpec> chemicals;
  QSSAConfig qssa;
  FixPluginsConfig fixes;
  HDF5Config hdf5;

  struct CheckpointConfig {
    std::string file;
    std::string step;
  };
  CheckpointConfig checkpoint;

  // Closed midstream restart artifacts (separate from live hdf5 analysis trail).
  // When enabled, Simulation writes restart/step_NNNNNN.h5 every interval_steps.
  struct RestartConfig {
    bool enabled = false;
    std::string directory = "restart";
    Int interval_steps = 0;  // 0 = disabled even if enabled=true
  };
  RestartConfig restart;

  struct InitialStrain {
    Int type = 0;
    Int count = 0;
    Real mu_max = 5.0e-4;
    std::vector<std::string> plasmids;
    bool conjugative = false;
    uint16_t cdi_type = 0;
    uint16_t cdi_immunity = 0;
    std::map<std::string, Real, std::less<>> receptor_expression;
  };
  struct PlasmidOverride {
    std::optional<Real> retardation;
    std::optional<Real> diff_coeff;
    std::optional<Real> burst_size;
  };
  std::map<std::string, PlasmidOverride, std::less<>> plasmid_overrides;
  std::vector<InitialStrain> initial_strains;
  InitialPopulationConfig initial_population;
  std::vector<std::string> enabled_fixes;
  // Exact Fix names skipped by FixRegistry; audit labels belong below.
  std::vector<std::string> disabled_fixes;
  std::vector<std::string> disabled_mechanisms;
  uint64_t seed = 42;
  GpuConfig gpu;
  ChemicalEnvironmentConfig chem_env;
  CellBiologyConfig cell_bio;
  QuorumSensingConfig quorum_sensing;
  ImmigrationConfig immigration;
  bool profile_steps = false;

  // Spec 5 §4 — Dysbiosis safety net. When > 0, the run halts after a
  // persistent, increasing, non-decelerating density excursion above this
  // threshold. See docs/OPERATING_ENVELOPE.md. 0 disables the check.
  Real dysbiosis_threshold = 1.0e8;   // cells/mL, 0 = disabled
  Real dysbiosis_sampling_interval = 300.0;  // simulated seconds
  Int dysbiosis_sample_count = 7;            // 30-minute trajectory window
};

class InputParser {
 public:
  // Parse a simple JSON-like config file
  static SimulationConfig parse(const std::string& filename);

  // Create default config with standard gut parameters
  static SimulationConfig default_config();

  // Register optional chemical species and apply feature-flag side effects
  static void finalize_config(SimulationConfig& cfg);

  // Apply a single flat config key (used by JSON and legacy line parsers).
  // Returns true if the key was recognized and applied, false otherwise.
  static bool apply_flat_key(SimulationConfig& cfg,
                             std::string_view key,
                             const std::string& val);

  // Apply the configured strictness policy to an unrecognized key.
  static void handle_unknown_config_key(std::string_view key);

 private:
  static std::string trim(std::string_view s);
  static Real parse_real(const std::string& key, const std::string& val);
  static Int parse_int(const std::string& key, const std::string& val);
};

}  // namespace gutibm

#endif  // GUTIBM_INPUT_PARSER_H
