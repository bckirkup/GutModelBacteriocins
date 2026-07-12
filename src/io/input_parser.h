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
#include "qssa_solver.h"
#include "fix_metabolism.h"
#include "fix_receptor.h"
#include "fix_bacteriocin.h"
#include "fix_conjugation.h"
#include "fix_mutation.h"
#include "fix_mechanics.h"
#include "hdf5_writer.h"
#include "hdf5_reader.h"
#include "gpu_config.h"

#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
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
};

struct CellBiologyConfig {
  FurConfig fur;
  CdiConfig cdi;
  MotilityConfig motility;
};

struct FixPluginsConfig {
  MetabolismConfig metabolism;
  ReceptorConfig receptor;
  BacteriocinConfig bacteriocin;
  ConjugationConfig conjugation;
  MutationConfig mutation;
  MechanicsConfig mechanics;
};

struct SimulationConfig {
  TimeControlConfig time;
  AdaptiveTimestepConfig adaptive_dt;

  DomainConfig domain;
  AdvectionConfig advection;
  VBFConfig vbf;
  std::vector<ChemicalSpec> chemicals;
  QSSAConfig qssa;
  FixPluginsConfig fixes;
  HDF5Config hdf5;

  struct CheckpointConfig {
    std::string file;
    std::string step;
  };
  CheckpointConfig checkpoint;

  struct InitialStrain {
    Int type = 0;
    Int count = 0;
    Real mu_max = 5.0e-4;
    std::vector<std::string> plasmids;
    bool conjugative = false;
    uint16_t cdi_type = 0;
    uint16_t cdi_immunity = 0;
  };
  std::vector<InitialStrain> initial_strains;
  std::vector<std::string> enabled_fixes;
  uint64_t seed = 42;
  GpuConfig gpu;
  ChemicalEnvironmentConfig chem_env;
  CellBiologyConfig cell_bio;
  bool profile_steps = false;

  // Spec 5 §4 — Dysbiosis safety net. When > 0, the run halts if global agent
  // density (cells/mL) exceeds this threshold, keeping the model in the
  // homeostatic regime it is calibrated for. 0 disables the check.
  Real dysbiosis_threshold = 0.0;   // cells/mL, 0 = disabled
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
                             const std::string& key,
                             const std::string& val);

 private:
  static std::string trim(std::string_view s);
  static Real parse_real(const std::string& key, const std::string& val);
  static Int parse_int(const std::string& key, const std::string& val);
};

}  // namespace gutibm

#endif  // GUTIBM_INPUT_PARSER_H
