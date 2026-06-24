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
#include "qssa_solver.h"
#include "fix_metabolism.h"
#include "fix_receptor.h"
#include "fix_bacteriocin.h"
#include "fix_conjugation.h"
#include "fix_mutation.h"
#include "fix_mechanics.h"
#include "hdf5_writer.h"
#include "hdf5_reader.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>

namespace gutibm {

struct SimulationConfig {
  // Time control
  Real total_time        = 86400.0;    // 24 h default
  Real bio_dt            = 60.0;       // 1 min biological timestep
  Real output_interval   = 3600.0;     // hourly output

  // Adaptive timestep (CFL-like)
  bool adaptive_dt_enabled = false;
  Real dt_min            = 1.0;        // minimum timestep (s)
  Real dt_max            = 300.0;      // maximum timestep (s)
  Real dt_safety         = 0.8;        // CFL safety factor
  Real dt_growth_limit   = 0.1;        // max mu*dt allowed

  // Domain
  DomainConfig domain;

  // Fields
  AdvectionConfig advection;
  VBFConfig vbf;
  std::vector<ChemicalSpec> chemicals;

  // QSSA
  QSSAConfig qssa;

  // Fixes
  MetabolismConfig metabolism;
  ReceptorConfig receptor;
  BacteriocinConfig bacteriocin;
  ConjugationConfig conjugation;
  MutationConfig mutation;
  MechanicsConfig mechanics;

  // Output
  HDF5Config hdf5;

  // Checkpoint restart (non-empty file → resume instead of initial_strains)
  struct CheckpointConfig {
    std::string file;  // HDF5 snapshot to load
    std::string step;  // step group name; empty = latest
  };
  CheckpointConfig checkpoint;

  // Initial population
  struct InitialStrain {
    Int type;
    Int count;
    Real mu_max;
    std::vector<std::string> plasmids;
    bool conjugative;
  };
  std::vector<InitialStrain> initial_strains;

  // Fix plugins to instantiate (empty = all registered defaults in order)
  std::vector<std::string> enabled_fixes;

  // Random seed
  uint64_t seed = 42;
};

class InputParser {
 public:
  // Parse a simple JSON-like config file
  static SimulationConfig parse(const std::string& filename);

  // Create default config with standard gut parameters
  static SimulationConfig default_config();

 private:
  static std::string trim(const std::string& s);
  static Real parse_real(const std::string& val);
  static Int parse_int(const std::string& val);
};

}  // namespace gutibm

#endif  // GUTIBM_INPUT_PARSER_H
