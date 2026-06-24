/* -----------------------------------------------------------------------
   GutIBM – Fix plugin registry (NUFEB/LAMMPS-style extensibility)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_REGISTRY_H
#define GUTIBM_FIX_REGISTRY_H

#include "fix.h"
#include "input_parser.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace gutibm {

class Simulation;

class FixRegistry {
 public:
  using Factory = std::function<std::unique_ptr<Fix>(Simulation&, const SimulationConfig&)>;

  // Register a Fix factory under a unique name (e.g. "metabolism").
  static void register_fix(const std::string& name, Factory factory);

  // Built-in Fix modules in default execution order.
  static void register_defaults();

  // Instantiate registered fixes (order = registration order, or cfg.enabled_fixes).
  static std::vector<std::unique_ptr<Fix>> create_all(Simulation& sim,
                                                      const SimulationConfig& cfg);

  // Default fix names in registration order.
  static std::vector<std::string> default_fix_names();

  // Names of currently registered fixes (for diagnostics/tests).
  static std::vector<std::string> registered_names();

 private:
  static std::vector<std::pair<std::string, Factory>>& entries();
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_REGISTRY_H
