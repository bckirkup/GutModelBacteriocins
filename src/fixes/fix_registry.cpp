/* -----------------------------------------------------------------------
   GutIBM – Fix plugin registry implementation
   ----------------------------------------------------------------------- */

#include "fix_registry.h"
#include "fix_metabolism.h"
#include "fix_bacteriocin.h"
#include "fix_receptor.h"
#include "fix_conjugation.h"
#include "fix_mutation.h"
#include "fix_mechanics.h"
#include "simulation.h"

namespace gutibm {

std::vector<std::pair<std::string, FixRegistry::Factory>>& FixRegistry::entries() {
  static std::vector<std::pair<std::string, Factory>> reg;
  return reg;
}

void FixRegistry::register_fix(const std::string& name, Factory factory) {
  auto& reg = entries();
  for (const auto& e : reg) {
    if (e.first == name) return;  // first registration wins
  }
  reg.emplace_back(name, std::move(factory));
}

void FixRegistry::register_defaults() {
  auto& reg = entries();
  if (!reg.empty()) return;

  register_fix("metabolism",
    [](Simulation& sim, const SimulationConfig& cfg) {
      return std::make_unique<FixMetabolism>(sim, cfg.metabolism);
    });
  register_fix("bacteriocin",
    [](Simulation& sim, const SimulationConfig& cfg) {
      return std::make_unique<FixBacteriocin>(sim, cfg.bacteriocin);
    });
  register_fix("receptor",
    [](Simulation& sim, const SimulationConfig& cfg) {
      return std::make_unique<FixReceptor>(sim, cfg.receptor);
    });
  register_fix("conjugation",
    [](Simulation& sim, const SimulationConfig& cfg) {
      return std::make_unique<FixConjugation>(sim, cfg.conjugation);
    });
  register_fix("mutation",
    [](Simulation& sim, const SimulationConfig& cfg) {
      return std::make_unique<FixMutation>(sim, cfg.mutation);
    });
  register_fix("mechanics",
    [](Simulation& sim, const SimulationConfig& cfg) {
      return std::make_unique<FixMechanics>(sim, cfg.mechanics);
    });
}

std::vector<std::unique_ptr<Fix>> FixRegistry::create_all(Simulation& sim,
                                                          const SimulationConfig& cfg) {
  register_defaults();
  std::vector<std::unique_ptr<Fix>> fixes;
  fixes.reserve(entries().size());
  for (auto& [name, factory] : entries()) {
    (void)name;
    fixes.push_back(factory(sim, cfg));
  }
  return fixes;
}

std::vector<std::string> FixRegistry::registered_names() {
  register_defaults();
  std::vector<std::string> names;
  names.reserve(entries().size());
  for (const auto& e : entries()) {
    names.push_back(e.first);
  }
  return names;
}

}  // namespace gutibm
