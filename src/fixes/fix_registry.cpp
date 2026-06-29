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

#include <iostream>
#include <unordered_map>

namespace gutibm {

std::vector<std::pair<std::string, FixRegistry::Factory>>& FixRegistry::entries() {
  static std::vector<std::pair<std::string, Factory>> reg;
  return reg;
}

void FixRegistry::register_fix(const std::string& name, Factory factory) {
  auto& reg = entries();
  for (const auto& entry : reg) {
    if (entry.first == name) return;  // first registration wins
  }
  reg.emplace_back(name, std::move(factory));
}

void FixRegistry::register_defaults() {
  if (const auto& reg = entries(); !reg.empty()) return;

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

  std::vector<std::string> order = cfg.enabled_fixes;
  if (order.empty()) {
    order = default_fix_names();
  }

  std::unordered_map<std::string, Factory> factories;
  factories.reserve(entries().size());
  for (auto& [name, factory] : entries()) {
    factories.emplace(name, factory);
  }

  std::vector<std::unique_ptr<Fix>> fixes;
  fixes.reserve(order.size());
  for (const auto& name : order) {
    auto it = factories.find(name);
    if (it == factories.end()) {
      std::cerr << "Warning: unknown fix '" << name << "' — skipping\n";
      continue;
    }
    fixes.push_back(it->second(sim, cfg));
  }
  return fixes;
}

std::vector<std::string> FixRegistry::default_fix_names() {
  register_defaults();
  std::vector<std::string> names;
  names.reserve(entries().size());
  for (const auto& [name, factory] : entries()) {
    names.push_back(name);
  }
  return names;
}

std::vector<std::string> FixRegistry::registered_names() {
  register_defaults();
  std::vector<std::string> names;
  names.reserve(entries().size());
  for (const auto& [name, factory] : entries()) {
    names.push_back(name);
  }
  return names;
}

}  // namespace gutibm
