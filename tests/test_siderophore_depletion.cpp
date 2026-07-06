/* -----------------------------------------------------------------------
   GutIBM – Siderophore secretion test (Spec 4 §3)
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "fix_metabolism.h"
#include "species_names.h"

#include <cassert>
#include <iostream>

using namespace gutibm;

int main() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {30e-6, 30e-6, 15e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.hdf5.enabled = false;
  cfg.seed = 444;
  cfg.chem_env.siderophore.enabled = true;
  cfg.chem_env.siderophore.secretion_rate = 1.0e-10;
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 1;
  strain.mu_max = 5e-4;
  cfg.initial_strains = {strain};
  InputParser::finalize_config(cfg);

  bool has_species = false;
  for (const auto& spec : cfg.chemicals) {
    if (spec.name == species::SIDEROPHORE) has_species = true;
  }
  assert(has_species);

  Simulation sim;
  sim.init(cfg);
  assert(sim.config().chem_env.siderophore.enabled);

  const Int cell = sim.agents()[0].grid_cell;
  assert(cell >= 0);
  const Int i_sid = sim.chemical_field().find(species::SIDEROPHORE);
  assert(i_sid >= 0);

  FixMetabolism metab(sim, sim.config().fixes.metabolism);
  metab.compute(60.0);
  assert(sim.chemical_field().reac(i_sid, cell) > 0.0);

  std::cout << "All siderophore depletion tests passed.\n";
  return 0;
}
