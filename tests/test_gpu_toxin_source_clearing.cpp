/* -----------------------------------------------------------------------
   GutIBM – GPU toxin source clearing

   A source-materialized toxin field must be cleared on both host and device
   when the next QSSA solve has no sources for that species.
   ----------------------------------------------------------------------- */

#include "gpu_test_support.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

namespace {

SimulationConfig make_config() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {50.0e-6, 50.0e-6, 25.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.hdf5.enabled = false;
  cfg.gpu.enabled = true;
  cfg.gpu.device_id = 0;
  cfg.initial_strains.clear();
  cfg.qssa.toxin_lumping = "per_receptor";
  return cfg;
}

ToxinBurstSource make_burst(const Simulation& sim) {
  ToxinBurstSource burst;
  burst.pos = {
      0.5 * (sim.domain().lo()[0] + sim.domain().hi()[0]),
      0.5 * (sim.domain().lo()[1] + sim.domain().hi()[1]),
      0.5 * (sim.domain().lo()[2] + sim.domain().hi()[2])};
  burst.params.diff_coeff = 4.0e-11;
  burst.params.retardation = 1.0;
  burst.params.pI = 9.0;
  burst.params.source_rate = 1.0e-12;
  burst.params.decay_rate = 0.0;
  burst.release_tau = 300.0;
  burst.target = ReceptorType::BtuB;
  return burst;
}

bool all_zero(const ChemicalField& chem, Int species_idx) {
  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    if (std::abs(chem.conc(species_idx, cell)) > 1.0e-30) {
      return false;
    }
  }
  return true;
}

}  // namespace

int main() {
  std::cout << "=== GPU Toxin Source Clearing ===\n";
  if (const int gpu_status = test::require_gpu("gpu_toxin_source_clearing");
      gpu_status != 0) {
    return gpu_status;
  }

  Simulation sim;
  const SimulationConfig cfg = make_config();
  sim.init(cfg);

  const Int toxin = sim.chemical_field().find(species::BACTERIOCIN_BTUB);
  assert(toxin >= 0);

  const std::vector<ToxinBurstSource> bursts{make_burst(sim)};
  sim.qssa().solve_bacteriocin_field(
      sim.agents(), bursts, 0.0, cfg.chem_env.protease, sim.advection(),
      sim.chemical_field(), toxin, ReceptorType::BtuB, &sim.chem_gpu(), true);
  sim.chem_gpu().sync_species_concentrations_to_host(
      sim.chemical_field(), toxin);
  bool materialized = false;
  for (Int cell = 0; cell < sim.chemical_field().ncells(); ++cell) {
    if (sim.chemical_field().conc(toxin, cell) > 0.0) {
      materialized = true;
      break;
    }
  }
  assert(materialized);

  sim.qssa().solve_bacteriocin_field(
      sim.agents(), {}, 60.0, cfg.chem_env.protease, sim.advection(),
      sim.chemical_field(), toxin, ReceptorType::BtuB, &sim.chem_gpu(), true);
  assert(all_zero(sim.chemical_field(), toxin));

  sim.chem_gpu().sync_species_concentrations_to_host(
      sim.chemical_field(), toxin);
  assert(all_zero(sim.chemical_field(), toxin));

  std::cout << "  host and device toxin fields cleared: PASSED\n";
  return 0;
}
