#include "chemical_field.h"
#include "domain.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

namespace {

Domain make_domain() {
  DomainConfig config;
  config.lo = {0.0, 0.0, 0.0};
  config.hi = {20.0e-6, 20.0e-6, 20.0e-6};
  config.grid_dx = 5.0e-6;
  config.hash_cell_size = 10.0e-6;
  Domain domain;
  domain.init(config);
  return domain;
}

void test_delivery_axis_cap_defers_excess_mass() {
  constexpr Real initial = 1.0e-3;
  constexpr Real dt = 60.0;
  Domain domain = make_domain();
  ChemicalSpec spec;
  spec.name = "oxygen";
  spec.diff_coeff = 1.0e-20;
  spec.initial_conc = initial;
  spec.diffusion_enabled = true;
  spec.delivery_enabled = true;
  spec.epithelial_boundary_mode = EpithelialBoundaryMode::Flux;

  ChemicalField chemical;
  chemical.init(domain, {spec});
  const Int species_index = chemical.find("oxygen");
  const Real cell_volume = domain.cell_volume();
  const Int target = domain.cell_index(1, 1, 1);
  chemical.add_prescribed_sink_global(
      species_index, target, 6.0 * initial * cell_volume);

  chemical.apply_diffusion(domain, dt);
  chemical.flux_accounting().commit_boundary_and_reaction_step();

  const Real deferred =
      chemical.flux_accounting().delivery_axis_deferred_mass_for_step(
          species_index);
  assert(std::isfinite(deferred));
  assert(deferred > cell_volume * initial);

  Real minimum = chemical.conc_global(species_index, 0);
  for (Int cell = 1; cell < chemical.global_ncells(); ++cell) {
    minimum = std::min(minimum, chemical.conc_global(species_index, cell));
  }
  assert(minimum >= 0.0);
  std::cout << "  test_delivery_axis_cap_defers_excess_mass: PASSED\n";
}

}  // namespace

int main() {
  test_delivery_axis_cap_defers_excess_mass();
  return 0;
}
