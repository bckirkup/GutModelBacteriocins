#include "chemical_field.h"
#include "diffusion_gpu.h"
#include "domain.h"

#include <cassert>
#include <iostream>

using namespace gutibm;

namespace {

Domain make_domain(Int nx, Int ny, Int nz, Int grid_halo_width = 1) {
  DomainConfig cfg;
  cfg.lo = {0.0, 0.0, 0.0};
  cfg.hi = {nx * 5.0e-6, ny * 5.0e-6, nz * 5.0e-6};
  cfg.grid_dx = 5.0e-6;
  cfg.hash_cell_size = 10.0e-6;
  cfg.grid_halo_width = grid_halo_width;
  Domain domain;
  domain.init(cfg);
  return domain;
}

Domain make_domain(Int nz) {
  return make_domain(1, 1, nz);
}

ChemicalSpec diffusing_species(const char* name,
                               EpithelialBoundaryMode mode) {
  ChemicalSpec spec;
  spec.name = name;
  spec.diff_coeff = 1.0e-9;
  spec.retardation = 1.0;
  spec.diffusion_enabled = true;
  spec.epithelial_boundary_mode = mode;
  return spec;
}

ChemicalSpec delivery_species(EpithelialBoundaryMode mode) {
  ChemicalSpec spec = diffusing_species("delivery", mode);
  spec.delivery_enabled = true;
  spec.epithelial_transfer_coeff =
      mode == EpithelialBoundaryMode::Robin ? 1.0e-5 : 0.0;
  return spec;
}

}  // namespace

int main() {
  const Domain tall_domain = make_domain(1025);
  assert(tall_domain.nz() == 1025);
  assert(diffusion_line_lengths_within(
      tall_domain, EpithelialBoundaryMode::Dirichlet, 1024));
  assert(!diffusion_line_lengths_within(
      tall_domain, EpithelialBoundaryMode::Robin, 1024));

  ChemicalField mixed;
  mixed.init(tall_domain, {
      diffusing_species("dirichlet", EpithelialBoundaryMode::Dirichlet),
      diffusing_species("robin", EpithelialBoundaryMode::Robin),
  });
  assert(!diffusion_all_species_within(tall_domain, mixed, 1024));

  const Domain supported_domain = make_domain(1024);
  assert(supported_domain.nz() == 1024);
  ChemicalField supported;
  supported.init(supported_domain, {
      diffusing_species("dirichlet", EpithelialBoundaryMode::Dirichlet),
      diffusing_species("robin", EpithelialBoundaryMode::Robin),
  });
  assert(diffusion_all_species_within(supported_domain, supported, 1024));

  const Domain route_domain = make_domain(4, 5, 6);
  ChemicalField route_field;
  route_field.init(
      route_domain, {delivery_species(EpithelialBoundaryMode::Robin)});
  assert(delivery_route_b_eligible(route_domain, route_field));
  assert(!delivery_route_b_eligible(route_domain, route_field, 4));

  const Domain exact_cap = make_domain(1, 1, 512);
  ChemicalField exact_cap_field;
  exact_cap_field.init(
      exact_cap, {delivery_species(EpithelialBoundaryMode::Robin)});
  assert(delivery_route_b_eligible(exact_cap, exact_cap_field));

  const Domain over_cap = make_domain(1, 1, 513);
  ChemicalField over_cap_field;
  over_cap_field.init(
      over_cap, {delivery_species(EpithelialBoundaryMode::Robin)});
  assert(!delivery_route_b_eligible(over_cap, over_cap_field));

  ChemicalField slab_field;
  const Domain slab_domain = make_domain(4, 5, 6, 2);
  slab_field.init(
      slab_domain, {delivery_species(EpithelialBoundaryMode::Robin)},
      "slab");
  assert(!delivery_route_b_eligible(slab_domain, slab_field));

  ChemicalSpec disabled_robin =
      diffusing_species("disabled_robin", EpithelialBoundaryMode::Robin);
  disabled_robin.diffusion_enabled = false;
  ChemicalField mixed_disabled;
  mixed_disabled.init(tall_domain, {
      diffusing_species("dirichlet", EpithelialBoundaryMode::Dirichlet),
      disabled_robin,
  });
  assert(diffusion_all_species_within(tall_domain, mixed_disabled, 1024));

  std::cout << "GPU diffusion species-mask predicates passed.\n";
  return 0;
}
