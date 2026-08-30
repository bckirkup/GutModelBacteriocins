#include "chemical_field.h"
#include "diffusion_gpu.h"
#include "domain.h"

#include <cassert>
#include <iostream>

using namespace gutibm;

namespace {

Domain make_domain(Int nz) {
  DomainConfig cfg;
  cfg.lo = {0.0, 0.0, 0.0};
  cfg.hi = {5.0e-6, 5.0e-6, nz * 5.0e-6};
  cfg.grid_dx = 5.0e-6;
  cfg.hash_cell_size = 10.0e-6;
  Domain domain;
  domain.init(cfg);
  return domain;
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
