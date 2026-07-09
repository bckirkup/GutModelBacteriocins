#include "chemical_field.h"
#include "domain.h"
#include "input_parser.h"
#include "species_names.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string_view>
#include <vector>

using namespace gutibm;

namespace {

Domain make_domain(Int nx, Int ny, Int nz) {
  DomainConfig cfg;
  cfg.lo = {0.0, 0.0, 0.0};
  cfg.hi = {nx * 5.0e-6, ny * 5.0e-6, nz * 5.0e-6};
  cfg.grid_dx = 5.0e-6;
  cfg.hash_cell_size = 10.0e-6;
  Domain domain;
  domain.init(cfg);
  return domain;
}

ChemicalSpec diffusing_species(Real diffusion, Real initial, Real boundary) {
  ChemicalSpec spec;
  spec.name = species::OXYGEN;
  spec.diff_coeff = diffusion;
  spec.retardation = 1.0;
  spec.initial_conc = initial;
  spec.boundary_conc = boundary;
  spec.diffusion_enabled = true;
  return spec;
}

void test_uniform_field_is_fixed_point() {
  Domain domain = make_domain(5, 4, 6);
  ChemicalField chem;
  chem.init(domain, {diffusing_species(2.1e-9, 0.25, 0.25)});

  chem.apply_diffusion(domain, 300.0);

  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    assert(std::abs(chem.conc(0, cell) - 0.25) < 1.0e-10);
  }
  std::cout << "  test_uniform_field_is_fixed_point: PASSED\n";
}

void test_point_source_golden_profile() {
  Domain domain = make_domain(5, 1, 2);
  ChemicalField chem;
  chem.init(domain, {diffusing_species(1.0e-12, 0.0, 0.0)});
  chem.conc(0, domain.cell_index(2, 0, 1)) = 1.0;

  chem.apply_diffusion(domain, 2.5);

  constexpr std::array<Real, 5> expected = {
      0.00586510263929619,
      0.0645161290322581,
      0.768328445747801,
      0.0645161290322581,
      0.00586510263929619,
  };
  for (Int ix = 0; ix < domain.nx(); ++ix) {
    const Real actual = chem.conc(0, domain.cell_index(ix, 0, 1));
    assert(std::abs(actual - expected[static_cast<size_t>(ix)]) < 1.0e-12);
  }
  std::cout << "  test_point_source_golden_profile: PASSED\n";
}

ChemicalField diffuse_point_source(const Domain& domain, Real diffusion,
                                   bool enabled) {
  ChemicalSpec spec = diffusing_species(diffusion, 0.0, 0.0);
  spec.diffusion_enabled = enabled;
  ChemicalField chem;
  chem.init(domain, {spec});
  chem.conc(0, domain.cell_index(2, 0, 1)) = 1.0;
  chem.apply_diffusion(domain, 2.5);
  return chem;
}

void test_diffusion_enable_and_coefficient_sensitivity() {
  Domain domain = make_domain(5, 1, 2);
  const ChemicalField disabled = diffuse_point_source(domain, 1.0e-12, false);
  const ChemicalField slow = diffuse_point_source(domain, 5.0e-13, true);
  const ChemicalField fast = diffuse_point_source(domain, 2.0e-12, true);
  const Int center = domain.cell_index(2, 0, 1);
  const Int neighbor = domain.cell_index(1, 0, 1);

  assert(disabled.conc(0, center) == 1.0);
  assert(disabled.conc(0, neighbor) == 0.0);
  assert(fast.conc(0, center) < slow.conc(0, center));
  assert(fast.conc(0, neighbor) > slow.conc(0, neighbor));
  std::cout << "  test_diffusion_enable_and_coefficient_sensitivity: PASSED\n";
}

void test_dirichlet_neumann_boundary_gradient() {
  Domain domain = make_domain(4, 4, 8);
  ChemicalField chem;
  chem.init(domain, {diffusing_species(2.1e-9, 0.0, 1.0)});

  chem.apply_diffusion(domain, 60.0);

  Real previous = chem.conc(0, domain.cell_index(0, 0, 0));
  assert(std::abs(previous - 1.0) < 1.0e-15);
  for (Int iz = 1; iz < domain.nz(); ++iz) {
    const Real current = chem.conc(0, domain.cell_index(0, 0, iz));
    assert(std::isfinite(current));
    assert(current > 0.0);
    assert(current <= previous);
    previous = current;
  }
  assert(previous < 1.0);
  std::cout << "  test_dirichlet_neumann_boundary_gradient: PASSED\n";
}

void test_configured_z_gradient_is_background_fixed_point() {
  Domain domain = make_domain(4, 3, 8);
  ChemicalSpec spec = diffusing_species(5.0e-10, 5.0e-3, 5.0e-3);
  spec.z_gradient_enabled = true;
  spec.z_gradient_lambda = 25.0e-6;
  ChemicalField chem;
  chem.init(domain, {spec});

  chem.apply_diffusion(domain, 60.0);

  assert(std::abs(chem.conc(0, domain.cell_index(0, 0, 0))
                  - spec.boundary_conc) < 1.0e-15);
  for (Int iz = 1; iz < domain.nz() - 1; ++iz) {
    const Real z_rel = (iz + 0.5) * domain.dx();
    const Real expected = spec.initial_conc
        * std::exp(-z_rel / spec.z_gradient_lambda);
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.nx(); ++ix) {
        const Real actual = chem.conc(0, domain.cell_index(ix, iy, iz));
        assert(std::abs(actual - expected) < 1.0e-12);
      }
    }
  }
  const Real top = chem.conc(0, domain.cell_index(0, 0, domain.nz() - 1));
  const Real below = chem.conc(0, domain.cell_index(0, 0, domain.nz() - 2));
  assert(std::abs(top - below) < 1.0e-15);

  chem.conc(0, domain.cell_index(1, 1, 1)) = 0.0;
  chem.apply_diffusion(domain, 300.0);
  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    assert(chem.conc(0, cell) >= 0.0);
  }
  std::cout << "  test_configured_z_gradient_is_background_fixed_point: PASSED\n";
}

void test_large_timestep_is_positive_and_finite() {
  Domain domain = make_domain(6, 5, 7);
  ChemicalField chem;
  chem.init(domain, {diffusing_species(2.1e-9, 0.0, 0.0)});
  for (Int iz = 1; iz < domain.nz(); ++iz) {
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.nx(); ++ix) {
        chem.conc(0, domain.cell_index(ix, iy, iz)) =
            ((ix + iy + iz) % 2 == 0) ? 1.0 : 0.0;
      }
    }
  }

  chem.apply_diffusion(domain, 300.0);

  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    const Real value = chem.conc(0, cell);
    assert(std::isfinite(value));
    assert(value >= 0.0);
    assert(value <= 1.0 + 1.0e-12);
  }
  std::cout << "  test_large_timestep_is_positive_and_finite: PASSED\n";
}

void test_default_species_configuration() {
  const SimulationConfig cfg = InputParser::default_config();
  const auto diffusion_enabled = [&cfg](std::string_view name) {
    const auto it = std::ranges::find_if(
        cfg.chemicals, [name](const ChemicalSpec& spec) { return spec.name == name; });
    assert(it != cfg.chemicals.end());
    return it->diffusion_enabled;
  };

  assert(diffusion_enabled(species::CARBON));
  assert(diffusion_enabled(species::IRON));
  assert(diffusion_enabled(species::B12));
  assert(diffusion_enabled(species::ACETATE));
  assert(diffusion_enabled(species::ETHANOLAMINE));
  assert(!diffusion_enabled(species::BACTERIOCIN_BTUB));
  assert(!diffusion_enabled(species::BACTERIOCIN_FEPA));
  assert(!diffusion_enabled(species::BACTERIOCIN_CIRA));
  assert(!diffusion_enabled(species::BACTERIOCIN_FHUA));
  std::cout << "  test_default_species_configuration: PASSED\n";
}

}  // namespace

int main() {
  std::cout << "=== Nutrient Diffusion Tests ===\n";
  test_uniform_field_is_fixed_point();
  test_point_source_golden_profile();
  test_diffusion_enable_and_coefficient_sensitivity();
  test_dirichlet_neumann_boundary_gradient();
  test_configured_z_gradient_is_background_fixed_point();
  test_large_timestep_is_positive_and_finite();
  test_default_species_configuration();
  std::cout << "All nutrient diffusion tests passed.\n";
  return 0;
}
