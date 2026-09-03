#include "chemical_field.h"
#include "domain.h"

#include <cassert>
#include <iostream>

using namespace gutibm;

namespace {

Domain make_domain() {
  DomainConfig config;
  config.lo = {0.0, 0.0, 0.0};
  config.hi = {10.0e-6, 10.0e-6, 10.0e-6};
  config.grid_dx = 5.0e-6;
  config.grid_halo_width = 2;
  config.hash_cell_size = 10.0e-6;
  Domain domain;
  domain.init(config);
  return domain;
}

ChemicalSpec diffusing_spec() {
  ChemicalSpec spec;
  spec.name = "diffusing";
  spec.initial_conc = 1.0;
  spec.boundary_conc = 1.0;
  spec.diffusion_enabled = true;
  spec.diff_coeff = 1.0e-10;
  return spec;
}

ChemicalSpec nondiffusing_spec() {
  ChemicalSpec spec;
  spec.name = "nondiffusing";
  spec.initial_conc = 2.0;
  spec.diffusion_enabled = false;
  return spec;
}

void expect_all_dirty(const ChemicalField& field, bool expected) {
  for (Int spec = 0; spec < field.num_species(); ++spec) {
    assert(field.host_conc_dirty(spec) == expected);
  }
}

void test_dirty_state_and_marked_writers() {
  const Domain domain = make_domain();
  ChemicalField field;
  field.init(domain, {diffusing_spec(), nondiffusing_spec()});
  expect_all_dirty(field, true);

  field.clear_host_conc_dirty(0);
  field.clear_host_conc_dirty(1);
  expect_all_dirty(field, false);

  field.mark_host_conc_dirty(1);
  assert(field.host_conc_dirty(1));
  assert(!field.host_conc_dirty(0));
  field.clear_host_conc_dirty(1);
  assert(!field.host_conc_dirty(1));

  field.mark_all_host_conc_dirty();
  expect_all_dirty(field, true);
  field.clear_host_conc_dirty(0);
  field.clear_host_conc_dirty(1);

  field.apply_diffusion(domain, 1.0);
  assert(field.host_conc_dirty(0));
  assert(!field.host_conc_dirty(1));

  field.clear_host_conc_dirty(0);
  field.apply_boundaries(domain);
  assert(field.host_conc_dirty(0));
  assert(field.host_conc_dirty(1));

  field.clear_host_conc_dirty(0);
  field.clear_host_conc_dirty(1);
  auto& row = field.mutable_species_concentration(0);
  row.front() = 3.0;
  assert(!field.host_conc_dirty(0));
  field.mark_host_conc_dirty(0);
  assert(field.host_conc_dirty(0));
}

void test_slab_halo_exchange_marks_rows() {
  const Domain domain = make_domain();
  ChemicalField field;
  field.init(domain, {diffusing_spec(), nondiffusing_spec()}, "slab");
  field.clear_host_conc_dirty(0);
  field.clear_host_conc_dirty(1);
  field.exchange_concentration_halos();
  expect_all_dirty(field, true);
}

void test_reinitialization_starts_dirty() {
  const Domain domain = make_domain();
  ChemicalField field;
  field.init(domain, {diffusing_spec(), nondiffusing_spec()});
  field.clear_host_conc_dirty(0);
  field.clear_host_conc_dirty(1);
  field.init(domain, {diffusing_spec()});
  expect_all_dirty(field, true);
}

}  // namespace

int main() {
  test_dirty_state_and_marked_writers();
  test_slab_halo_exchange_marks_rows();
  test_reinitialization_starts_dirty();
  std::cout << "ChemicalField residency bookkeeping tests passed.\n";
  return 0;
}
