#include "chemical_field.h"
#include "domain.h"

#include <cassert>
#include <iostream>

using namespace gutibm;

namespace {

ChemicalField make_field() {
  DomainConfig config;
  config.hi = {20.0e-6, 20.0e-6, 20.0e-6};
  config.grid_dx = 5.0e-6;
  Domain domain;
  domain.init(config);

  ChemicalSpec carbon;
  carbon.name = "carbon";
  ChemicalSpec oxygen;
  oxygen.name = "oxygen";

  ChemicalField field;
  field.init(domain, {carbon, oxygen});
  return field;
}

void test_reaction_dirty_bookkeeping() {
  ChemicalField field = make_field();
  field.zero_reactions();
  for (Int spec = 0; spec < field.num_species(); ++spec) {
    assert(!field.host_reac_dirty(spec));
  }

  field.mark_host_reac_dirty(0);
  assert(field.host_reac_dirty(0));
  assert(!field.host_reac_dirty(1));
  field.clear_host_reac_dirty(0);
  assert(!field.host_reac_dirty(0));

  field.mark_all_host_reac_dirty();
  assert(field.host_reac_dirty(0));
  assert(field.host_reac_dirty(1));
  field.clear_host_reac_dirty(0);
  field.clear_host_reac_dirty(1);

  (void)field.mutable_species_reaction(1);
  assert(field.host_reac_dirty(1));
  field.clear_host_reac_dirty(1);

  field.reac(0, 0) = 1.0;
  assert(!field.host_reac_dirty(0));
}

}  // namespace

int main() {
  test_reaction_dirty_bookkeeping();
  std::cout << "ChemicalField reaction residency bookkeeping tests passed.\n";
  return 0;
}
