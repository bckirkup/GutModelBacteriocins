/* -----------------------------------------------------------------------
   GutIBM – Agent and genome tests
   ----------------------------------------------------------------------- */

#include "agent.h"
#include "plasmid.h"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace gutibm;

void test_agent_creation() {
  Agent a = Agent::create_default(1, 1, {100e-6, 100e-6, 50e-6}, 5e-4);

  assert(a.tag == 1);
  assert(a.type == 1);
  assert(std::abs(a.radius - CELL_RADIUS_DEFAULT) < 1e-12);
  assert(a.mu_max > 0.0);
  assert(a.state == PhenoState::NORMAL);
  assert(a.receptor_expr[0] == 1.0);  // full expression
  assert(a.genome.lineage_id == 1);
  assert(a.genome.bi_loci.empty());

  std::cout << "  test_agent_creation: PASSED\n";
}

void test_agent_pool() {
  AgentPool pool;
  pool.reserve(100);

  for (int i = 0; i < 10; ++i) {
    Agent a = Agent::create_default(pool.next_tag(), 1,
                                     {i * 10e-6, 0, 50e-6}, 5e-4);
    pool.push_back(std::move(a));
  }
  assert(pool.size() == 10);

  pool.remove(5);
  assert(pool.size() == 9);

  std::cout << "  test_agent_pool: PASSED\n";
}

void test_plasmid_library() {
  auto col_e1 = PlasmidLibrary::colicin_E1();
  assert(col_e1.target == ReceptorType::BtuB);
  assert(col_e1.pI > 8.5);
  assert(col_e1.bclass == BacteriocinClass::LETHAL_CORE);

  auto col_b = PlasmidLibrary::colicin_B();
  assert(col_b.target == ReceptorType::FepA);
  assert(col_b.pI < 6.0);
  assert(col_b.bclass == BacteriocinClass::LETHAL_HALO);

  // ColE2 correction: secreted as complex, net acidic
  auto col_e2 = PlasmidLibrary::colicin_E2();
  assert(col_e2.pI < 7.0);  // corrected per VADI review

  const auto& lib = PlasmidLibrary::entries();
  assert(lib.size() >= 5);

  std::cout << "  test_plasmid_library: PASSED\n";
}

void test_genome_operations() {
  Agent a = Agent::create_default(1, 1, {0, 0, 0}, 5e-4);

  // Add BI clusters
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_B());
  assert(a.genome.bi_loci.size() == 2);

  // Receptor downregulation
  a.receptor_expr[static_cast<int>(ReceptorType::BtuB)] = 0.1;
  assert(a.receptor_expr[static_cast<int>(ReceptorType::BtuB)] < 0.5);

  std::cout << "  test_genome_operations: PASSED\n";
}

void test_sphere_volume() {
  Real r = 0.5e-6;
  Real vol = sphere_volume(r);
  Real expected = (4.0/3.0) * PI * r * r * r;
  assert(std::abs(vol - expected) < 1e-30);

  std::cout << "  test_sphere_volume: PASSED\n";
}

int main() {
  std::cout << "=== Agent Tests ===\n";
  test_agent_creation();
  test_agent_pool();
  test_plasmid_library();
  test_genome_operations();
  test_sphere_volume();
  std::cout << "All agent tests passed.\n";
  return 0;
}
