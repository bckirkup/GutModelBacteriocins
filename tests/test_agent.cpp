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
  // ColE1: pI 9.0 → LETHAL_CORE
  auto col_e1 = PlasmidLibrary::colicin_E1();
  assert(col_e1.target == ReceptorType::BtuB);
  assert(std::abs(col_e1.pI - 9.0) < 1e-6);
  assert(col_e1.bclass == BacteriocinClass::LETHAL_CORE);

  // ColE2: pI 6.5 → LETHAL_HALO (secreted as acidic Im2 complex)
  auto col_e2 = PlasmidLibrary::colicin_E2();
  assert(col_e2.target == ReceptorType::BtuB);
  assert(std::abs(col_e2.pI - 6.5) < 1e-6);
  assert(col_e2.bclass == BacteriocinClass::LETHAL_HALO);

  // ColB: pI 5.4 → LETHAL_HALO
  auto col_b = PlasmidLibrary::colicin_B();
  assert(col_b.target == ReceptorType::FepA);
  assert(std::abs(col_b.pI - 5.4) < 1e-6);
  assert(col_b.bclass == BacteriocinClass::LETHAL_HALO);

  // ColIa: pI 7.2 → NEUTRAL
  auto col_ia = PlasmidLibrary::colicin_Ia();
  assert(col_ia.target == ReceptorType::CirA);
  assert(std::abs(col_ia.pI - 7.2) < 1e-6);
  assert(col_ia.bclass == BacteriocinClass::NEUTRAL);

  // ColM: pI 9.3 → LETHAL_CORE
  auto col_m = PlasmidLibrary::colicin_M();
  assert(col_m.target == ReceptorType::FhuA);
  assert(std::abs(col_m.pI - 9.3) < 1e-6);
  assert(col_m.bclass == BacteriocinClass::LETHAL_CORE);

  // MccV: pI 5.0 → LETHAL_HALO
  auto mcc_v = PlasmidLibrary::microcin_V();
  assert(mcc_v.target == ReceptorType::CirA);
  assert(std::abs(mcc_v.pI - 5.0) < 1e-6);
  assert(mcc_v.bclass == BacteriocinClass::LETHAL_HALO);

  const auto& lib = PlasmidLibrary::entries();
  assert(lib.size() == 6);

  std::cout << "  test_plasmid_library: PASSED\n";
}

void test_classify_by_pI() {
  // LETHAL_CORE: pI > 8.5
  assert(classify_by_pI(9.0) == BacteriocinClass::LETHAL_CORE);
  assert(classify_by_pI(8.6) == BacteriocinClass::LETHAL_CORE);
  assert(classify_by_pI(12.0) == BacteriocinClass::LETHAL_CORE);

  // LETHAL_HALO: pI < 7.0
  assert(classify_by_pI(6.9) == BacteriocinClass::LETHAL_HALO);
  assert(classify_by_pI(6.5) == BacteriocinClass::LETHAL_HALO);
  assert(classify_by_pI(5.0) == BacteriocinClass::LETHAL_HALO);
  assert(classify_by_pI(3.0) == BacteriocinClass::LETHAL_HALO);

  // NEUTRAL: 7.0 <= pI <= 8.5
  assert(classify_by_pI(7.0) == BacteriocinClass::NEUTRAL);
  assert(classify_by_pI(7.5) == BacteriocinClass::NEUTRAL);
  assert(classify_by_pI(8.0) == BacteriocinClass::NEUTRAL);
  assert(classify_by_pI(8.5) == BacteriocinClass::NEUTRAL);

  // Boundary checks
  assert(classify_by_pI(7.0) == BacteriocinClass::NEUTRAL);   // exactly 7.0 → NEUTRAL
  assert(classify_by_pI(8.5) == BacteriocinClass::NEUTRAL);   // exactly 8.5 → NEUTRAL

  std::cout << "  test_classify_by_pI: PASSED\n";
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

void test_partial_resistance_fields() {
  Agent a = Agent::create_default(1, 1, {0, 0, 0}, 5e-4);

  // Default: all affinities should be 1.0 (wild-type)
  for (int r = 0; r < NUM_RECEPTORS; ++r) {
    assert(a.genome.toxin_affinity[r] == 1.0);
    assert(a.genome.ligand_affinity[r] == 1.0);
  }

  // Simulate a partial resistance mutation on BtuB
  int btuB = static_cast<int>(ReceptorType::BtuB);
  a.genome.toxin_affinity[btuB] = 0.05;   // 20x reduced toxin binding
  a.genome.ligand_affinity[btuB] = 0.8;   // within 2x of wild-type

  // Receptor expression should remain unchanged (distinct from downregulation)
  assert(a.receptor_expr[btuB] == 1.0);
  assert(a.genome.receptor_expression[btuB] == 1.0);

  // Verify toxin affinity is heavily reduced
  assert(a.genome.toxin_affinity[btuB] < 0.11);
  assert(a.genome.toxin_affinity[btuB] > 0.0);

  // Verify ligand affinity stays near wild-type
  assert(a.genome.ligand_affinity[btuB] >= 0.5);
  assert(a.genome.ligand_affinity[btuB] <= 1.0);

  // Other receptors should be unaffected
  int fepA = static_cast<int>(ReceptorType::FepA);
  assert(a.genome.toxin_affinity[fepA] == 1.0);
  assert(a.genome.ligand_affinity[fepA] == 1.0);

  std::cout << "  test_partial_resistance_fields: PASSED\n";
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
  test_classify_by_pI();
  test_genome_operations();
  test_partial_resistance_fields();
  test_sphere_volume();
  std::cout << "All agent tests passed.\n";
  return 0;
}
