/* -----------------------------------------------------------------------
   GutIBM – MPI agent transfer round-trip tests
   ----------------------------------------------------------------------- */

#include "agent_transfer.h"
#include "plasmid.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

void test_round_trip_preserves_state() {
  Agent a = Agent::create_default(42, 3, {10e-6, 20e-6, 30e-6}, 4e-4);
  a.in_crypt = true;
  a.mu_realized = 1.2e-5;
  a.genome.toxin_affinity[static_cast<int>(ReceptorType::BtuB)] = 0.05;
  a.genome.ligand_affinity[static_cast<int>(ReceptorType::FepA)] = 0.75;
  a.genome.has_conjugative_plasmid = true;

  BICluster bi = PlasmidLibrary::colicin_E1();
  bi.immunity_binding_affinity = 0.12;
  a.genome.bi_loci.push_back(bi);

  BICluster bi2 = PlasmidLibrary::colicin_B();
  bi2.immunity_binding_affinity = 1.0;
  a.genome.bi_loci.push_back(bi2);

  std::vector<Agent> in = {a};
  std::vector<char> buf;
  agent_transfer_serialize(in, buf);
  std::vector<Agent> out = agent_transfer_deserialize(buf);

  assert(out.size() == 1);
  const Agent& b = out[0];

  assert(b.tag == a.tag);
  assert(b.in_crypt == true);
  assert(std::abs(b.mu_realized - a.mu_realized) < 1e-20);
  assert(b.genome.has_conjugative_plasmid == true);
  assert(std::abs(b.genome.toxin_affinity[static_cast<int>(ReceptorType::BtuB)] - 0.05) < 1e-12);
  assert(std::abs(b.genome.ligand_affinity[static_cast<int>(ReceptorType::FepA)] - 0.75) < 1e-12);
  assert(b.genome.bi_loci.size() == 2);
  assert(std::abs(b.genome.bi_loci[0].immunity_binding_affinity - 0.12) < 1e-12);
  assert(std::abs(b.genome.bi_loci[1].immunity_binding_affinity - 1.0) < 1e-12);

  std::cout << "  test_round_trip_preserves_state: PASSED\n";
}

void test_round_trip_multiple_agents() {
  std::vector<Agent> in;
  for (int i = 0; i < 5; ++i) {
    Agent a = Agent::create_default(i + 1, 1, {i * 5e-6, 0, 10e-6}, 5e-4);
    a.in_crypt = (i % 2 == 0);
    in.push_back(a);
  }

  std::vector<char> buf;
  agent_transfer_serialize(in, buf);
  auto out = agent_transfer_deserialize(buf);
  assert(out.size() == in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    assert(out[i].tag == in[i].tag);
    assert(out[i].in_crypt == in[i].in_crypt);
  }

  std::cout << "  test_round_trip_multiple_agents: PASSED\n";
}

int main() {
  std::cout << "=== Agent Transfer Tests ===\n";
  test_round_trip_preserves_state();
  test_round_trip_multiple_agents();
  std::cout << "All agent transfer tests passed.\n";
  return 0;
}
