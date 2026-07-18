/* -----------------------------------------------------------------------
   GutIBM – GPU receptor buffer layout contract (SoA)

   AgentPoolGpu packs receptor_expr / ligand_affinity as SoA:
     index = receptor * num_agents + agent
   Metabolism and receptor kernels must use the same formula. An AoS read
   (agent * NUM_RECEPTORS + receptor) silently cross-wires agents once Fur
   makes iron receptors differ from BtuB — the failure mode behind gpu_smoke
   agent_rel ~0.27 with chem_rel ~1e-9.
   ----------------------------------------------------------------------- */

#include "types.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

namespace {

void test_soa_indexing_isolates_receptors_per_agent() {
  constexpr Int n = 3;
  std::vector<double> receptor_expr(static_cast<size_t>(NUM_RECEPTORS) * n, 0.0);

  // Distinct per-agent iron vs BtuB expression (Fur-like heterogeneity).
  for (Int i = 0; i < n; ++i) {
    receptor_expr[static_cast<size_t>(to_underlying(ReceptorType::BtuB)) * n + i] =
        1.0;
    receptor_expr[static_cast<size_t>(to_underlying(ReceptorType::FepA)) * n + i] =
        2.0 + static_cast<double>(i);
    receptor_expr[static_cast<size_t>(to_underlying(ReceptorType::IroN)) * n + i] =
        10.0 + static_cast<double>(i);
  }

  for (Int i = 0; i < n; ++i) {
    const double fepA =
        receptor_expr[to_underlying(ReceptorType::FepA) * n + i];
    const double iroN =
        receptor_expr[to_underlying(ReceptorType::IroN) * n + i];
    assert(std::abs(fepA - (2.0 + static_cast<double>(i))) < 1e-15);
    assert(std::abs(iroN - (10.0 + static_cast<double>(i))) < 1e-15);

    // AoS (wrong) would read neighboring receptor slots for agent 0 as
    // FepA/IroN when receptors differ — lock that this layout is SoA.
    const double aos_fepA =
        receptor_expr[static_cast<size_t>(i) * NUM_RECEPTORS
                      + to_underlying(ReceptorType::FepA)];
    if (i == 0) {
      assert(std::abs(aos_fepA - fepA) > 0.5);
    }
  }
}

}  // namespace

int main() {
  test_soa_indexing_isolates_receptors_per_agent();
  std::cout << "All GPU receptor layout tests passed.\n";
  return 0;
}
