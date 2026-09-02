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
#include "gpu_profile.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;
using enum ReceptorType;

namespace {

void test_soa_indexing_isolates_receptors_per_agent() {
  constexpr Int n = 3;
  std::vector receptor_expr(static_cast<size_t>(NUM_RECEPTORS) * n, 0.0);

  // Distinct per-agent iron vs BtuB expression (Fur-like heterogeneity).
  for (Int i = 0; i < n; ++i) {
    receptor_expr[static_cast<size_t>(to_underlying(BtuB)) * n + i] = 1.0;
    receptor_expr[static_cast<size_t>(to_underlying(FepA)) * n + i] =
        2.0 + static_cast<double>(i);
    receptor_expr[static_cast<size_t>(to_underlying(IroN)) * n + i] =
        10.0 + static_cast<double>(i);
  }

  for (Int i = 0; i < n; ++i) {
    const double fepA = receptor_expr[to_underlying(FepA) * n + i];
    const double iroN = receptor_expr[to_underlying(IroN) * n + i];
    assert(std::abs(fepA - (2.0 + static_cast<double>(i))) < 1e-15);
    assert(std::abs(iroN - (10.0 + static_cast<double>(i))) < 1e-15);

    // AoS (wrong) would read neighboring receptor slots for agent 0 as
    // FepA/IroN when receptors differ — lock that this layout is SoA.
    const double aos_fepA =
        receptor_expr[static_cast<size_t>(i) * NUM_RECEPTORS + to_underlying(FepA)];
    if (i == 0) {
      assert(std::abs(aos_fepA - fepA) > 0.5);
    }
  }
}

void test_transfer_profile_counters() {
  gpu_transfer_profile_set_enabled(true);
  gpu_transfer_record_h2d(0.25, 128);
  gpu_transfer_record_h2d(0.75, 256);
  gpu_transfer_record_d2h(0.5, 512);

  const GpuTransferProfile profile = gpu_transfer_profile_snapshot();
  assert(std::abs(profile.h2d_s - 1.0) < 1.0e-15);
  assert(std::abs(profile.d2h_s - 0.5) < 1.0e-15);
  assert(profile.h2d_bytes == 384);
  assert(profile.d2h_bytes == 512);
  assert(profile.h2d_calls == 2);
  assert(profile.d2h_calls == 1);

  gpu_transfer_profile_reset();
  const GpuTransferProfile reset = gpu_transfer_profile_snapshot();
  assert(reset.h2d_s == 0.0);
  assert(reset.d2h_s == 0.0);
  assert(reset.h2d_bytes == 0);
  assert(reset.d2h_bytes == 0);
  assert(reset.h2d_calls == 0);
  assert(reset.d2h_calls == 0);

  gpu_transfer_record_h2d(0.125, 64);
  gpu_transfer_record_d2h(0.25, 96);
  const GpuTransferProfile round_trip = gpu_transfer_profile_snapshot();
  assert(std::abs(round_trip.h2d_s - 0.125) < 1.0e-15);
  assert(std::abs(round_trip.d2h_s - 0.25) < 1.0e-15);
  assert(round_trip.h2d_bytes == 64);
  assert(round_trip.d2h_bytes == 96);
  assert(round_trip.h2d_calls == 1);
  assert(round_trip.d2h_calls == 1);
  gpu_transfer_profile_set_enabled(false);
}

}  // namespace

int main() {
  test_soa_indexing_isolates_receptors_per_agent();
  test_transfer_profile_counters();
  std::cout << "All GPU receptor layout tests passed.\n";
  return 0;
}
