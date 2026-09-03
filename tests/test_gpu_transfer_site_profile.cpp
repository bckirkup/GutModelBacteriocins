/* -----------------------------------------------------------------------
   GutIBM - host-side GPU transfer site attribution contract
   ----------------------------------------------------------------------- */

#include "gpu_profile.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace gutibm;

namespace {

GpuTransferSiteProfile find_site(
    const std::vector<GpuTransferSiteProfile>& profiles,
    const std::string& label) {
  for (const auto& profile : profiles) {
    if (profile.label == label) return profile;
  }
  assert(false);
  return {};
}

void test_site_attribution_and_invariants() {
  gpu_transfer_profile_reset();
  gpu_transfer_profile_set_enabled(true);

  gpu_transfer_record_h2d(1.0, 100);
  {
    GpuTransferSite site_a("a");
    gpu_transfer_record_h2d(2.0, 200);
    {
      GpuTransferSite site_b("b");
      gpu_transfer_record_d2h(3.0, 300);
    }
    gpu_transfer_record_h2d(4.0, 150);
  }
  {
    GpuTransferSite site_a("a");
    gpu_transfer_record_d2h(5.0, 50);
  }

  const GpuTransferProfile total = gpu_transfer_profile_snapshot();
  const auto profiles = gpu_transfer_site_profiles();
  assert(profiles.size() == 3);
  assert(profiles.front().label == "a");

  const auto unattributed = find_site(profiles, "unattributed");
  const auto site_a = find_site(profiles, "a");
  const auto site_b = find_site(profiles, "b");
  assert(unattributed.h2d_bytes == 100);
  assert(unattributed.h2d_calls == 1);
  assert(site_a.h2d_bytes == 350);
  assert(site_a.h2d_calls == 2);
  assert(site_a.d2h_bytes == 50);
  assert(site_a.d2h_calls == 1);
  assert(site_b.d2h_bytes == 300);
  assert(site_b.d2h_calls == 1);
  assert(std::abs(site_a.h2d_s - 6.0) < 1.0e-15);
  assert(std::abs(site_a.d2h_s - 5.0) < 1.0e-15);

  unsigned long long h2d_bytes = 0;
  unsigned long long d2h_bytes = 0;
  unsigned long long h2d_calls = 0;
  unsigned long long d2h_calls = 0;
  for (const auto& profile : profiles) {
    h2d_bytes += profile.h2d_bytes;
    d2h_bytes += profile.d2h_bytes;
    h2d_calls += profile.h2d_calls;
    d2h_calls += profile.d2h_calls;
  }
  assert(h2d_bytes == total.h2d_bytes);
  assert(d2h_bytes == total.d2h_bytes);
  assert(h2d_calls == total.h2d_calls);
  assert(d2h_calls == total.d2h_calls);

  gpu_transfer_profile_reset();
  assert(gpu_transfer_site_profiles().empty());
  assert(gpu_transfer_profile_snapshot().h2d_bytes == 0);
  assert(gpu_transfer_profile_snapshot().d2h_bytes == 0);

  gpu_transfer_profile_set_enabled(false);
  gpu_transfer_record_h2d(6.0, 600);
  gpu_transfer_record_d2h(7.0, 700);
  assert(gpu_transfer_site_profiles().empty());
  const GpuTransferProfile disabled = gpu_transfer_profile_snapshot();
  assert(disabled.h2d_bytes == 0);
  assert(disabled.d2h_bytes == 0);
  assert(disabled.h2d_calls == 0);
  assert(disabled.d2h_calls == 0);
}

}  // namespace

int main() {
  test_site_attribution_and_invariants();
  gpu_transfer_profile_set_enabled(false);
  std::cout << "GPU transfer site profile tests passed.\n";
  return 0;
}
