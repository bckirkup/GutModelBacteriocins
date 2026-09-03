#include "gpu_profile.h"

#include <algorithm>
#include <mutex>
#include <ranges>
#include <unordered_map>

namespace gutibm {

namespace {
struct GpuTransferState {
  bool enabled = false;
  double h2d_s = 0.0;
  double d2h_s = 0.0;
  unsigned long long h2d_bytes = 0;
  unsigned long long d2h_bytes = 0;
  unsigned long long h2d_calls = 0;
  unsigned long long d2h_calls = 0;
  double slab_x_roundtrip_s = 0.0;
};

thread_local const char* current_transfer_site = "unattributed";

GpuTransferState& state() {
  static GpuTransferState s;
  return s;
}

std::unordered_map<std::string, GpuTransferSiteProfile>& site_profiles() {
  static std::unordered_map<std::string, GpuTransferSiteProfile> profiles;
  return profiles;
}

std::mutex& site_profiles_mutex() {
  static std::mutex mutex;
  return mutex;
}
}  // namespace

GpuTransferSite::GpuTransferSite(const char* label)
    : previous_(current_transfer_site) {
  current_transfer_site = label;
}

GpuTransferSite::~GpuTransferSite() {
  current_transfer_site = previous_;
}

bool gpu_transfer_profiling_enabled() { return state().enabled; }

void gpu_transfer_profile_set_enabled(bool enabled) {
  state().enabled = enabled;
  if (enabled) gpu_transfer_profile_reset();
}

void gpu_transfer_profile_reset() {
  auto& s = state();
  s.h2d_s = 0.0;
  s.d2h_s = 0.0;
  s.h2d_bytes = 0;
  s.d2h_bytes = 0;
  s.h2d_calls = 0;
  s.d2h_calls = 0;
  s.slab_x_roundtrip_s = 0.0;
  std::scoped_lock lock(site_profiles_mutex());
  site_profiles().clear();
}

void gpu_transfer_record_h2d(double seconds, unsigned long long bytes) {
  if (state().enabled) {
    auto& s = state();
    s.h2d_s += seconds;
    s.h2d_bytes += bytes;
    ++s.h2d_calls;
    std::scoped_lock lock(site_profiles_mutex());
    auto& site = site_profiles()[current_transfer_site];
    site.label = current_transfer_site;
    site.h2d_s += seconds;
    site.h2d_bytes += bytes;
    ++site.h2d_calls;
  }
}

void gpu_transfer_record_d2h(double seconds, unsigned long long bytes) {
  if (state().enabled) {
    auto& s = state();
    s.d2h_s += seconds;
    s.d2h_bytes += bytes;
    ++s.d2h_calls;
    std::scoped_lock lock(site_profiles_mutex());
    auto& site = site_profiles()[current_transfer_site];
    site.label = current_transfer_site;
    site.d2h_s += seconds;
    site.d2h_bytes += bytes;
    ++site.d2h_calls;
  }
}

GpuTransferProfile gpu_transfer_profile_snapshot() {
  const auto& s = state();
  return {s.h2d_s, s.d2h_s, s.h2d_bytes, s.d2h_bytes, s.h2d_calls,
          s.d2h_calls, s.slab_x_roundtrip_s};
}

std::vector<GpuTransferSiteProfile> gpu_transfer_site_profiles() {
  std::scoped_lock lock(site_profiles_mutex());
  std::vector<GpuTransferSiteProfile> result;
  result.reserve(site_profiles().size());
  for (const auto& [label, profile] : site_profiles()) {
    (void)label;
    result.push_back(profile);
  }
  std::ranges::sort(result, [](const auto& lhs, const auto& rhs) {
    const auto lhs_bytes = lhs.h2d_bytes + lhs.d2h_bytes;
    const auto rhs_bytes = rhs.h2d_bytes + rhs.d2h_bytes;
    if (lhs_bytes != rhs_bytes) return lhs_bytes > rhs_bytes;
    return lhs.label < rhs.label;
  });
  return result;
}

void gpu_transfer_record_slab_x_roundtrip(double seconds) {
  if (state().enabled) state().slab_x_roundtrip_s += seconds;
}

}  // namespace gutibm
