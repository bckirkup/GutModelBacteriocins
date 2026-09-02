#include "gpu_profile.h"

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

GpuTransferState& state() {
  static GpuTransferState s;
  return s;
}
}  // namespace

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
}

void gpu_transfer_record_h2d(double seconds, unsigned long long bytes) {
  if (state().enabled) {
    auto& s = state();
    s.h2d_s += seconds;
    s.h2d_bytes += bytes;
    ++s.h2d_calls;
  }
}

void gpu_transfer_record_d2h(double seconds, unsigned long long bytes) {
  if (state().enabled) {
    auto& s = state();
    s.d2h_s += seconds;
    s.d2h_bytes += bytes;
    ++s.d2h_calls;
  }
}

GpuTransferProfile gpu_transfer_profile_snapshot() {
  const auto& s = state();
  return {s.h2d_s, s.d2h_s, s.h2d_bytes, s.d2h_bytes, s.h2d_calls,
          s.d2h_calls, s.slab_x_roundtrip_s};
}

void gpu_transfer_record_slab_x_roundtrip(double seconds) {
  if (state().enabled) state().slab_x_roundtrip_s += seconds;
}

}  // namespace gutibm
