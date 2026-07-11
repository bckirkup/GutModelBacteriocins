#include "gpu_profile.h"

namespace gutibm {

namespace {
struct GpuTransferState {
  bool enabled = false;
  double h2d_s = 0.0;
  double d2h_s = 0.0;
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
}

void gpu_transfer_record_h2d(double seconds) {
  if (state().enabled) state().h2d_s += seconds;
}

void gpu_transfer_record_d2h(double seconds) {
  if (state().enabled) state().d2h_s += seconds;
}

GpuTransferProfile gpu_transfer_profile_snapshot() {
  const auto& s = state();
  return {s.h2d_s, s.d2h_s};
}

}  // namespace gutibm
