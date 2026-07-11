/* -----------------------------------------------------------------------
   GutIBM – Per-step wall-clock profiling accumulators
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_STEP_PROFILER_H
#define GUTIBM_STEP_PROFILER_H

#include <chrono>

namespace gutibm {

struct StepProfile {
  double ghost_exchange_s = 0.0;
  double spatial_hash_s = 0.0;
  double biology_s = 0.0;
  double chemistry_s = 0.0;
  double physics_s = 0.0;
  double mpi_migrate_s = 0.0;
  double cleanup_s = 0.0;
  double gpu_h2d_s = 0.0;
  double gpu_d2h_s = 0.0;
  double mpi_reaction_reduce_s = 0.0;
  double hdf5_s = 0.0;
  int step_count = 0;

  void reset() { *this = StepProfile{}; }

  double total_s() const {
    return ghost_exchange_s + spatial_hash_s + biology_s + chemistry_s +
           physics_s + mpi_migrate_s + cleanup_s + gpu_h2d_s + gpu_d2h_s +
           mpi_reaction_reduce_s + hdf5_s;
  }
};

class StepProfiler {
 public:
  using Clock = std::chrono::steady_clock;

  explicit StepProfiler(bool enabled) : enabled_(enabled) {}

  void start() {
    if (!enabled_) return;
    t_ = Clock::now();
  }

  void lap(double& accumulator) {
    if (!enabled_) return;
    const auto now = Clock::now();
    accumulator += std::chrono::duration<double>(now - t_).count();
    t_ = now;
  }

 private:
  bool enabled_ = false;
  Clock::time_point t_{};
};

}  // namespace gutibm

#endif  // GUTIBM_STEP_PROFILER_H
