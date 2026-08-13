/* -----------------------------------------------------------------------
   GutIBM – Accelerating dysbiosis guard
   ----------------------------------------------------------------------- */

#include "dysbiosis_guard.h"

namespace gutibm {

DysbiosisGuard::DysbiosisGuard(Real threshold,
                               Real sampling_interval,
                               Int sample_count) {
  configure(threshold, sampling_interval, sample_count);
}

void DysbiosisGuard::configure(Real threshold,
                               Real sampling_interval,
                               Int sample_count) {
  threshold_ = threshold;
  sampling_interval_ = sampling_interval;
  sample_count_ = sample_count;
}

void DysbiosisGuard::reset(Real current_time) {
  halted_ = false;
  halt_density_cells_per_mL_ = 0.0;
  density_rate_cells_per_mL_per_s_ = 0.0;
  density_history_.clear();
  next_sample_time_ = current_time;
}

bool DysbiosisGuard::observe(Real current_time, Real density_cells_per_mL) {
  if (threshold_ <= 0.0 || sampling_interval_ <= 0.0 || sample_count_ < 3) {
    return false;
  }

  if (current_time < next_sample_time_) return halted_;

  density_history_.push_back(density_cells_per_mL);
  if (density_history_.size() > static_cast<size_t>(sample_count_)) {
    density_history_.erase(density_history_.begin());
  }
  next_sample_time_ = current_time + sampling_interval_;

  if (!is_accelerating_window()) return false;

  halted_ = true;
  halt_density_cells_per_mL_ = density_cells_per_mL;
  const size_t last = density_history_.size() - 1;
  density_rate_cells_per_mL_per_s_ =
      (density_history_[last] - density_history_[last - 1]) /
      sampling_interval_;
  return true;
}

bool DysbiosisGuard::is_accelerating_window() const {
  if (density_history_.size() < static_cast<size_t>(sample_count_)) {
    return false;
  }

  const size_t first =
      density_history_.size() - static_cast<size_t>(sample_count_);
  for (size_t i = first; i < density_history_.size(); ++i) {
    if (density_history_[i] <= threshold_) return false;
    if (i == first) continue;

    const Real increment =
        density_history_[i] - density_history_[i - 1];
    if (increment <= 0.0) return false;
    if (i == first + 1) continue;

    const Real previous_increment =
        density_history_[i - 1] - density_history_[i - 2];
    if (increment < previous_increment) return false;
  }
  return true;
}

}  // namespace gutibm
