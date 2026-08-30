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
  const size_t first =
      density_history_.size() - static_cast<size_t>(sample_count_);
  // Report the average net density rise across the accepted window.
  density_rate_cells_per_mL_per_s_ =
      (density_history_.back() - density_history_[first]) /
      (static_cast<Real>(sample_count_ - 1) * sampling_interval_);
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
  }

  const auto increment_count =
      static_cast<size_t>(sample_count_ - 1);
  const size_t first_half_count = increment_count / 2;
  const size_t second_half_count = increment_count - first_half_count;
  Real first_half_sum = 0.0;
  Real second_half_sum = 0.0;
  for (size_t i = 0; i < increment_count; ++i) {
    const Real increment =
        density_history_[first + i + 1] - density_history_[first + i];
    if (i < first_half_count) {
      first_half_sum += increment;
    } else {
      second_half_sum += increment;
    }
  }

  const Real first_half_mean =
      first_half_sum / static_cast<Real>(first_half_count);
  const Real second_half_mean =
      second_half_sum / static_cast<Real>(second_half_count);
  return density_history_.back() > density_history_[first] &&
         first_half_mean > 0.0 && second_half_mean > 0.0 &&
         second_half_mean >= first_half_mean;
}

}  // namespace gutibm
