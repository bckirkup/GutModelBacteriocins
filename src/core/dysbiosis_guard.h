/* -----------------------------------------------------------------------
   GutIBM – Accelerating dysbiosis guard
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_DYSBIOSIS_GUARD_H
#define GUTIBM_DYSBIOSIS_GUARD_H

#include "types.h"

#include <vector>

namespace gutibm {

class DysbiosisGuard {
 public:
  DysbiosisGuard() = default;
  DysbiosisGuard(Real threshold, Real sampling_interval, Int sample_count);

  void configure(Real threshold, Real sampling_interval, Int sample_count);
  void reset(Real current_time);
  bool observe(Real current_time, Real density_cells_per_mL);

  bool halted() const { return halted_; }
  Real halt_density_cells_per_mL() const { return halt_density_cells_per_mL_; }
  Real density_rate_cells_per_mL_per_s() const {
    return density_rate_cells_per_mL_per_s_;
  }
  const std::vector<Real>& density_history() const { return density_history_; }

 private:
  bool is_accelerating_window() const;

  Real threshold_ = 0.0;
  Real sampling_interval_ = 0.0;
  Int sample_count_ = 0;
  bool halted_ = false;
  Real halt_density_cells_per_mL_ = 0.0;
  Real density_rate_cells_per_mL_per_s_ = 0.0;
  std::vector<Real> density_history_;
  Real next_sample_time_ = 0.0;
};

}  // namespace gutibm

#endif  // GUTIBM_DYSBIOSIS_GUARD_H
