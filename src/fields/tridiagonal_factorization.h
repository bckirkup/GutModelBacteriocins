#pragma once

#include "types.h"
#include <vector>

namespace gutibm {

class TridiagonalFactorization {
 public:
  void factorize(const std::vector<Real>& lower,
                 const std::vector<Real>& diagonal,
                 const std::vector<Real>& upper) {
    diagonal_ = diagonal;
    upper_ = upper;
    multipliers_.assign(lower.size(), 0.0);
    for (size_t i = 1; i < diagonal_.size(); ++i) {
      const Real multiplier = lower[i - 1] / diagonal_[i - 1];
      multipliers_[i - 1] = multiplier;
      diagonal_[i] -= multiplier * upper_[i - 1];
    }
  }

  void solve_in_place(std::vector<Real>& values) const {
    for (size_t i = 1; i < values.size(); ++i) {
      values[i] -= multipliers_[i - 1] * values[i - 1];
    }
    values.back() /= diagonal_.back();
    for (size_t i = values.size() - 1; i > 0; --i) {
      values[i - 1] =
          (values[i - 1] - upper_[i - 1] * values[i]) / diagonal_[i - 1];
    }
  }

 private:
  std::vector<Real> diagonal_;
  std::vector<Real> upper_;
  std::vector<Real> multipliers_;
};

}  // namespace gutibm
