/* -----------------------------------------------------------------------
   GutIBM – Deterministic (seeded) random number generation

   THREAD SAFETY: RNG wraps a single std::mt19937_64 and is NOT safe for
   concurrent calls from multiple threads — doing so would corrupt the engine
   state and destroy run-to-run reproducibility. Under OpenMP, use one RNG per
   thread (seeded distinctly, e.g. seed + thread_id) or draw all random variates
   serially. Fix modules follow the latter pattern: e.g. fix_receptor computes
   kill probabilities in parallel but applies the RNG draws / kills serially.
   (A per-thread RNG pool is tracked as future work — Spec 0 §9.)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_RANDOM_H
#define GUTIBM_RANDOM_H

#include "types.h"
#include <cmath>
#include <random>

namespace gutibm {

class RNG {
 public:
  explicit RNG(uint64_t seed = 42) : gen_(seed) {}

  void seed(uint64_t s) { gen_.seed(s); }

  Real uniform()     { return udist_(gen_); }
  Real uniform(Real lo, Real hi) {
    return lo + (hi - lo) * udist_(gen_);
  }

  Real gaussian(Real mean, Real stddev) {
    std::normal_distribution d(mean, stddev);
    return d(gen_);
  }

  Int randint(Int lo, Int hi) {
    std::uniform_int_distribution d(lo, hi);
    return d(gen_);
  }

  bool bernoulli(Real p) {
    return udist_(gen_) < p;
  }

  Real exponential(Real lambda) {
    std::exponential_distribution d(lambda);
    return d(gen_);
  }

  // Poisson-distributed random variate
  Int poisson(Real mean) {
    std::poisson_distribution d(mean);
    return d(gen_);
  }

  std::mt19937_64& engine() { return gen_; }  // NOSONAR cpp:S2245 — simulation PRNG, not cryptography

 private:
  std::mt19937_64 gen_;  // NOSONAR cpp:S2245 — reproducible IbM stochasticity via explicit seed
  std::uniform_real_distribution<Real> udist_{0.0, 1.0};
};

}  // namespace gutibm

#endif  // GUTIBM_RANDOM_H
