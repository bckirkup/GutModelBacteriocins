#include "immigration.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace gutibm {

namespace {

constexpr Int kCandidateBatchSize = 512;
constexpr Real kBoundaryEpsilon = 1.0e-15;

Vec3 random_position(const ImmigrationConfig& cfg, const Vec3& lo,
                     const Vec3& hi, RNG& rng) {
  Vec3 pos = {
      rng.uniform(lo[0], hi[0]),
      rng.uniform(lo[1], hi[1]),
      rng.uniform(lo[2], hi[2]),
  };
  if (cfg.placement == "z_slab") {
    const Real upper = std::max(cfg.z_min, cfg.z_max - kBoundaryEpsilon);
    pos[2] = rng.uniform(cfg.z_min, upper);
  }
  pos[2] = std::clamp(pos[2], lo[2], hi[2] - kBoundaryEpsilon);
  return pos;
}

std::vector<Vec3> candidate_batch(const ImmigrationConfig& cfg, const Vec3& lo,
                                  const Vec3& hi, RNG& rng) {
  std::vector<Vec3> candidates;
  candidates.reserve(kCandidateBatchSize);
  for (Int i = 0; i < kCandidateBatchSize; ++i) {
    candidates.push_back(random_position(cfg, lo, hi, rng));
  }
  return candidates;
}

}  // namespace

Int immigration_event_count(const ImmigrationConfig& cfg, Int relative_step,
                            Real dt, RNG& rng) {
  if (!cfg.enabled) return 0;
  if (cfg.schedule == "pulse") return relative_step == cfg.step ? 1 : 0;
  return rng.poisson(std::max<Real>(0.0, cfg.rate * dt));
}

std::vector<Vec3> immigration_positions(
    const ImmigrationConfig& cfg, const Vec3& lo, const Vec3& hi, RNG& rng,
    bool has_live_agents, const ImmigrationDistanceReducer& reduce_distances) {
  std::vector<Vec3> result;
  result.reserve(static_cast<size_t>(std::max<Int>(0, cfg.count)));
  for (Int immigrant = 0; immigrant < cfg.count; ++immigrant) {
    if (cfg.placement != "at_distance") {
      result.push_back(random_position(cfg, lo, hi, rng));
      continue;
    }
    if (!has_live_agents) {
      std::cerr << "Warning: immigration at_distance has no live biomass; "
                   "falling back to uniform placement\n";
      const std::vector<Vec3> candidates = candidate_batch(cfg, lo, hi, rng);
      std::vector<Real> unused_distances(kCandidateBatchSize,
                                         std::numeric_limits<Real>::max());
      reduce_distances(candidates, unused_distances);
      result.push_back(random_position(cfg, lo, hi, rng));
      continue;
    }

    Vec3 selected{};
    bool found = false;
    for (Int attempt = 0; attempt < 2; ++attempt) {
      const std::vector<Vec3> candidates = candidate_batch(cfg, lo, hi, rng);
      std::vector<Real> distances_sq(kCandidateBatchSize,
                                     std::numeric_limits<Real>::max());
      reduce_distances(candidates, distances_sq);
      Int best = 0;
      Real best_error = std::numeric_limits<Real>::max();
      for (Int i = 0; i < kCandidateBatchSize; ++i) {
        const Real error = std::abs(std::sqrt(distances_sq[static_cast<size_t>(i)]) -
                                    cfg.distance);
        if (error < best_error) {
          best = i;
          best_error = error;
        }
      }
      selected = candidates[static_cast<size_t>(best)];
      if (best_error <= cfg.distance_tolerance) {
        found = true;
        break;
      }
    }
    if (!found) {
      std::cerr << "Warning: immigration candidate did not meet distance "
                   "tolerance; skipping cell\n";
      continue;
    }
    result.push_back(selected);
  }
  return result;
}

}  // namespace gutibm
