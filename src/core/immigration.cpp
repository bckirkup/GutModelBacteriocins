#include "immigration.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <format>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {

namespace {

constexpr Int kCandidateBatchSize = 512;
constexpr Real kBoundaryEpsilon = 1.0e-15;

Vec3 random_unit_direction(RNG& rng) {
  Vec3 direction;
  Real norm_sq = 0.0;
  do {
    for (Real& value : direction) value = rng.gaussian(0.0, 1.0);
    norm_sq = direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2];
  } while (norm_sq <= 0.0);
  const Real inverse_norm = 1.0 / std::sqrt(norm_sq);
  for (Real& value : direction) value *= inverse_norm;
  return direction;
}

constexpr std::array<Vec3, 6> kAxisDirections = {
    Vec3{1.0, 0.0, 0.0}, Vec3{-1.0, 0.0, 0.0},
    Vec3{0.0, 1.0, 0.0}, Vec3{0.0, -1.0, 0.0},
    Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, -1.0}};

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

std::vector<Vec3> shell_candidate_batch(
    const ImmigrationConfig& cfg, const Vec3& lo, const Vec3& hi, RNG& rng,
    const std::vector<Vec3>& anchors,
    const ImmigrationPositionProjector& project_position) {
  std::vector<Vec3> candidates;
  candidates.reserve(kCandidateBatchSize);
  const Int deterministic_count = std::min<Int>(
      static_cast<Int>(anchors.size()), kCandidateBatchSize);
  for (Int i = 0; i < deterministic_count; ++i) {
    Vec3 candidate = anchors[static_cast<size_t>(i)];
    const Vec3& direction =
        kAxisDirections[static_cast<size_t>(i % kAxisDirections.size())];
    for (Int axis = 0; axis < 3; ++axis) {
      candidate[axis] += cfg.distance * direction[axis];
    }
    project_position(candidate);
    candidate[2] = std::clamp(candidate[2], lo[2], hi[2] - kBoundaryEpsilon);
    candidates.push_back(candidate);
  }
  for (Int i = deterministic_count; i < kCandidateBatchSize; ++i) {
    const Int anchor_index =
        rng.randint(0, static_cast<Int>(anchors.size()) - 1);
    Vec3 candidate = anchors[static_cast<size_t>(anchor_index)];
    const Vec3 direction = random_unit_direction(rng);
    for (Int axis = 0; axis < 3; ++axis) {
      candidate[axis] += cfg.distance * direction[axis];
    }
    project_position(candidate);
    candidate[2] = std::clamp(candidate[2], lo[2], hi[2] - kBoundaryEpsilon);
    candidates.push_back(candidate);
  }
  return candidates;
}

std::pair<Int, Real> best_distance_candidate(
    const std::vector<Real>& distances_sq, Real target_distance) {
  Int best = 0;
  Real best_error = std::numeric_limits<Real>::max();
  for (Int i = 0; i < kCandidateBatchSize; ++i) {
    const Real error = std::abs(
        std::sqrt(distances_sq[static_cast<size_t>(i)]) - target_distance);
    if (error < best_error) {
      best = i;
      best_error = error;
    }
  }
  return {best, best_error};
}

}  // namespace

void ImmigrationEngine::validate(const ImmigrationConfig& immigration,
                                 Int initial_strain_count, const Vec3& lo,
                                 const Vec3& hi) const {
  if (!immigration.enabled) return;
  if (immigration.count < 0) {
    throw ConfigError("immigration.count must be non-negative");
  }
  if (immigration.strain_index < 0 ||
      immigration.strain_index >= initial_strain_count) {
    throw ConfigError("immigration.strain_index is outside initial_strains");
  }
  if (immigration.placement == "at_distance" &&
      immigration.distance <= 0.0) {
    throw ConfigError("immigration.distance must be positive for at_distance");
  }
  if (immigration.placement == "z_slab" &&
      immigration.z_min >= immigration.z_max) {
    throw ConfigError("immigration.z_min must be less than z_max");
  }
  if (immigration.placement == "z_slab" &&
      (immigration.z_min < lo[2] || immigration.z_max > hi[2])) {
    throw ConfigError(std::format(
        "immigration z_slab band [{}, {}] lies outside domain z [{}, {}]",
        immigration.z_min, immigration.z_max, lo[2], hi[2]));
  }
  if (immigration.distance_tolerance < 0.0 ||
      immigration.rate < 0.0) {
    throw ConfigError("immigration tolerance and rate must be non-negative");
  }
}

std::vector<Vec3> ImmigrationEngine::anchors(
    const ImmigrationConfig& immigration, const AgentPool& agents,
    const Domain& domain, Int global_live_count) const {
  if (global_live_count <= 0) return {};
  if (immigration.distance_reference == "centroid") {
    return {centroid_anchor(agents, domain, global_live_count)};
  }
  return support_anchors(agents, domain);
}

Vec3 ImmigrationEngine::centroid_anchor(const AgentPool& agents,
                                         const Domain& domain,
                                         Int global_live_count) const {
#ifndef GUTIBM_MPI
  (void)domain;
#endif
  Vec3 local_sum = {0.0, 0.0, 0.0};
  for (const Agent& agent : agents) {
    if (agent.state == PhenoState::DEAD) continue;
    for (Int axis = 0; axis < 3; ++axis) local_sum[axis] += agent.x[axis];
  }
  Vec3 global_sum = local_sum;
#ifdef GUTIBM_MPI
  if (domain.nprocs() > 1) {
    MPI_Allreduce(local_sum.data(), global_sum.data(), 3, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  for (Real& value : global_sum) {
    value /= static_cast<Real>(global_live_count);
  }
  return global_sum;
}

std::vector<Vec3> ImmigrationEngine::support_anchors(
    const AgentPool& agents, const Domain& domain) const {
#ifndef GUTIBM_MPI
  (void)domain;
#endif
  constexpr std::array<Vec3, 6> directions = {
      Vec3{1.0, 0.0, 0.0}, Vec3{-1.0, 0.0, 0.0},
      Vec3{0.0, 1.0, 0.0}, Vec3{0.0, -1.0, 0.0},
      Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, -1.0}};
  std::array<Vec3, directions.size()> local_support{};
  std::array<Real, directions.size()> local_best;
  local_best.fill(-std::numeric_limits<Real>::max());
  for (const Agent& agent : agents) {
    if (agent.state == PhenoState::DEAD) continue;
    for (size_t i = 0; i < directions.size(); ++i) {
      const Real projection = agent.x[0] * directions[i][0]
                            + agent.x[1] * directions[i][1]
                            + agent.x[2] * directions[i][2];
      if (projection > local_best[i]) {
        local_best[i] = projection;
        local_support[i] = agent.x;
      }
    }
  }
  std::vector<Vec3> sampled(local_support.begin(), local_support.end());
#ifdef GUTIBM_MPI
  if (const Int nprocs = domain.nprocs(); nprocs > 1) {
    const auto local_values = static_cast<Int>(sampled.size() * 3);
    std::vector<Int> counts(static_cast<size_t>(nprocs));
    MPI_Allgather(&local_values, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    std::vector displacements(static_cast<size_t>(nprocs), Int{0});
    Int total_values = 0;
    for (Int rank = 0; rank < nprocs; ++rank) {
      displacements[static_cast<size_t>(rank)] = total_values;
      total_values += counts[static_cast<size_t>(rank)];
    }
    std::vector<Real> gathered_values(static_cast<size_t>(total_values));
    MPI_Allgatherv(sampled.empty() ? nullptr : sampled.front().data(),
                   local_values, MPI_DOUBLE, gathered_values.data(),
                   counts.data(), displacements.data(), MPI_DOUBLE,
                   MPI_COMM_WORLD);
    std::vector<Vec3> anchors;
    anchors.reserve(static_cast<size_t>(total_values / 3));
    for (Int i = 0; i < total_values; i += 3) {
      anchors.emplace_back(Vec3{
          {gathered_values[static_cast<size_t>(i)],
           gathered_values[static_cast<size_t>(i + 1)],
           gathered_values[static_cast<size_t>(i + 2)]}});
    }
    return anchors;
  }
#endif
  return sampled;
}

void ImmigrationEngine::reduce_distances(
    const ImmigrationConfig& immigration,
    const std::vector<Vec3>& candidates, std::vector<Real>& distances_sq,
    const AgentPool& agents, const Domain& domain) const {
  if (immigration.distance_reference == "centroid") {
    Vec3 centroid = {0.0, 0.0, 0.0};
    Int local_count = 0;
    for (const Agent& agent : agents) {
      if (agent.state == PhenoState::DEAD) continue;
      for (Int axis = 0; axis < 3; ++axis) centroid[axis] += agent.x[axis];
      ++local_count;
    }
#ifdef GUTIBM_MPI
    if (domain.nprocs() > 1) {
      Vec3 global_sum = {0.0, 0.0, 0.0};
      Int global_count = 0;
      MPI_Allreduce(centroid.data(), global_sum.data(), 3, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);
      centroid = global_sum;
      local_count = global_count;
    }
#endif
    if (local_count > 0) {
      for (Real& value : centroid) value /= static_cast<Real>(local_count);
      for (size_t i = 0; i < candidates.size(); ++i) {
        distances_sq[i] = domain.min_image_dist_sq(candidates[i], centroid);
      }
    }
  } else {
    for (size_t i = 0; i < candidates.size(); ++i) {
      Real nearest = std::numeric_limits<Real>::max();
      for (const Agent& agent : agents) {
        if (agent.state == PhenoState::DEAD) continue;
        nearest = std::min(
            nearest, domain.min_image_dist_sq(candidates[i], agent.x));
      }
      distances_sq[i] = nearest;
    }
  }
#ifdef GUTIBM_MPI
  if (domain.nprocs() > 1) {
    MPI_Allreduce(MPI_IN_PLACE, distances_sq.data(),
                  static_cast<int>(distances_sq.size()), MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
  }
#endif
}

// MPI contract: rng_ is replicated and never observes rank-local state, so
// every rank agrees on event counts and candidates. When an event fires, every
// rank also participates in the replicated anchor Allgatherv before the
// distance reduction. Both collectives are unconditional for that event. Only
// the owning rank constructs a cell; AgentPool::next_tag() supplies that rank's
// stride stream, making IDs globally unique without migration duplicates.
void ImmigrationEngine::inject(const ImmigrationConfig& immigration,
                               Int current_step, Real dt,
                               const AgentPool& agents, const Domain& domain,
                               const ImmigrationAgentFactory& create_agent) {
  if (!immigration.enabled) return;
  const Int relative_step = current_step - start_step_;
  const Int event_count =
      immigration_event_count(immigration, relative_step, dt, rng_);
  if (event_count == 0) return;

  const bool log_warnings = domain.rank() == 0;
  for (Int event = 0; event < event_count; ++event) {
    Int local_live_count = 0;
    for (const Agent& agent : agents) {
      if (agent.state != PhenoState::DEAD) ++local_live_count;
    }
    Int global_live_count = local_live_count;
#ifdef GUTIBM_MPI
    if (domain.nprocs() > 1) {
      MPI_Allreduce(&local_live_count, &global_live_count, 1, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    const bool has_live_agents = global_live_count > 0;
    const std::vector<Vec3> immigration_anchors =
        anchors(immigration, agents, domain, global_live_count);
    const auto reducer = [this, &immigration, &agents, &domain](
                             const std::vector<Vec3>& candidates,
                             std::vector<Real>& distances_sq) {
      reduce_distances(immigration, candidates, distances_sq, agents, domain);
    };
    const std::vector<Vec3> positions = immigration_positions(
        immigration, domain.lo(), domain.hi(), rng_, immigration_anchors,
        has_live_agents, log_warnings, reducer,
        [&domain](Vec3& position) { domain.apply_pbc(position); });
    if (immigration.placement == "at_distance" &&
        static_cast<Int>(positions.size()) != immigration.count) {
      throw SimulationError(
          "immigration at_distance could not place all requested immigrants");
    }
    for (const Vec3& position : positions) {
      if (domain.is_local(position)) create_agent(position);
    }
  }
}

Int immigration_event_count(const ImmigrationConfig& cfg, Int relative_step,
                            Real dt, RNG& rng) {
  if (!cfg.enabled) return 0;
  if (cfg.schedule == "pulse") return relative_step == cfg.step ? 1 : 0;
  return rng.poisson(std::max<Real>(0.0, cfg.rate * dt));
}

std::vector<Vec3> immigration_positions(
    const ImmigrationConfig& cfg, const Vec3& lo, const Vec3& hi, RNG& rng,
    const std::vector<Vec3>& anchors, bool has_live_agents, bool log_warnings,
    const ImmigrationDistanceReducer& reduce_distances,
    const ImmigrationPositionProjector& project_position) {
  std::vector<Vec3> result;
  result.reserve(static_cast<size_t>(std::max<Int>(0, cfg.count)));
  for (Int immigrant = 0; immigrant < cfg.count; ++immigrant) {
    if (cfg.placement != "at_distance") {
      result.push_back(random_position(cfg, lo, hi, rng));
      continue;
    }
    if (!has_live_agents) {
      if (log_warnings) {
        std::cerr << "Warning: immigration at_distance has no live biomass; "
                     "falling back to uniform placement\n";
      }
      const std::vector<Vec3> candidates = candidate_batch(cfg, lo, hi, rng);
      std::vector unused_distances(kCandidateBatchSize,
                                   std::numeric_limits<Real>::max());
      reduce_distances(candidates, unused_distances);
      result.push_back(random_position(cfg, lo, hi, rng));
      continue;
    }

    Vec3 selected{};
    bool found = false;
    for (Int attempt = 0; attempt < 2; ++attempt) {
      const std::vector<Vec3> candidates = shell_candidate_batch(
          cfg, lo, hi, rng, anchors, project_position);
      std::vector distances_sq(kCandidateBatchSize,
                               std::numeric_limits<Real>::max());
      reduce_distances(candidates, distances_sq);
      const auto [best, best_error] =
          best_distance_candidate(distances_sq, cfg.distance);
      selected = candidates[static_cast<size_t>(best)];
      if (best_error <= cfg.distance_tolerance) {
        found = true;
        break;
      }
    }
    if (!found) {
      if (log_warnings) {
        std::cerr << "Warning: immigration candidate did not meet distance "
                     "tolerance; skipping cell\n";
      }
      continue;
    }
    result.push_back(selected);
  }
  return result;
}

}  // namespace gutibm
