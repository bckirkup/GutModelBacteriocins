#ifndef GUTIBM_IMMIGRATION_H
#define GUTIBM_IMMIGRATION_H

#include "immigration_config.h"
#include "random.h"

#include <functional>
#include <vector>

namespace gutibm {

using ImmigrationDistanceReducer =
    std::function<void(const std::vector<Vec3>&, std::vector<Real>&)>;

Int immigration_event_count(const ImmigrationConfig& cfg, Int relative_step,
                            Real dt, RNG& rng);

std::vector<Vec3> immigration_positions(
    const ImmigrationConfig& cfg, const Vec3& lo, const Vec3& hi, RNG& rng,
    bool has_live_agents, bool log_warnings,
    const ImmigrationDistanceReducer& reduce_distances);

}  // namespace gutibm

#endif  // GUTIBM_IMMIGRATION_H
