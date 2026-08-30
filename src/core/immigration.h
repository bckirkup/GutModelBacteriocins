#ifndef GUTIBM_IMMIGRATION_H
#define GUTIBM_IMMIGRATION_H

#include "immigration_config.h"
#include "random.h"
#include "agent.h"
#include "domain.h"
#include "error.h"

#include <functional>
#include <vector>

namespace gutibm {

using ImmigrationDistanceReducer =
    std::function<void(const std::vector<Vec3>&, std::vector<Real>&)>;
using ImmigrationPositionProjector = std::function<void(Vec3&)>;
using ImmigrationAgentFactory = std::function<void(const Vec3&)>;

class ImmigrationEngine {
 public:
  void seed(uint64_t seed) { rng_.seed(seed); }
  void set_start_step(Int step) { start_step_ = step; }

  void validate(const ImmigrationConfig& cfg, Int initial_strain_count,
                const Vec3& lo, const Vec3& hi) const;

  void inject(const ImmigrationConfig& cfg, Int current_step, Real dt,
              const AgentPool& agents, const Domain& domain,
              const ImmigrationAgentFactory& create_agent);

 private:
  std::vector<Vec3> anchors(const ImmigrationConfig& cfg,
                            const AgentPool& agents, const Domain& domain,
                            Int global_live_count) const;
  Vec3 centroid_anchor(const AgentPool& agents, const Domain& domain,
                       Int global_live_count) const;
  std::vector<Vec3> support_anchors(const AgentPool& agents,
                                    const Domain& domain) const;
  void reduce_distances(const ImmigrationConfig& cfg,
                        const std::vector<Vec3>& candidates,
                        std::vector<Real>& distances_sq,
                        const AgentPool& agents, const Domain& domain) const;
  RNG rng_;
  Int start_step_ = 0;
};

Int immigration_event_count(const ImmigrationConfig& cfg, Int relative_step,
                            Real dt, RNG& rng);

std::vector<Vec3> immigration_positions(
    const ImmigrationConfig& cfg, const Vec3& lo, const Vec3& hi, RNG& rng,
    const std::vector<Vec3>& anchors, bool has_live_agents, bool log_warnings,
    const ImmigrationDistanceReducer& reduce_distances,
    const ImmigrationPositionProjector& project_position);

}  // namespace gutibm

#endif  // GUTIBM_IMMIGRATION_H
