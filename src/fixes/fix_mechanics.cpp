/* -----------------------------------------------------------------------
   GutIBM – Soft-sphere mechanical repulsion implementation
   ----------------------------------------------------------------------- */

#include "fix_mechanics.h"
#include "simulation.h"
#include "mechanics_gpu.h"
#include "spatial_hash_gpu.h"
#include "dispatch.h"

#include <cmath>
#include <algorithm>
#include <vector>

namespace gutibm {

namespace {

bool participates_in_mechanics(const Agent& a, Real sim_time,
                               const SimulationConfig& cfg) {
  if (a.state != PhenoState::DEAD) return true;
  return cfg.cell_bio.cdi.enabled && a.timers.death_time >= 0.0
      && (sim_time - a.timers.death_time) < cfg.cell_bio.cdi.corpse_persistence;
}

bool is_active_corpse(const Agent& a) {
  return a.state == PhenoState::DEAD;
}

struct PairGeometry {
  Vec3 n;
  Real d;
  Real sum_r;
};

bool compute_pair_geometry(const Vec3& delta, Real sum_r, PairGeometry& geom) {
  Real d2 = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
  if (d2 <= 0.0) return false;
  geom.d = std::sqrt(d2);
  geom.sum_r = sum_r;
  geom.n = {delta[0] / geom.d, delta[1] / geom.d, delta[2] / geom.d};
  return true;
}

constexpr Real kMinimumViscosity = 1.0e-30;
constexpr Real kMinimumRadius = 1.0e-30;
constexpr Real kPi = 3.14159265358979323846;

Real mobility(const Agent& agent, Real viscosity) {
  const Real radius = std::max(agent.radius, kMinimumRadius);
  const Real effective_viscosity = std::max(viscosity, kMinimumViscosity);
  return 1.0 / (6.0 * kPi * effective_viscosity * radius);
}

void apply_pair_displacement(Vec3& displacement_i, Vec3& displacement_j,
                             const Agent& ai, const Agent& aj, const Vec3& n,
                             Real force_mag, Real dt, Real viscosity,
                             Real direction) {
  const Real mobility_i = mobility(ai, viscosity);
  const Real mobility_j = mobility(aj, viscosity);
  const Real mobility_sum = mobility_i + mobility_j;
  const Real relative_displacement = force_mag * dt * mobility_sum;
  const Real push_i = relative_displacement * mobility_i / mobility_sum;
  const Real push_j = relative_displacement * mobility_j / mobility_sum;

  displacement_i[0] -= direction * n[0] * push_i;
  displacement_i[1] -= direction * n[1] * push_i;
  displacement_i[2] -= direction * n[2] * push_i;
  displacement_j[0] += direction * n[0] * push_j;
  displacement_j[1] += direction * n[1] * push_j;
  displacement_j[2] += direction * n[2] * push_j;
}

void apply_repulsion(const Agent& ai, const Agent& aj,
                     const PairGeometry& geom, const MechanicsConfig& cfg,
                     Real dt, Real viscosity, Vec3& displacement_i,
                     Vec3& displacement_j) {
  Real overlap = geom.sum_r - geom.d;
  if (overlap <= 0.0) return;

  Real force_mag = cfg.hertzian_enabled
      ? cfg.hertz_k * std::pow(overlap, 1.5)
      : cfg.hertz_k * overlap;
  apply_pair_displacement(displacement_i, displacement_j, ai, aj, geom.n,
                          force_mag, dt, viscosity, 1.0);
}

void apply_adhesion(const Agent& ai, const Agent& aj,
                    const PairGeometry& geom, const MechanicsConfig& cfg,
                    Real dt, Real viscosity, Vec3& displacement_i,
                    Vec3& displacement_j) {
  if (!cfg.adhesion_enabled) return;
  if (is_active_corpse(ai) || is_active_corpse(aj)) return;

  Real gap = geom.d - geom.sum_r;
  if (gap <= 0.0 || gap >= cfg.adhesion_range) return;

  Real adhesion_frac = 1.0 - (gap / cfg.adhesion_range);
  Real adhesion_force = cfg.adhesion_strength * adhesion_frac;
  apply_pair_displacement(displacement_i, displacement_j, ai, aj, geom.n,
                          adhesion_force, dt, viscosity, -1.0);
}

void resolve_agent_pair(const Agent& ai, const Agent& aj, const Domain& domain,
                        const MechanicsConfig& cfg, Real dt, Real viscosity,
                        Vec3& displacement_i, Vec3& displacement_j) {
  Vec3 delta = domain.min_image_delta(ai.x, aj.x);
  PairGeometry geom{};
  if (!compute_pair_geometry(delta, ai.radius + aj.radius, geom)) return;

  apply_repulsion(ai, aj, geom, cfg, dt, viscosity, displacement_i,
                  displacement_j);
  apply_adhesion(ai, aj, geom, cfg, dt, viscosity, displacement_i,
                 displacement_j);
}

void apply_domain_constraints(Vec3& position, const Domain& domain) {
  const Vec3 size = domain.size();
  const auto& periodic = domain.config().periodic;
  for (Int axis = 0; axis < 2; ++axis) {
    if (!periodic[static_cast<size_t>(axis)]) {
      position[static_cast<size_t>(axis)] = std::clamp(
          position[static_cast<size_t>(axis)], domain.lo()[axis],
          std::nextafter(domain.hi()[axis], domain.lo()[axis]));
      continue;
    }
    while (position[static_cast<size_t>(axis)] < domain.lo()[axis]) {
      position[static_cast<size_t>(axis)] += size[axis];
    }
    while (position[static_cast<size_t>(axis)] >= domain.hi()[axis]) {
      position[static_cast<size_t>(axis)] -= size[axis];
    }
  }
  position[2] = std::max(position[2], domain.lo()[2]);
}

Int apply_displacements(Simulation& sim, std::vector<Vec3>& displacements) {
  const auto& domain = sim.domain();
  Int clamp_count = 0;
  auto& agents = sim.agents();
  for (Int i = 0; i < agents.size(); ++i) {
    Vec3& displacement = displacements[static_cast<size_t>(i)];
    const Real norm = std::sqrt(
        displacement[0] * displacement[0]
        + displacement[1] * displacement[1]
        + displacement[2] * displacement[2]);
    const Real limit =
        kMechanicsMaxDisplacementRadiusFraction * agents[i].radius;
    if (norm > limit && norm > 0.0) {
      const Real scale = limit / norm;
      displacement[0] *= scale;
      displacement[1] *= scale;
      displacement[2] *= scale;
      if (!agents[i].flags.is_ghost) ++clamp_count;
    }
    agents[i].x[0] += displacement[0];
    agents[i].x[1] += displacement[1];
    agents[i].x[2] += displacement[2];
    apply_domain_constraints(agents[i].x, domain);
  }
  return clamp_count;
}

bool try_gpu_mechanics(Simulation& sim, const MechanicsConfig& cfg, Real dt) {
  if (!sim.gpu_active()) return false;

  auto& agents = sim.agents();
  const Int n = agents.size();
  if (n <= 0) return false;

  auto& ag = sim.agents_gpu();
  ag.sync_from_host(agents);

  SpatialHashGpu hash;
  const auto& dom = sim.domain();
  if (!gpu_build_spatial_hash(
          ag, n, dom.lo(), dom.hi(), dom.spatial_hash().cell_size(), hash)) {
    return false;
  }

  Int clamp_count = 0;
  if (!gpu_run_mechanics(ag, n, hash, dom, cfg, dt,
                         sim.config().vbf.viscosity, clamp_count)) {
    return false;
  }

  ag.sync_positions_to_host(agents);
  sim.mechanics_step_stats().displacement_clamps += clamp_count;
  return true;
}

}  // namespace

void FixMechanics::compute(Real dt) {
  if (try_gpu_mechanics(sim_, cfg_, dt)) return;

  auto& agents = sim_.agents();
  const auto& hash = sim_.domain().spatial_hash();
  const Real sim_time = sim_.time();
  const auto& sim_cfg = sim_.config();
  const Real viscosity = sim_cfg.vbf.viscosity;
  std::vector<Vec3> displacements(static_cast<size_t>(agents.size()));

  for (Int i = 0; i < agents.size(); ++i) {
    Agent& ai = agents[i];
    if (!participates_in_mechanics(ai, sim_time, sim_cfg)) continue;

    auto neighbors = hash.query_neighbors(ai.x);
    for (Int j : neighbors) {
      if (j <= i) continue;
      Agent& aj = agents[j];
      if (!participates_in_mechanics(aj, sim_time, sim_cfg)) continue;
      resolve_agent_pair(ai, aj, sim_.domain(), cfg_, dt, viscosity,
                         displacements[static_cast<size_t>(i)],
                         displacements[static_cast<size_t>(j)]);
    }
  }
  sim_.mechanics_step_stats().displacement_clamps +=
      apply_displacements(sim_, displacements);
}

}  // namespace gutibm
