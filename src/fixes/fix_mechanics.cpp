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

void apply_pair_displacement(Agent& ai, Agent& aj, const Vec3& n,
                             Real force_mag, Real dt) {
  Real mi = std::max(ai.mass, 1.0e-30);
  Real mj = std::max(aj.mass, 1.0e-30);
  Real inv_mi = 1.0 / mi;
  Real inv_mj = 1.0 / mj;
  Real inv_sum = 1.0 / (inv_mi + inv_mj);

  Real push_i = force_mag * dt * inv_mi * inv_mi * inv_sum;
  Real push_j = force_mag * dt * inv_mj * inv_mj * inv_sum;

  ai.x[0] -= n[0] * push_i;
  ai.x[1] -= n[1] * push_i;
  ai.x[2] -= n[2] * push_i;
  aj.x[0] += n[0] * push_j;
  aj.x[1] += n[1] * push_j;
  aj.x[2] += n[2] * push_j;
}

void apply_repulsion(Agent& ai, Agent& aj, const PairGeometry& geom,
                     const MechanicsConfig& cfg, Real dt) {
  Real overlap = geom.sum_r - geom.d;
  if (overlap <= 0.0) return;

  Real force_mag = cfg.hertzian_enabled
      ? cfg.hertz_k * std::pow(overlap, 1.5)
      : cfg.hertz_k * overlap;
  apply_pair_displacement(ai, aj, geom.n, force_mag, dt);
}

void apply_adhesion(Agent& ai, Agent& aj, const PairGeometry& geom,
                    const MechanicsConfig& cfg, Real dt) {
  if (!cfg.adhesion_enabled) return;
  if (is_active_corpse(ai) || is_active_corpse(aj)) return;

  Real gap = geom.d - geom.sum_r;
  if (gap <= 0.0 || gap >= cfg.adhesion_range) return;

  Real adhesion_frac = 1.0 - (gap / cfg.adhesion_range);
  Real adhesion_force = cfg.adhesion_strength * adhesion_frac;

  Real mi = std::max(ai.mass, 1.0e-30);
  Real mj = std::max(aj.mass, 1.0e-30);
  Real inv_mi = 1.0 / mi;
  Real inv_mj = 1.0 / mj;
  Real inv_sum = 1.0 / (inv_mi + inv_mj);

  Real pull_i = adhesion_force * dt * inv_mi * inv_mi * inv_sum;
  Real pull_j = adhesion_force * dt * inv_mj * inv_mj * inv_sum;

  ai.x[0] += geom.n[0] * pull_i;
  ai.x[1] += geom.n[1] * pull_i;
  ai.x[2] += geom.n[2] * pull_i;
  aj.x[0] -= geom.n[0] * pull_j;
  aj.x[1] -= geom.n[1] * pull_j;
  aj.x[2] -= geom.n[2] * pull_j;
}

void resolve_agent_pair(Agent& ai, Agent& aj, const Domain& domain,
                        const MechanicsConfig& cfg, Real dt) {
  Vec3 delta = domain.min_image_delta(ai.x, aj.x);
  PairGeometry geom{};
  if (!compute_pair_geometry(delta, ai.radius + aj.radius, geom)) return;

  apply_repulsion(ai, aj, geom, cfg, dt);
  apply_adhesion(ai, aj, geom, cfg, dt);
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

  if (!gpu_run_mechanics(ag, n, hash, dom, cfg, dt)) {
    return false;
  }

  ag.sync_positions_to_host(agents);
  return true;
}

}  // namespace

void FixMechanics::compute(Real dt) {
  if (try_gpu_mechanics(sim_, cfg_, dt)) return;

  auto& agents = sim_.agents();
  const auto& hash = sim_.domain().spatial_hash();
  const Real sim_time = sim_.time();
  const auto& sim_cfg = sim_.config();

  for (Int i = 0; i < agents.size(); ++i) {
    Agent& ai = agents[i];
    if (!participates_in_mechanics(ai, sim_time, sim_cfg)) continue;

    auto neighbors = hash.query_neighbors(ai.x);
    for (Int j : neighbors) {
      if (j <= i) continue;
      Agent& aj = agents[j];
      if (!participates_in_mechanics(aj, sim_time, sim_cfg)) continue;
      resolve_agent_pair(ai, aj, sim_.domain(), cfg_, dt);
    }
  }
}

}  // namespace gutibm
