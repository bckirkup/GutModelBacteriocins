/* -----------------------------------------------------------------------
   GutIBM – Soft-sphere mechanical repulsion implementation
   ----------------------------------------------------------------------- */

#include "fix_mechanics.h"
#include "simulation.h"

#include <cmath>
#include <algorithm>

namespace gutibm {

void FixMechanics::compute(Real dt) {
  auto& agents = sim_.agents();
  const auto& hash = sim_.domain().spatial_hash();

  for (Int i = 0; i < agents.size(); ++i) {
    Agent& ai = agents[i];
    if (ai.state == PhenoState::DEAD) continue;

    auto neighbors = hash.query_neighbors(ai.x);
    for (Int j : neighbors) {
      if (j <= i) continue;
      Agent& aj = agents[j];
      if (aj.state == PhenoState::DEAD) continue;

      Vec3 delta = sim_.domain().min_image_delta(ai.x, aj.x);
      Real d2 = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
      Real sum_r = ai.radius + aj.radius;

      if (d2 <= 0.0) continue;

      Real d = std::sqrt(d2);

      // Unit normal from i to j
      Vec3 n = {delta[0]/d, delta[1]/d, delta[2]/d};

      // --- Repulsive force (overlap > 0) ---
      if (Real overlap = sum_r - d; overlap > 0.0) {
        Real force_mag;
        if (cfg_.hertzian_enabled) {
          // Hertzian contact: F = k * overlap^(3/2)
          force_mag = cfg_.hertz_k * std::pow(overlap, 1.5);
        } else {
          // Fallback linear spring (legacy behavior)
          force_mag = cfg_.hertz_k * overlap;
        }

        // Apply equal and opposite displacement impulse
        // Using reduced mass for proper two-body dynamics
        Real mi = std::max(ai.mass, 1.0e-30);
        Real mj = std::max(aj.mass, 1.0e-30);
        Real inv_mi = 1.0 / mi;
        Real inv_mj = 1.0 / mj;
        Real inv_sum = 1.0 / (inv_mi + inv_mj);

        Real push_i = force_mag * dt * inv_mi * inv_mi * inv_sum;
        Real push_j = force_mag * dt * inv_mj * inv_mj * inv_sum;

        // Push i away (opposite to n), push j along n
        ai.x[0] -= n[0] * push_i;
        ai.x[1] -= n[1] * push_i;
        ai.x[2] -= n[2] * push_i;
        aj.x[0] += n[0] * push_j;
        aj.x[1] += n[1] * push_j;
        aj.x[2] += n[2] * push_j;
      }

      // --- Adhesion (EPS-mediated, attractive at short gaps) ---
      if (cfg_.adhesion_enabled) {
        Real gap = d - sum_r;  // positive when cells not overlapping
        if (gap > 0.0 && gap < cfg_.adhesion_range) {
          // Linear decay of adhesion with distance
          Real adhesion_frac = 1.0 - (gap / cfg_.adhesion_range);
          Real adhesion_force = cfg_.adhesion_strength * adhesion_frac;

          Real mi = std::max(ai.mass, 1.0e-30);
          Real mj = std::max(aj.mass, 1.0e-30);
          Real inv_mi = 1.0 / mi;
          Real inv_mj = 1.0 / mj;
          Real inv_sum = 1.0 / (inv_mi + inv_mj);

          Real pull_i = adhesion_force * dt * inv_mi * inv_mi * inv_sum;
          Real pull_j = adhesion_force * dt * inv_mj * inv_mj * inv_sum;

          // Pull i toward j (along n), pull j toward i (opposite n)
          ai.x[0] += n[0] * pull_i;
          ai.x[1] += n[1] * pull_i;
          ai.x[2] += n[2] * pull_i;
          aj.x[0] -= n[0] * pull_j;
          aj.x[1] -= n[1] * pull_j;
          aj.x[2] -= n[2] * pull_j;
        }
      }
    }
  }
}

}  // namespace gutibm
