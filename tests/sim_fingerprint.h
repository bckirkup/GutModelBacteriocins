/* -----------------------------------------------------------------------
   GutIBM – Deterministic simulation fingerprint for cross-build parity
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_TEST_SIM_FINGERPRINT_H
#define GUTIBM_TEST_SIM_FINGERPRINT_H

#include "simulation.h"
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>

namespace gutibm {
namespace test_util {

inline uint64_t hash_combine(uint64_t h, uint64_t v) {
  return h ^ (v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
}

inline int64_t quantize(Real v, Real scale = 1e12) {
  return static_cast<int64_t>(std::llround(v * scale));
}

// Deterministic fingerprint of post-run simulation state (agents + grid).
inline uint64_t simulation_fingerprint(const Simulation& sim) {
  uint64_t h = 0;
  h = hash_combine(h, static_cast<uint64_t>(sim.step_count()));
  h = hash_combine(h, static_cast<uint64_t>(quantize(sim.time())));

  const auto& agents = sim.agents();
  h = hash_combine(h, static_cast<uint64_t>(agents.size()));
  for (const Agent& a : agents) {
    if (a.state == PhenoState::DEAD) continue;
    h = hash_combine(h, static_cast<uint64_t>(a.tag));
    h = hash_combine(h, static_cast<uint64_t>(a.type));
    h = hash_combine(h, static_cast<uint64_t>(quantize(a.x[0])));
    h = hash_combine(h, static_cast<uint64_t>(quantize(a.x[1])));
    h = hash_combine(h, static_cast<uint64_t>(quantize(a.x[2])));
    h = hash_combine(h, static_cast<uint64_t>(quantize(a.biomass)));
    h = hash_combine(h, static_cast<uint64_t>(quantize(a.mu_realized)));
  }

  const auto& chem = sim.chemical_field();
  for (const auto& row : chem.conc_data()) {
    Real sum = 0.0;
    for (Real val : row) {
      sum += val;
    }
    h = hash_combine(h, static_cast<uint64_t>(quantize(sum)));
  }

  return h;
}

inline std::string fingerprint_hex(const Simulation& sim) {
  std::ostringstream oss;
  oss << std::hex << std::setw(16) << std::setfill('0')
      << simulation_fingerprint(sim);
  return oss.str();
}

}  // namespace test_util
}  // namespace gutibm

#endif  // GUTIBM_TEST_SIM_FINGERPRINT_H
