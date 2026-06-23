/* -----------------------------------------------------------------------
   GutIBM – Bacterial agent data structure
   Each agent represents one discrete Enterobacteriaceae cell.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_AGENT_H
#define GUTIBM_AGENT_H

#include "types.h"
#include <cstdint>
#include <vector>

namespace gutibm {

struct Agent {
  // ── Identity ──────────────────────────────────────────────────────────
  TagID  tag;              // globally unique identifier
  Int    type;             // species / phylogroup index (1-based)
  Int    owner_rank;       // MPI rank owning this agent

  // ── Spatial ───────────────────────────────────────────────────────────
  Vec3   x;                // position (m)
  Vec3   v;                // velocity (m/s), set by advection + mechanics
  Real   radius;           // cell radius (m)
  Real   mass;             // cell mass (kg)
  Real   outer_radius;     // EPS/capsule outer radius

  // ── Metabolism ────────────────────────────────────────────────────────
  Real   mu_max;           // intrinsic max specific growth rate (1/s)
  Real   mu_realized;      // current realized growth rate after penalties
  Real   biomass;          // dry biomass (kg)
  Real   maintenance;      // maintenance energy coefficient

  // ── Receptor state ────────────────────────────────────────────────────
  //  expression level 0.0 (fully downregulated) → 1.0 (wild-type)
  std::array<Real, NUM_RECEPTORS> receptor_expr;

  // Km values (modified by receptor expression)
  Real   km_iron;          // iron half-saturation (mol/m^3)
  Real   km_b12;           // B12 half-saturation (mol/m^3)
  Real   km_carbon;        // carbon half-saturation (mol/m^3)

  // ── Phenotype ─────────────────────────────────────────────────────────
  PhenoState state;

  // ── Genome ────────────────────────────────────────────────────────────
  Genome genome;

  // ── Stochastic timers ─────────────────────────────────────────────────
  Real   age;              // time since last division (s)
  Real   sos_timer;        // SOS response countdown (s), <0 = inactive

  // ── Crypt state ──────────────────────────────────────────────────────
  bool   in_crypt = false; // true when agent resides in a crypt refugium

  // ── Grid coupling ─────────────────────────────────────────────────────
  Int    grid_cell;        // index into the chemical field grid

  // ── Factory ───────────────────────────────────────────────────────────
  static Agent create_default(TagID id, Int type, Vec3 pos, Real mu_max_val);
};

// Agent population container
class AgentPool {
 public:
  AgentPool() = default;

  Int  size() const { return static_cast<Int>(agents_.size()); }
  void reserve(Int n) { agents_.reserve(n); }
  void push_back(Agent a) { agents_.push_back(std::move(a)); }
  void remove(Int idx);

  Agent&       operator[](Int i)       { return agents_[i]; }
  const Agent& operator[](Int i) const { return agents_[i]; }

  auto begin()        { return agents_.begin(); }
  auto end()          { return agents_.end(); }
  auto begin()  const { return agents_.begin(); }
  auto end()    const { return agents_.end(); }

  TagID next_tag() { return next_tag_++; }
  void  set_next_tag(TagID t) { next_tag_ = t; }

 private:
  std::vector<Agent> agents_;
  TagID next_tag_ = 1;
};

}  // namespace gutibm

#endif  // GUTIBM_AGENT_H
