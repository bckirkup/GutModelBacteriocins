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

struct AgentIdentity {
  TagID tag = 0;
  Int type = 0;
  Int owner_rank = 0;
};

struct AgentKm {
  Real km_iron = 0;
  Real km_b12 = 0;
  Real km_carbon = 0;
};

struct AgentTimers {
  Real age = 0;
  Real sos_timer = 0;
  Real death_time = -1.0;
};

struct AgentFlags {
  bool in_crypt = false;
  bool just_divided = false;
  bool microcin_penalty_applied = false;
};

struct Agent {
  // ── Identity ──────────────────────────────────────────────────────────
  AgentIdentity identity;

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
  std::array<Real, NUM_RECEPTORS> receptor_expr_base;
  std::array<Real, NUM_RECEPTORS> receptor_expr;

  AgentKm km;

  // ── Phenotype ─────────────────────────────────────────────────────────
  PhenoState state;

  // ── Genome ────────────────────────────────────────────────────────────
  Genome genome;

  AgentTimers timers;

  // ── Motility (Spec 3) ───────────────────────────────────────────────
  struct MotilityState {
    Vec3 swim_direction = {1.0, 0.0, 0.0};
    Vec3 step_displacement = {0.0, 0.0, 0.0};
    Real swim_speed = 0.0;
    Real run_timer = 0.0;
    bool is_stopped = false;
    Real stop_timer = 0.0;
    Real prev_carbon = 0.0;
    Real prev_oxygen = 0.0;
  };
  MotilityState motility;

  AgentFlags flags;

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

  TagID next_tag() {
    const TagID id = next_tag_;
    next_tag_ += tag_stride_;
    return id;
  }
  void set_next_tag(TagID t) { next_tag_ = t; }

  // Globally unique tags: id = rank + 1 + k * nprocs (k = 0, 1, 2, ...)
  static TagID first_tag_for_rank(Int rank, Int nprocs);
  static TagID tag_stride(Int nprocs);
  static TagID next_tag_after_max(TagID max_seen, Int rank, Int nprocs);
  void configure_tags(TagID first_tag, TagID stride);

 private:
  std::vector<Agent> agents_;
  TagID next_tag_ = 1;
  TagID tag_stride_ = 1;
};

}  // namespace gutibm

#endif  // GUTIBM_AGENT_H
