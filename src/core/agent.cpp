/* -----------------------------------------------------------------------
   GutIBM – Agent implementation
   ----------------------------------------------------------------------- */

#include "agent.h"
#include <algorithm>

namespace gutibm {

TagID AgentPool::first_tag_for_rank(Int rank, [[maybe_unused]] Int nprocs) {
  return static_cast<TagID>(rank) + 1;
}

TagID AgentPool::tag_stride(Int nprocs) {
  return static_cast<TagID>(std::max(1, nprocs));
}

TagID AgentPool::next_tag_after_max(TagID max_seen, Int rank, Int nprocs) {
  const TagID stride = tag_stride(nprocs);
  TagID candidate = std::max(max_seen + 1, first_tag_for_rank(rank, nprocs));
  const auto rank_offset = static_cast<TagID>(rank);
  if (const auto offset = (candidate - 1) % stride; offset != rank_offset) {
    candidate += (rank_offset + stride - offset) % stride;
  }
  return candidate;
}

void AgentPool::configure_tags(TagID first_tag, TagID stride) {
  next_tag_ = first_tag;
  tag_stride_ = std::max<TagID>(1, stride);
}

Agent Agent::create_default(TagID id, Int type, Vec3 pos, Real mu_max_val) {
  Agent a{};
  a.identity.tag        = id;
  a.identity.type       = type;
  a.identity.owner_rank = 0;
  a.x             = pos;
  a.v             = {0.0, 0.0, 0.0};
  a.radius        = CELL_RADIUS_DEFAULT;
  a.mass          = sphere_mass(CELL_RADIUS_DEFAULT, CELL_DENSITY_DEFAULT);
  a.outer_radius  = CELL_RADIUS_DEFAULT * 1.05;
  a.mu_max        = mu_max_val;
  a.mu_realized   = mu_max_val;
  a.biomass       = a.mass;
  a.maintenance   = 0.0;

  a.receptor_expr_base.fill(1.0);
  a.receptor_expr = a.receptor_expr_base;

  a.km.km_iron   = 1.0e-6;
  a.km.km_b12    = 1.0e-9;
  a.km.km_carbon = 5.0e-3;

  a.state = PhenoState::NORMAL;

  a.genome.lineage_id   = id;
  a.genome.parent_id    = 0;
  a.genome.generation   = 0;
  a.genome.has_conjugative_plasmid = false;
  a.genome.mutations    = 0;
  a.genome.plasmid_cost_amelioration = 0.0;
  a.genome.receptor_expression.fill(1.0);
  a.genome.toxin_affinity.fill(1.0);
  a.genome.ligand_affinity.fill(1.0);

  a.timers.age       = 0.0;
  a.timers.sos_timer = -1.0;
  a.flags.in_crypt   = false;
  a.grid_cell = -1;

  return a;
}

void AgentPool::remove(Int idx) {
  if (idx < 0 || idx >= size()) return;
  if (idx < size() - 1) {
    agents_[idx] = std::move(agents_.back());
  }
  agents_.pop_back();
}

}  // namespace gutibm
