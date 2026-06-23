/* -----------------------------------------------------------------------
   GutIBM – Agent implementation
   ----------------------------------------------------------------------- */

#include "agent.h"
#include <algorithm>

namespace gutibm {

Agent Agent::create_default(TagID id, Int type, Vec3 pos, Real mu_max_val) {
  Agent a{};
  a.tag           = id;
  a.type          = type;
  a.owner_rank    = 0;
  a.x             = pos;
  a.v             = {0.0, 0.0, 0.0};
  a.radius        = CELL_RADIUS_DEFAULT;
  a.mass          = sphere_mass(CELL_RADIUS_DEFAULT, CELL_DENSITY_DEFAULT);
  a.outer_radius  = CELL_RADIUS_DEFAULT * 1.05;
  a.mu_max        = mu_max_val;
  a.mu_realized   = mu_max_val;
  a.biomass       = a.mass;
  a.maintenance   = 0.0;

  a.receptor_expr.fill(1.0);  // wild-type: all receptors fully expressed

  a.km_iron   = 1.0e-6;   // 1 uM baseline
  a.km_b12    = 1.0e-9;   // 1 nM baseline
  a.km_carbon = 5.0e-3;   // 5 mM baseline

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

  a.age       = 0.0;
  a.sos_timer = -1.0;
  a.in_crypt  = false;
  a.grid_cell = -1;

  return a;
}

void AgentPool::remove(Int idx) {
  if (idx < 0 || idx >= size()) return;
  // swap with last element for O(1) removal
  if (idx < size() - 1) {
    agents_[idx] = std::move(agents_.back());
  }
  agents_.pop_back();
}

}  // namespace gutibm
