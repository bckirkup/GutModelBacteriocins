/* -----------------------------------------------------------------------
   GutIBM – Metabolism fix implementation
   ----------------------------------------------------------------------- */

#include "fix_metabolism.h"
#include "simulation.h"
#include <cmath>
#include <algorithm>

namespace gutibm {

FixMetabolism::FixMetabolism(Simulation& sim, const MetabolismConfig& cfg)
    : Fix("metabolism", sim), cfg_(cfg) {}

void FixMetabolism::init() {}

void FixMetabolism::compute(Real dt) {
  auto& agents = sim_.agents();
  for (Int i = 0; i < agents.size(); ++i) {
    Agent& a = agents[i];
    if (a.state == PhenoState::DEAD) continue;
    compute_growth_rate(a);
    grow_agent(a, dt);
  }
}

void FixMetabolism::post_step(Real dt) {
  // Division and death checks after growth
  auto& agents = sim_.agents();
  std::vector<Agent> new_agents;

  for (Int i = agents.size() - 1; i >= 0; --i) {
    Agent& a = agents[i];
    if (a.state == PhenoState::DEAD) continue;
    check_death(a);
  }

  // Division pass (separate to avoid invalidating indices)
  Int n = agents.size();
  for (Int i = 0; i < n; ++i) {
    Agent& a = agents[i];
    if (a.state == PhenoState::DEAD) continue;

    Real initial_mass = sphere_mass(CELL_RADIUS_DEFAULT, CELL_DENSITY_DEFAULT);
    if (a.biomass >= cfg_.division_threshold * initial_mass) {
      // Create daughter cell
      Agent daughter = a;
      daughter.tag = agents.next_tag();
      daughter.genome.parent_id = a.tag;
      daughter.genome.generation = a.genome.generation + 1;
      daughter.genome.lineage_id = a.genome.lineage_id;

      // Split biomass equally
      a.biomass *= 0.5;
      daughter.biomass = a.biomass;

      // Update radii from biomass
      Real vol = a.biomass / CELL_DENSITY_DEFAULT;
      a.radius = std::cbrt(3.0 * vol / (4.0 * PI));
      a.mass   = a.biomass;
      daughter.radius = a.radius;
      daughter.mass   = a.biomass;

      // Offset daughter position by one cell diameter
      auto& rng = sim_.rng();
      Real theta = rng.uniform(0.0, 2.0 * PI);
      Real phi   = rng.uniform(0.0, PI);
      Real offset = 2.0 * a.radius;
      daughter.x[0] = a.x[0] + offset * std::sin(phi) * std::cos(theta);
      daughter.x[1] = a.x[1] + offset * std::sin(phi) * std::sin(theta);
      daughter.x[2] = a.x[2] + offset * std::cos(phi);

      sim_.domain().apply_pbc(daughter.x);

      daughter.age = 0.0;
      a.age = 0.0;

      new_agents.push_back(std::move(daughter));
    }
  }

  // Add new agents
  for (auto& na : new_agents) {
    agents.push_back(std::move(na));
  }
}

void FixMetabolism::compute_growth_rate(Agent& agent) {
  auto& chem = sim_.chemical_field();
  Int cell = agent.grid_cell;
  if (cell < 0) {
    agent.mu_realized = 0.0;
    return;
  }

  // Get local concentrations
  Int i_carbon = chem.find("carbon");
  Int i_iron   = chem.find("iron");
  Int i_b12    = chem.find("b12");

  Real S_carbon = (i_carbon >= 0) ? chem.conc(i_carbon, cell) : 1.0;
  Real S_iron   = (i_iron >= 0)   ? chem.conc(i_iron, cell)   : 1.0;
  Real S_b12    = (i_b12 >= 0)    ? chem.conc(i_b12, cell)    : 1.0;

  // Receptor-modified Km values
  // When receptor expression drops, effective Km increases (worse affinity)
  Real expr_fepA = agent.receptor_expr[static_cast<int>(ReceptorType::FepA)];
  Real expr_iroN = agent.receptor_expr[static_cast<int>(ReceptorType::IroN)];
  Real expr_iutA = agent.receptor_expr[static_cast<int>(ReceptorType::IutA)];
  Real expr_fiu  = agent.receptor_expr[static_cast<int>(ReceptorType::Fiu)];
  Real expr_btuB = agent.receptor_expr[static_cast<int>(ReceptorType::BtuB)];

  // Prevent division by zero
  expr_btuB = std::max(expr_btuB, 0.01);

  // Graded iron uptake: each receptor contributes proportional to expression and affinity
  Real iron_uptake = 0.0;
  iron_uptake += expr_fepA * S_iron / (cfg_.km_iron_primary + S_iron);
  iron_uptake += expr_iroN * S_iron / (cfg_.km_iron_iroN + S_iron);
  iron_uptake += expr_iutA * S_iron / (cfg_.km_iron_iutA + S_iron);
  iron_uptake += expr_fiu  * S_iron / (cfg_.km_iron_fiu  + S_iron);
  // Normalize: wild-type with all receptors at 1.0 should give ~same as before
  Real monod_iron = iron_uptake / (1.0 + expr_iroN + expr_iutA + expr_fiu);

  Real Km_b12  = agent.km_b12  / expr_btuB;
  Real Km_carb = agent.km_carbon;

  // Triple Monod kinetics (uncoupled)
  Real monod_carbon = S_carbon / (Km_carb + S_carbon);
  Real monod_b12    = S_b12    / (Km_b12  + S_b12);

  Real mu = agent.mu_max * monod_carbon * monod_iron * monod_b12;

  // Metabolic penalties for receptor downregulation
  // BtuB loss → MetE pathway required (proteome cost + ethanolamine loss)
  // Acetate inhibits MetE, scaling the penalty with local [acetate]
  if (expr_btuB < 0.5) {
    Real metE_eff = cfg_.metE_penalty;
    Int i_acetate = chem.find("acetate");
    if (i_acetate >= 0) {
      Real acetate_conc = chem.conc(i_acetate, cell);
      Real acetate_factor = 1.0
          + (cfg_.metE_acetate_max_factor - 1.0)
            * acetate_conc / (cfg_.metE_acetate_km + acetate_conc);
      metE_eff *= acetate_factor;
    }
    mu *= (1.0 - metE_eff - cfg_.eut_penalty);
  }

  // Plasmid maintenance cost (reduced by compensatory mutations per VADI §79)
  if (!agent.genome.bi_loci.empty()) {
    Real per_locus = std::max(0.0, 0.02 - agent.genome.plasmid_cost_amelioration);
    Real plasmid_cost = per_locus * agent.genome.bi_loci.size();
    plasmid_cost = std::min(plasmid_cost, 0.10);  // cap at 10%
    mu *= (1.0 - plasmid_cost);
  }

  // Subtract maintenance
  mu -= cfg_.maintenance_rate;

  agent.mu_realized = mu;
}

void FixMetabolism::grow_agent(Agent& agent, Real dt) {
  // Biomass increase
  Real d_biomass = agent.mu_realized * agent.biomass * dt;
  agent.biomass += d_biomass;
  agent.biomass = std::max(agent.biomass, 1.0e-20);

  // Update radius from biomass
  Real vol = agent.biomass / CELL_DENSITY_DEFAULT;
  agent.radius = std::cbrt(3.0 * vol / (4.0 * PI));
  agent.mass   = agent.biomass;
  agent.age   += dt;

  // Nutrient consumption from grid
  auto& chem = sim_.chemical_field();
  Int cell = agent.grid_cell;
  if (cell < 0 || d_biomass <= 0.0) return;

  Int i_carbon = chem.find("carbon");
  Int i_iron   = chem.find("iron");
  Int i_b12    = chem.find("b12");

  Real cell_vol = sim_.domain().dx() * sim_.domain().dx() * sim_.domain().dx();

  if (i_carbon >= 0 && cell_vol > 0.0) {
    chem.reac(i_carbon, cell) -= d_biomass * cfg_.yield_carbon / cell_vol;
  }
  if (i_iron >= 0 && cell_vol > 0.0) {
    chem.reac(i_iron, cell) -= d_biomass * cfg_.yield_iron / cell_vol;
  }
  if (i_b12 >= 0 && cell_vol > 0.0) {
    chem.reac(i_b12, cell) -= d_biomass * cfg_.yield_b12 / cell_vol;
  }
}

void FixMetabolism::check_death(Agent& agent) {
  if (agent.mu_realized < cfg_.death_threshold && agent.age > 3600.0) {
    agent.state = PhenoState::DEAD;
  }
}

}  // namespace gutibm
