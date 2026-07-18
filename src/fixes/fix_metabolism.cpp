/* -----------------------------------------------------------------------
   GutIBM – Metabolism fix implementation
   ----------------------------------------------------------------------- */

#include "fix_metabolism.h"
#include "species_names.h"
#include "simulation.h"
#include "receptor_utils.h"
#include <cmath>
#include <algorithm>
#ifdef GUTIBM_OPENMP
#include <omp.h>
#include <utility>
#endif

namespace gutibm {

FixMetabolism::FixMetabolism(Simulation& sim, const MetabolismConfig& cfg)
    : Fix("metabolism", sim), cfg_(cfg) {}

void FixMetabolism::init() { /* no-op: parameters set via cfg_ at construction */ }

namespace {

void apply_fur_receptor_expr(Simulation& sim) {
  const auto& fur_cfg = sim.config().cell_bio.fur;
  if (!fur_cfg.enabled) return;

  auto& agents = sim.agents();
  auto& chem = sim.chemical_field();
  Int i_iron = chem.find(species::IRON);

  for (Agent& agent : agents) {
    if (agent.state == PhenoState::DEAD) continue;
    Int cell = agent.grid_cell;
    if (cell < 0) continue;

    Real S_iron = (i_iron >= 0) ? chem.conc(i_iron, cell) : 0.0;
    const Real fur_factor = 1.0 + fur_cfg.upregulation_max * fur_cfg.Km
        / (fur_cfg.Km + S_iron);
    for (int r = 0; r < NUM_RECEPTORS; ++r) {
      if (!is_iron_receptor(r)) {
        agent.receptor_expr[r] = agent.receptor_expr_base[r];
        continue;
      }
      agent.receptor_expr[r] = std::min(
          agent.receptor_expr_base[r] * fur_factor, fur_cfg.receptor_max);
    }
  }
}

bool try_gpu_metabolism(Simulation& sim, const MetabolismConfig& cfg, Real dt) {
  if (!sim.gpu_active()) return false;

  apply_fur_receptor_expr(sim);

  auto& agents = sim.agents();
  auto& ag = sim.agents_gpu();
  auto& cg = sim.chem_gpu();
  ag.sync_from_host(agents);
  const auto& chem = sim.chemical_field();
  Int i_carbon = chem.find(species::CARBON);
  Int i_iron = chem.find(species::IRON);
  Int i_b12 = chem.find(species::B12);
  Int i_acetate = chem.find(species::ACETATE);
  Int i_eut = chem.find(species::ETHANOLAMINE);
  Int i_o2 = chem.find(species::OXYGEN);
  const auto& o2cfg = sim.config().chem_env.oxygen;
  if (!ag.run_metabolism(
          sim.domain(), cfg,
          {
            i_carbon >= 0 ? cg.conc_device(i_carbon) : nullptr,
            i_iron >= 0 ? cg.conc_device(i_iron) : nullptr,
            i_b12 >= 0 ? cg.conc_device(i_b12) : nullptr,
            i_acetate >= 0 ? cg.conc_device(i_acetate) : nullptr,
            i_eut >= 0 ? cg.conc_device(i_eut) : nullptr,
            o2cfg.enabled && i_o2 >= 0 ? cg.conc_device(i_o2) : nullptr,
            i_carbon >= 0 ? cg.reac_device(i_carbon) : nullptr,
            i_iron >= 0 ? cg.reac_device(i_iron) : nullptr,
            i_b12 >= 0 ? cg.reac_device(i_b12) : nullptr,
            o2cfg.enabled ? 1 : 0,
            o2cfg.boost_max,
            o2cfg.Km,
          },
          dt)) {
    return false;
  }
  ag.sync_to_host(agents);
  return true;
}

}  // namespace

void FixMetabolism::compute(Real dt) {
  if (try_gpu_metabolism(sim_, cfg_, dt)) {
    // Division stays in compute (same as CPU path) so fix_bacteriocin in this
    // biology pass can observe just_divided during the division timestep.
    perform_divisions();
    return;
  }

  auto& agents = sim_.agents();
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Agent& a : agents) {
    if (a.state == PhenoState::DEAD) continue;
    compute_growth_rate(a);
    grow_agent(a, dt);
  }

  // Division must run in compute (not post_step) so fix_bacteriocin in the same
  // biology pass can observe just_divided during the division timestep.
  perform_divisions();
}

void FixMetabolism::post_step(Real /*dt*/) {
  auto& agents = sim_.agents();

  for (Int i = agents.size() - 1; i >= 0; --i) {
    Agent& a = agents[i];
    if (a.state == PhenoState::DEAD) continue;
    check_death(a);
  }
}

void FixMetabolism::perform_divisions() {
  auto& agents = sim_.agents();
  std::vector<Agent> new_agents;

  for (Agent& a : agents) {
    if (a.state == PhenoState::DEAD) continue;

    Real initial_mass = sphere_mass(CELL_RADIUS_DEFAULT, CELL_DENSITY_DEFAULT);
    if (a.biomass >= cfg_.division_threshold * initial_mass) {
      // Create daughter cell
      Agent daughter = a;
      daughter.identity.tag = agents.next_tag();
      daughter.genome.parent_id = a.identity.tag;
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

      daughter.timers.age = 0.0;
      a.timers.age = 0.0;
      daughter.receptor_expr_base = a.receptor_expr_base;
      daughter.genome.receptor_expression = a.genome.receptor_expression;
      daughter.motility = a.motility;

      a.flags.just_divided = true;
      daughter.flags.just_divided = true;

      sim_.step_events().divisions++;
      new_agents.push_back(std::move(daughter));
    }
  }

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
  Int i_carbon = chem.find(species::CARBON);
  Int i_iron   = chem.find(species::IRON);
  Int i_b12    = chem.find(species::B12);

  Real S_carbon = (i_carbon >= 0) ? chem.conc(i_carbon, cell) : 1.0;
  Real S_iron   = (i_iron >= 0)   ? chem.conc(i_iron, cell)   : 1.0;
  Real S_b12    = (i_b12 >= 0)    ? chem.conc(i_b12, cell)    : 1.0;

  if (const auto& fur_cfg = sim_.config().cell_bio.fur; fur_cfg.enabled) {
    const Real fur_factor = 1.0 + fur_cfg.upregulation_max * fur_cfg.Km
        / (fur_cfg.Km + S_iron);
    for (int r = 0; r < NUM_RECEPTORS; ++r) {
      if (!is_iron_receptor(r)) {
        agent.receptor_expr[r] = agent.receptor_expr_base[r];
        continue;
      }
      agent.receptor_expr[r] = std::min(
          agent.receptor_expr_base[r] * fur_factor, fur_cfg.receptor_max);
    }
  }

  // Receptor-modified Km values
  // When receptor expression drops, effective Km increases (worse affinity)
  // Receptor expressions for graded iron uptake + partial resistance
  int ri_fepA = to_underlying(ReceptorType::FepA);
  int ri_btuB = to_underlying(ReceptorType::BtuB);
  int ri_iroN = to_underlying(ReceptorType::IroN);
  int ri_iutA = to_underlying(ReceptorType::IutA);
  int ri_fiu  = to_underlying(ReceptorType::Fiu);
  Real expr_fepA = agent.receptor_expr[ri_fepA];
  Real expr_iroN = agent.receptor_expr[ri_iroN];
  Real expr_iutA = agent.receptor_expr[ri_iutA];
  Real expr_fiu  = agent.receptor_expr[ri_fiu];
  Real expr_btuB = agent.receptor_expr[ri_btuB];

  // Prevent division by zero
  expr_btuB = std::max(expr_btuB, 0.01);

  // Partial resistance ligand_affinity modulates effective Km
  Real lig_aff_fepA = std::max(agent.genome.ligand_affinity[ri_fepA], 0.01);
  Real lig_aff_btuB = std::max(agent.genome.ligand_affinity[ri_btuB], 0.01);
  Real lig_aff_iroN = std::max(agent.genome.ligand_affinity[ri_iroN], 0.01);
  Real lig_aff_iutA = std::max(agent.genome.ligand_affinity[ri_iutA], 0.01);
  Real lig_aff_fiu  = std::max(agent.genome.ligand_affinity[ri_fiu],  0.01);

  // Graded iron uptake: each receptor contributes proportional to expression,
  // ligand affinity, and its own Km
  Real iron_uptake = 0.0;
  iron_uptake += expr_fepA * lig_aff_fepA * S_iron / (cfg_.km_iron_primary + S_iron);
  iron_uptake += expr_iroN * lig_aff_iroN * S_iron / (cfg_.km_iron_iroN + S_iron);
  iron_uptake += expr_iutA * lig_aff_iutA * S_iron / (cfg_.km_iron_iutA + S_iron);
  iron_uptake += expr_fiu  * lig_aff_fiu  * S_iron / (cfg_.km_iron_fiu  + S_iron);
  Real monod_iron = iron_uptake / (1.0 + expr_iroN + expr_iutA + expr_fiu);

  Real Km_b12  = agent.km.km_b12  / (expr_btuB * lig_aff_btuB);
  Real Km_carb = agent.km.km_carbon;

  // Triple Monod kinetics (uncoupled)
  Real monod_carbon = S_carbon / (Km_carb + S_carbon);
  Real monod_b12    = S_b12    / (Km_b12  + S_b12);

  Real mu = agent.mu_max * monod_carbon * monod_iron * monod_b12;

  if (const auto& o2cfg = sim_.config().chem_env.oxygen; o2cfg.enabled) {
    if (Int i_o2 = chem.find(species::OXYGEN); i_o2 >= 0) {
      const Real s_o2 = chem.conc(i_o2, cell);
      const Real monod_o2_boost =
          1.0 + o2cfg.boost_max * s_o2 / (o2cfg.Km + s_o2);
      mu *= monod_o2_boost;
    }
  }

  // Metabolic penalties for receptor downregulation
  // BtuB loss → MetE pathway required (proteome cost)
  // Acetate inhibits MetE, scaling the penalty with local [acetate]
  // + concentration-dependent ethanolamine utilization loss
  if (expr_btuB < 0.5) {
    Real metE_eff = cfg_.metE_penalty;
    if (Int i_acetate = chem.find(species::ACETATE); i_acetate >= 0) {
      Real acetate_conc = chem.conc(i_acetate, cell);
      Real acetate_factor = 1.0
          + (cfg_.metE_acetate_max_factor - 1.0)
            * acetate_conc / (cfg_.metE_acetate_km + acetate_conc);
      metE_eff *= acetate_factor;
    }
    Real eut_conc = 0.0;
    if (Int i_eut = chem.find(species::ETHANOLAMINE); i_eut >= 0)
      eut_conc = chem.conc(i_eut, cell);
    Real eut_effect = cfg_.eut_max_penalty * eut_conc / (cfg_.eut_km + eut_conc);
    mu *= (1.0 - metE_eff - eut_effect);
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
  agent.timers.age   += dt;

  // Nutrient consumption and siderophore coupling from grid
  auto& chem = sim_.chemical_field();
  Int cell = agent.grid_cell;
  if (cell < 0) return;

  Int i_carbon = chem.find(species::CARBON);
  Int i_iron   = chem.find(species::IRON);
  Int i_acetate = chem.find(species::ACETATE);

  Real cell_vol = sim_.domain().dx() * sim_.domain().dx() * sim_.domain().dx();

  const auto& sid_cfg = sim_.config().chem_env.siderophore;
  if (sid_cfg.enabled && cell_vol > 0.0) {
    Int i_sid = chem.find(species::SIDEROPHORE);
    if (i_sid >= 0) {
      const Real s_iron = (i_iron >= 0) ? chem.conc(i_iron, cell) : 0.0;
      Real fur_Km = 1.0e-6;
      if (sim_.config().cell_bio.fur.enabled) {
        fur_Km = sim_.config().cell_bio.fur.Km;
      }
      const Real fur_activity = 1.0 - s_iron / (fur_Km + s_iron);
      const Real sid_rate = sid_cfg.secretion_rate * std::max(0.0, fur_activity)
          * agent.biomass / cell_vol;
      chem.reac(i_sid, cell) += sid_rate;

      const Real s_sid = chem.conc(i_sid, cell);
      const Real chelation = sid_cfg.chelation_rate * s_sid * s_iron;
      if (i_iron >= 0) {
        chem.reac(i_iron, cell) -= chelation;
      }
      chem.reac(i_sid, cell) -= chelation;
      if (i_iron >= 0) {
        chem.reac(i_iron, cell) += sid_rate * sid_cfg.recapture_fraction;
      }

      const Real expr_fepA = agent.receptor_expr[to_underlying(ReceptorType::FepA)];
      if (i_iron >= 0 && expr_fepA > 0.0) {
        const Real reimport = expr_fepA * s_sid / (sid_cfg.Km_reimport + s_sid)
            * agent.biomass / cell_vol;
        chem.reac(i_iron, cell) += reimport;
      }
    }
  }

  if (d_biomass <= 0.0 || dt <= 0.0) return;

  if (i_carbon >= 0 && cell_vol > 0.0) {
    Real delta_c = d_biomass * cfg_.yield_carbon / (cell_vol * dt);
    #ifdef GUTIBM_OPENMP
    #pragma omp atomic
    #endif
    chem.reac(i_carbon, cell) -= delta_c;
  }
  if (i_iron >= 0 && cell_vol > 0.0) {
    Real delta_fe = d_biomass * cfg_.yield_iron / (cell_vol * dt);
    #ifdef GUTIBM_OPENMP
    #pragma omp atomic
    #endif
    chem.reac(i_iron, cell) -= delta_fe;
  }
  // Spec 6 §3 — the B12/corrinoid field is NOT depleted. It represents the
  // total bioavailable corrinoid pool (~1 uM), produced by the anaerobic
  // majority at rates far exceeding E. coli demand; the pool is effectively
  // constant relative to the modeled population. (cfg_.yield_b12 retained for
  // config compatibility but no longer removes corrinoid from the field.)

  const auto& acfg = sim_.config().chem_env.acetate;
  if (acfg.enabled && i_acetate >= 0 && cell_vol > 0.0) {
    const Real acetate_conc = chem.conc(i_acetate, cell);
    if (agent.mu_realized > acfg.overflow_threshold) {
      #ifdef GUTIBM_OPENMP
      #pragma omp atomic
      #endif
      chem.reac(i_acetate, cell) +=
          acfg.overflow_rate * agent.biomass / cell_vol;
    }
    const Real scavenge = acfg.scavenge_rate * acetate_conc
        / (acfg.scavenge_Km + acetate_conc) * agent.biomass / cell_vol;
    #ifdef GUTIBM_OPENMP
    #pragma omp atomic
    #endif
    chem.reac(i_acetate, cell) -= scavenge;
  }
}

void FixMetabolism::check_death(Agent& agent) {
  if (agent.mu_realized < cfg_.death_threshold && agent.timers.age > 3600.0) {
    agent.state = PhenoState::DEAD;
    sim_.step_events().starvation_deaths++;
  }
}

}  // namespace gutibm
