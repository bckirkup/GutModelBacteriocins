/* -----------------------------------------------------------------------
   GutIBM – Metabolism fix implementation
   ----------------------------------------------------------------------- */

#include "fix_metabolism.h"
#include "species_names.h"
#include "simulation.h"
#include "receptor_utils.h"
#include "carbon_maintenance.h"
#include "metabolic_mode.h"
#include "error.h"
#include "delivery_support.h"
#include <cassert>
#include <cmath>
#include <algorithm>
#include <vector>
#ifdef GUTIBM_OPENMP
#include <omp.h>
#include <utility>
#endif

namespace gutibm {

FixMetabolism::FixMetabolism(Simulation& sim, const MetabolismConfig& cfg)
    : Fix("metabolism", sim), cfg_(cfg) {}

void FixMetabolism::init() {
  ensure_delivery_support_stencil();
}

Real FixMetabolism::delivery_concentration(
    const Agent& agent, Int species_index) const {
  const auto& chem = sim_.chemical_field();
  const auto& domain = sim_.domain();
  if (cfg_.delivery_far_field_radius <= 0.0) {
    return chem.total_conc_global(species_index, agent.grid_cell, domain);
  }

  const auto& support = delivery_support_cells(agent);
  Real concentration_sum = 0.0;
  for (const Int cell : support) {
    concentration_sum += chem.total_conc_global(
        species_index, cell, domain);
  }
  return support.empty()
      ? 0.0
      : concentration_sum / static_cast<Real>(support.size());
}

std::vector<Int> FixMetabolism::enumerate_delivery_support_cells(
    const Agent& agent) const {
  const auto& domain = sim_.domain();
  const Real radius = cfg_.delivery_far_field_radius;
  if (radius <= 0.0) return {agent.grid_cell};
  ensure_delivery_support_stencil();
  std::vector<Int> support;
  enumerate_physical_delivery_ball(
      domain, agent.x, delivery_support_stencil_, support);
  if (support.empty() && agent.grid_cell >= 0) {
    support.push_back(agent.grid_cell);
  }
  return support;
}

void FixMetabolism::ensure_delivery_support_stencil() const {
  if (!delivery_support_stencil_.matches(
          sim_.domain(), cfg_.delivery_far_field_radius)) {
    delivery_support_stencil_ = make_delivery_support_stencil(
        sim_.domain(), cfg_.delivery_far_field_radius);
  }
}

const std::vector<Int>& FixMetabolism::delivery_support_cells(
    const Agent& agent) const {
  auto cached = delivery_support_cache_.find(agent.identity.tag);
  if (cached != delivery_support_cache_.end()) {
    return cached->second;
  }
  // The cache is populated before the OpenMP biology loop. Keep an
  // unexpected miss safe if a caller is added to that loop later.
#ifdef GUTIBM_OPENMP
#pragma omp critical(delivery_support_cache_insert)
#endif
  {
    cached = delivery_support_cache_.find(agent.identity.tag);
    if (cached == delivery_support_cache_.end()) {
      cached = delivery_support_cache_.emplace(
          agent.identity.tag, enumerate_delivery_support_cells(agent)).first;
    }
  }
  return cached->second;
}

void FixMetabolism::prepare_delivery_support_cache() {
  delivery_support_cache_.clear();
  if (cfg_.uptake_limit_mode != UptakeLimitMode::Delivery
      || cfg_.delivery_far_field_radius <= 0.0) {
    return;
  }
  delivery_support_cache_.reserve(sim_.agents().size());
  for (const Agent& agent : sim_.agents()) {
    if (agent.state == PhenoState::DEAD || agent.flags.is_ghost
        || agent.grid_cell < 0) {
      continue;
    }
    delivery_support_cache_.emplace(
        agent.identity.tag, enumerate_delivery_support_cells(agent));
  }
}

void FixMetabolism::add_delivery_mass(
    Int species_index, const Agent& agent, Real amount) const {
  if (amount <= 0.0) return;
  auto& chem = sim_.chemical_field();
  const auto& support = delivery_support_cells(agent);
  if (support.empty()) return;
  const Real per_cell = amount / static_cast<Real>(support.size());
  Real deposited = 0.0;
  for (size_t i = 0; i < support.size(); ++i) {
    const Real share = i + 1 == support.size()
        ? amount - deposited : per_cell;
    chem.add_prescribed_sink_global(species_index, support[i], share);
    deposited += share;
  }
}

Real FixMetabolism::delivery_field_funding(
    Int species_index, const Agent& agent, Real amount,
    const std::vector<Real>& requested_by_cell) const {
  if (amount <= 0.0) return 0.0;
  const auto& chem = sim_.chemical_field();
  const auto& support = delivery_support_cells(agent);
  if (support.empty()) return 0.0;
  const Real per_cell = amount / static_cast<Real>(support.size());
  Real deposited = 0.0;
  Real funding = 0.0;
  for (size_t i = 0; i < support.size(); ++i) {
    const Real share = i + 1 == support.size()
        ? amount - deposited : per_cell;
    const Int cell = support[i];
    const Real requested = requested_by_cell[static_cast<size_t>(cell)];
    if (requested > 0.0) {
      funding += chem.prescribed_sink_global(species_index, cell)
          * (share / requested);
    }
    deposited += share;
  }
  return funding;
}

namespace {

struct DeliveryFunding {
  Real maintenance_paid = 0.0;
  Real growth_funded = 0.0;
  Real growth_fraction = 0.0;
};

void apply_delivery_shrinkage(Agent& agent, Real dt) {
  agent.biomass = std::max(agent.biomass + agent.pending_biomass, 1.0e-20);
  const Real volume = agent.biomass / CELL_DENSITY_DEFAULT;
  agent.radius = std::cbrt(3.0 * volume / (4.0 * PI));
  agent.mass = agent.biomass;
  agent.timers.age += dt;
}

DeliveryFunding calculate_delivery_funding(
    Real growth_demand, Real maintenance_demand, Real funded_total) {
  const Real total_demand = growth_demand + maintenance_demand;
  const Real share = std::min(total_demand, std::max(0.0, funded_total));
  DeliveryFunding funding;
  funding.maintenance_paid = std::min(maintenance_demand, share);
  funding.growth_funded = std::min(
      growth_demand,
      std::max(0.0, share - funding.maintenance_paid));
  funding.growth_fraction = growth_demand > 0.0
      ? funding.growth_funded / growth_demand : 0.0;
  return funding;
}

void record_delivery_funding(
    ChemicalField& chem, Int species_index, Real growth_demand,
    Real maint_demand, const DeliveryFunding& funding) {
  chem.flux_accounting().add_maintenance(
      species_index, funding.maintenance_paid);
  if (maint_demand > funding.maintenance_paid) {
    chem.flux_accounting().add_maintenance_shortfall(
        species_index, maint_demand - funding.maintenance_paid);
    chem.flux_accounting().add_maintenance_limited_agents(
        species_index, 1.0);
  }
  chem.flux_accounting().add_agent_uptake(
      species_index, funding.growth_funded);
  if (growth_demand > funding.growth_funded) {
    chem.flux_accounting().add_uptake_shortfall(
        species_index, growth_demand - funding.growth_funded);
    chem.flux_accounting().add_uptake_limited(species_index, 1.0);
  }
}

void apply_delivery_funding(
    Agent& agent, const DeliveryFunding& funding, Real dt) {
  if (agent.mu_realized > 0.0) {
    agent.mu_realized *= funding.growth_fraction;
  }
  agent.biomass += agent.pending_biomass * funding.growth_fraction;
  const Real volume = std::max(agent.biomass, 1.0e-20)
      / CELL_DENSITY_DEFAULT;
  agent.biomass = std::max(agent.biomass, 1.0e-20);
  agent.radius = std::cbrt(3.0 * volume / (4.0 * PI));
  agent.mass = agent.biomass;
  agent.timers.age += dt;
}

Real implicit_ferric_enterobactin_reimport(
    Real ferric_after_production, Real vmax, Real km, Real dt) {
  if (ferric_after_production <= 0.0 || vmax <= 0.0 || dt <= 0.0) {
    return 0.0;
  }
  const Real linear_term = km + dt * vmax - ferric_after_production;
  const Real discriminant = linear_term * linear_term
      + 4.0 * ferric_after_production * km;
  const Real root = std::sqrt(discriminant);
  const Real ferric_after_reimport = linear_term >= 0.0
      ? 2.0 * ferric_after_production * km / (linear_term + root)
      : 0.5 * (-linear_term + root);
  return (ferric_after_production - ferric_after_reimport) / dt;
}

bool try_gpu_metabolism(Simulation& sim, const MetabolismConfig& cfg, Real dt) {
  if (!sim.gpu_active()) return false;
  if (cfg.uptake_limit_mode == UptakeLimitMode::Delivery) return false;

  auto& agents = sim.agents();
  Int local_agent_count = 0;
  bool ghost_seen = false;
  for (const Agent& agent : agents) {
    if (agent.flags.is_ghost) {
      ghost_seen = true;
    } else {
      if (ghost_seen) return false;
      ++local_agent_count;
    }
  }
  if (local_agent_count <= 0) return false;
  auto& ag = sim.agents_gpu();
  auto& cg = sim.chem_gpu();
  cg.reset_agent_uptake();
  cg.reset_uptake_limit_totals();
  ag.sync_from_host(agents);
  const auto& chem = sim.chemical_field();
  Int i_carbon = chem.find(species::CARBON);
  if (cfg.carbon_maintenance_rate > 0.0
      && cfg.uptake_limit_mode == UptakeLimitMode::Voxel) {
    cg.prepare_maintenance(chem, i_carbon, sim.domain().cell_volume());
  }
  Int i_iron = chem.find(species::IRON);
  Int i_b12 = chem.find(species::B12);
  Int i_acetate = chem.find(species::ACETATE);
  Int i_eut = chem.find(species::ETHANOLAMINE);
  Int i_o2 = chem.find(species::OXYGEN);
  const auto& o2cfg = sim.config().chem_env.oxygen;
  const auto& fur_cfg = sim.config().cell_bio.fur;
  const auto& acetate_cfg = sim.config().chem_env.acetate;
  GpuMetabolismBuffers buffers;
  buffers.d_conc_carbon = i_carbon >= 0 ? cg.conc_device(i_carbon) : nullptr;
  buffers.d_conc_iron = i_iron >= 0 ? cg.conc_device(i_iron) : nullptr;
  buffers.d_conc_b12 = i_b12 >= 0 ? cg.conc_device(i_b12) : nullptr;
  buffers.d_conc_acetate =
      i_acetate >= 0 ? cg.conc_device(i_acetate) : nullptr;
  buffers.d_conc_eut = i_eut >= 0 ? cg.conc_device(i_eut) : nullptr;
  buffers.d_conc_oxygen =
      o2cfg.enabled && i_o2 >= 0 ? cg.conc_device(i_o2) : nullptr;
  buffers.d_reac_carbon = i_carbon >= 0 ? cg.reac_device(i_carbon) : nullptr;
  buffers.d_reac_iron = i_iron >= 0 ? cg.reac_device(i_iron) : nullptr;
  buffers.d_reac_b12 = i_b12 >= 0 ? cg.reac_device(i_b12) : nullptr;
  buffers.d_reac_acetate =
      i_acetate >= 0 ? cg.reac_device(i_acetate) : nullptr;
  buffers.iron_uptake_enabled = cfg.iron_uptake_enabled ? 1 : 0;
  buffers.b12_uptake_enabled = cfg.b12_uptake_enabled ? 1 : 0;
  buffers.eut_enabled = cfg.eut_enabled ? 1 : 0;
  buffers.fur_enabled = fur_cfg.enabled ? 1 : 0;
  buffers.fur_Km = fur_cfg.Km;
  buffers.fur_upregulation_max = fur_cfg.upregulation_max;
  buffers.fur_receptor_max = fur_cfg.receptor_max;
  buffers.acetate_enabled = acetate_cfg.enabled ? 1 : 0;
  buffers.acetate_overflow_threshold = acetate_cfg.overflow_threshold;
  buffers.acetate_overflow_rate = acetate_cfg.overflow_rate;
  buffers.acetate_scavenge_rate = acetate_cfg.scavenge_rate;
  buffers.acetate_scavenge_Km = acetate_cfg.scavenge_Km;
  buffers.o2_enabled = o2cfg.enabled ? 1 : 0;
  buffers.o2_boost_max = o2cfg.boost_max;
  buffers.o2_Km = o2cfg.Km;
  buffers.metabolic_switch_enabled = o2cfg.metabolic_switch_enabled ? 1 : 0;
  buffers.mu_crit = o2cfg.mu_crit;
  buffers.aerobic_mu_factor = o2cfg.aerobic_mu_factor;
  buffers.anaerobic_mu_factor = o2cfg.anaerobic_mu_factor;
  buffers.aerobic_carbon_cost_factor =
      o2cfg.aerobic_carbon_cost_factor;
  buffers.anaerobic_carbon_cost_factor =
      o2cfg.anaerobic_carbon_cost_factor;
  buffers.tau_metabolic_switch = o2cfg.tau_metabolic_switch;
  buffers.ferm_acid_yield = o2cfg.ferm_acid_yield;
  buffers.anaerobic_maintenance_factor =
      o2cfg.anaerobic_maintenance_factor;
  buffers.acid_inhibition_enabled =
      cfg.acid_inhibition_enabled ? 1 : 0;
  buffers.acid_inhibition_max = cfg.acid_inhibition_max;
  buffers.acid_inhibition_Ki = cfg.acid_inhibition_Ki;
  buffers.acetate_pKa = cfg.acetate_pKa;
  buffers.environment_pH =
      sim.config().fixes.bacteriocin.mucin_charge.ph;
  buffers.global_nx = sim.domain().nx();
  buffers.global_ny = sim.domain().ny();
  buffers.storage_nx = cg.storage_nx();
  buffers.owned_global_x_begin =
      cg.slab_mode() ? cg.owned_x_begin() : 0;
  buffers.owned_global_x_end =
      cg.slab_mode() ? cg.owned_x_end() : sim.domain().nx();
  buffers.owned_storage_x_begin = cg.owned_storage_x_begin();
  buffers.receptor_count = NUM_RECEPTORS;
  if (i_carbon >= 0) {
    const ChemicalSpec& carbon_spec = chem.spec(i_carbon);
    buffers.effective_diffusivity_carbon = uptake::effective_diffusivity(
        carbon_spec.diff_coeff, carbon_spec.retardation);
  }
  if (i_iron >= 0) {
    const ChemicalSpec& iron_spec = chem.spec(i_iron);
    buffers.effective_diffusivity_iron = uptake::effective_diffusivity(
        iron_spec.diff_coeff, iron_spec.retardation);
  }
  buffers.d_uptake_limit_totals = cg.uptake_limit_totals_device();
  buffers.d_maintenance_available = cg.maintenance_available_device();
  if (!ag.run_metabolism(
          sim.domain(), cfg, buffers, cg.agent_uptake_device(), dt,
          local_agent_count)) {
    return false;
  }
  ag.sync_to_host(agents);
  ag.sync_receptor_expression_to_host(agents);
  cg.accumulate_reactions_to_host(sim.chemical_field());
  cg.download_agent_uptake(sim.chemical_field());
  cg.download_uptake_limit_totals(sim.chemical_field());
  return true;
}

}  // namespace

void FixMetabolism::compute(Real dt) {
  if (try_gpu_metabolism(sim_, cfg_, dt)) {
    // Division stays in compute (same as CPU path) so fix_bacteriocin in this
    // biology pass can observe just_divided during the division timestep.
    apply_siderophore_chemistry(dt);
    perform_divisions();
    return;
  }

  prepare_delivery_support_cache();
  auto& agents = sim_.agents();
  prepare_carbon_maintenance();
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Agent& agent : agents) compute_agent(agent, dt);
  apply_siderophore_chemistry(dt);

  // Division must run in compute (not post_step) so fix_bacteriocin in the same
  // biology pass can observe just_divided during the division timestep.
  perform_divisions();
}

void FixMetabolism::compute_agent(Agent& agent, Real dt) {
  if (agent.state == PhenoState::DEAD) return;
  compute_growth_rate(agent, dt);
  if (cfg_.uptake_limit_mode == UptakeLimitMode::Delivery
      && !agent.flags.is_ghost
      && agent.pending_growth_chemistry_biomass > 0.0) {
    apply_growth_chemistry(agent, agent.pending_growth_chemistry_biomass, dt);
    agent.pending_growth_chemistry_biomass = 0.0;
  }
  agent.pending_growth_carbon = 0.0;
  agent.pending_maintenance_carbon = 0.0;
  agent.pending_carbon_ceiling = 0.0;
  agent.pending_carbon_funding = 0.0;
  agent.pending_biomass = 0.0;
  agent.pending_oxygen_growth = 0.0;
  agent.pending_oxygen_maintenance = 0.0;
  agent.pending_oxygen_ceiling = 0.0;
  agent.pending_oxygen_funding = 0.0;
  if (agent.flags.is_ghost) {
    const Real demanded_biomass = agent.mu_realized * agent.biomass * dt;
    agent.mu_realized *= uptake_limit_fraction(
        agent, demanded_biomass, dt, false);
    return;
  }
  charge_carbon_maintenance(agent, dt);
  if (cfg_.uptake_limit_mode == UptakeLimitMode::Delivery) {
    prepare_delivery_uptake(agent, dt);
    if (sim_.config().chem_env.oxygen.delivery_uptake_enabled) {
      prepare_delivery_oxygen(agent, dt);
    }
    return;
  }
  grow_agent(agent, dt);
}

void FixMetabolism::prepare_delivery_uptake(Agent& agent, Real dt) {
  const Real d_biomass = agent.mu_realized * agent.biomass * dt;
  const Real growth_demand = std::max(d_biomass, 0.0)
      * realized_carbon_cost(agent);
  auto& chem = sim_.chemical_field();
  const Int carbon = chem.find(species::CARBON);
  if (carbon < 0 || agent.grid_cell < 0) return;
  const Real cell_volume = sim_.domain().cell_volume();
  const Real concentration = delivery_concentration(agent, carbon);
  const ChemicalSpec& carbon_spec = chem.spec(carbon);
  const Real effective_diffusivity = uptake::effective_diffusivity(
      carbon_spec.diff_coeff, carbon_spec.retardation);
  const Real demand_total = growth_demand
      + agent.pending_maintenance_carbon;
  const Real ceiling = std::min(
      demand_total,
      uptake::allowed_uptake_mol(
          to_underlying(UptakeLimitMode::Sherwood), concentration,
          effective_diffusivity, agent.radius, cell_volume, dt));
  agent.pending_growth_carbon = growth_demand;
  agent.pending_carbon_funding = std::max(0.0, ceiling);
  agent.pending_carbon_ceiling = std::max(0.0, ceiling);
  agent.pending_biomass = d_biomass;
  chem.flux_accounting().add_uptake_demand(carbon, growth_demand);
  // This zero-order prescribed-mass sink is not the retired legacy sink:
  // legacy uptake was applied and clipped to voxel content before diffusion,
  // whereas this draw is inside implicit diffusion and is capped analytically.
  // It has no conductance or k/3 splitting bias; the solve supplies its
  // diffusive neighbourhood rather than restricting uptake to the agent voxel.
  add_delivery_mass(carbon, agent, ceiling);
}

void FixMetabolism::prepare_delivery_oxygen(Agent& agent, Real dt) {
  auto& chem = sim_.chemical_field();
  const Int oxygen = chem.find(species::OXYGEN);
  if (oxygen < 0 || agent.grid_cell < 0 || dt <= 0.0) return;
  const auto& oxygen_cfg = sim_.config().chem_env.oxygen;
  const Real respiratory_fraction = oxygen_cfg.metabolic_switch_enabled
      ? 1.0 - agent.realized_fermentation_fraction : 1.0;
  const Real growth_demand = oxygen_cfg.q_consumption
      * std::max(agent.mu_realized, 0.0)
      * std::clamp(respiratory_fraction, 0.0, 1.0) * dt;
  const Real maintenance_demand = oxygen_cfg.q_maintenance * dt;
  const Real demand_total = growth_demand + maintenance_demand;
  const Real cell_volume = sim_.domain().cell_volume();
  const Real concentration = delivery_concentration(agent, oxygen);
  const ChemicalSpec& oxygen_spec = chem.spec(oxygen);
  const Real effective_diffusivity = uptake::effective_diffusivity(
      oxygen_spec.diff_coeff, oxygen_spec.retardation);
  const Real ceiling = std::min(
      demand_total,
      uptake::allowed_uptake_mol(
          to_underlying(UptakeLimitMode::Sherwood), concentration,
          effective_diffusivity, agent.radius, cell_volume, dt));
  agent.pending_oxygen_growth = growth_demand;
  agent.pending_oxygen_maintenance = maintenance_demand;
  agent.pending_oxygen_funding = std::max(0.0, ceiling);
  agent.pending_oxygen_ceiling = std::max(0.0, ceiling);
  chem.flux_accounting().add_uptake_demand(oxygen, growth_demand);
  // This zero-order prescribed-mass sink is not the retired legacy sink:
  // it is inside implicit diffusion, not clipped to one voxel, and is capped
  // at analytic diffusive delivery rather than unbounded demand. It has no
  // conductance or k/3 splitting bias; the solve supplies the neighbourhood.
  add_delivery_mass(oxygen, agent, ceiling);
}

void FixMetabolism::commit_delivery_uptake(Real dt) {
  if (cfg_.uptake_limit_mode != UptakeLimitMode::Delivery || dt <= 0.0) {
    return;
  }
  auto& chem = sim_.chemical_field();
  const Int carbon = chem.find(species::CARBON);
  commit_delivery_carbon(dt, carbon);
  const Int oxygen = chem.find(species::OXYGEN);
  if (oxygen >= 0
      && sim_.config().chem_env.oxygen.delivery_uptake_enabled) {
    commit_delivery_oxygen(dt, oxygen);
  }
}

void FixMetabolism::commit_delivery_carbon(Real dt, Int carbon) {
  auto& chem = sim_.chemical_field();
  std::vector carbon_demand_by_cell(
      static_cast<size_t>(chem.global_ncells()), 0.0);
  for (const Agent& agent : sim_.agents()) {
    if (agent.state == PhenoState::DEAD || agent.flags.is_ghost
        || agent.grid_cell < 0) continue;
    if (agent.pending_carbon_funding <= 0.0) continue;
    if (cfg_.delivery_far_field_radius <= 0.0) {
      carbon_demand_by_cell[static_cast<size_t>(agent.grid_cell)] +=
          agent.pending_carbon_funding;
      continue;
    }
    const auto& support = delivery_support_cells(agent);
    const Real per_cell = support.empty() ? 0.0
        : agent.pending_carbon_funding / static_cast<Real>(support.size());
    Real deposited = 0.0;
    for (size_t i = 0; i < support.size(); ++i) {
      const Real share = i + 1 == support.size()
          ? agent.pending_carbon_funding - deposited : per_cell;
      const Int cell = support[i];
      carbon_demand_by_cell[static_cast<size_t>(cell)] += share;
      deposited += share;
    }
  }
  if (!chem.slab_mode()) {
    // Replicated fields reduce prescribed draws globally, so pro-rata
    // denominators must include every rank's non-ghost agents. Slab fields
    // retain rank-local denominators because each rank owns its agent cells.
    chem.sum_values_across_ranks(carbon_demand_by_cell);
  }
  for (Agent& agent : sim_.agents()) {
    if (agent.state == PhenoState::DEAD || agent.flags.is_ghost
        || agent.grid_cell < 0) {
      agent.pending_carbon_funding = 0.0;
      continue;
    }
    const Real total_demand = agent.pending_growth_carbon
        + agent.pending_maintenance_carbon;
    if (agent.pending_biomass < 0.0 && total_demand <= 0.0) {
      apply_delivery_shrinkage(agent, dt);
      continue;
    }
    if (total_demand <= 0.0) {
      agent.timers.age += dt;
      continue;
    }
    Real field_funding = 0.0;
    if (cfg_.delivery_far_field_radius <= 0.0) {
      const Real cell_requested = agent.grid_cell >= 0
          ? carbon_demand_by_cell[static_cast<size_t>(agent.grid_cell)] : 0.0;
      field_funding = cell_requested > 0.0 && carbon >= 0
          ? chem.prescribed_sink_global(carbon, agent.grid_cell)
              * (agent.pending_carbon_funding / cell_requested)
          : 0.0;
    } else {
      field_funding = delivery_field_funding(
          carbon, agent, agent.pending_carbon_funding,
          carbon_demand_by_cell);
    }
    agent.pending_carbon_funding = field_funding;
    const DeliveryFunding funding = calculate_delivery_funding(
        agent.pending_growth_carbon, agent.pending_maintenance_carbon,
        field_funding);
    record_delivery_funding(
        chem, carbon, agent.pending_growth_carbon,
        agent.pending_maintenance_carbon, funding);
    if (agent.pending_biomass < 0.0) {
      apply_delivery_shrinkage(agent, dt);
    } else {
      apply_delivery_funding(agent, funding, dt);
    }
    if (agent.pending_biomass > 0.0) {
      agent.pending_growth_chemistry_biomass +=
          agent.pending_biomass * funding.growth_fraction;
    }
  }
}

void FixMetabolism::commit_delivery_oxygen(Real dt, Int oxygen) {
  auto& chem = sim_.chemical_field();
  std::vector oxygen_demand_by_cell(
      static_cast<size_t>(chem.global_ncells()), 0.0);
  for (const Agent& agent : sim_.agents()) {
    if (agent.state == PhenoState::DEAD || agent.flags.is_ghost
        || agent.grid_cell < 0) continue;
    if (agent.pending_oxygen_funding <= 0.0) continue;
    if (cfg_.delivery_far_field_radius <= 0.0) {
      oxygen_demand_by_cell[static_cast<size_t>(agent.grid_cell)] +=
          agent.pending_oxygen_funding;
      continue;
    }
    const auto& support = delivery_support_cells(agent);
    const Real per_cell = support.empty() ? 0.0
        : agent.pending_oxygen_funding / static_cast<Real>(support.size());
    Real deposited = 0.0;
    for (size_t i = 0; i < support.size(); ++i) {
      const Real share = i + 1 == support.size()
          ? agent.pending_oxygen_funding - deposited : per_cell;
      const Int cell = support[i];
      oxygen_demand_by_cell[static_cast<size_t>(cell)] += share;
      deposited += share;
    }
  }
  if (!chem.slab_mode()) {
    // Replicated fields reduce prescribed draws globally, so pro-rata
    // denominators must include every rank's non-ghost agents. Slab fields
    // retain rank-local denominators because each rank owns its agent cells.
    chem.sum_values_across_ranks(oxygen_demand_by_cell);
  }
  for (Agent& agent : sim_.agents()) {
    if (agent.state == PhenoState::DEAD || agent.flags.is_ghost
        || agent.grid_cell < 0) {
      agent.pending_oxygen_funding = 0.0;
      agent.respired_oxygen_rate = 0.0;
      continue;
    }
    const Real growth_demand = agent.pending_oxygen_growth;
    const Real maintenance_demand = agent.pending_oxygen_maintenance;
    if (const Real total_demand =
            growth_demand + maintenance_demand; total_demand <= 0.0) {
      agent.respired_oxygen_rate = 0.0;
      continue;
    }
    Real field_funding = 0.0;
    if (cfg_.delivery_far_field_radius <= 0.0) {
      const Real cell_requested = agent.grid_cell >= 0
          ? oxygen_demand_by_cell[static_cast<size_t>(agent.grid_cell)] : 0.0;
      field_funding = cell_requested > 0.0
          ? chem.prescribed_sink_global(oxygen, agent.grid_cell)
              * (agent.pending_oxygen_funding / cell_requested)
          : 0.0;
    } else {
      field_funding = delivery_field_funding(
          oxygen, agent, agent.pending_oxygen_funding,
          oxygen_demand_by_cell);
    }
    agent.pending_oxygen_funding = field_funding;
    const DeliveryFunding funding = calculate_delivery_funding(
        growth_demand, maintenance_demand,
        field_funding);
    record_delivery_funding(
        chem, oxygen, growth_demand, maintenance_demand, funding);
    const auto& oxygen_cfg = sim_.config().chem_env.oxygen;
    if (oxygen_cfg.respiration_driver_mode == RespirationDriver::Funded
        && oxygen_cfg.metabolic_switch_enabled && growth_demand > 0.0) {
      const Real instantaneous = 1.0 - funding.growth_fraction;
      agent.realized_fermentation_fraction = metabolic_mode::relax(
          agent.realized_fermentation_fraction, instantaneous,
          dt, oxygen_cfg.tau_metabolic_switch);
    }
    agent.respired_oxygen_rate =
        (funding.growth_funded + funding.maintenance_paid) / dt;
  }
}

void FixMetabolism::post_chemistry(Real dt) {
  if (cfg_.uptake_limit_mode != UptakeLimitMode::Delivery) return;
  commit_delivery_uptake(dt);
  auto& chem = sim_.chemical_field();
  chem.sum_agent_uptake_across_ranks();
  chem.flux_accounting().commit_agent_uptake_step();
}

void FixMetabolism::charge_carbon_maintenance(Agent& agent, Real dt) {
  if (cfg_.carbon_maintenance_rate <= 0.0 || dt <= 0.0) return;
  const Int cell = agent.grid_cell;
  if (cell < 0) return;
  auto& chem = sim_.chemical_field();
  const Int carbon = chem.find(species::CARBON);
  const Real cell_volume = sim_.domain().cell_volume();
  if (carbon < 0 || cell_volume <= 0.0) return;
  const auto& oxygen = sim_.config().chem_env.oxygen;
  const Real maintenance_factor = oxygen.metabolic_switch_enabled
      ? metabolic_mode::interpolate(1.0, oxygen.anaerobic_maintenance_factor,
                                    agent.realized_fermentation_fraction)
      : 1.0;
  const Real requested_amount = carbon_maintenance::requested(
      cfg_.carbon_maintenance_rate * maintenance_factor, agent.biomass, dt);
  if (cfg_.uptake_limit_mode == UptakeLimitMode::Delivery) {
    agent.pending_maintenance_carbon = requested_amount;
    return;
  }
  Real draw = 0.0;
  if (cfg_.uptake_limit_mode == UptakeLimitMode::Voxel) {
    const Int storage_cell = chem.global_to_storage_cell(cell);
    if (storage_cell >= 0
        && static_cast<size_t>(storage_cell)
               < carbon_maintenance_available_.size()) {
      #ifdef GUTIBM_OPENMP
      #pragma omp critical(gutibm_carbon_maintenance)
      #endif
      {
        Real& available = carbon_maintenance_available_[
            static_cast<size_t>(storage_cell)];
        draw = std::min(requested_amount, available);
        available -= draw;
      }
    }
  } else {
    const ChemicalSpec& spec = chem.spec(carbon);
    const Real allowed = uptake::allowed_uptake_mol(
        to_underlying(cfg_.uptake_limit_mode), chem.conc_global(carbon, cell),
        uptake::effective_diffusivity(spec.diff_coeff, spec.retardation),
        agent.radius, cell_volume, dt);
    draw = allowed < 0.0
        ? requested_amount : std::min(requested_amount, allowed);
  }
  chem.flux_accounting().add_maintenance(carbon, draw);
  if (requested_amount > draw) {
    chem.flux_accounting().add_maintenance_shortfall(
        carbon, requested_amount - draw);
    chem.flux_accounting().add_maintenance_limited_agents(carbon, 1.0);
  }
  #ifdef GUTIBM_OPENMP
  #pragma omp atomic
  #endif
  chem.reac_global(carbon, cell) -= draw / (cell_volume * dt);
}

void FixMetabolism::prepare_carbon_maintenance() {
  carbon_maintenance_available_.clear();
  if (cfg_.carbon_maintenance_rate <= 0.0
      || cfg_.uptake_limit_mode != UptakeLimitMode::Voxel) return;
  const auto& chem = sim_.chemical_field();
  const Int carbon = chem.find(species::CARBON);
  const Real cell_volume = sim_.domain().cell_volume();
  if (carbon < 0 || cell_volume <= 0.0) return;
  carbon_maintenance_available_.resize(
      static_cast<size_t>(chem.ncells()), 0.0);
  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    carbon_maintenance_available_[static_cast<size_t>(cell)] =
        carbon_maintenance::available(chem.conc(carbon, cell), cell_volume);
  }
}

void FixMetabolism::apply_siderophore_chemistry(Real dt) {
  if (const auto& sid_cfg = sim_.config().chem_env.siderophore;
      !sid_cfg.enabled) return;

  const auto& chem = sim_.chemical_field();
  const Int i_sid = chem.find(species::SIDEROPHORE);
  const Int i_iron = chem.find(species::IRON);
  const Int i_ferric_enterobactin =
      chem.find(species::FERRIC_ENTEROBACTIN);
  if (i_sid < 0 || i_ferric_enterobactin < 0) return;

  const Real cell_volume = sim_.domain().cell_volume();
  if (cell_volume <= 0.0) return;

  const Int num_cells = chem.ncells();
  if (biomass_by_cell_.size() != static_cast<size_t>(num_cells)) {
    biomass_by_cell_.assign(static_cast<size_t>(num_cells), 0.0);
    fepA_biomass_by_cell_.assign(static_cast<size_t>(num_cells), 0.0);
    occupancy_by_cell_.assign(static_cast<size_t>(num_cells), 0);
    chelation_by_cell_.assign(static_cast<size_t>(num_cells), 0.0);
    touched_cells_.clear();
  } else {
    for (const Int cell : touched_cells_) {
      const auto index = static_cast<size_t>(cell);
      biomass_by_cell_[index] = 0.0;
      fepA_biomass_by_cell_[index] = 0.0;
      occupancy_by_cell_[index] = 0;
    }
    touched_cells_.clear();
  }
  for (const Agent& agent : sim_.agents()) {
    if (agent.state == PhenoState::DEAD || agent.flags.is_ghost) continue;
    const Int global_cell = agent.grid_cell;
    if (global_cell < 0) continue;
    const Int cell = chem.global_to_storage_cell(global_cell);
    if (cell < 0 || cell >= num_cells) continue;
    const auto index = static_cast<size_t>(cell);
    if (occupancy_by_cell_[index] == 0) touched_cells_.push_back(cell);
    biomass_by_cell_[index] += agent.biomass;
    fepA_biomass_by_cell_[index] +=
        agent.receptor_expr[to_underlying(ReceptorType::FepA)]
        * agent.biomass;
    occupancy_by_cell_[index] += 1;
  }

  if (i_iron >= 0) {
    apply_siderophore_chelation(
        i_sid, i_iron, i_ferric_enterobactin, num_cells);
  }

  apply_siderophore_reimport(i_sid, i_iron, i_ferric_enterobactin, num_cells,
                             cell_volume, dt);
}

void FixMetabolism::apply_siderophore_chelation(
    Int i_sid, Int i_iron, Int i_ferric_enterobactin, Int num_cells) {
  const auto& sid_cfg = sim_.config().chem_env.siderophore;
  auto& chem = sim_.chemical_field();
  const Int x_begin = chem.slab_mode()
      ? chem.owned_storage_x_begin() : sim_.domain().local_grid_x_begin();
  const Int x_end = chem.slab_mode()
      ? chem.owned_storage_x_end() : sim_.domain().local_grid_x_end();
  const Int storage_nx = chem.slab_mode()
      ? chem.storage_nx() : sim_.domain().nx();
  for (Int iz = 0; iz < chem.global_nz(); ++iz) {
    for (Int iy = 0; iy < chem.global_ny(); ++iy) {
      for (Int ix = x_begin; ix < x_end; ++ix) {
        const Int cell = (iz * chem.global_ny() + iy) * storage_nx + ix;
        if (cell < 0 || cell >= num_cells) continue;
        const Real s_sid = chem.conc(i_sid, cell);
        const Real s_iron = chem.conc(i_iron, cell);
        const Real chelation = sid_cfg.chelation_rate * s_sid * s_iron;
        chelation_by_cell_[static_cast<size_t>(cell)] = chelation;
        chem.reac(i_iron, cell) -= chelation;
        chem.reac(i_sid, cell) -= chelation;
        chem.reac(i_ferric_enterobactin, cell) += chelation;
      }
    }
  }
}

void FixMetabolism::apply_siderophore_reimport(
    Int i_sid, Int i_iron, Int i_ferric_enterobactin, Int num_cells,
    Real cell_volume, Real dt) {
  const auto& sid_cfg = sim_.config().chem_env.siderophore;
  auto& chem = sim_.chemical_field();
  const Int x_begin = chem.slab_mode()
      ? chem.owned_storage_x_begin() : sim_.domain().local_grid_x_begin();
  const Int x_end = chem.slab_mode()
      ? chem.owned_storage_x_end() : sim_.domain().local_grid_x_end();
  const Int storage_nx = chem.slab_mode()
      ? chem.storage_nx() : sim_.domain().nx();
  for (Int iz = 0; iz < chem.global_nz(); ++iz) {
    for (Int iy = 0; iy < chem.global_ny(); ++iy) {
      for (Int ix = x_begin; ix < x_end; ++ix) {
        const Int cell = (iz * chem.global_ny() + iy) * storage_nx + ix;
        if (cell < 0 || cell >= num_cells) continue;
        const auto index = static_cast<size_t>(cell);
        if (occupancy_by_cell_[index] == 0) continue;

        const Real s_iron = (i_iron >= 0) ? chem.conc(i_iron, cell) : 0.0;
        Real fur_Km = 1.0e-6;
        if (sim_.config().cell_bio.fur.enabled) {
          fur_Km = sim_.config().cell_bio.fur.Km;
        }
        const Real fur_activity = 1.0 - s_iron / (fur_Km + s_iron);
        const Real biomass_density = biomass_by_cell_[index] / cell_volume;
        const Real sid_rate = sid_cfg.secretion_rate
            * std::max(0.0, fur_activity) * biomass_density;
        chem.reac(i_sid, cell) += sid_rate;

        if (i_iron < 0) continue;
        const Real s_ferric_enterobactin =
            chem.conc(i_ferric_enterobactin, cell);
        const Real vmax = sid_cfg.Vmax_reimport
            * fepA_biomass_by_cell_[index] / cell_volume;
        const Real ferric_after_production = s_ferric_enterobactin
            + chelation_by_cell_[index] * dt;
        const Real reimport = implicit_ferric_enterobactin_reimport(
            ferric_after_production, vmax, sid_cfg.Km_reimport, dt);
        chem.reac(i_ferric_enterobactin, cell) -= reimport;
        chem.reac(i_iron, cell) += reimport;
      }
    }
  }
}

void FixMetabolism::perform_divisions() {
  auto& agents = sim_.agents();
  std::vector<Agent> new_agents;

  for (Agent& a : agents) {
    if (a.state == PhenoState::DEAD || a.flags.is_ghost) continue;

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

Real FixMetabolism::realized_carbon_cost(const Agent& agent) const {
  const auto& oxygen = sim_.config().chem_env.oxygen;
  if (!oxygen.metabolic_switch_enabled) return cfg_.yield_carbon;
  return cfg_.yield_carbon * metabolic_mode::interpolate(
      oxygen.aerobic_carbon_cost_factor,
      oxygen.anaerobic_carbon_cost_factor,
      agent.realized_fermentation_fraction);
}

void FixMetabolism::compute_growth_rate(Agent& agent, Real dt) {
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

  Real S_carbon = (i_carbon >= 0)
      ? chem.conc_global(i_carbon, cell) : 1.0;
  Real S_iron   = (cfg_.iron_uptake_enabled && i_iron >= 0)
      ? chem.conc_global(i_iron, cell) : 1.0;
  Real S_b12    = (cfg_.b12_uptake_enabled && i_b12 >= 0)
      ? chem.conc_global(i_b12, cell) : 1.0;

  if (const auto& fur_cfg = sim_.config().cell_bio.fur;
      fur_cfg.enabled && cfg_.iron_uptake_enabled) {
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
  Real monod_iron = 1.0;
  if (cfg_.iron_uptake_enabled) {
    Real iron_uptake = 0.0;
    iron_uptake += expr_fepA * lig_aff_fepA * S_iron
        / (cfg_.km_iron_primary + S_iron);
    iron_uptake += expr_iroN * lig_aff_iroN * S_iron
        / (cfg_.km_iron_iroN + S_iron);
    iron_uptake += expr_iutA * lig_aff_iutA * S_iron
        / (cfg_.km_iron_iutA + S_iron);
    iron_uptake += expr_fiu * lig_aff_fiu * S_iron
        / (cfg_.km_iron_fiu + S_iron);
    monod_iron = iron_uptake / (1.0 + expr_iroN + expr_iutA + expr_fiu);
  }

  Real Km_b12  = agent.km.km_b12  / (expr_btuB * lig_aff_btuB);
  Real Km_carb = agent.km.km_carbon;

  // Triple Monod kinetics (uncoupled)
  Real monod_carbon = S_carbon / (Km_carb + S_carbon);
  Real monod_b12 = cfg_.b12_uptake_enabled
      ? S_b12 / (Km_b12 + S_b12) : 1.0;

  Real mu = agent.mu_max * monod_carbon * monod_iron * monod_b12;

  if (const auto& o2cfg = sim_.config().chem_env.oxygen; o2cfg.enabled) {
    if (Int i_o2 = chem.find(species::OXYGEN); i_o2 >= 0) {
      const Real s_o2 = chem.conc_global(i_o2, cell);
      if (o2cfg.metabolic_switch_enabled) {
        const Real availability =
            metabolic_mode::oxygen_availability(s_o2, o2cfg.Km);
        if (o2cfg.respiration_driver_mode == RespirationDriver::Ambient) {
          const Real instantaneous = metabolic_mode::fermentation_fraction(
              availability, mu, o2cfg.mu_crit);
          agent.realized_fermentation_fraction = metabolic_mode::relax(
              agent.realized_fermentation_fraction, instantaneous,
              dt, o2cfg.tau_metabolic_switch);
        }
        mu *= metabolic_mode::interpolate(o2cfg.aerobic_mu_factor,
                                          o2cfg.anaerobic_mu_factor,
                                          agent.realized_fermentation_fraction);
      } else {
        const Real availability =
            metabolic_mode::oxygen_availability(s_o2, o2cfg.Km);
        const Real monod_o2_boost =
            1.0 + o2cfg.boost_max * availability;
        mu *= monod_o2_boost;
      }
    }
  }

  if (cfg_.acid_inhibition_enabled) {
    const Int i_acetate = chem.find(species::ACETATE);
    if (i_acetate >= 0) {
      const Real acetate = chem.conc_global(i_acetate, cell);
      const Real inhibition = metabolic_mode::acid_inhibition(
          acetate, sim_.config().fixes.bacteriocin.mucin_charge.ph,
          cfg_.acetate_pKa, cfg_.acid_inhibition_Ki,
          cfg_.acid_inhibition_max);
      mu *= 1.0 - inhibition;
    }
  }

  // Metabolic penalties for receptor downregulation
  // BtuB loss → MetE pathway required (proteome cost)
  // Acetate inhibits MetE, scaling the penalty with local [acetate]
  // + concentration-dependent ethanolamine utilization loss
  if (expr_btuB < 0.5) {
    Real metE_eff = cfg_.metE_penalty;
    if (Int i_acetate = chem.find(species::ACETATE); i_acetate >= 0) {
      Real acetate_conc = chem.conc_global(i_acetate, cell);
      Real acetate_factor = 1.0
          + (cfg_.metE_acetate_max_factor - 1.0)
            * acetate_conc / (cfg_.metE_acetate_km + acetate_conc);
      metE_eff *= acetate_factor;
    }
    Real eut_conc = 0.0;
    if (cfg_.eut_enabled) {
      if (Int i_eut = chem.find(species::ETHANOLAMINE); i_eut >= 0) {
      eut_conc = chem.conc_global(i_eut, cell);
      }
    }
    Real eut_effect = cfg_.eut_enabled
        ? cfg_.eut_max_penalty * eut_conc / (cfg_.eut_km + eut_conc) : 0.0;
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

Real FixMetabolism::uptake_limit_fraction(
    const Agent& agent, Real d_biomass, Real dt, bool record_diagnostics) {
  const Int cell = agent.grid_cell;
  if (cell < 0 || d_biomass <= 0.0 || dt <= 0.0) return 1.0;

  auto& chem = sim_.chemical_field();
  const Real cell_vol = sim_.domain().cell_volume();
  const int mode = to_underlying(cfg_.uptake_limit_mode);
  Real fraction = 1.0;

  const auto limit_species = [&](Int index, Real yield) {
    if (index < 0) return;
    const Real demanded = d_biomass * yield;
    if (demanded <= 0.0) return;
    if (record_diagnostics) {
      chem.flux_accounting().add_uptake_demand(index, demanded);
    }
    const ChemicalSpec& spec = chem.spec(index);
    const Real allowed = uptake::allowed_uptake_mol(
        mode, chem.conc_global(index, cell),
        uptake::effective_diffusivity(spec.diff_coeff, spec.retardation),
        agent.radius, cell_vol, dt);
    const Real species_fraction = uptake::limit_fraction(allowed, demanded);
    if (species_fraction < 1.0) {
      if (record_diagnostics) {
        chem.flux_accounting().add_uptake_limited(index, 1.0);
      }
      fraction = std::min(fraction, species_fraction);
    }
  };

  limit_species(chem.find(species::CARBON), realized_carbon_cost(agent));
  if (cfg_.iron_uptake_enabled) {
    limit_species(chem.find(species::IRON), cfg_.yield_iron);
  }
  return fraction;
}

void FixMetabolism::grow_agent(Agent& agent, Real dt) {
  // Biomass increase, funded by the uptake the agent can actually acquire
  Real d_biomass = agent.mu_realized * agent.biomass * dt;
  const Real funded = uptake_limit_fraction(agent, d_biomass, dt, true);
  if (funded < 1.0) {
    d_biomass *= funded;
    agent.mu_realized *= funded;
  }
  agent.biomass += d_biomass;
  agent.biomass = std::max(agent.biomass, 1.0e-20);

  // Update radius from biomass
  Real vol = agent.biomass / CELL_DENSITY_DEFAULT;
  agent.radius = std::cbrt(3.0 * vol / (4.0 * PI));
  agent.mass   = agent.biomass;
  agent.timers.age   += dt;
  apply_growth_chemistry(agent, d_biomass, dt);
}

void FixMetabolism::apply_growth_chemistry(
    const Agent& agent, Real d_biomass, Real dt) {
  // Nutrient consumption and siderophore coupling from grid
  auto& chem = sim_.chemical_field();
  Int cell = agent.grid_cell;
  if (cell < 0) return;

  Int i_carbon = chem.find(species::CARBON);
  Int i_iron   = chem.find(species::IRON);
  Int i_acetate = chem.find(species::ACETATE);

  Real cell_vol = sim_.domain().cell_volume();

  if (d_biomass <= 0.0 || dt <= 0.0) return;

  if (i_carbon >= 0 && cell_vol > 0.0) {
    const Real carbon_yield = realized_carbon_cost(agent);
    Real delta_c = d_biomass * carbon_yield / (cell_vol * dt);
    if (cfg_.uptake_limit_mode != UptakeLimitMode::Delivery) {
      chem.flux_accounting().add_agent_uptake(
          i_carbon, d_biomass * carbon_yield);
      #ifdef GUTIBM_OPENMP
      #pragma omp atomic
      #endif
      chem.reac_global(i_carbon, cell) -= delta_c;
    }
  }
  if (cfg_.iron_uptake_enabled && i_iron >= 0 && cell_vol > 0.0) {
    Real delta_fe = d_biomass * cfg_.yield_iron / (cell_vol * dt);
    chem.flux_accounting().add_agent_uptake(
        i_iron, d_biomass * cfg_.yield_iron);
    #ifdef GUTIBM_OPENMP
    #pragma omp atomic
    #endif
    chem.reac_global(i_iron, cell) -= delta_fe;
  }
  // Spec 6 §3 — the B12/corrinoid field is NOT depleted. It represents the
  // total bioavailable corrinoid pool (~1 uM), produced by the anaerobic
  // majority at rates far exceeding E. coli demand; the pool is effectively
  // constant relative to the modeled population. (cfg_.yield_b12 retained for
  // config compatibility but no longer removes corrinoid from the field.)

  const auto& acfg = sim_.config().chem_env.acetate;
  if (acfg.enabled && i_acetate >= 0 && cell_vol > 0.0) {
    const Real acetate_conc = chem.conc_global(i_acetate, cell);
    if (!sim_.config().chem_env.oxygen.metabolic_switch_enabled
        && agent.mu_realized > acfg.overflow_threshold) {
      #ifdef GUTIBM_OPENMP
      #pragma omp atomic
      #endif
      chem.reac_global(i_acetate, cell) +=
          acfg.overflow_rate * agent.biomass / cell_vol;
    }
    if (sim_.config().chem_env.oxygen.metabolic_switch_enabled) {
      const Real acid =
          sim_.config().chem_env.oxygen.ferm_acid_yield
          * agent.realized_fermentation_fraction * d_biomass
          * realized_carbon_cost(agent);
      if (acid > 0.0) {
        #ifdef GUTIBM_OPENMP
        #pragma omp atomic
        #endif
        chem.reac_global(i_acetate, cell) += acid / (cell_vol * dt);
      }
    }
    const Real scavenge = acfg.scavenge_rate * acetate_conc
        / (acfg.scavenge_Km + acetate_conc) * agent.biomass / cell_vol;
    #ifdef GUTIBM_OPENMP
    #pragma omp atomic
    #endif
    chem.reac_global(i_acetate, cell) -= scavenge;
  }
}

}  // namespace gutibm
