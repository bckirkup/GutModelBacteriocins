#include "agent_pool_gpu.h"
#include "agent.h"
#include "chemical_field.h"
#include "chemical_field_gpu.h"
#include "chem_environment_config.h"
#include "domain.h"
#include "fix_metabolism.h"
#include "gpu_kernels.h"
#include "dispatch.h"
#include "receptor_utils.h"
#include "species_names.h"

#include <cassert>

namespace gutibm {

namespace {

#ifdef GUTIBM_CUDA
void record_device_delivery_ledger(
    ChemicalField& chem, Int species_index, const std::vector<double>& ledger) {
  if (species_index < 0 || ledger.size() < 6) return;
  auto& flux = chem.flux_accounting();
  flux.add_maintenance(species_index, ledger[0]);
  if (ledger[1] > 0.0) flux.add_maintenance_shortfall(species_index, ledger[1]);
  if (ledger[2] > 0.0) {
    flux.add_maintenance_limited_agents(species_index, ledger[2]);
  }
  flux.add_agent_uptake(species_index, ledger[3]);
  if (ledger[4] > 0.0) flux.add_uptake_shortfall(species_index, ledger[4]);
  if (ledger[5] > 0.0) flux.add_uptake_limited(species_index, ledger[5]);
}
#endif

}  // namespace

void AgentPoolGpu::resize(Int n) {
  size_ = n;
  if (n <= 0) return;
  const auto count = static_cast<size_t>(n);
  d_x_.allocate(count);
  d_y_.allocate(count);
  d_z_.allocate(count);
  d_death_time_.allocate(count);
  d_grid_cell_.allocate(count);
  d_state_.allocate(count);
  d_is_ghost_.allocate(count);
  d_has_prophage_immunity_.allocate(count);
  d_has_inducible_prophage_.allocate(count);
  d_prophage_induced_.allocate(count);
  d_mu_realized_.allocate(count);
  d_fermentation_fraction_.allocate(count);
  d_biomass_.allocate(count);
  d_radius_.allocate(count);
  d_mass_.allocate(count);
  d_age_.allocate(count);
  d_mu_max_.allocate(count);
  d_km_b12_.allocate(count);
  d_km_carbon_.allocate(count);
  d_receptor_expr_.allocate(static_cast<size_t>(NUM_RECEPTORS) * count);
  d_receptor_expr_base_.allocate(static_cast<size_t>(NUM_RECEPTORS) * count);
  d_ligand_affinity_.allocate(static_cast<size_t>(NUM_RECEPTORS) * count);
  d_iron_receptor_.allocate(static_cast<size_t>(NUM_RECEPTORS));
  d_bi_loci_count_.allocate(count);
  d_plasmid_amelioration_.allocate(count);
  d_pending_growth_carbon_.allocate(count);
  d_pending_maintenance_carbon_.allocate(count);
  d_pending_carbon_ceiling_.allocate(count);
  d_pending_carbon_funding_.allocate(count);
  d_pending_biomass_.allocate(count);
  d_pending_growth_chemistry_biomass_.allocate(count);
  d_pending_oxygen_growth_.allocate(count);
  d_pending_oxygen_maintenance_.allocate(count);
  d_pending_oxygen_ceiling_.allocate(count);
  d_pending_oxygen_funding_.allocate(count);
  d_respired_oxygen_rate_.allocate(count);
}

void AgentPoolGpu::sync_from_host(const AgentPool& pool) {
  Int n = pool.size();
  resize(n);
  if (n <= 0) return;

  std::vector<double> x(n);
  std::vector<double> y(n);
  std::vector<double> z(n);
  std::vector<double> death_time(n);
  std::vector<int> grid_cell(n);
  std::vector<int> state(n);
  std::vector<int> is_ghost(n);
  std::vector<int> has_prophage_immunity(n);
  std::vector<int> has_inducible_prophage(n);
  std::vector<int> prophage_induced(n);
  std::vector<int> bi_loci_count(n);
  std::vector<double> mu_realized(n);
  std::vector<double> fermentation_fraction(n);
  std::vector<double> biomass(n);
  std::vector<double> radius(n);
  std::vector<double> mass(n);
  std::vector<double> age(n);
  std::vector<double> mu_max(n);
  std::vector<double> km_b12(n);
  std::vector<double> km_carbon(n);
  std::vector<double> receptor_expr(static_cast<size_t>(NUM_RECEPTORS) * n);
  std::vector<double> receptor_expr_base(
      static_cast<size_t>(NUM_RECEPTORS) * n);
  std::vector<double> ligand_affinity(static_cast<size_t>(NUM_RECEPTORS) * n);
  std::vector<int> iron_receptor(static_cast<size_t>(NUM_RECEPTORS));
  std::vector<double> plasmid_amelioration(n);
  std::vector<double> pending_growth_carbon(n);
  std::vector<double> pending_maintenance_carbon(n);
  std::vector<double> pending_carbon_ceiling(n);
  std::vector<double> pending_carbon_funding(n);
  std::vector<double> pending_biomass(n);
  std::vector<double> pending_growth_chemistry_biomass(n);
  std::vector<double> pending_oxygen_growth(n);
  std::vector<double> pending_oxygen_maintenance(n);
  std::vector<double> pending_oxygen_ceiling(n);
  std::vector<double> pending_oxygen_funding(n);
  std::vector<double> respired_oxygen_rate(n);

  for (Int i = 0; i < n; ++i) {
    const Agent& a = pool[i];
    x[i] = a.x[0];
    y[i] = a.x[1];
    z[i] = a.x[2];
    death_time[i] = a.timers.death_time;
    grid_cell[i] = a.grid_cell;
    state[i] = static_cast<int>(to_underlying(a.state));
    is_ghost[i] = a.flags.is_ghost ? 1 : 0;
    has_prophage_immunity[i] = a.flags.has_prophage_immunity ? 1 : 0;
    has_inducible_prophage[i] = a.flags.has_inducible_prophage ? 1 : 0;
    prophage_induced[i] = a.flags.prophage_induced ? 1 : 0;
    mu_realized[i] = a.mu_realized;
    fermentation_fraction[i] = a.realized_fermentation_fraction;
    biomass[i] = a.biomass;
    radius[i] = a.radius;
    mass[i] = a.mass;
    age[i] = a.timers.age;
    mu_max[i] = a.mu_max;
    km_b12[i] = a.km.km_b12;
    km_carbon[i] = a.km.km_carbon;
    bi_loci_count[i] = static_cast<int>(a.genome.bi_loci.size());
    plasmid_amelioration[i] = a.genome.plasmid_cost_amelioration;
    pending_growth_carbon[i] = a.pending_growth_carbon;
    pending_maintenance_carbon[i] = a.pending_maintenance_carbon;
    pending_carbon_ceiling[i] = a.pending_carbon_ceiling;
    pending_carbon_funding[i] = a.pending_carbon_funding;
    pending_biomass[i] = a.pending_biomass;
    pending_growth_chemistry_biomass[i] = a.pending_growth_chemistry_biomass;
    pending_oxygen_growth[i] = a.pending_oxygen_growth;
    pending_oxygen_maintenance[i] = a.pending_oxygen_maintenance;
    pending_oxygen_ceiling[i] = a.pending_oxygen_ceiling;
    pending_oxygen_funding[i] = a.pending_oxygen_funding;
    respired_oxygen_rate[i] = a.respired_oxygen_rate;
    for (int r = 0; r < NUM_RECEPTORS; ++r) {
      // SoA layout: [receptor][agent]. GPU metabolism/receptor kernels index as
      // receptor_expr[r * n + i] — do not pack as AoS (i * NUM_RECEPTORS + r).
      receptor_expr[static_cast<size_t>(r) * n + i] = a.receptor_expr[r];
      receptor_expr_base[static_cast<size_t>(r) * n + i] =
          a.receptor_expr_base[r];
      ligand_affinity[static_cast<size_t>(r) * n + i] = a.genome.ligand_affinity[r];
    }
  }
  for (int r = 0; r < NUM_RECEPTORS; ++r) {
    iron_receptor[static_cast<size_t>(r)] = is_iron_receptor(r) ? 1 : 0;
  }

  d_x_.upload(x);
  d_y_.upload(y);
  d_z_.upload(z);
  d_death_time_.upload(death_time);
  d_grid_cell_.upload(grid_cell);
  d_state_.upload(state);
  d_is_ghost_.upload(is_ghost);
  d_has_prophage_immunity_.upload(has_prophage_immunity);
  d_has_inducible_prophage_.upload(has_inducible_prophage);
  d_prophage_induced_.upload(prophage_induced);
  d_mu_realized_.upload(mu_realized);
  d_fermentation_fraction_.upload(fermentation_fraction);
  d_biomass_.upload(biomass);
  d_radius_.upload(radius);
  d_mass_.upload(mass);
  d_age_.upload(age);
  d_mu_max_.upload(mu_max);
  d_km_b12_.upload(km_b12);
  d_km_carbon_.upload(km_carbon);
  d_receptor_expr_.upload(receptor_expr);
  d_receptor_expr_base_.upload(receptor_expr_base);
  d_ligand_affinity_.upload(ligand_affinity);
  d_iron_receptor_.upload(iron_receptor);
  d_bi_loci_count_.upload(bi_loci_count);
  d_plasmid_amelioration_.upload(plasmid_amelioration);
  d_pending_growth_carbon_.upload(pending_growth_carbon);
  d_pending_maintenance_carbon_.upload(pending_maintenance_carbon);
  d_pending_carbon_ceiling_.upload(pending_carbon_ceiling);
  d_pending_carbon_funding_.upload(pending_carbon_funding);
  d_pending_biomass_.upload(pending_biomass);
  d_pending_growth_chemistry_biomass_.upload(pending_growth_chemistry_biomass);
  d_pending_oxygen_growth_.upload(pending_oxygen_growth);
  d_pending_oxygen_maintenance_.upload(pending_oxygen_maintenance);
  d_pending_oxygen_ceiling_.upload(pending_oxygen_ceiling);
  d_pending_oxygen_funding_.upload(pending_oxygen_funding);
  d_respired_oxygen_rate_.upload(respired_oxygen_rate);
}

void AgentPoolGpu::sync_to_host(AgentPool& pool) const {
  Int n = pool.size();
  if (n <= 0 || n != size_) return;

  std::vector<double> mu_realized(n);
  std::vector<double> fermentation_fraction(n);
  std::vector<double> biomass(n);
  std::vector<double> radius(n);
  std::vector<double> mass(n);
  std::vector<double> age(n);
  std::vector<int> grid_cell(n);
  d_mu_realized_.download(mu_realized);
  d_fermentation_fraction_.download(fermentation_fraction);
  d_biomass_.download(biomass);
  d_radius_.download(radius);
  d_mass_.download(mass);
  d_age_.download(age);
  d_grid_cell_.download(grid_cell);

  for (Int i = 0; i < n; ++i) {
    Agent& a = pool[i];
    a.mu_realized = mu_realized[i];
    a.realized_fermentation_fraction = fermentation_fraction[i];
    a.biomass = biomass[i];
    a.radius = radius[i];
    a.mass = mass[i];
    a.timers.age = age[i];
    a.grid_cell = grid_cell[i];
  }
}

void AgentPoolGpu::sync_delivery_to_host(AgentPool& pool) const {
  Int n = pool.size();
  if (n <= 0 || n != size_) return;

  std::vector<double> pending_growth_carbon(n);
  std::vector<double> pending_maintenance_carbon(n);
  std::vector<double> pending_carbon_ceiling(n);
  std::vector<double> pending_carbon_funding(n);
  std::vector<double> pending_biomass(n);
  std::vector<double> pending_growth_chemistry_biomass(n);
  std::vector<double> pending_oxygen_growth(n);
  std::vector<double> pending_oxygen_maintenance(n);
  std::vector<double> pending_oxygen_ceiling(n);
  std::vector<double> pending_oxygen_funding(n);
  std::vector<double> respired_oxygen_rate(n);

  d_pending_growth_carbon_.download(pending_growth_carbon);
  d_pending_maintenance_carbon_.download(pending_maintenance_carbon);
  d_pending_carbon_ceiling_.download(pending_carbon_ceiling);
  d_pending_carbon_funding_.download(pending_carbon_funding);
  d_pending_biomass_.download(pending_biomass);
  d_pending_growth_chemistry_biomass_.download(
      pending_growth_chemistry_biomass);
  d_pending_oxygen_growth_.download(pending_oxygen_growth);
  d_pending_oxygen_maintenance_.download(pending_oxygen_maintenance);
  d_pending_oxygen_ceiling_.download(pending_oxygen_ceiling);
  d_pending_oxygen_funding_.download(pending_oxygen_funding);
  d_respired_oxygen_rate_.download(respired_oxygen_rate);

  for (Int i = 0; i < n; ++i) {
    Agent& a = pool[i];
    a.pending_growth_carbon = pending_growth_carbon[i];
    a.pending_maintenance_carbon = pending_maintenance_carbon[i];
    a.pending_carbon_ceiling = pending_carbon_ceiling[i];
    a.pending_carbon_funding = pending_carbon_funding[i];
    a.pending_biomass = pending_biomass[i];
    a.pending_growth_chemistry_biomass = pending_growth_chemistry_biomass[i];
    a.pending_oxygen_growth = pending_oxygen_growth[i];
    a.pending_oxygen_maintenance = pending_oxygen_maintenance[i];
    a.pending_oxygen_ceiling = pending_oxygen_ceiling[i];
    a.pending_oxygen_funding = pending_oxygen_funding[i];
    a.respired_oxygen_rate = respired_oxygen_rate[i];
  }
}

void AgentPoolGpu::sync_receptor_expression_to_host(AgentPool& pool) const {
  Int n = pool.size();
  if (n <= 0 || n != size_) return;

  std::vector<double> receptor_expr(
      static_cast<size_t>(NUM_RECEPTORS) * n);
  d_receptor_expr_.download(receptor_expr);
  for (Int i = 0; i < n; ++i) {
    Agent& agent = pool[i];
    for (int r = 0; r < NUM_RECEPTORS; ++r) {
      agent.receptor_expr[r] =
          receptor_expr[static_cast<size_t>(r) * n + i];
    }
  }
}

void AgentPoolGpu::sync_positions_to_host(AgentPool& pool) const {
  Int n = pool.size();
  if (n <= 0 || n != size_) return;

  std::vector<double> x(n);
  std::vector<double> y(n);
  std::vector<double> z(n);
  d_x_.download(x);
  d_y_.download(y);
  d_z_.download(z);

  for (Int i = 0; i < n; ++i) {
    pool[i].x[0] = x[i];
    pool[i].x[1] = y[i];
    pool[i].x[2] = z[i];
  }
}

bool AgentPoolGpu::run_metabolism(
    const Domain& domain, const MetabolismConfig& cfg,
    const GpuMetabolismBuffers& buffers, double* uptake_totals, double dt,
    Int num_agents) {

#ifndef GUTIBM_CUDA
  (void)domain;
  (void)cfg;
  (void)buffers;
  (void)uptake_totals;
  (void)dt;
  (void)num_agents;
  return false;
#else
  if (!gpu_runtime_enabled() || size_ <= 0) return false;
  assert(buffers.receptor_count == NUM_RECEPTORS);

  gpu::launch_metabolism_kernel(
      buffers.d_conc_carbon, buffers.d_conc_iron, buffers.d_conc_b12,
      buffers.d_conc_acetate, buffers.d_conc_eut,
      buffers.d_conc_oxygen,
      buffers.d_reac_carbon, buffers.d_reac_iron, buffers.d_reac_b12,
      buffers.d_reac_acetate,
      d_mu_realized_.data(), d_biomass_.data(), d_radius_.data(),
      d_mass_.data(), d_age_.data(),
      d_fermentation_fraction_.data(),
      d_grid_cell_.data(), d_state_.data(),
      d_mu_max_.data(), d_km_b12_.data(), d_km_carbon_.data(),
      d_receptor_expr_.data(), d_receptor_expr_base_.data(),
      d_ligand_affinity_.data(), d_iron_receptor_.data(),
      d_bi_loci_count_.data(), d_plasmid_amelioration_.data(),
      size_, num_agents, size_, buffers.receptor_count, dt,
      domain.cell_volume(), CELL_DENSITY_DEFAULT,
      cfg.km_iron_primary, cfg.km_iron_iroN, cfg.km_iron_iutA, cfg.km_iron_fiu,
      cfg.maintenance_rate, cfg.metE_penalty, cfg.metE_acetate_max_factor,
      cfg.carbon_maintenance_rate,
      cfg.metE_acetate_km, cfg.eut_max_penalty, cfg.eut_km,
      cfg.yield_carbon, cfg.yield_iron, cfg.yield_b12,
      buffers.iron_uptake_enabled,
      buffers.b12_uptake_enabled,
      buffers.eut_enabled,
      buffers.fur_enabled,
      buffers.fur_Km, buffers.fur_upregulation_max, buffers.fur_receptor_max,
      buffers.acetate_enabled,
      buffers.acetate_overflow_threshold, buffers.acetate_overflow_rate,
      buffers.acetate_scavenge_rate, buffers.acetate_scavenge_Km,
      buffers.o2_enabled, buffers.o2_boost_max, buffers.o2_Km,
      buffers.metabolic_switch_enabled, buffers.respiration_driver_mode,
      buffers.mu_crit,
      buffers.aerobic_mu_factor, buffers.anaerobic_mu_factor,
      buffers.aerobic_carbon_cost_factor,
      buffers.anaerobic_carbon_cost_factor,
      buffers.tau_metabolic_switch, buffers.ferm_acid_yield,
      buffers.anaerobic_maintenance_factor,
      buffers.acid_inhibition_enabled, buffers.acid_inhibition_max,
      buffers.acid_inhibition_Ki, buffers.acetate_pKa,
      buffers.environment_pH,
      uptake_totals,
      buffers.d_maintenance_available,
      to_underlying(cfg.uptake_limit_mode),
      buffers.effective_diffusivity_carbon,
      buffers.effective_diffusivity_iron,
      buffers.d_uptake_limit_totals,
      domain.nx(), domain.ny(), buffers.storage_nx,
      buffers.owned_global_x_begin, buffers.owned_global_x_end,
      buffers.owned_storage_x_begin, gpu_compute_stream());

  gpu_sync_compute();
  gpu_check_error("metabolism_kernel");
  ++metabolism_gpu_steps_;
  return true;
#endif
}

bool AgentPoolGpu::run_delivery_prepare(
    const Domain& domain, const MetabolismConfig& cfg,
    const GpuMetabolismBuffers& buffers,
    const GpuDeliverySupportCsr& support,
    Int ncells, double dt, Int local_agent_count,
    GpuDeliveryPrepareOutput& output) {
#ifndef GUTIBM_CUDA
  (void)domain;
  (void)cfg;
  (void)buffers;
  (void)support;
  (void)ncells;
  (void)dt;
  (void)local_agent_count;
  (void)output;
  return false;
#else
  if (!gpu_runtime_enabled() || size_ <= 0 || ncells <= 0) return false;
  assert(buffers.receptor_count == NUM_RECEPTORS);

  const auto cell_count = static_cast<size_t>(ncells);
  DeviceBuffer<double> d_prescribed_carbon_local;
  double* prescribed_carbon_ptr = buffers.d_prescribed_carbon;
  if (prescribed_carbon_ptr == nullptr) {
    d_prescribed_carbon_local.allocate(cell_count);
    d_prescribed_carbon_local.zero();
    prescribed_carbon_ptr = d_prescribed_carbon_local.data();
  }
  DeviceBuffer<double> d_prescribed_oxygen_local;
  double* prescribed_oxygen_ptr = nullptr;
  if (buffers.o2_delivery_enabled != 0) {
    prescribed_oxygen_ptr = buffers.d_prescribed_oxygen;
    if (prescribed_oxygen_ptr == nullptr) {
      d_prescribed_oxygen_local.allocate(cell_count);
      d_prescribed_oxygen_local.zero();
      prescribed_oxygen_ptr = d_prescribed_oxygen_local.data();
    }
  }
  DeviceBuffer<double> d_totals(3);
  d_totals.zero();

  DeviceBuffer<int> d_support_offsets;
  DeviceBuffer<int> d_support_cells;
  if (support.active()) {
    d_support_offsets.upload(support.offsets);
    if (!support.cells.empty()) {
      d_support_cells.upload(support.cells);
    }
  }

  gpu::DeliveryPrepareParams params;
  params.num_agents = size_;
  params.local_agent_count = local_agent_count;
  params.agent_stride = size_;
  params.receptor_count = buffers.receptor_count;
  params.dt = dt;
  params.cell_volume = domain.cell_volume();
  params.km_iron_primary = cfg.km_iron_primary;
  params.km_iron_iroN = cfg.km_iron_iroN;
  params.km_iron_iutA = cfg.km_iron_iutA;
  params.km_iron_fiu = cfg.km_iron_fiu;
  params.maintenance_rate = cfg.maintenance_rate;
  params.metE_penalty = cfg.metE_penalty;
  params.metE_acetate_max_factor = cfg.metE_acetate_max_factor;
  params.carbon_maintenance_rate = cfg.carbon_maintenance_rate;
  params.metE_acetate_km = cfg.metE_acetate_km;
  params.eut_max_penalty = cfg.eut_max_penalty;
  params.eut_km = cfg.eut_km;
  params.yield_carbon = cfg.yield_carbon;
  params.yield_iron = cfg.yield_iron;
  params.iron_uptake_enabled = buffers.iron_uptake_enabled;
  params.b12_uptake_enabled = buffers.b12_uptake_enabled;
  params.eut_enabled = buffers.eut_enabled;
  params.fur_enabled = buffers.fur_enabled;
  params.fur_Km = buffers.fur_Km;
  params.fur_upregulation_max = buffers.fur_upregulation_max;
  params.fur_receptor_max = buffers.fur_receptor_max;
  params.acetate_enabled = buffers.acetate_enabled;
  params.acetate_overflow_threshold = buffers.acetate_overflow_threshold;
  params.acetate_overflow_rate = buffers.acetate_overflow_rate;
  params.acetate_scavenge_rate = buffers.acetate_scavenge_rate;
  params.acetate_scavenge_Km = buffers.acetate_scavenge_Km;
  params.o2_enabled = buffers.o2_enabled;
  params.o2_delivery_enabled = buffers.o2_delivery_enabled;
  params.o2_boost_max = buffers.o2_boost_max;
  params.o2_Km = buffers.o2_Km;
  params.o2_q_consumption = buffers.o2_q_consumption;
  params.o2_q_maintenance = buffers.o2_q_maintenance;
  params.metabolic_switch_enabled = buffers.metabolic_switch_enabled;
  params.respiration_driver_mode = buffers.respiration_driver_mode;
  params.mu_crit = buffers.mu_crit;
  params.aerobic_mu_factor = buffers.aerobic_mu_factor;
  params.anaerobic_mu_factor = buffers.anaerobic_mu_factor;
  params.aerobic_carbon_cost_factor = buffers.aerobic_carbon_cost_factor;
  params.anaerobic_carbon_cost_factor = buffers.anaerobic_carbon_cost_factor;
  params.tau_metabolic_switch = buffers.tau_metabolic_switch;
  params.ferm_acid_yield = buffers.ferm_acid_yield;
  params.anaerobic_maintenance_factor = buffers.anaerobic_maintenance_factor;
  params.acid_inhibition_enabled = buffers.acid_inhibition_enabled;
  params.acid_inhibition_max = buffers.acid_inhibition_max;
  params.acid_inhibition_Ki = buffers.acid_inhibition_Ki;
  params.acetate_pKa = buffers.acetate_pKa;
  params.environment_pH = buffers.environment_pH;
  params.effective_diffusivity_carbon = buffers.effective_diffusivity_carbon;
  params.effective_diffusivity_oxygen = buffers.effective_diffusivity_oxygen;
  params.global_nx = domain.nx();
  params.global_ny = domain.ny();
  params.storage_nx = buffers.storage_nx;
  params.owned_global_x_begin = buffers.owned_global_x_begin;
  params.owned_global_x_end = buffers.owned_global_x_end;
  params.owned_storage_x_begin = buffers.owned_storage_x_begin;

  gpu::DeliveryPrepareDeviceBuffers dev;
  dev.conc_carbon = buffers.d_conc_carbon;
  dev.conc_iron = buffers.d_conc_iron;
  dev.conc_b12 = buffers.d_conc_b12;
  dev.conc_acetate = buffers.d_conc_acetate;
  dev.conc_eut = buffers.d_conc_eut;
  dev.conc_oxygen = buffers.d_conc_oxygen;
  dev.reac_iron = buffers.d_reac_iron;
  dev.reac_acetate = buffers.d_reac_acetate;
  dev.prescribed_carbon = prescribed_carbon_ptr;
  dev.prescribed_oxygen = prescribed_oxygen_ptr;
  dev.mu_realized = d_mu_realized_.data();
  dev.biomass = d_biomass_.data();
  dev.radius = d_radius_.data();
  dev.fermentation_fraction = d_fermentation_fraction_.data();
  dev.pending_growth_carbon = d_pending_growth_carbon_.data();
  dev.pending_maintenance_carbon = d_pending_maintenance_carbon_.data();
  dev.pending_carbon_ceiling = d_pending_carbon_ceiling_.data();
  dev.pending_carbon_funding = d_pending_carbon_funding_.data();
  dev.pending_biomass = d_pending_biomass_.data();
  dev.pending_growth_chemistry_biomass =
      d_pending_growth_chemistry_biomass_.data();
  dev.pending_oxygen_growth = d_pending_oxygen_growth_.data();
  dev.pending_oxygen_maintenance = d_pending_oxygen_maintenance_.data();
  dev.pending_oxygen_ceiling = d_pending_oxygen_ceiling_.data();
  dev.pending_oxygen_funding = d_pending_oxygen_funding_.data();
  dev.grid_cell = d_grid_cell_.data();
  dev.state = d_state_.data();
  dev.mu_max = d_mu_max_.data();
  dev.km_b12 = d_km_b12_.data();
  dev.km_carbon = d_km_carbon_.data();
  dev.receptor_expr = d_receptor_expr_.data();
  dev.receptor_expr_base = d_receptor_expr_base_.data();
  dev.ligand_affinity = d_ligand_affinity_.data();
  dev.iron_receptor = d_iron_receptor_.data();
  dev.bi_loci_count = d_bi_loci_count_.data();
  dev.plasmid_amelioration = d_plasmid_amelioration_.data();
  dev.support_offsets = support.active() ? d_support_offsets.data() : nullptr;
  dev.support_cells =
      (support.active() && !support.cells.empty()) ? d_support_cells.data()
                                                   : nullptr;
  dev.iron_uptake_total = d_totals.data();
  dev.carbon_demand_total = d_totals.data() + 1;
  dev.oxygen_demand_total = d_totals.data() + 2;

  gpu::launch_delivery_prepare_kernel(params, dev, gpu_compute_stream());
  gpu_sync_compute();
  gpu_check_error("delivery_prepare_kernel");

  if (buffers.d_prescribed_carbon == nullptr) {
    d_prescribed_carbon_local.download(output.prescribed_carbon);
  }
  if (buffers.o2_delivery_enabled != 0) {
    if (buffers.d_prescribed_oxygen == nullptr) {
      d_prescribed_oxygen_local.download(output.prescribed_oxygen);
    }
  } else {
    output.prescribed_oxygen.clear();
  }
  std::vector totals_host(3, 0.0);
  d_totals.download(totals_host);
  output.iron_uptake = totals_host[0];
  output.carbon_demand = totals_host[1];
  output.oxygen_demand = totals_host[2];

  ++metabolism_gpu_steps_;
  return true;
#endif
}

bool AgentPoolGpu::run_delivery_commit(const GpuDeliveryCommitContext& ctx) {
#ifndef GUTIBM_CUDA
  (void)ctx;
  return false;
#else
  if (!gpu_runtime_enabled() || size_ <= 0 || ctx.dt <= 0.0) return false;
  const Int ncells = ctx.chem.ncells();
  if (ncells <= 0) return false;

  const Int carbon = ctx.chem.find(species::CARBON);
  const Int oxygen = ctx.chem.find(species::OXYGEN);
  const bool o2_delivery =
      (oxygen >= 0 && ctx.oxygen.delivery_uptake_enabled);

  DeviceBuffer<int> d_support_offsets;
  DeviceBuffer<int> d_support_cells;
  if (ctx.support.active()) {
    d_support_offsets.upload(ctx.support.offsets);
    if (!ctx.support.cells.empty()) {
      d_support_cells.upload(ctx.support.cells);
    }
  }
  const int* support_offsets_ptr =
      ctx.support.active() ? d_support_offsets.data() : nullptr;
  const int* support_cells_ptr =
      (ctx.support.active() && !ctx.support.cells.empty())
          ? d_support_cells.data()
          : nullptr;

  const auto cell_count = static_cast<size_t>(ncells);
  DeviceBuffer<double> d_demand_carbon;
  DeviceBuffer<double> d_prescribed_carbon;
  const double* prescribed_carbon_ptr = nullptr;
  if (carbon >= 0) {
    d_demand_carbon.allocate(cell_count);
    d_demand_carbon.zero();
    gpu::launch_delivery_demand_by_cell_kernel(
        d_pending_carbon_funding_.data(), d_grid_cell_.data(), d_state_.data(),
        d_is_ghost_.data(), support_offsets_ptr, support_cells_ptr,
        d_demand_carbon.data(), size_, ctx.domain.nx(), ctx.domain.ny(),
        ctx.storage_nx, ctx.owned_global_x_begin, ctx.owned_global_x_end,
        ctx.owned_storage_x_begin, gpu_compute_stream());
    gpu_sync_compute();
    if (!ctx.chem.slab_mode() && ctx.domain.nprocs() > 1) {
      std::vector<double> demand_host;
      d_demand_carbon.download(demand_host);
      ctx.chem.sum_values_across_ranks(demand_host);
      d_demand_carbon.upload(demand_host);
    }
    prescribed_carbon_ptr =
        ctx.chem_gpu != nullptr
            ? ctx.chem_gpu->post_diffusion_prescribed_device(carbon)
            : nullptr;
    if (prescribed_carbon_ptr == nullptr) {
      d_prescribed_carbon.upload(ctx.chem.prescribed_sink_storage(carbon));
      prescribed_carbon_ptr = d_prescribed_carbon.data();
    }
  }

  DeviceBuffer<double> d_demand_oxygen;
  DeviceBuffer<double> d_prescribed_oxygen;
  const double* prescribed_oxygen_ptr = nullptr;
  if (o2_delivery) {
    d_demand_oxygen.allocate(cell_count);
    d_demand_oxygen.zero();
    gpu::launch_delivery_demand_by_cell_kernel(
        d_pending_oxygen_funding_.data(), d_grid_cell_.data(), d_state_.data(),
        d_is_ghost_.data(), support_offsets_ptr, support_cells_ptr,
        d_demand_oxygen.data(), size_, ctx.domain.nx(), ctx.domain.ny(),
        ctx.storage_nx, ctx.owned_global_x_begin, ctx.owned_global_x_end,
        ctx.owned_storage_x_begin, gpu_compute_stream());
    gpu_sync_compute();
    if (!ctx.chem.slab_mode() && ctx.domain.nprocs() > 1) {
      std::vector<double> demand_host;
      d_demand_oxygen.download(demand_host);
      ctx.chem.sum_values_across_ranks(demand_host);
      d_demand_oxygen.upload(demand_host);
    }
    prescribed_oxygen_ptr =
        ctx.chem_gpu != nullptr
            ? ctx.chem_gpu->post_diffusion_prescribed_device(oxygen)
            : nullptr;
    if (prescribed_oxygen_ptr == nullptr) {
      d_prescribed_oxygen.upload(ctx.chem.prescribed_sink_storage(oxygen));
      prescribed_oxygen_ptr = d_prescribed_oxygen.data();
    }
  }

  DeviceBuffer<double> d_carbon_ledger(6);
  d_carbon_ledger.zero();
  DeviceBuffer<double> d_oxygen_ledger(6);
  d_oxygen_ledger.zero();

  gpu::DeliveryCommitParams params;
  params.num_agents = size_;
  params.dt = ctx.dt;
  params.cell_density = CELL_DENSITY_DEFAULT;
  params.o2_delivery_enabled = o2_delivery ? 1 : 0;
  params.metabolic_switch_enabled =
      ctx.oxygen.metabolic_switch_enabled ? 1 : 0;
  params.respiration_driver_mode =
      ctx.oxygen.respiration_driver_mode == RespirationDriver::Funded ? 1 : 0;
  params.tau_metabolic_switch = ctx.oxygen.tau_metabolic_switch;
  params.global_nx = ctx.domain.nx();
  params.global_ny = ctx.domain.ny();
  params.storage_nx = ctx.storage_nx;
  params.owned_global_x_begin = ctx.owned_global_x_begin;
  params.owned_global_x_end = ctx.owned_global_x_end;
  params.owned_storage_x_begin = ctx.owned_storage_x_begin;

  gpu::DeliveryCommitDeviceBuffers dev;
  dev.prescribed_carbon = prescribed_carbon_ptr;
  dev.demand_carbon_by_cell = carbon >= 0 ? d_demand_carbon.data() : nullptr;
  dev.prescribed_oxygen = prescribed_oxygen_ptr;
  dev.demand_oxygen_by_cell = o2_delivery ? d_demand_oxygen.data() : nullptr;
  dev.mu_realized = d_mu_realized_.data();
  dev.biomass = d_biomass_.data();
  dev.radius = d_radius_.data();
  dev.mass = d_mass_.data();
  dev.age = d_age_.data();
  dev.fermentation_fraction = d_fermentation_fraction_.data();
  dev.respired_oxygen_rate = d_respired_oxygen_rate_.data();
  dev.pending_growth_carbon = d_pending_growth_carbon_.data();
  dev.pending_maintenance_carbon = d_pending_maintenance_carbon_.data();
  dev.pending_carbon_funding = d_pending_carbon_funding_.data();
  dev.pending_biomass = d_pending_biomass_.data();
  dev.pending_growth_chemistry_biomass =
      d_pending_growth_chemistry_biomass_.data();
  dev.pending_oxygen_growth = d_pending_oxygen_growth_.data();
  dev.pending_oxygen_maintenance = d_pending_oxygen_maintenance_.data();
  dev.pending_oxygen_funding = d_pending_oxygen_funding_.data();
  dev.grid_cell = d_grid_cell_.data();
  dev.state = d_state_.data();
  dev.is_ghost = d_is_ghost_.data();
  dev.support_offsets = support_offsets_ptr;
  dev.support_cells = support_cells_ptr;
  dev.carbon_ledger = d_carbon_ledger.data();
  dev.oxygen_ledger = d_oxygen_ledger.data();

  gpu::launch_delivery_commit_kernel(params, dev, gpu_compute_stream());
  gpu_sync_compute();
  gpu_check_error("delivery_commit_kernel");

  std::vector<double> carbon_ledger;
  d_carbon_ledger.download(carbon_ledger);
  record_device_delivery_ledger(ctx.chem, carbon, carbon_ledger);
  if (o2_delivery) {
    std::vector<double> oxygen_ledger;
    d_oxygen_ledger.download(oxygen_ledger);
    record_device_delivery_ledger(ctx.chem, oxygen, oxygen_ledger);
  }

  ++delivery_commit_gpu_steps_;
  return true;
#endif
}

}  // namespace gutibm
