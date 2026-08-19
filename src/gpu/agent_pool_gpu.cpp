#include "agent_pool_gpu.h"
#include "agent.h"
#include "domain.h"
#include "fix_metabolism.h"
#include "gpu_kernels.h"
#include "dispatch.h"
#include "receptor_utils.h"

#include <cassert>

namespace gutibm {

void AgentPoolGpu::resize(Int n) {
  size_ = n;
  if (n <= 0) return;
  d_x_.allocate(static_cast<size_t>(n));
  d_y_.allocate(static_cast<size_t>(n));
  d_z_.allocate(static_cast<size_t>(n));
  d_grid_cell_.allocate(static_cast<size_t>(n));
  d_state_.allocate(static_cast<size_t>(n));
  d_mu_realized_.allocate(static_cast<size_t>(n));
  d_biomass_.allocate(static_cast<size_t>(n));
  d_radius_.allocate(static_cast<size_t>(n));
  d_mass_.allocate(static_cast<size_t>(n));
  d_age_.allocate(static_cast<size_t>(n));
  d_mu_max_.allocate(static_cast<size_t>(n));
  d_km_b12_.allocate(static_cast<size_t>(n));
  d_km_carbon_.allocate(static_cast<size_t>(n));
  d_receptor_expr_.allocate(static_cast<size_t>(NUM_RECEPTORS) * static_cast<size_t>(n));
  d_receptor_expr_base_.allocate(
      static_cast<size_t>(NUM_RECEPTORS) * static_cast<size_t>(n));
  d_ligand_affinity_.allocate(static_cast<size_t>(NUM_RECEPTORS) * static_cast<size_t>(n));
  d_iron_receptor_.allocate(static_cast<size_t>(NUM_RECEPTORS));
  d_bi_loci_count_.allocate(static_cast<size_t>(n));
  d_plasmid_amelioration_.allocate(static_cast<size_t>(n));
}

void AgentPoolGpu::sync_from_host(const AgentPool& pool) {
  Int n = pool.size();
  resize(n);
  if (n <= 0) return;

  std::vector<double> x(n);
  std::vector<double> y(n);
  std::vector<double> z(n);
  std::vector<int> grid_cell(n);
  std::vector<int> state(n);
  std::vector<int> bi_loci_count(n);
  std::vector<double> mu_realized(n);
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

  for (Int i = 0; i < n; ++i) {
    const Agent& a = pool[i];
    x[i] = a.x[0];
    y[i] = a.x[1];
    z[i] = a.x[2];
    grid_cell[i] = a.grid_cell;
    state[i] = static_cast<int>(to_underlying(a.state));
    mu_realized[i] = a.mu_realized;
    biomass[i] = a.biomass;
    radius[i] = a.radius;
    mass[i] = a.mass;
    age[i] = a.timers.age;
    mu_max[i] = a.mu_max;
    km_b12[i] = a.km.km_b12;
    km_carbon[i] = a.km.km_carbon;
    bi_loci_count[i] = static_cast<int>(a.genome.bi_loci.size());
    plasmid_amelioration[i] = a.genome.plasmid_cost_amelioration;
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
  d_grid_cell_.upload(grid_cell);
  d_state_.upload(state);
  d_mu_realized_.upload(mu_realized);
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
}

void AgentPoolGpu::sync_to_host(AgentPool& pool) const {
  Int n = pool.size();
  if (n <= 0 || n != size_) return;

  std::vector<double> mu_realized(n);
  std::vector<double> biomass(n);
  std::vector<double> radius(n);
  std::vector<double> mass(n);
  std::vector<double> age(n);
  std::vector<int> grid_cell(n);
  d_mu_realized_.download(mu_realized);
  d_biomass_.download(biomass);
  d_radius_.download(radius);
  d_mass_.download(mass);
  d_age_.download(age);
  d_grid_cell_.download(grid_cell);

  for (Int i = 0; i < n; ++i) {
    Agent& a = pool[i];
    a.mu_realized = mu_realized[i];
    a.biomass = biomass[i];
    a.radius = radius[i];
    a.mass = mass[i];
    a.timers.age = age[i];
    a.grid_cell = grid_cell[i];
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
      d_grid_cell_.data(), d_state_.data(),
      d_mu_max_.data(), d_km_b12_.data(), d_km_carbon_.data(),
      d_receptor_expr_.data(), d_receptor_expr_base_.data(),
      d_ligand_affinity_.data(), d_iron_receptor_.data(),
      d_bi_loci_count_.data(), d_plasmid_amelioration_.data(),
      size_, num_agents, size_, buffers.receptor_count, dt,
      domain.cell_volume(), CELL_DENSITY_DEFAULT,
      cfg.km_iron_primary, cfg.km_iron_iroN, cfg.km_iron_iutA, cfg.km_iron_fiu,
      cfg.maintenance_rate, cfg.metE_penalty, cfg.metE_acetate_max_factor,
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
      uptake_totals,
      domain.nx(), domain.ny(), buffers.storage_nx,
      buffers.owned_global_x_begin, buffers.owned_global_x_end,
      buffers.owned_storage_x_begin, gpu_compute_stream());

  gpu_sync_compute();
  gpu_check_error("metabolism_kernel");
  ++metabolism_gpu_steps_;
  return true;
#endif
}

}  // namespace gutibm
