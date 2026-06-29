#include "agent_pool_gpu.h"
#include "agent.h"
#include "domain.h"
#include "fix_metabolism.h"
#include "gpu_kernels.h"
#include "dispatch.h"

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
  d_ligand_affinity_.allocate(static_cast<size_t>(NUM_RECEPTORS) * static_cast<size_t>(n));
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
  std::vector<double> ligand_affinity(static_cast<size_t>(NUM_RECEPTORS) * n);
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
    age[i] = a.age;
    mu_max[i] = a.mu_max;
    km_b12[i] = a.km_b12;
    km_carbon[i] = a.km_carbon;
    bi_loci_count[i] = static_cast<int>(a.genome.bi_loci.size());
    plasmid_amelioration[i] = a.genome.plasmid_cost_amelioration;
    for (int r = 0; r < NUM_RECEPTORS; ++r) {
      receptor_expr[static_cast<size_t>(r) * n + i] = a.receptor_expr[r];
      ligand_affinity[static_cast<size_t>(r) * n + i] = a.genome.ligand_affinity[r];
    }
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
  d_ligand_affinity_.upload(ligand_affinity);
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
    a.age = age[i];
    a.grid_cell = grid_cell[i];
  }
}

bool AgentPoolGpu::run_metabolism(
    const Domain& domain, const MetabolismConfig& cfg,
    const double* d_conc_carbon, const double* d_conc_iron, const double* d_conc_b12,
    const double* d_conc_acetate, const double* d_conc_eut,
    double* d_reac_carbon, double* d_reac_iron, double* d_reac_b12,
    double dt) {

#ifndef GUTIBM_CUDA
  (void)domain;
  (void)cfg;
  (void)d_conc_carbon;
  (void)d_conc_iron;
  (void)d_conc_b12;
  (void)d_conc_acetate;
  (void)d_conc_eut;
  (void)d_reac_carbon;
  (void)d_reac_iron;
  (void)d_reac_b12;
  (void)dt;
  return false;
#else
  if (!gpu_runtime_enabled() || size_ <= 0) return false;

  gpu::launch_metabolism_kernel(
      d_conc_carbon, d_conc_iron, d_conc_b12, d_conc_acetate, d_conc_eut,
      d_reac_carbon, d_reac_iron, d_reac_b12,
      d_mu_realized_.data(), d_biomass_.data(), d_radius_.data(),
      d_mass_.data(), d_age_.data(),
      d_grid_cell_.data(), d_state_.data(),
      d_mu_max_.data(), d_km_b12_.data(), d_km_carbon_.data(),
      d_receptor_expr_.data(), d_ligand_affinity_.data(),
      d_bi_loci_count_.data(), d_plasmid_amelioration_.data(),
      size_, dt, domain.dx(), CELL_DENSITY_DEFAULT,
      cfg.km_iron_primary, cfg.km_iron_iroN, cfg.km_iron_iutA, cfg.km_iron_fiu,
      cfg.maintenance_rate, cfg.metE_penalty, cfg.metE_acetate_max_factor,
      cfg.metE_acetate_km, cfg.eut_max_penalty, cfg.eut_km,
      cfg.yield_carbon, cfg.yield_iron, cfg.yield_b12, nullptr);

  cudaDeviceSynchronize();
  gpu_check_error("metabolism_kernel");
  return true;
#endif
}

}  // namespace gutibm
