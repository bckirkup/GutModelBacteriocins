#ifndef GUTIBM_AGENT_POOL_GPU_H
#define GUTIBM_AGENT_POOL_GPU_H

#include "types.h"
#include "device_memory.h"
#include <vector>

namespace gutibm {

class Agent;
class AgentPool;
class Domain;
struct MetabolismConfig;

struct GpuMetabolismBuffers {
  const double* d_conc_carbon = nullptr;
  const double* d_conc_iron = nullptr;
  const double* d_conc_b12 = nullptr;
  const double* d_conc_acetate = nullptr;
  const double* d_conc_eut = nullptr;
  const double* d_conc_oxygen = nullptr;
  double* d_reac_carbon = nullptr;
  double* d_reac_iron = nullptr;
  double* d_reac_b12 = nullptr;
  double* d_reac_acetate = nullptr;
  int iron_uptake_enabled = 1;
  int b12_uptake_enabled = 1;
  int eut_enabled = 1;
  int fur_enabled = 0;
  double fur_Km = 0.0;
  double fur_upregulation_max = 0.0;
  double fur_receptor_max = 0.0;
  int acetate_enabled = 0;
  double acetate_overflow_threshold = 0.0;
  double acetate_overflow_rate = 0.0;
  double acetate_scavenge_rate = 0.0;
  double acetate_scavenge_Km = 0.0;
  int o2_enabled = 0;
  double o2_boost_max = 0.0;
  double o2_Km = 0.0;
  int metabolic_switch_enabled = 0;
  double mu_crit = 0.0;
  double aerobic_mu_factor = 1.0;
  double anaerobic_mu_factor = 1.0;
  double aerobic_carbon_cost_factor = 1.0;
  double anaerobic_carbon_cost_factor = 1.0;
  double tau_metabolic_switch = 0.0;
  double ferm_acid_yield = 0.0;
  double anaerobic_maintenance_factor = 1.0;
  int acid_inhibition_enabled = 0;
  double acid_inhibition_max = 0.0;
  double acid_inhibition_Ki = 0.0;
  double acetate_pKa = 0.0;
  double environment_pH = 7.0;
  int global_nx = 0;
  int global_ny = 0;
  int storage_nx = 0;
  int owned_global_x_begin = 0;
  int owned_global_x_end = 0;
  int owned_storage_x_begin = 0;
  int receptor_count = 0;
  double effective_diffusivity_carbon = 0.0;
  double effective_diffusivity_iron = 0.0;
  double* d_uptake_limit_totals = nullptr;
  double* d_maintenance_available = nullptr;
};

class AgentPoolGpu {
 public:
  void resize(Int n);
  void sync_from_host(const AgentPool& pool);
  void sync_to_host(AgentPool& pool) const;
  void sync_receptor_expression_to_host(AgentPool& pool) const;

  Int size() const { return size_; }

  double* x() { return d_x_.data(); }
  double* y() { return d_y_.data(); }
  double* z() { return d_z_.data(); }
  const double* x() const { return d_x_.data(); }
  const double* y() const { return d_y_.data(); }
  const double* z() const { return d_z_.data(); }
  int*    grid_cell() { return d_grid_cell_.data(); }
  const int* grid_cell() const { return d_grid_cell_.data(); }
  int*    state() { return d_state_.data(); }
  const int* state() const { return d_state_.data(); }
  double* mu_realized() { return d_mu_realized_.data(); }
  const double* mu_realized() const { return d_mu_realized_.data(); }
  double* fermentation_fraction() { return d_fermentation_fraction_.data(); }
  const double* fermentation_fraction() const {
    return d_fermentation_fraction_.data();
  }
  double* biomass() { return d_biomass_.data(); }
  const double* biomass() const { return d_biomass_.data(); }
  double* radius() { return d_radius_.data(); }
  double* mass() { return d_mass_.data(); }
  double* age() { return d_age_.data(); }
  double* mu_max() { return d_mu_max_.data(); }
  double* km_b12() { return d_km_b12_.data(); }
  double* km_carbon() { return d_km_carbon_.data(); }
  double* receptor_expr() { return d_receptor_expr_.data(); }
  const double* receptor_expr() const { return d_receptor_expr_.data(); }
  const double* receptor_expr_base() const { return d_receptor_expr_base_.data(); }
  double* ligand_affinity() { return d_ligand_affinity_.data(); }
  const double* ligand_affinity() const { return d_ligand_affinity_.data(); }
  const int* iron_receptor() const { return d_iron_receptor_.data(); }
  int*    bi_loci_count() { return d_bi_loci_count_.data(); }
  double* plasmid_amelioration() { return d_plasmid_amelioration_.data(); }

  bool run_metabolism(const Domain& domain, const MetabolismConfig& cfg,
                      const GpuMetabolismBuffers& buffers, double* uptake_totals,
                      double dt, Int num_agents);
  Int metabolism_gpu_steps() const { return metabolism_gpu_steps_; }

  void sync_positions_to_host(AgentPool& pool) const;

 private:
  Int size_ = 0;
  Int metabolism_gpu_steps_ = 0;
  DeviceBuffer<double> d_x_;
  DeviceBuffer<double> d_y_;
  DeviceBuffer<double> d_z_;
  DeviceBuffer<int> d_grid_cell_;
  DeviceBuffer<int> d_state_;
  DeviceBuffer<int> d_bi_loci_count_;
  DeviceBuffer<double> d_mu_realized_;
  DeviceBuffer<double> d_fermentation_fraction_;
  DeviceBuffer<double> d_biomass_;
  DeviceBuffer<double> d_radius_;
  DeviceBuffer<double> d_mass_;
  DeviceBuffer<double> d_age_;
  DeviceBuffer<double> d_mu_max_;
  DeviceBuffer<double> d_km_b12_;
  DeviceBuffer<double> d_km_carbon_;
  DeviceBuffer<double> d_receptor_expr_;
  DeviceBuffer<double> d_receptor_expr_base_;
  DeviceBuffer<double> d_ligand_affinity_;
  DeviceBuffer<int> d_iron_receptor_;
  DeviceBuffer<double> d_plasmid_amelioration_;
};

}  // namespace gutibm

#endif  // GUTIBM_AGENT_POOL_GPU_H
