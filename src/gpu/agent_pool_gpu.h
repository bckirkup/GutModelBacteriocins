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
  double* d_reac_carbon = nullptr;
  double* d_reac_iron = nullptr;
  double* d_reac_b12 = nullptr;
};

class AgentPoolGpu {
 public:
  void resize(Int n);
  void sync_from_host(const AgentPool& pool);
  void sync_to_host(AgentPool& pool) const;

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
  double* biomass() { return d_biomass_.data(); }
  const double* biomass() const { return d_biomass_.data(); }
  double* radius() { return d_radius_.data(); }
  double* mass() { return d_mass_.data(); }
  double* age() { return d_age_.data(); }
  double* mu_max() { return d_mu_max_.data(); }
  double* km_b12() { return d_km_b12_.data(); }
  double* km_carbon() { return d_km_carbon_.data(); }
  double* receptor_expr() { return d_receptor_expr_.data(); }
  double* ligand_affinity() { return d_ligand_affinity_.data(); }
  int*    bi_loci_count() { return d_bi_loci_count_.data(); }
  double* plasmid_amelioration() { return d_plasmid_amelioration_.data(); }

  bool run_metabolism(const Domain& domain, const MetabolismConfig& cfg,
                      const GpuMetabolismBuffers& buffers, double dt);

 private:
  Int size_ = 0;
  DeviceBuffer<double> d_x_;
  DeviceBuffer<double> d_y_;
  DeviceBuffer<double> d_z_;
  DeviceBuffer<int> d_grid_cell_;
  DeviceBuffer<int> d_state_;
  DeviceBuffer<int> d_bi_loci_count_;
  DeviceBuffer<double> d_mu_realized_;
  DeviceBuffer<double> d_biomass_;
  DeviceBuffer<double> d_radius_;
  DeviceBuffer<double> d_mass_;
  DeviceBuffer<double> d_age_;
  DeviceBuffer<double> d_mu_max_;
  DeviceBuffer<double> d_km_b12_;
  DeviceBuffer<double> d_km_carbon_;
  DeviceBuffer<double> d_receptor_expr_;
  DeviceBuffer<double> d_ligand_affinity_;
  DeviceBuffer<double> d_plasmid_amelioration_;
};

}  // namespace gutibm

#endif  // GUTIBM_AGENT_POOL_GPU_H
