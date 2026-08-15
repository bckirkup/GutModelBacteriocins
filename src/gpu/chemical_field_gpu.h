#ifndef GUTIBM_CHEMICAL_FIELD_GPU_H
#define GUTIBM_CHEMICAL_FIELD_GPU_H

#include "types.h"
#include "device_memory.h"
#include <memory>
#include <vector>

namespace gutibm {

class ChemicalField;
class Domain;

class ChemicalFieldGpu {
 public:
  void init(ChemicalField& field);
  void sync_to_device(const ChemicalField& field);
  void sync_to_host(ChemicalField& field);
  void sync_concentrations_to_device(const ChemicalField& field);
  void sync_reactions_to_device(const ChemicalField& field);
  void sync_concentrations_to_host(ChemicalField& field);
  void sync_reactions_to_host(ChemicalField& field);
  void sync_species_concentrations_to_host(ChemicalField& field, Int spec);
  void sync_species_concentrations_to_device(const ChemicalField& field, Int spec);
  void zero_reactions_on_device();

  bool apply_reactions(double dt, const Domain& domain);
  void reset_reaction_clip();
  void download_reaction_clip(ChemicalField& field);
  bool apply_diffusion(const Domain& domain, ChemicalField& field, Real dt);
  bool apply_boundaries(const Domain& domain, ChemicalField& field);
  void reset_agent_uptake();
  void download_agent_uptake(ChemicalField& field);
  void reset_vbf_totals();
  double* vbf_totals_device();
  void download_vbf_totals(std::vector<double>& values);
  double* agent_uptake_device() {
    return active_ ? d_agent_uptake_.data() : nullptr;
  }

  bool try_sum_reactions_on_device(ChemicalField& field);

  double* conc_device(Int spec);
  const double* conc_device(Int spec) const;
  double* reac_device(Int spec);
  const double* reac_device(Int spec) const;

  bool active() const { return active_; }
  Int storage_nx() const { return storage_nx_; }
  Int global_nx() const { return global_nx_; }
  Int global_ny() const { return global_ny_; }
  Int global_nz() const { return global_nz_; }
  Int owned_x_begin() const { return owned_x_begin_; }
  Int owned_x_end() const { return owned_x_end_; }
  Int halo_width() const { return halo_width_; }
  Int owned_storage_x_begin() const { return slab_mode_ ? halo_width_ : 0; }
  Int owned_storage_x_end() const {
    return slab_mode_ ? halo_width_ + (owned_x_end_ - owned_x_begin_)
                      : global_nx_;
  }
  bool slab_mode() const { return slab_mode_; }

 private:
  bool active_ = false;
  Int nspec_ = 0;
  Int ncells_ = 0;
  Int global_nx_ = 0;
  Int global_ny_ = 0;
  Int global_nz_ = 0;
  Int owned_x_begin_ = 0;
  Int owned_x_end_ = 0;
  Int halo_width_ = 0;
  Int storage_nx_ = 0;
  bool slab_mode_ = false;
  std::vector<DeviceBuffer<double>> d_conc_;
  std::vector<DeviceBuffer<double>> d_reac_;
  DeviceBuffer<double> d_boundary_conc_;
  DeviceBuffer<double> d_boundary_injected_;
  DeviceBuffer<double> d_agent_uptake_;
  DeviceBuffer<double> d_vbf_totals_;
  DeviceBuffer<double> d_reaction_clip_;
};

}  // namespace gutibm

#endif  // GUTIBM_CHEMICAL_FIELD_GPU_H
