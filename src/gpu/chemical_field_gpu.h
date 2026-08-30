#ifndef GUTIBM_CHEMICAL_FIELD_GPU_H
#define GUTIBM_CHEMICAL_FIELD_GPU_H

#include "types.h"
#include "device_memory.h"
#include <memory>
#include <vector>

namespace gutibm {

class ChemicalField;
class Domain;
struct ChemicalSpec;

class ChemicalFieldGpu {
 public:
  void init(ChemicalField& field);
  void sync_to_device(const ChemicalField& field);
  void sync_to_host(ChemicalField& field);
  void sync_concentrations_to_device(const ChemicalField& field);
  void sync_reactions_to_device(const ChemicalField& field);
  void sync_concentrations_to_host(ChemicalField& field);
  void sync_reactions_to_host(ChemicalField& field);
  void accumulate_reactions_to_host(ChemicalField& field);
  void sync_species_concentrations_to_host(ChemicalField& field, Int spec);
  void sync_species_concentrations_to_device(const ChemicalField& field, Int spec);
  void zero_reactions_on_device();

  bool apply_reactions(double dt, const Domain& domain);
  void reset_reaction_clip();
  void download_reaction_clip(ChemicalField& field) const;
  bool apply_diffusion(const Domain& domain, ChemicalField& field, Real dt);
  bool apply_delivery_species(const Domain& domain, const ChemicalSpec& spec,
                              Int species, Real dt, bool prescribed_active);
  void reset_delivery_boundary(Int species);
  Real download_delivery_boundary(Int species) const;
  bool delivery_has_negative(Int spec);
  Real delivery_negative_fraction(Int spec);
  void prepare_delivery_species(
      Int spec, const std::vector<Real>& sink,
      const std::vector<Real>& prescribed);
  void snapshot_delivery_species(Int spec);
  void restore_delivery_species(Int spec);
  void restore_delivery_original(Int spec);
  void download_delivery_species(Int spec, std::vector<Real>& concentration,
                                 std::vector<Real>& realized) const;
  void upload_delivery_concentration(Int spec,
                                     const std::vector<Real>& concentration);
  void upload_delivery_prescribed(const std::vector<Real>& prescribed);
  Real download_delivery_gradient_source() const;
  bool apply_boundaries(const Domain& domain, ChemicalField& field);
  void reset_agent_uptake();
  void download_agent_uptake(ChemicalField& field) const;
  void prepare_maintenance(const ChemicalField& field, Int carbon,
                           Real cell_volume);
  void reset_uptake_limit_totals();
  void download_uptake_limit_totals(ChemicalField& field) const;
  double* uptake_limit_totals_device() {
    return active_ ? d_uptake_limit_totals_.data() : nullptr;
  }
  void reset_vbf_totals();
  double* vbf_totals_device();
  void download_vbf_totals(std::vector<double>& values) const;
  void reset_agent_counts();
  int* agent_counts_device();
  double* agent_uptake_device() {
    return active_ ? d_agent_uptake_.data() : nullptr;
  }
  double* maintenance_available_device() {
    return active_ ? d_maintenance_available_.data() : nullptr;
  }

  bool try_sum_reactions_on_device(ChemicalField& field);

  double* conc_device(Int spec);
  const double* conc_device(Int spec) const;
  double* reac_device(Int spec);
  const double* reac_device(Int spec) const;

  bool active() const { return active_; }
  Int storage_nx() const { return storage_nx_; }
  Int ncells() const { return ncells_; }
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
  bool diffusion_fallback_warning_emitted_ = false;
  std::vector<DeviceBuffer<double>> d_conc_;
  std::vector<DeviceBuffer<double>> d_reac_;
  DeviceBuffer<double> d_boundary_conc_;
  DeviceBuffer<double> d_boundary_injected_;
  DeviceBuffer<double> d_agent_uptake_;
  DeviceBuffer<double> d_maintenance_available_;
  DeviceBuffer<double> d_uptake_limit_totals_;
  DeviceBuffer<double> d_vbf_totals_;
  DeviceBuffer<double> d_reaction_clip_;
  DeviceBuffer<int> d_agent_counts_;
  DeviceBuffer<double> d_delivery_sink_;
  DeviceBuffer<double> d_delivery_prescribed_;
  DeviceBuffer<double> d_delivery_realized_;
  DeviceBuffer<double> d_delivery_concentration_backup_;
  DeviceBuffer<double> d_delivery_realized_backup_;
  DeviceBuffer<double> d_delivery_gradient_source_;
  DeviceBuffer<unsigned long long> d_delivery_negative_count_;
  unsigned long long delivery_negative_count_host_ = 0;
  Int delivery_negative_count_spec_ = -1;
};

}  // namespace gutibm

#endif  // GUTIBM_CHEMICAL_FIELD_GPU_H
