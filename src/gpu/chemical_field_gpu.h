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
  bool apply_diffusion(const Domain& domain, const ChemicalField& field, Real dt);
  bool apply_boundaries(const Domain& domain, const ChemicalField& field);

  bool try_sum_reactions_on_device(ChemicalField& field);

  double* conc_device(Int spec);
  const double* conc_device(Int spec) const;
  double* reac_device(Int spec);
  const double* reac_device(Int spec) const;

  bool active() const { return active_; }

 private:
  bool active_ = false;
  Int nspec_ = 0;
  Int ncells_ = 0;
  std::vector<DeviceBuffer<double>> d_conc_;
  std::vector<DeviceBuffer<double>> d_reac_;
  DeviceBuffer<double> d_boundary_conc_;
};

}  // namespace gutibm

#endif  // GUTIBM_CHEMICAL_FIELD_GPU_H
