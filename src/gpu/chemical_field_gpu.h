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
  void sync_to_device(ChemicalField& field);
  void sync_to_host(ChemicalField& field);
  void zero_reactions_on_device();

  bool apply_reactions(double dt, const Domain& domain);
  bool apply_boundaries(const Domain& domain, const ChemicalField& field);

  double* conc_device(Int spec);
  double* reac_device(Int spec);

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
