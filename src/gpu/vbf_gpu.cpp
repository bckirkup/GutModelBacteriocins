#include "vbf_gpu.h"
#include "chemical_field.h"
#include "chemical_field_gpu.h"
#include "chem_environment_config.h"
#include "dispatch.h"
#include "domain.h"
#include "gpu_kernels.h"
#include "species_names.h"
#include "vbf.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

namespace gutibm {

bool gpu_apply_vbf_coupling(ChemicalFieldGpu& chem_gpu,
                            const ChemicalField& chem,
                            const Domain& domain,
                            const VBF& vbf,
                            const OxygenConfig& oxygen,
                            const AcetateConfig& acetate,
                            const MucinConfig& mucin) {
#ifndef GUTIBM_CUDA
  (void)chem_gpu;
  (void)chem;
  (void)domain;
  (void)vbf;
  (void)oxygen;
  (void)acetate;
  (void)mucin;
  return false;
#else
  if (!gpu_runtime_enabled() || !chem_gpu.active()) return false;

  const VBFConfig& cfg = vbf.config();
  gpu::VbfLaunchParams params;
  params.nx = domain.nx();
  params.ny = domain.ny();
  params.nz = domain.nz();
  params.dx = domain.dx();
  params.nutrient_sink = cfg.nutrient_sink;
  params.carbon_sink_vmax = cfg.carbon_sink_vmax;
  params.carbon_sink_km = cfg.carbon_sink_km;
  params.use_dynamic_mucin = cfg.use_dynamic_mucin ? 1 : 0;
  params.mucin_z_gradient_enabled = cfg.mucin_z_gradient_enabled ? 1 : 0;
  params.mucin_z_gradient_lambda = cfg.mucin_z_gradient_lambda;
  params.mucin_liberation = cfg.mucin_liberation;
  params.vbf_density = cfg.density;
  params.oxygen_enabled = oxygen.enabled ? 1 : 0;
  params.oxygen_vbf_sink = oxygen.vbf_sink;
  params.acetate_enabled = acetate.enabled ? 1 : 0;
  params.acetate_vbf_production = acetate.vbf_production;
  params.acetate_vbf_consumption = acetate.vbf_consumption;
  params.acetate_epithelial_uptake = acetate.epithelial_uptake;
  params.mucin_enabled = mucin.enabled ? 1 : 0;
  params.mucin_secretion_rate = mucin.secretion_rate;
  params.mucin_Km_degradation = mucin.Km_degradation;
  params.mucin_k_liberation = mucin.k_liberation;

  const Int i_carbon = chem.find(species::CARBON);
  const Int i_iron = chem.find(species::IRON);
  const Int i_oxygen = chem.find(species::OXYGEN);
  const Int i_acetate = chem.find(species::ACETATE);
  const Int i_mucin = chem.find(species::MUCIN);

  gpu::launch_vbf_coupling_kernel(
      domain.ncells(), params,
      i_carbon >= 0 ? chem_gpu.reac_device(i_carbon) : nullptr,
      i_carbon >= 0 ? chem_gpu.conc_device(i_carbon) : nullptr,
      i_iron >= 0 ? chem_gpu.reac_device(i_iron) : nullptr,
      i_iron >= 0 ? chem_gpu.conc_device(i_iron) : nullptr,
      i_oxygen >= 0 ? chem_gpu.reac_device(i_oxygen) : nullptr,
      i_oxygen >= 0 ? chem_gpu.conc_device(i_oxygen) : nullptr,
      i_acetate >= 0 ? chem_gpu.reac_device(i_acetate) : nullptr,
      i_mucin >= 0 ? chem_gpu.reac_device(i_mucin) : nullptr,
      i_mucin >= 0 ? chem_gpu.conc_device(i_mucin) : nullptr,
      gpu_compute_stream());

  gpu_sync_compute();
  gpu_check_error("gpu_apply_vbf_coupling");
  return true;
#endif
}

}  // namespace gutibm
