#include "chemical_field.h"
#include "chem_environment_config.h"
#include "domain.h"
#include "chemical_field_gpu.h"
#include "device.h"
#include "dispatch.h"
#include "species_names.h"
#include "vbf.h"
#include "vbf_gpu.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

int main() {
  DomainConfig domain_cfg;
  domain_cfg.hi = {10e-6, 10e-6, 10e-6};
  domain_cfg.grid_dx = 5e-6;
  Domain domain;
  domain.init(domain_cfg);

  ChemicalSpec carbon;
  carbon.name = species::CARBON;
  carbon.initial_conc = 1.0e-4;
  carbon.diffusion_enabled = false;
  ChemicalField chem;
  chem.init(domain, {carbon});
  const Int carbon_index = chem.find(species::CARBON);
  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    chem.conc(carbon_index, cell) = 1.0e-4;
  }

  VBFConfig vbf_cfg;
  vbf_cfg.mucin_liberation = 0.0;
  vbf_cfg.carbon_sink_vmax = 2.0e-7;
  vbf_cfg.carbon_sink_km = 1.0e-4;
  vbf_cfg.nutrient_sink = 0.0;
  VBF vbf;
  vbf.init(vbf_cfg, domain);
  OxygenConfig oxygen;
  oxygen.enabled = false;
  AcetateConfig acetate;
  acetate.enabled = false;
  MucinConfig mucin;
  mucin.enabled = false;
  constexpr Real dt = 300.0;
  VbfFluxTotals totals;
  vbf.apply_nutrient_coupling(
      chem, domain, dt, oxygen, acetate, mucin, &totals);

  const Real cell_volume = domain.dx() * domain.dx() * domain.dx();
  Real applied_amount = 0.0;
  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    applied_amount += -chem.reac(carbon_index, cell) * cell_volume * dt;
  }
  assert(applied_amount > 0.0);
  assert(std::abs(applied_amount - totals.carbon_sink)
         < 1.0e-12 * applied_amount);
  const Real initial_concentration = 1.0e-4;
  const Real pre_update_sink_rate = vbf_cfg.carbon_sink_vmax
      * initial_concentration
      / (vbf_cfg.carbon_sink_km + initial_concentration);
  const Real pre_update_amount = pre_update_sink_rate
      * static_cast<Real>(chem.ncells()) * cell_volume * dt;
  const Real post_update_concentration = initial_concentration
      - pre_update_sink_rate * dt;
  const Real post_update_sink_rate = vbf_cfg.carbon_sink_vmax
      * post_update_concentration
      / (vbf_cfg.carbon_sink_km + post_update_concentration);
  const Real post_update_amount = post_update_sink_rate
      * static_cast<Real>(chem.ncells()) * cell_volume * dt;
  constexpr Real kMinimumDistinguishableFraction = 0.05;
  assert(std::abs(pre_update_amount - post_update_amount)
         > kMinimumDistinguishableFraction * applied_amount);
  assert(std::abs(pre_update_amount - applied_amount)
         < 1.0e-12 * applied_amount);
#ifdef GUTIBM_CUDA
  GpuConfig gpu_cfg;
  gpu_cfg.enabled = true;
  gpu_cfg.device_id = 0;
  gpu_set_config(gpu_cfg);
  if (gpu_init_for_rank(0, 1)) {
    ChemicalField gpu_reference;
    gpu_reference.init(domain, {carbon});
    for (Int cell = 0; cell < gpu_reference.ncells(); ++cell) {
      gpu_reference.conc(carbon_index, cell) = 1.0e-4;
    }
    ChemicalFieldGpu gpu_field;
    gpu_field.init(gpu_reference);
    gpu_field.sync_to_device(gpu_reference);
    gpu_field.zero_reactions_on_device();
    VbfFluxTotals gpu_totals;
    assert(gpu_apply_vbf_coupling(
        gpu_field, gpu_reference, domain, vbf, oxygen, acetate, mucin,
        gpu_totals, dt));
    gpu_field.sync_reactions_to_host(gpu_reference);
    Real gpu_applied_amount = 0.0;
    for (Int cell = 0; cell < gpu_reference.ncells(); ++cell) {
      gpu_applied_amount += -gpu_reference.reac(carbon_index, cell)
          * cell_volume * dt;
    }
    assert(std::abs(gpu_applied_amount - gpu_totals.carbon_sink)
           < 1.0e-12 * gpu_applied_amount);
    assert(std::abs(gpu_totals.carbon_sink - totals.carbon_sink)
           < 1.0e-12 * totals.carbon_sink);
  }
#endif
  std::cout << "PASS: VBF reported carbon sink equals applied reaction integral\n";
  return 0;
}
