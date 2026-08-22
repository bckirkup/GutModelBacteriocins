#include "chemical_field.h"
#include "chem_environment_config.h"
#include "domain.h"
#include "chemical_field_gpu.h"
#include "device.h"
#include "dispatch.h"
#include "species_names.h"
#include "vbf.h"
#include "vbf_gpu.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

using namespace gutibm;

namespace {

struct SinkMeasurement {
  Real amount = 0.0;
  Real updated_concentration = 0.0;
};

SinkMeasurement measure_sink(const Domain& domain, const ChemicalSpec& carbon,
                             const VBFConfig& vbf_cfg, Real dt,
                             const std::vector<Int>& agent_counts = {}) {
  ChemicalField chem;
  chem.init(domain, {carbon});
  const Int carbon_index = chem.find(species::CARBON);
  VBF vbf;
  vbf.init(vbf_cfg, domain);
  OxygenConfig oxygen;
  AcetateConfig acetate;
  MucinConfig mucin;
  VbfFluxTotals totals;
  vbf.apply_nutrient_coupling(
      chem, domain, dt, oxygen, acetate, mucin, &totals, agent_counts);

  const Real cell_volume = domain.cell_volume();
  const Real updated = chem.conc(carbon_index, 0)
      + chem.reac(carbon_index, 0) * dt;
  return {
      -chem.reac(carbon_index, 0) * cell_volume * dt,
      updated,
  };
}

void test_implicit_sink_mass_closure(const Domain& domain,
                                     const ChemicalSpec& carbon) {
  VBFConfig vbf_cfg;
  vbf_cfg.mucin_liberation = 5.0e-5;
  vbf_cfg.carbon_sink_vmax = 5.0e-3;
  vbf_cfg.carbon_sink_km = 1.0e-4;
  vbf_cfg.nutrient_sink = 0.0;
  VBF vbf;
  vbf.init(vbf_cfg, domain);
  OxygenConfig oxygen;
  AcetateConfig acetate;
  MucinConfig mucin;
  ChemicalField chem;
  ChemicalSpec initial = carbon;
  initial.initial_conc = 2.0e-4;
  chem.init(domain, {initial});
  const Int carbon_index = chem.find(species::CARBON);
  const Real cell_volume = domain.cell_volume();
  constexpr Real dt = 60.0;
  constexpr int steps = 5;
  const Real before = chem.conc(carbon_index, 0)
      * static_cast<Real>(chem.ncells()) * cell_volume;
  Real source = 0.0;
  Real sink = 0.0;
  for (int step = 0; step < steps; ++step) {
    chem.zero_reactions();
    VbfFluxTotals totals;
    vbf.apply_nutrient_coupling(
        chem, domain, dt, oxygen, acetate, mucin, &totals);
    source += totals.carbon_source;
    sink += totals.carbon_sink;
    for (Int cell = 0; cell < chem.ncells(); ++cell) {
      const Real updated =
          chem.conc(carbon_index, cell) + chem.reac(carbon_index, cell) * dt;
      if (updated < 0.0) {
        chem.flux_accounting().add_reaction_clip(
            carbon_index, -updated * cell_volume);
      }
      chem.conc(carbon_index, cell) = std::max(updated, 0.0);
    }
    chem.flux_accounting().commit_boundary_and_reaction_step();
  }
  const Real after = chem.conc(carbon_index, 0)
      * static_cast<Real>(chem.ncells()) * cell_volume;
  const Real clip = chem.flux_accounting().reaction_clip_interval[
      static_cast<size_t>(carbon_index)];
  const Real scale = std::abs(before) + std::abs(after)
      + std::abs(source) + std::abs(sink) + std::abs(clip);
  assert(std::abs((after - before) - (source - sink + clip))
         < 1.0e-12 * scale);
}

}  // namespace

int main() {
  DomainConfig domain_cfg;
  domain_cfg.hi = {10e-6, 10e-6, 10e-6};
  domain_cfg.grid_dx = 5e-6;
  domain_cfg.chemistry_stride = {2, 2, 1};
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

  const Real cell_volume = domain.cell_volume();
  Real applied_amount = 0.0;
  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    applied_amount += -chem.reac(carbon_index, cell) * cell_volume * dt;
  }
  assert(applied_amount > 0.0);
  assert(std::abs(applied_amount - totals.carbon_sink)
         < 1.0e-12 * applied_amount);
  assert(applied_amount < 1.0e-4 * static_cast<Real>(chem.ncells())
         * cell_volume);

  const std::array<Real, 4> vmax_values = {
      2.0e-7, 2.0e-6, 2.0e-4, 2.0e-2};
  std::array<SinkMeasurement, 4> measurements{};
  for (size_t i = 0; i < vmax_values.size(); ++i) {
    VBFConfig sweep_cfg = vbf_cfg;
    sweep_cfg.carbon_sink_vmax = vmax_values[i];
    measurements[i] = measure_sink(domain, carbon, sweep_cfg, dt);
    assert(measurements[i].amount <= carbon.initial_conc * cell_volume);
    assert(measurements[i].updated_concentration > 0.0);
    if (i > 0) {
      assert(measurements[i].amount > measurements[i - 1].amount);
    }
  }
  test_implicit_sink_mass_closure(domain, carbon);
  std::vector<Int> empty_counts(static_cast<size_t>(domain.ncells()), 0);
  std::vector<Int> dense_counts = empty_counts;
  dense_counts[0] = 1;
  const std::array<Real, 4> coupling_values = {
      0.0, 1.0e-21, 1.0e-20, 1.0e-19};
  std::array<Real, 4> dense_empty_gaps{};
  for (size_t i = 0; i < coupling_values.size(); ++i) {
    VBFConfig coupling_cfg = vbf_cfg;
    coupling_cfg.agent_carbon_coupling = coupling_values[i];
    const auto empty = measure_sink(
        domain, carbon, coupling_cfg, dt, empty_counts);
    const auto dense = measure_sink(
        domain, carbon, coupling_cfg, dt, dense_counts);
    dense_empty_gaps[i] =
        empty.updated_concentration - dense.updated_concentration;
    assert(i == 0
           || dense.updated_concentration < empty.updated_concentration);
    if (i > 0) {
      assert(dense_empty_gaps[i] >= dense_empty_gaps[i - 1]);
    }
  }
  std::cout << "PASS: agent-density VBF coupling is monotone\n";
  NutrientFluxAccounting flux;
  flux.init(1);
  flux.agent_uptake_interval[0] = 2.0;
  flux.refresh_nutrient_blocking_fraction();
  assert(std::abs(flux.nutrient_blocking_fraction[0] - 1.0) < 1.0e-15);
  flux.vbf_sink_interval[0] = 2.0;
  flux.refresh_nutrient_blocking_fraction();
  assert(std::abs(flux.nutrient_blocking_fraction[0] - 0.5) < 1.0e-15);
  flux.agent_uptake_interval[0] = 0.0;
  flux.vbf_sink_interval[0] = 0.0;
  flux.refresh_nutrient_blocking_fraction();
  assert(std::abs(flux.nutrient_blocking_fraction[0]) < 1.0e-15);
  std::cout << "PASS: nutrient blocking fraction tracks agent/VBF uptake\n";
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
