/* -----------------------------------------------------------------------
   GutIBM – ChemicalFieldGpu facade parity (Spec 9 PR2)
   Compares CPU ChemicalField paths with ChemicalFieldGpu device facades.
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "chemical_field_gpu.h"
#include "dispatch.h"
#include "domain.h"
#include "species_names.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

namespace {

Domain make_domain(Int nx, Int ny, Int nz) {
  DomainConfig cfg;
  cfg.lo = {0.0, 0.0, 0.0};
  cfg.hi = {nx * 5.0e-6, ny * 5.0e-6, nz * 5.0e-6};
  cfg.grid_dx = 5.0e-6;
  cfg.hash_cell_size = 10.0e-6;
  Domain domain;
  domain.init(cfg);
  return domain;
}

ChemicalSpec diffusing_species(Real diffusion, Real initial, Real boundary) {
  ChemicalSpec spec;
  spec.name = species::OXYGEN;
  spec.diff_coeff = diffusion;
  spec.retardation = 1.0;
  spec.initial_conc = initial;
  spec.boundary_conc = boundary;
  spec.diffusion_enabled = true;
  return spec;
}

ChemicalSpec static_species(Real initial, Real boundary) {
  ChemicalSpec spec;
  spec.name = species::CARBON;
  spec.diff_coeff = 0.0;
  spec.retardation = 1.0;
  spec.initial_conc = initial;
  spec.boundary_conc = boundary;
  spec.diffusion_enabled = false;
  return spec;
}

Real max_abs_diff(const std::vector<Real>& a, const std::vector<Real>& b) {
  Real max_diff = 0.0;
  const size_t n = std::min(a.size(), b.size());
  for (size_t i = 0; i < n; ++i) {
    max_diff = std::max(max_diff, std::abs(a[i] - b[i]));
  }
  return max_diff;
}

void perturb_interior(ChemicalField& field, Int spec, const Domain& domain) {
  field.conc(spec, domain.cell_index(1, 1, 2)) += 0.05;
  field.conc(spec, domain.cell_index(2, 2, 3)) -= 0.03;
}

}  // namespace

void test_gpu_facade_diffusion_matches_cpu() {
  Domain domain = make_domain(5, 4, 6);
  ChemicalField chem_cpu;
  chem_cpu.init(domain, {diffusing_species(2.1e-9, 0.25, 0.25)});
  perturb_interior(chem_cpu, 0, domain);

  ChemicalField chem_gpu_ref = chem_cpu;
  ChemicalFieldGpu chem_gpu;
  chem_gpu.init(chem_gpu_ref);

  chem_cpu.apply_diffusion(domain, 300.0);
  assert(chem_gpu.apply_diffusion(domain, chem_gpu_ref, 300.0));
  chem_gpu.sync_concentrations_to_host(chem_gpu_ref);

  const Real max_diff = max_abs_diff(
      chem_cpu.conc_data().front(), chem_gpu_ref.conc_data().front());
  assert(max_diff < 1.0e-10);

  std::cout << "  test_gpu_facade_diffusion_matches_cpu: PASSED"
            << " (max_diff=" << max_diff << ")\n";
}

void test_gpu_facade_z_gradient_diffusion() {
  Domain domain = make_domain(4, 3, 8);
  ChemicalSpec spec = diffusing_species(5.0e-10, 5.0e-3, 5.0e-3);
  spec.z_gradient_enabled = true;
  spec.z_gradient_lambda = 25.0e-6;

  ChemicalField chem_cpu;
  chem_cpu.init(domain, {spec});

  ChemicalField chem_gpu_ref = chem_cpu;
  ChemicalFieldGpu chem_gpu;
  chem_gpu.init(chem_gpu_ref);

  chem_cpu.apply_diffusion(domain, 60.0);
  assert(chem_gpu.apply_diffusion(domain, chem_gpu_ref, 60.0));
  chem_gpu.sync_concentrations_to_host(chem_gpu_ref);

  const Real max_diff = max_abs_diff(
      chem_cpu.conc_data().front(), chem_gpu_ref.conc_data().front());
  assert(max_diff < 1.0e-10);

  std::cout << "  test_gpu_facade_z_gradient_diffusion: PASSED"
            << " (max_diff=" << max_diff << ")\n";
}

void test_gpu_facade_boundaries_diffusing_vs_static() {
  Domain domain = make_domain(4, 4, 6);
  ChemicalField chem_cpu;
  chem_cpu.init(domain, {
      diffusing_species(2.1e-9, 0.1, 1.0),
      static_species(0.2, 0.5),
  });
  chem_cpu.conc(0, domain.cell_index(1, 1, 4)) = 0.4;
  chem_cpu.conc(1, domain.cell_index(2, 2, 4)) = 0.7;

  const Int top_diffusing = domain.cell_index(0, 0, domain.nz() - 1);
  const Int below_diffusing = domain.cell_index(0, 0, domain.nz() - 2);
  const Int top_static = domain.cell_index(1, 1, domain.nz() - 1);
  const Int below_static = domain.cell_index(1, 1, domain.nz() - 2);
  chem_cpu.conc(0, top_diffusing) = 0.01;
  chem_cpu.conc(0, below_diffusing) = 0.99;
  chem_cpu.conc(1, top_static) = 0.01;
  chem_cpu.conc(1, below_static) = 0.99;

  ChemicalField chem_gpu_ref = chem_cpu;
  ChemicalFieldGpu chem_gpu;
  chem_gpu.init(chem_gpu_ref);

  chem_cpu.apply_boundaries(domain);
  assert(chem_gpu.apply_boundaries(domain, chem_gpu_ref));
  chem_gpu.sync_concentrations_to_host(chem_gpu_ref);

  for (Int s = 0; s < chem_cpu.num_species(); ++s) {
    const Real max_diff = max_abs_diff(
        chem_cpu.conc_data()[static_cast<size_t>(s)],
        chem_gpu_ref.conc_data()[static_cast<size_t>(s)]);
    assert(max_diff < 1.0e-10);
  }

  assert(std::abs(chem_cpu.conc(0, top_diffusing) - 0.01) < 1.0e-10);
  assert(std::abs(chem_cpu.conc(0, below_diffusing) - 0.99) < 1.0e-10);
  assert(std::abs(chem_cpu.conc(1, top_static) - chem_cpu.conc(1, below_static))
         < 1.0e-10);
  assert(std::abs(chem_cpu.conc(1, below_static) - 0.99) < 1.0e-10);

  std::cout << "  test_gpu_facade_boundaries_diffusing_vs_static: PASSED\n";
}

int main() {
  std::cout << "=== GPU ChemicalField Facade Tests ===\n";

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  std::cout << "All GPU chemical-field facade tests passed.\n";
  return 0;
#else
  GpuConfig gcfg;
  gcfg.enabled = true;
  gcfg.device_id = 0;
  gpu_set_config(gcfg);

  if (!gpu_init_for_rank(0, 1)) {
    std::cout << "  SKIPPED (no CUDA device)\n";
    std::cout << "All GPU chemical-field facade tests passed.\n";
    return 0;
  }

  test_gpu_facade_diffusion_matches_cpu();
  test_gpu_facade_z_gradient_diffusion();
  test_gpu_facade_boundaries_diffusing_vs_static();

  std::cout << "All GPU chemical-field facade tests passed.\n";
  return 0;
#endif
}
