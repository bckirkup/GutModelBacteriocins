/* -----------------------------------------------------------------------
   GutIBM – GPU nutrient diffusion parity (Spec 9 PR1)
   Compares CPU ChemicalField::apply_diffusion with gpu_apply_species_diffusion.
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "diffusion_gpu.h"
#include "dispatch.h"
#include "domain.h"
#include "gpu_diagnostic_format.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"
#include "gpu_test_support.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace gutibm;
using gutibm::gpu_diagnostic::format_real;

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

ChemicalSpec delivery_species(EpithelialBoundaryMode mode) {
  ChemicalSpec spec = diffusing_species(2.1e-9, 0.0, 1.0);
  spec.epithelial_boundary_mode = mode;
  spec.z_gradient_enabled = false;
  spec.epithelial_transfer_coeff =
      mode == EpithelialBoundaryMode::Robin ? 2.0e-5 : 0.0;
  spec.epithelial_flux =
      mode == EpithelialBoundaryMode::Flux ? 1.0e-10 : 0.0;
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

void copy_conc(const ChemicalField& field, std::vector<Real>& out) {
  out = field.conc_data().front();
}

}  // namespace

void test_gpu_uniform_field_fixed_point() {
  Domain domain = make_domain(5, 4, 6);
  ChemicalField chem_cpu;
  chem_cpu.init(domain, {diffusing_species(2.1e-9, 0.25, 0.25)});

  std::vector<Real> conc_gpu;
  copy_conc(chem_cpu, conc_gpu);

  chem_cpu.apply_diffusion(domain, 300.0);
  assert(gpu_apply_species_diffusion(
      domain, chem_cpu.spec(0), conc_gpu, 300.0));

  const Real max_diff = max_abs_diff(chem_cpu.conc_data().front(), conc_gpu);
  if (!(max_diff < 1.0e-10)) {
    const auto& cpu_field = chem_cpu.conc_data().front();
    const auto [cpu_min_it, cpu_max_it] =
        std::ranges::minmax_element(cpu_field);
    Real cpu_sum = 0.0;
    for (const Real value : cpu_field) cpu_sum += value;
    const Real cpu_mean = cpu_sum / static_cast<Real>(cpu_field.size());
    const Real cpu_scale = std::max(
        {std::abs(*cpu_min_it), std::abs(cpu_mean), std::abs(*cpu_max_it),
         1.0e-30});
    std::cerr
        << "[gpu_diag][gpu_diffusion][uniform_field] measured_max_diff="
        << format_real(max_diff) << " reference=0 abs_diff="
        << format_real(max_diff) << " rel_diff="
        << format_real(max_diff / cpu_scale)
        << " tolerance=1.0e-10 cpu_min=" << format_real(*cpu_min_it)
        << " cpu_mean=" << format_real(cpu_mean)
        << " cpu_max=" << format_real(*cpu_max_it) << "\n";
  }
  assert(max_diff < 1.0e-10);

  std::cout << "  test_gpu_uniform_field_fixed_point: PASSED"
            << " (max_diff=" << max_diff << ")\n";
}

void test_gpu_point_source_invariants() {
  Domain domain = make_domain(5, 1, 2);
  ChemicalField chem_cpu;
  chem_cpu.init(domain, {diffusing_species(1.0e-12, 0.0, 0.0)});
  chem_cpu.conc(0, domain.cell_index(2, 0, 1)) = 1.0;

  std::vector<Real> conc_gpu;
  copy_conc(chem_cpu, conc_gpu);

  chem_cpu.apply_diffusion(domain, 2.5);
  assert(gpu_apply_species_diffusion(
      domain, chem_cpu.spec(0), conc_gpu, 2.5));

  const Real max_diff = max_abs_diff(chem_cpu.conc_data().front(), conc_gpu);
  assert(max_diff < 1.0e-10);
  std::cout << "  test_gpu_point_source_invariants: PASSED"
            << " (max_diff=" << max_diff << ")\n";
}

void test_gpu_singleton_periodic_axes() {
  Domain domain = make_domain(1, 1, 3);
  ChemicalField chem_cpu;
  chem_cpu.init(domain, {diffusing_species(1.0e-12, 0.0, 0.0)});
  chem_cpu.conc(0, domain.cell_index(0, 0, 1)) = 1.0;

  std::vector<Real> conc_gpu;
  copy_conc(chem_cpu, conc_gpu);

  chem_cpu.apply_diffusion(domain, 2.5);
  assert(gpu_apply_species_diffusion(
      domain, chem_cpu.spec(0), conc_gpu, 2.5));

  const Real max_diff = max_abs_diff(chem_cpu.conc_data().front(), conc_gpu);
  assert(max_diff < 1.0e-10);
  for (const Real value : conc_gpu) {
    assert(std::isfinite(value));
  }

  std::cout << "  test_gpu_singleton_periodic_axes: PASSED"
            << " (max_diff=" << max_diff << ")\n";
}

void test_gpu_dirichlet_neumann_boundary() {
  Domain domain = make_domain(4, 4, 8);
  ChemicalField chem_cpu;
  chem_cpu.init(domain, {diffusing_species(2.1e-9, 0.0, 1.0)});

  std::vector<Real> conc_gpu;
  copy_conc(chem_cpu, conc_gpu);

  chem_cpu.apply_diffusion(domain, 60.0);
  assert(gpu_apply_species_diffusion(
      domain, chem_cpu.spec(0), conc_gpu, 60.0));

  const Real max_diff = max_abs_diff(chem_cpu.conc_data().front(), conc_gpu);
  assert(max_diff < 1.0e-10);

  std::cout << "  test_gpu_dirichlet_neumann_boundary: PASSED"
            << " (max_diff=" << max_diff << ")\n";
}

void test_gpu_z_gradient_background_fixed_point() {
  Domain domain = make_domain(4, 3, 8);
  ChemicalSpec spec = diffusing_species(5.0e-10, 5.0e-3, 5.0e-3);
  spec.z_gradient_enabled = true;
  spec.z_gradient_lambda = 25.0e-6;

  ChemicalField chem_cpu;
  chem_cpu.init(domain, {spec});

  std::vector<Real> conc_gpu;
  copy_conc(chem_cpu, conc_gpu);

  chem_cpu.apply_diffusion(domain, 60.0);
  assert(gpu_apply_species_diffusion(domain, spec, conc_gpu, 60.0));

  const Real max_diff = max_abs_diff(chem_cpu.conc_data().front(), conc_gpu);
  assert(max_diff < 1.0e-10);

  std::cout << "  test_gpu_z_gradient_background_fixed_point: PASSED"
            << " (max_diff=" << max_diff << ")\n";
}

void test_gpu_delivery_boundary_modes() {
  Domain domain = make_domain(4, 3, 8);
  for (const auto mode : {EpithelialBoundaryMode::Robin,
                          EpithelialBoundaryMode::Flux}) {
    const ChemicalSpec spec = delivery_species(mode);
    ChemicalField chem_cpu;
    chem_cpu.init(domain, {spec});
    for (Int cell = 0; cell < domain.ncells(); ++cell) {
      chem_cpu.conc(0, cell) = 0.01 * static_cast<Real>(cell + 1);
    }
    std::vector<Real> conc_gpu;
    copy_conc(chem_cpu, conc_gpu);
    chem_cpu.apply_diffusion(domain, 60.0);
    assert(gpu_apply_species_diffusion(domain, spec, conc_gpu, 60.0));
    const Real max_diff = max_abs_diff(
        chem_cpu.conc_data().front(), conc_gpu);
    assert(max_diff < 1.0e-10);
  }
  std::cout << "  test_gpu_delivery_boundary_modes: PASSED\n";
}

#ifdef GUTIBM_CUDA
SimulationConfig mixed_boundary_config() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.time.total_time = 60.0;
  cfg.time.bio_dt = 60.0;
  cfg.hdf5.enabled = false;
  cfg.gpu.enabled = true;
  cfg.initial_strains.clear();
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {5.0e-6, 5.0e-6, 1025 * 5.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.chem_env.oxygen.enabled = false;
  cfg.chem_env.acetate.enabled = false;
  cfg.chem_env.mucin.enabled = false;
  cfg.chem_env.siderophore.enabled = false;
  cfg.chem_env.ferrichrome.enabled = false;
  for (ChemicalSpec& spec : cfg.chemicals) {
    spec.diffusion_enabled = false;
  }

  const auto carbon = std::find_if(
      cfg.chemicals.begin(), cfg.chemicals.end(),
      [](const ChemicalSpec& spec) { return spec.name == species::CARBON; });
  const auto oxygen = std::find_if(
      cfg.chemicals.begin(), cfg.chemicals.end(),
      [](const ChemicalSpec& spec) { return spec.name == species::OXYGEN; });
  assert(carbon != cfg.chemicals.end());
  assert(oxygen != cfg.chemicals.end());
  carbon->diffusion_enabled = true;
  carbon->diff_coeff = 2.1e-9;
  carbon->retardation = 1.0;
  carbon->initial_conc = 0.25;
  carbon->boundary_conc = 0.25;
  carbon->epithelial_boundary_mode = EpithelialBoundaryMode::Dirichlet;
  oxygen->diffusion_enabled = true;
  oxygen->diff_coeff = 2.1e-9;
  oxygen->retardation = 1.0;
  oxygen->initial_conc = 0.15;
  oxygen->boundary_conc = 0.15;
  oxygen->epithelial_boundary_mode = EpithelialBoundaryMode::Robin;
  oxygen->epithelial_transfer_coeff = 2.0e-5;
  return cfg;
}

void test_mixed_boundary_falls_back_to_cpu() {
  SimulationConfig cpu_config = mixed_boundary_config();
  cpu_config.gpu.enabled = false;
  Simulation cpu;
  cpu.init(cpu_config);
  cpu.step(cpu_config.time.bio_dt);

  SimulationConfig gpu_config = mixed_boundary_config();
  Simulation gpu;
  gpu.init(gpu_config);
  assert(gpu.gpu_active());
  assert(gpu.domain().nz() == 1025);
  const Int carbon = gpu.chemical_field().find(species::CARBON);
  const Int oxygen = gpu.chemical_field().find(species::OXYGEN);
  assert(carbon >= 0);
  assert(oxygen >= 0);
  assert(gpu.chemical_field().spec(carbon).epithelial_boundary_mode
         != gpu.chemical_field().spec(oxygen).epithelial_boundary_mode);
  gpu.step(gpu_config.time.bio_dt);

  assert(std::string(gpu.chemistry_placement()) == "host");
  const ChemicalField& cpu_field = cpu.chemical_field();
  const ChemicalField& gpu_field = gpu.chemical_field();
  assert(cpu_field.num_species() == gpu_field.num_species());
  for (Int s = 0; s < cpu_field.num_species(); ++s) {
    const Real max_diff = max_abs_diff(
        cpu_field.conc_data()[static_cast<size_t>(s)],
        gpu_field.conc_data()[static_cast<size_t>(s)]);
    assert(max_diff < 1.0e-10);
  }
  std::cout << "  test_mixed_boundary_falls_back_to_cpu: PASSED\n";
}
#endif

int main() {
  std::cout << "=== GPU Diffusion Parity Tests ===\n";
  if (const int gpu_status = test::require_gpu("gpu_diffusion");
      gpu_status != 0) {
    return gpu_status;
  }

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  std::cout << "All GPU diffusion tests passed.\n";
  return 0;
#else
  GpuConfig gcfg;
  gcfg.enabled = true;
  gcfg.device_id = 0;
  gpu_set_config(gcfg);

  if (!gpu_init_for_rank(0, 1)) return 1;

  test_gpu_uniform_field_fixed_point();
  test_gpu_point_source_invariants();
  test_gpu_singleton_periodic_axes();
  test_gpu_dirichlet_neumann_boundary();
  test_gpu_z_gradient_background_fixed_point();
  test_gpu_delivery_boundary_modes();
  test_mixed_boundary_falls_back_to_cpu();

  std::cout << "All GPU diffusion tests passed.\n";
  return 0;
#endif
}
