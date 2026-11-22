/* -----------------------------------------------------------------------
   GutIBM – ChemicalFieldGpu facade parity (Spec 9 PR2)
   Compares CPU ChemicalField paths with ChemicalFieldGpu device facades.
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "chemical_field_gpu.h"
#include "dispatch.h"
#include "domain.h"
#include "gpu_profile.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"
#include "gpu_test_support.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
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

Real max_rel_diff(const std::vector<Real>& host,
                  const std::vector<Real>& device) {
  Real max_err = 0.0;
  const size_t n = std::min(host.size(), device.size());
  for (size_t i = 0; i < n; ++i) {
    const Real scale = std::max({1.0e-15, std::abs(host[i]), std::abs(device[i])});
    max_err = std::max(max_err, std::abs(host[i] - device[i]) / scale);
  }
  return max_err;
}

void perturb_interior(ChemicalField& field, Int spec, const Domain& domain) {
  field.conc(spec, domain.cell_index(1, 1, 2)) += 0.05;
  field.conc(spec, domain.cell_index(2, 2, 3)) -= 0.03;
}

SimulationConfig make_ten_step_residency_config(
    const std::string& uptake_mode, bool gpu_enabled) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.seed = 480;
  cfg.time.bio_dt = 60.0;
  cfg.time.total_time = 600.0;
  cfg.time.output_interval = 600.0;
  cfg.hdf5.enabled = false;
  cfg.profile_steps = false;
  cfg.gpu.enabled = gpu_enabled;
  cfg.gpu.device_id = 0;

  // 24x24x12 grid keeps 10-step CPU+GPU parity fast while large enough that
  // full-grid concentration fields dominate small scalar/agent transfers.
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {120.0e-6, 120.0e-6, 60.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;

  cfg.fixes.metabolism.uptake_limit = uptake_mode;
  cfg.chem_env.siderophore.enabled = false;
  cfg.chem_env.oxygen.k_ROS = 0.0;
  cfg.qssa.use_fmm = false;

  // Eliminate stochastic discrete branching so CPU and GPU stay locked on the
  // same trajectory and concentration parity can be asserted to <= 1e-12.
  cfg.fixes.bacteriocin.sos_lysis_prob = 0.0;
  cfg.fixes.bacteriocin.sos_basal_rate = 0.0;
  cfg.disabled_fixes = {"conjugation", "mutation"};

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain producer;
  producer.type = 1;
  producer.count = 16;
  producer.mu_max = 5.5e-4;
  producer.plasmids = {"ColE1", "ColB"};
  producer.conjugative = false;
  cfg.initial_strains.push_back(producer);

  SimulationConfig::InitialStrain sensitive;
  sensitive.type = 2;
  sensitive.count = 16;
  sensitive.mu_max = 5.0e-4;
  sensitive.plasmids = {};
  sensitive.conjugative = false;
  cfg.initial_strains.push_back(sensitive);

  return cfg;
}

uint64_t site_bytes_for(const std::vector<GpuTransferSiteProfile>& sites,
                        std::string_view label,
                        bool h2d,
                        bool d2h) {
  uint64_t total = 0;
  for (const auto& site : sites) {
    if (site.label != label) continue;
    if (h2d) total += site.h2d_bytes;
    if (d2h) total += site.d2h_bytes;
  }
  return total;
}

void run_ten_step_residency_mode(const std::string& uptake_mode) {
  constexpr int kTotalSteps = 10;
  const SimulationConfig cfg_cpu =
      make_ten_step_residency_config(uptake_mode, false);
  const SimulationConfig cfg_gpu =
      make_ten_step_residency_config(uptake_mode, true);

  Simulation sim_cpu;
  sim_cpu.init(cfg_cpu);

  gpu_transfer_profile_set_enabled(false);
  gpu_transfer_profile_reset();
  Simulation sim_gpu;
  sim_gpu.init(cfg_gpu);
  assert(sim_gpu.gpu_active());

  gpu_transfer_profile_set_enabled(true);
  gpu_transfer_profile_reset();

  uint64_t total_h2d_10 = 0;
  uint64_t total_d2h_10 = 0;
  uint64_t mid_step_conc_upload = 0;
  uint64_t mid_step_diff_download = 0;
  uint64_t mid_step_qssa_sync = 0;
  uint64_t mid_step_chem_reactions = 0;

  for (int step_idx = 0; step_idx < kTotalSteps; ++step_idx) {
    gpu_transfer_profile_reset();
    sim_cpu.step(cfg_cpu.time.bio_dt);
    sim_gpu.step(cfg_gpu.time.bio_dt);

    const GpuProfile xfer = gpu_transfer_profile_snapshot();
    const auto sites = gpu_transfer_site_profiles();
    total_h2d_10 += xfer.h2d_bytes;
    total_d2h_10 += xfer.d2h_bytes;

    // Steps 2..9 (0-based indices 1..8): full-grid concentrations must remain
    // strictly device-resident with zero full-field H2D/D2H transfers.
    if (step_idx >= 1 && step_idx <= 8) {
      mid_step_conc_upload +=
          site_bytes_for(sites, "step_conc_upload", true, false);
      mid_step_diff_download +=
          site_bytes_for(sites, "chem_diffusion_result", false, true);
      mid_step_qssa_sync +=
          site_bytes_for(sites, "qssa_toxin_host_sync", false, true);
      mid_step_chem_reactions +=
          site_bytes_for(sites, "chem_reactions", true, true);
    }
  }

  assert(mid_step_conc_upload == 0ULL);
  assert(mid_step_diff_download == 0ULL);
  assert(mid_step_qssa_sync == 0ULL);
  if (uptake_mode != "delivery") {
    assert(mid_step_chem_reactions == 0ULL);
    assert(std::string(sim_gpu.chemistry_placement()) == "device");
  } else {
    assert(std::string(sim_gpu.chemistry_placement()) == "device_delivery");
  }

  // Compare total 10-step PCIe volume against the historical un-resident
  // baseline (~5 full-field grid transfers per step -> ~57 GB at 500x500x50).
  const auto& chem_gpu = sim_gpu.chemical_field();
  const auto& chem_cpu = sim_cpu.chemical_field();
  const uint64_t full_field_bytes =
      static_cast<uint64_t>(chem_gpu.num_species())
      * static_cast<uint64_t>(chem_gpu.ncells()) * sizeof(Real);
  const uint64_t historical_baseline_bytes =
      5ULL * full_field_bytes * static_cast<uint64_t>(kTotalSteps);
  const uint64_t total_pcie_bytes = total_h2d_10 + total_d2h_10;
  assert(total_pcie_bytes * 10ULL < historical_baseline_bytes);

  // Final step-10 read materializes device concentrations lazily on host and
  // verifies <= 1e-12 relative parity across every chemical species.
  Real worst_rel_err = 0.0;
  assert(chem_cpu.num_species() == chem_gpu.num_species());
  for (Int s = 0; s < chem_cpu.num_species(); ++s) {
    const Real rel_err = max_rel_diff(
        chem_cpu.species_concentration(s), chem_gpu.species_concentration(s));
    worst_rel_err = std::max(worst_rel_err, rel_err);
    assert(rel_err <= 1.0e-12);
  }
  assert(sim_cpu.agents().size() == sim_gpu.agents().size());

  gpu_transfer_profile_set_enabled(false);
  const double pct_of_baseline =
      100.0 * static_cast<double>(total_pcie_bytes)
      / static_cast<double>(historical_baseline_bytes);
  std::cout << "  test_ten_step_device_resident[" << uptake_mode << "]: PASSED"
            << " (max_rel_err=" << worst_rel_err
            << ", pcie_bytes=" << total_pcie_bytes
            << ", pct_of_baseline=" << pct_of_baseline << "%"
            << ", pinned_staging_active="
            << (gpu_transfer_pinned_staging_active() ? 1 : 0) << ")\n";
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

void test_host_reaction_survives_gpu_round_trip() {
  Domain domain = make_domain(4, 4, 4);
  ChemicalField host;
  host.init(domain, {static_species(0.2, 0.5)});
  const Int reaction_cell = domain.cell_index(1, 1, 1);
  host.reac(0, reaction_cell) = -0.125;
  host.mark_host_reac_dirty(0);

  ChemicalFieldGpu chem_gpu;
  chem_gpu.init(host);
  chem_gpu.sync_reactions_to_device(host);
  host.reac(0, reaction_cell) = 0.0;
  chem_gpu.sync_reactions_to_host(host);

  const Real expected = -0.125;
  assert(std::memcmp(&host.reac(0, reaction_cell), &expected,
                     sizeof(expected)) == 0);
  std::cout << "  test_host_reaction_survives_gpu_round_trip: PASSED\n";
}

void test_ten_step_device_resident_parity_and_pcie_budget() {
  run_ten_step_residency_mode("none");
  run_ten_step_residency_mode("sherwood");
  run_ten_step_residency_mode("delivery");
}

int main() {
  std::cout << "=== GPU ChemicalField Facade Tests ===\n";
  if (const int gpu_status = test::require_gpu("gpu_chemical_field");
      gpu_status != 0) {
    return gpu_status;
  }

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  std::cout << "All GPU chemical-field facade tests passed.\n";
  return 0;
#else
  GpuConfig gcfg;
  gcfg.enabled = true;
  gcfg.device_id = 0;
  gpu_set_config(gcfg);

  if (!gpu_init_for_rank(0, 1)) return 1;

  test_gpu_facade_diffusion_matches_cpu();
  test_gpu_facade_z_gradient_diffusion();
  test_gpu_facade_boundaries_diffusing_vs_static();
  test_host_reaction_survives_gpu_round_trip();
  test_ten_step_device_resident_parity_and_pcie_budget();

  std::cout << "All GPU chemical-field facade tests passed.\n";
  return 0;
#endif
}
