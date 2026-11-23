/* -----------------------------------------------------------------------
   GutIBM – Dynamic mucin stoichiometric limitation and z-gradient tests
   ----------------------------------------------------------------------- */

#include "chem_environment_config.h"
#include "chemical_field.h"
#include "domain.h"
#include "input_parser.h"
#include "simulation.h"
#include "species_names.h"
#include "vbf.h"
#include "vbf_mucin_liberation.h"

#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

using namespace gutibm;

namespace {

bool is_exact_zero(Real value) {
  return value >= 0.0 && !(value > 0.0);
}

bool bitwise_equal(Real a, Real b) {
  return std::memcmp(&a, &b, sizeof(Real)) == 0;
}

Real global_owned_sum(const ChemicalField& chem,
                      const Domain& domain,
                      Int spec) {
  Real local_sum = 0.0;
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    if (!chem.owns_global_cell(cell)) continue;
    const Real val = chem.conc_global(spec, cell);
    assert(std::isfinite(val));
    assert(val >= 0.0);
    local_sum += val;
  }
  Real total = local_sum * domain.cell_volume();
#ifdef GUTIBM_MPI
  if (chem.slab_mode() && domain.nprocs() > 1) {
    MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  return total;
}

Real global_max_value(Real local_value, const Domain& domain) {
  Real result = local_value;
#ifdef GUTIBM_MPI
  if (domain.nprocs() > 1) {
    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
  }
#else
  (void)domain;
#endif
  return result;
}

SimulationConfig make_closed_dynamic_mucin_config(
    const std::string& decomposition,
    bool mucin_z_gradient_enabled,
    bool mucin_diffusion_enabled,
    Real secretion_rate) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {40.0e-6, 40.0e-6, 80.0e-6};
  cfg.domain.grid_dx = 10.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.domain.ghost_width = 10.0e-6;
  cfg.domain.grid_halo_width = 1;
  cfg.chemistry_decomposition = decomposition;
  cfg.fixes.metabolism.delivery_far_field_radius = 0.0;

  cfg.time.bio_dt = 60.0;
  cfg.time.total_time = 20.0 * cfg.time.bio_dt;
  cfg.time.output_interval = cfg.time.total_time;
  cfg.hdf5.enabled = false;
  cfg.gpu.enabled = false;

  cfg.initial_strains.clear();
  cfg.enabled_fixes = {"metabolism"};
  cfg.chem_env.oxygen.enabled = false;
  cfg.chem_env.acetate.enabled = false;
  cfg.chem_env.siderophore.enabled = false;
  cfg.chem_env.ferrichrome.enabled = false;

  cfg.chem_env.mucin.enabled = true;
  cfg.chem_env.mucin.initial_conc = 1.0e-6;
  cfg.chem_env.mucin.secretion_rate = secretion_rate;
  cfg.chem_env.mucin.Km_degradation = 1.0e-3;
  // k_liberation * density = 1e-14 * 1e11 = 1e-3 mol/(m^3 s), so explicit
  // forward Euler at dt = 60 s would overshoot initial mucin (1e-6) by ~60x.
  cfg.chem_env.mucin.k_liberation = 1.0e-14;

  cfg.vbf.density = 1.0e11;
  cfg.vbf.use_dynamic_mucin = true;
  cfg.vbf.mucin_z_gradient_enabled = mucin_z_gradient_enabled;
  cfg.vbf.mucin_z_gradient_lambda = 25.0e-6;
  cfg.vbf.carbon_sink_vmax = 0.0;
  cfg.vbf.nutrient_sink = 0.0;
  cfg.qssa.drift_envelope_policy = "allow";

  cfg.carbon_epithelial_boundary = "flux";
  cfg.carbon_epithelial_flux = 0.0;
  for (auto& spec : cfg.chemicals) {
    if (spec.name == species::CARBON) {
      spec.initial_conc = 0.0;
      spec.boundary_conc = 0.0;
      spec.z_gradient_enabled = false;
      spec.diffusion_enabled = true;
    }
  }

  InputParser::finalize_config(cfg);

  for (auto& spec : cfg.chemicals) {
    if (spec.name == species::MUCIN) {
      spec.initial_conc = cfg.chem_env.mucin.initial_conc;
      spec.boundary_conc = cfg.chem_env.mucin.initial_conc;
      spec.z_gradient_enabled = false;
      spec.epithelial_boundary_mode = EpithelialBoundaryMode::Flux;
      spec.epithelial_flux = 0.0;
      spec.diffusion_enabled = mucin_diffusion_enabled;
    }
  }
  return cfg;
}

Real run_closed_domain_case(const std::string& decomposition,
                            bool mucin_z_gradient_enabled,
                            bool mucin_diffusion_enabled,
                            Real secretion_rate,
                            int steps) {
  SimulationConfig cfg = make_closed_dynamic_mucin_config(
      decomposition, mucin_z_gradient_enabled, mucin_diffusion_enabled,
      secretion_rate);
  Simulation sim;
  sim.init(cfg);

  const ChemicalField& chem = sim.chemical_field();
  const Domain& domain = sim.domain();
  const Int i_carbon = chem.find(species::CARBON);
  const Int i_mucin = chem.find(species::MUCIN);
  assert(i_carbon >= 0);
  assert(i_mucin >= 0);

  const Real initial_mucin_moles = global_owned_sum(chem, domain, i_mucin);
  const Real initial_carbon_moles = global_owned_sum(chem, domain, i_carbon);
  const Real initial_total_moles = initial_mucin_moles + initial_carbon_moles;
  assert(initial_mucin_moles > 0.0);

  const Real basal_layer_volume =
      static_cast<Real>(domain.nx() * domain.ny()) * domain.cell_volume();
  Real expected_secreted_moles = 0.0;
  Real max_rel_err = 0.0;

  for (int step = 0; step < steps; ++step) {
    sim.step(cfg.time.bio_dt);
    expected_secreted_moles +=
        secretion_rate * basal_layer_volume * cfg.time.bio_dt;

    const auto& flux = chem.flux_accounting();
    const Real mucin_clip = global_max_value(
        flux.reaction_clip_for_step(i_mucin), domain);
    const Real carbon_clip = global_max_value(
        flux.reaction_clip_for_step(i_carbon), domain);
    assert(is_exact_zero(mucin_clip));
    assert(is_exact_zero(carbon_clip));

    const Real mucin_moles = global_owned_sum(chem, domain, i_mucin);
    const Real carbon_moles = global_owned_sum(chem, domain, i_carbon);
    const Real total_moles = mucin_moles + carbon_moles;
    const Real expected_total = initial_total_moles + expected_secreted_moles;
    const Real rel_err =
        std::abs(total_moles - expected_total) / expected_total;
    max_rel_err = std::max(max_rel_err, rel_err);
    assert(rel_err <= 1.0e-12);
  }

  const Real final_mucin_moles = global_owned_sum(chem, domain, i_mucin);
  const Real final_carbon_moles = global_owned_sum(chem, domain, i_carbon);
  assert(final_carbon_moles > initial_carbon_moles);
  if (!(secretion_rate > 0.0)) {
    assert(final_mucin_moles < initial_mucin_moles);
    assert(final_mucin_moles > 0.0);
  }
  return max_rel_err;
}

void test_closed_domain_mass_conservation(int rank) {
  Real worst_err = 0.0;
  worst_err = std::max(
      worst_err, run_closed_domain_case("replicated", false, false, 0.0, 20));
  worst_err = std::max(
      worst_err, run_closed_domain_case("slab", false, false, 0.0, 20));
  worst_err = std::max(
      worst_err, run_closed_domain_case("replicated", true, true, 0.0, 20));
  worst_err = std::max(
      worst_err, run_closed_domain_case("slab", true, true, 0.0, 20));
  // Basal epithelial secretion at iz == 0 adds positive mucin flux cleanly.
  worst_err = std::max(
      worst_err, run_closed_domain_case("replicated", true, true, 2.0e-7, 10));

  if (rank == 0) {
    std::cout << "  test_closed_domain_mass_conservation: PASSED"
              << " (max_rel_err=" << worst_err << ")\n";
  }
}

void test_near_saturating_steady_limit_and_bounds(int rank) {
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {20.0e-6, 20.0e-6, 40.0e-6};
  dcfg.grid_dx = 5.0e-6;
  Domain domain;
  domain.init(dcfg);

  ChemicalSpec carbon_spec;
  carbon_spec.name = species::CARBON;
  carbon_spec.initial_conc = 0.0;
  carbon_spec.diffusion_enabled = false;

  ChemicalSpec mucin_spec;
  mucin_spec.name = species::MUCIN;
  mucin_spec.initial_conc = 1.0e-2;
  mucin_spec.boundary_conc = 1.0e-2;
  mucin_spec.diffusion_enabled = false;

  ChemicalField chem;
  chem.init(domain, {carbon_spec, mucin_spec});
  const Int i_carbon = chem.find(species::CARBON);
  const Int i_mucin = chem.find(species::MUCIN);

  VBFConfig vcfg;
  vcfg.density = 1.0e11;
  vcfg.use_dynamic_mucin = true;
  vcfg.mucin_z_gradient_enabled = true;
  vcfg.mucin_z_gradient_lambda = 25.0e-6;
  vcfg.carbon_sink_vmax = 0.0;
  vcfg.nutrient_sink = 0.0;

  MucinConfig mcfg;
  mcfg.enabled = true;
  mcfg.initial_conc = 1.0e-2;
  mcfg.secretion_rate = 0.0;
  mcfg.Km_degradation = 1.0e-3;
  mcfg.k_liberation = 5.0e-16;

  OxygenConfig oxygen;
  AcetateConfig acetate;
  VBF vbf;
  vbf.init(vcfg, domain);

  constexpr Real small_dt = 0.01;
  chem.zero_reactions();
  vbf.apply_nutrient_coupling(chem, domain, small_dt, oxygen, acetate, mcfg);

  for (Int iz = 0; iz < domain.nz(); ++iz) {
    const Real z_rel = (iz + 0.5) * domain.dx_z();
    const Real z_weight = std::exp(-z_rel / vcfg.mucin_z_gradient_lambda);
    const Real vmax = mcfg.k_liberation * vcfg.density * z_weight;
    const Real expected_explicit =
        vmax * mcfg.initial_conc / (mcfg.Km_degradation + mcfg.initial_conc);
    const Int cell = domain.cell_index(0, 0, iz);
    const Real r_carbon = chem.reac(i_carbon, cell);
    const Real r_mucin = chem.reac(i_mucin, cell);

    assert(r_carbon > 0.0);
    assert(std::abs(r_carbon + r_mucin) <= 1.0e-15 * r_carbon);
    const Real rel_err =
        std::abs(r_carbon - expected_explicit) / expected_explicit;
    assert(rel_err <= 1.0e-5);
  }

  // Direct edge-case bounds for implicit_mucin_liberation.
  assert(is_exact_zero(
      vbf::implicit_mucin_liberation(0.0, 5.0e-16, 1.0e11, 1.0, 1.0e-3, 60.0)));
  assert(is_exact_zero(vbf::implicit_mucin_liberation(
      -1.0e-6, 5.0e-16, 1.0e11, 1.0, 1.0e-3, 60.0)));
  assert(is_exact_zero(
      vbf::implicit_mucin_liberation(1.0e-2, 0.0, 1.0e11, 1.0, 1.0e-3, 60.0)));
  assert(is_exact_zero(
      vbf::implicit_mucin_liberation(1.0e-2, 5.0e-16, 1.0e11, 1.0, 1.0e-3, 0.0)));

  // Extreme stiffness with tiny and zero Km must never overshoot M / dt.
  for (const Real km : {0.0, 1.0e-18, 1.0e-9, 1.0e-3}) {
    for (const Real m0 : {1.0e-12, 1.0e-6, 1.0e-2}) {
      for (const Real dt : {0.01, 7.0, 60.0, 300.0}) {
        const Real rate = vbf::implicit_mucin_liberation(
            m0, 1.0e-6, 1.0e11, 1.0, km, dt);
        assert(std::isfinite(rate));
        assert(rate >= 0.0);
        assert(m0 - rate * dt >= 0.0);
      }
    }
  }

  if (rank == 0) {
    std::cout << "  test_near_saturating_steady_limit_and_bounds: PASSED\n";
  }
}

void verify_z_gradient_profile_for_mode(const std::string& decomposition) {
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {20.0e-6, 20.0e-6, 40.0e-6};
  dcfg.grid_dx = 5.0e-6;
  dcfg.grid_halo_width = 2;
  Domain domain;
  domain.init(dcfg);

  ChemicalSpec carbon_spec;
  carbon_spec.name = species::CARBON;
  carbon_spec.initial_conc = 0.0;
  carbon_spec.diffusion_enabled = false;

  ChemicalSpec mucin_spec;
  mucin_spec.name = species::MUCIN;
  mucin_spec.initial_conc = 1.0e-2;
  mucin_spec.boundary_conc = 1.0e-2;
  mucin_spec.diffusion_enabled = false;

  ChemicalField chem;
  chem.init(domain, {carbon_spec, mucin_spec}, decomposition);
  const Int i_carbon = chem.find(species::CARBON);
  const Int i_mucin = chem.find(species::MUCIN);

  VBFConfig vcfg;
  vcfg.density = 1.0e11;
  vcfg.use_dynamic_mucin = true;
  vcfg.mucin_z_gradient_enabled = true;
  vcfg.mucin_z_gradient_lambda = 25.0e-6;
  vcfg.carbon_sink_vmax = 0.0;
  vcfg.nutrient_sink = 0.0;

  MucinConfig mcfg;
  mcfg.enabled = true;
  mcfg.initial_conc = 1.0e-2;
  mcfg.secretion_rate = 0.0;
  mcfg.Km_degradation = 1.0e-3;
  mcfg.k_liberation = 5.0e-16;

  OxygenConfig oxygen;
  AcetateConfig acetate;
  VBF vbf;
  vbf.init(vcfg, domain);

  const Int owned_ix = domain.local_grid_x_begin();

  // Small dt recovers the exact exponential z-gradient ratio.
  constexpr Real small_dt = 0.01;
  chem.zero_reactions();
  vbf.apply_nutrient_coupling(chem, domain, small_dt, oxygen, acetate, mcfg);
  const Real r0_small = chem.reac_global(i_carbon, domain.cell_index(owned_ix, 0, 0));
  assert(r0_small > 0.0);
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    const Int gcell = domain.cell_index(owned_ix, 0, iz);
    const Real r_c = chem.reac_global(i_carbon, gcell);
    const Real r_m = chem.reac_global(i_mucin, gcell);
    assert(std::abs(r_c + r_m) <= 1.0e-15 * r_c);

    const Real ratio = r_c / r0_small;
    const Real expected_ratio = std::exp(
        -((static_cast<Real>(iz) + 0.5) - 0.5) * domain.dx_z()
        / vcfg.mucin_z_gradient_lambda);
    assert(std::abs(ratio - expected_ratio) / expected_ratio <= 1.0e-5);
  }

  // Full bio_dt = 60 s matches exact implicit formula at every z layer.
  constexpr Real bio_dt = 60.0;
  chem.zero_reactions();
  vbf.apply_nutrient_coupling(chem, domain, bio_dt, oxygen, acetate, mcfg);
  Real prev_rate = 0.0;
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    const Int gcell = domain.cell_index(owned_ix, 0, iz);
    const Real r_c = chem.reac_global(i_carbon, gcell);
    const Real r_m = chem.reac_global(i_mucin, gcell);
    const Real z_rel = (iz + 0.5) * domain.dx_z();
    const Real z_weight = std::exp(-z_rel / vcfg.mucin_z_gradient_lambda);
    const Real expected = vbf::implicit_mucin_liberation(
        mcfg.initial_conc, mcfg.k_liberation, vcfg.density,
        z_weight, mcfg.Km_degradation, bio_dt);
    assert(std::abs(r_c - expected) / expected <= 1.0e-14);
    assert(std::abs(r_c + r_m) <= 1.0e-15 * r_c);
    if (iz > 0) {
      assert(r_c < prev_rate);
    }
    prev_rate = r_c;
  }
}

void test_mucin_z_gradient_profile_parity(int rank) {
  verify_z_gradient_profile_for_mode("replicated");
  verify_z_gradient_profile_for_mode("slab");

  // Direct replicated vs slab cellwise parity and graded lambda sensitivity.
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {20.0e-6, 20.0e-6, 40.0e-6};
  dcfg.grid_dx = 5.0e-6;
  dcfg.grid_halo_width = 2;
  Domain domain;
  domain.init(dcfg);

  ChemicalSpec carbon_spec;
  carbon_spec.name = species::CARBON;
  carbon_spec.initial_conc = 0.0;
  carbon_spec.diffusion_enabled = false;

  ChemicalSpec mucin_spec;
  mucin_spec.name = species::MUCIN;
  mucin_spec.initial_conc = 1.0e-2;
  mucin_spec.boundary_conc = 1.0e-2;
  mucin_spec.diffusion_enabled = false;

  MucinConfig mcfg;
  mcfg.enabled = true;
  mcfg.initial_conc = 1.0e-2;
  mcfg.secretion_rate = 0.0;
  mcfg.Km_degradation = 1.0e-3;
  mcfg.k_liberation = 5.0e-16;
  OxygenConfig oxygen;
  AcetateConfig acetate;

  const std::array<Real, 3> lambdas = {10.0e-6, 25.0e-6, 100.0e-6};
  std::array<Real, 3> top_to_bottom_ratios{};
  for (size_t i = 0; i < lambdas.size(); ++i) {
    VBFConfig vcfg;
    vcfg.density = 1.0e11;
    vcfg.use_dynamic_mucin = true;
    vcfg.mucin_z_gradient_enabled = true;
    vcfg.mucin_z_gradient_lambda = lambdas[i];
    vcfg.carbon_sink_vmax = 0.0;
    vcfg.nutrient_sink = 0.0;

    VBF vbf;
    vbf.init(vcfg, domain);

    ChemicalField chem_rep;
    chem_rep.init(domain, {carbon_spec, mucin_spec}, "replicated");
    chem_rep.zero_reactions();
    vbf.apply_nutrient_coupling(chem_rep, domain, 60.0, oxygen, acetate, mcfg);

    ChemicalField chem_slab;
    chem_slab.init(domain, {carbon_spec, mucin_spec}, "slab");
    chem_slab.zero_reactions();
    vbf.apply_nutrient_coupling(chem_slab, domain, 60.0, oxygen, acetate, mcfg);

    const Int i_carbon = chem_rep.find(species::CARBON);
    const Int i_mucin = chem_rep.find(species::MUCIN);
    for (Int cell = 0; cell < chem_slab.global_ncells(); ++cell) {
      if (!chem_slab.owns_global_cell(cell)) continue;
      assert(bitwise_equal(chem_rep.reac_global(i_carbon, cell),
                           chem_slab.reac_global(i_carbon, cell)));
      assert(bitwise_equal(chem_rep.reac_global(i_mucin, cell),
                           chem_slab.reac_global(i_mucin, cell)));
    }

    const Int ix = domain.local_grid_x_begin();
    const Real r_bottom =
        chem_slab.reac_global(i_carbon, domain.cell_index(ix, 0, 0));
    const Real r_top = chem_slab.reac_global(
        i_carbon, domain.cell_index(ix, 0, domain.nz() - 1));
    top_to_bottom_ratios[i] = r_top / r_bottom;
    if (i > 0) {
      assert(top_to_bottom_ratios[i] > top_to_bottom_ratios[i - 1]);
    }
  }
  assert(top_to_bottom_ratios[2] - top_to_bottom_ratios[0] > 0.5);

  if (rank == 0) {
    std::cout << "  test_mucin_z_gradient_profile_parity: PASSED\n";
  }
}

void test_default_off_bitwise_guard(int rank) {
  const SimulationConfig defaults = InputParser::default_config();
  assert(!defaults.chem_env.mucin.enabled);
  assert(!defaults.vbf.use_dynamic_mucin);

  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {20.0e-6, 20.0e-6, 40.0e-6};
  dcfg.grid_dx = 5.0e-6;
  Domain domain;
  domain.init(dcfg);

  ChemicalSpec carbon_spec;
  carbon_spec.name = species::CARBON;
  carbon_spec.initial_conc = 1.0e-3;
  carbon_spec.diffusion_enabled = false;

  ChemicalSpec mucin_spec;
  mucin_spec.name = species::MUCIN;
  mucin_spec.initial_conc = 1.0e-2;
  mucin_spec.diffusion_enabled = false;

  OxygenConfig oxygen;
  AcetateConfig acetate;

  VBFConfig vcfg_static;
  vcfg_static.density = 1.0e11;
  vcfg_static.mucin_liberation = 5.0e-5;
  vcfg_static.mucin_z_gradient_enabled = true;
  vcfg_static.mucin_z_gradient_lambda = 25.0e-6;
  vcfg_static.use_dynamic_mucin = false;
  vcfg_static.carbon_sink_vmax = 5.5e-5;
  vcfg_static.carbon_sink_km = 1.0e-4;

  MucinConfig mcfg_off;
  mcfg_off.enabled = false;
  mcfg_off.secretion_rate = 0.0;

  ChemicalField chem_baseline;
  chem_baseline.init(domain, {carbon_spec, mucin_spec});
  chem_baseline.zero_reactions();
  VBF vbf_static;
  vbf_static.init(vcfg_static, domain);
  VbfFluxTotals totals_baseline;
  vbf_static.apply_nutrient_coupling(
      chem_baseline, domain, 60.0, oxygen, acetate, mcfg_off, &totals_baseline);

  // Case A: use_dynamic_mucin = false, mucin.enabled = true
  MucinConfig mcfg_on;
  mcfg_on.enabled = true;
  mcfg_on.secretion_rate = 0.0;
  ChemicalField chem_a;
  chem_a.init(domain, {carbon_spec, mucin_spec});
  chem_a.zero_reactions();
  VbfFluxTotals totals_a;
  vbf_static.apply_nutrient_coupling(
      chem_a, domain, 60.0, oxygen, acetate, mcfg_on, &totals_a);

  // Case B: use_dynamic_mucin = true, mucin.enabled = false
  VBFConfig vcfg_dyn_flag = vcfg_static;
  vcfg_dyn_flag.use_dynamic_mucin = true;
  VBF vbf_dyn_flag;
  vbf_dyn_flag.init(vcfg_dyn_flag, domain);
  ChemicalField chem_b;
  chem_b.init(domain, {carbon_spec, mucin_spec});
  chem_b.zero_reactions();
  VbfFluxTotals totals_b;
  vbf_dyn_flag.apply_nutrient_coupling(
      chem_b, domain, 60.0, oxygen, acetate, mcfg_off, &totals_b);

  const Int i_carbon = chem_baseline.find(species::CARBON);
  const Int i_mucin = chem_baseline.find(species::MUCIN);
  assert(bitwise_equal(totals_baseline.carbon_source, totals_a.carbon_source));
  assert(bitwise_equal(totals_baseline.carbon_source, totals_b.carbon_source));
  assert(bitwise_equal(totals_baseline.carbon_sink, totals_a.carbon_sink));
  assert(bitwise_equal(totals_baseline.carbon_sink, totals_b.carbon_sink));

  for (Int cell = 0; cell < chem_baseline.ncells(); ++cell) {
    assert(bitwise_equal(chem_baseline.reac(i_carbon, cell),
                         chem_a.reac(i_carbon, cell)));
    assert(bitwise_equal(chem_baseline.reac(i_carbon, cell),
                         chem_b.reac(i_carbon, cell)));
    assert(bitwise_equal(chem_baseline.reac(i_mucin, cell),
                         chem_a.reac(i_mucin, cell)));
    assert(bitwise_equal(chem_baseline.reac(i_mucin, cell),
                         chem_b.reac(i_mucin, cell)));
  }

  if (rank == 0) {
    std::cout << "  test_default_off_bitwise_guard: PASSED\n";
  }
}

}  // namespace

int main(int argc, char** argv) {
  int rank = 0;
#ifdef GUTIBM_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  (void)argc;
  (void)argv;
#endif

  if (rank == 0) {
    std::cout << "=== Dynamic Mucin Stoichiometry & z-Gradient Tests ===\n";
  }
  test_closed_domain_mass_conservation(rank);
  test_near_saturating_steady_limit_and_bounds(rank);
  test_mucin_z_gradient_profile_parity(rank);
  test_default_off_bitwise_guard(rank);
  if (rank == 0) {
    std::cout << "All dynamic mucin stoichiometry tests passed.\n";
  }

#ifdef GUTIBM_MPI
  MPI_Finalize();
#endif
  return 0;
}
