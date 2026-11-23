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
#include <limits>
#include <string>
#include <string_view>
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

Domain make_unit_domain() {
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {20.0e-6, 20.0e-6, 40.0e-6};
  dcfg.grid_dx = 5.0e-6;
  dcfg.grid_halo_width = 2;
  Domain domain;
  domain.init(dcfg);
  return domain;
}

std::vector<ChemicalSpec> make_carbon_mucin_specs(Real carbon_init,
                                                  Real mucin_init) {
  ChemicalSpec carbon_spec;
  carbon_spec.name = species::CARBON;
  carbon_spec.initial_conc = carbon_init;
  carbon_spec.diffusion_enabled = false;

  ChemicalSpec mucin_spec;
  mucin_spec.name = species::MUCIN;
  mucin_spec.initial_conc = mucin_init;
  mucin_spec.boundary_conc = mucin_init;
  mucin_spec.diffusion_enabled = false;
  return {carbon_spec, mucin_spec};
}

VBFConfig make_dynamic_vbf_config(Real lambda = 25.0e-6,
                                  bool use_dynamic = true) {
  VBFConfig vcfg;
  vcfg.density = 1.0e11;
  vcfg.use_dynamic_mucin = use_dynamic;
  vcfg.mucin_z_gradient_enabled = true;
  vcfg.mucin_z_gradient_lambda = lambda;
  vcfg.carbon_sink_vmax = 0.0;
  vcfg.nutrient_sink = 0.0;
  return vcfg;
}

MucinConfig make_dynamic_mucin_config(bool enabled = true,
                                      Real k_liberation = 5.0e-16) {
  MucinConfig mcfg;
  mcfg.enabled = enabled;
  mcfg.initial_conc = 1.0e-2;
  mcfg.secretion_rate = 0.0;
  mcfg.Km_degradation = 1.0e-3;
  mcfg.k_liberation = k_liberation;
  return mcfg;
}

ChemicalField run_vbf_step(const Domain& domain,
                           const VBFConfig& vcfg,
                           const MucinConfig& mcfg,
                           std::string_view decomposition,
                           Real carbon_init,
                           Real mucin_init,
                           Real dt,
                           VbfFluxTotals* totals = nullptr) {
  ChemicalField chem;
  chem.init(domain, make_carbon_mucin_specs(carbon_init, mucin_init),
            decomposition);
  chem.zero_reactions();
  VBF vbf;
  vbf.init(vcfg, domain);
  OxygenConfig oxygen;
  AcetateConfig acetate;
  vbf.apply_nutrient_coupling(chem, domain, dt, oxygen, acetate, mcfg, totals);
  return chem;
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
    std::string_view decomposition,
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
  cfg.chemistry_decomposition = std::string(decomposition);
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
  // Explicit forward Euler over 60 s would overshoot initial mucin by ~60x.
  cfg.chem_env.mucin.k_liberation = 1.0e-14;

  cfg.vbf = make_dynamic_vbf_config(25.0e-6, true);
  cfg.vbf.mucin_z_gradient_enabled = mucin_z_gradient_enabled;
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

Real run_closed_domain_case(std::string_view decomposition,
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
  // Uniform immobile mucin, depth-weighted immobile mucin (default polymer),
  // depth-weighted diffusing mucin, and basal epithelial secretion.
  worst_err = std::max(
      worst_err, run_closed_domain_case("replicated", false, false, 0.0, 20));
  worst_err = std::max(
      worst_err, run_closed_domain_case("slab", false, false, 0.0, 20));
  worst_err = std::max(
      worst_err, run_closed_domain_case("replicated", true, false, 0.0, 20));
  worst_err = std::max(
      worst_err, run_closed_domain_case("slab", true, false, 0.0, 20));
  worst_err = std::max(
      worst_err, run_closed_domain_case("replicated", true, true, 0.0, 20));
  worst_err = std::max(
      worst_err, run_closed_domain_case("slab", true, true, 0.0, 20));
  worst_err = std::max(
      worst_err, run_closed_domain_case("replicated", true, false, 2.0e-7, 10));

  if (rank == 0) {
    std::cout << "  test_closed_domain_mass_conservation: PASSED"
              << " (max_rel_err=" << worst_err << ")\n";
  }
}

void test_near_saturating_steady_limit_and_bounds(int rank) {
  const Domain domain = make_unit_domain();
  const VBFConfig vcfg = make_dynamic_vbf_config();
  const MucinConfig mcfg = make_dynamic_mucin_config();
  constexpr Real small_dt = 0.01;
  const ChemicalField chem = run_vbf_step(
      domain, vcfg, mcfg, "replicated", 0.0, mcfg.initial_conc, small_dt);

  const Int i_carbon = chem.find(species::CARBON);
  const Int i_mucin = chem.find(species::MUCIN);
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

  // Boundary inputs for the shared host/device helper.
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

  // Subnormal regression: rounding overshoot on tiny mucin concentrations
  // must terminate via nextafter and preserve mucin_conc - rate * dt >= 0.
  const Real denorm_min = std::numeric_limits<Real>::denorm_min();
  const Real subnormal_example = 7.7e-322;
  const Real ex_rate = vbf::implicit_mucin_liberation(
      subnormal_example, 5.0e-16, 1.0e11, std::exp(-3.0e-6 / 25.0e-6),
      1.0e-3, 60.0);
  assert(std::isfinite(ex_rate));
  assert(ex_rate >= 0.0);
  assert(subnormal_example - ex_rate * 60.0 >= 0.0);

  for (int k = 1; k <= 200; ++k) {
    const Real m_sub = static_cast<Real>(k) * denorm_min;
    const Real r_sub = vbf::implicit_mucin_liberation(
        m_sub, 5.0e-16, 1.0e11, 1.0, 1.0e-3, 60.0);
    assert(std::isfinite(r_sub));
    assert(r_sub >= 0.0);
    assert(m_sub - r_sub * 60.0 >= 0.0);
  }

  if (rank == 0) {
    std::cout << "  test_near_saturating_steady_limit_and_bounds: PASSED\n";
  }
}

void verify_z_gradient_profile_for_mode(std::string_view decomposition) {
  const Domain domain = make_unit_domain();
  const VBFConfig vcfg = make_dynamic_vbf_config();
  const MucinConfig mcfg = make_dynamic_mucin_config();
  const Int owned_ix = domain.local_grid_x_begin();

  // Small dt recovers the exact exponential z-gradient ratio.
  constexpr Real small_dt = 0.01;
  const ChemicalField chem_small = run_vbf_step(
      domain, vcfg, mcfg, decomposition, 0.0, mcfg.initial_conc, small_dt);
  const Int i_carbon = chem_small.find(species::CARBON);
  const Int i_mucin = chem_small.find(species::MUCIN);
  const Real r0_small =
      chem_small.reac_global(i_carbon, domain.cell_index(owned_ix, 0, 0));
  assert(r0_small > 0.0);
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    const Int gcell = domain.cell_index(owned_ix, 0, iz);
    const Real r_c = chem_small.reac_global(i_carbon, gcell);
    const Real r_m = chem_small.reac_global(i_mucin, gcell);
    assert(std::abs(r_c + r_m) <= 1.0e-15 * r_c);

    const Real ratio = r_c / r0_small;
    const Real expected_ratio = std::exp(
        -((static_cast<Real>(iz) + 0.5) - 0.5) * domain.dx_z()
        / vcfg.mucin_z_gradient_lambda);
    assert(std::abs(ratio - expected_ratio) / expected_ratio <= 1.0e-5);
  }

  // Full bio_dt matches the exact implicit formula at every z layer.
  constexpr Real bio_dt = 60.0;
  const ChemicalField chem_bio = run_vbf_step(
      domain, vcfg, mcfg, decomposition, 0.0, mcfg.initial_conc, bio_dt);
  Real prev_rate = 0.0;
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    const Int gcell = domain.cell_index(owned_ix, 0, iz);
    const Real r_c = chem_bio.reac_global(i_carbon, gcell);
    const Real r_m = chem_bio.reac_global(i_mucin, gcell);
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

  const Domain domain = make_unit_domain();
  const MucinConfig mcfg = make_dynamic_mucin_config();
  const std::array<Real, 3> lambdas = {10.0e-6, 25.0e-6, 100.0e-6};
  std::array<Real, 3> top_to_bottom_ratios{};
  for (size_t i = 0; i < lambdas.size(); ++i) {
    const VBFConfig vcfg = make_dynamic_vbf_config(lambdas[i], true);
    const ChemicalField chem_rep = run_vbf_step(
        domain, vcfg, mcfg, "replicated", 0.0, mcfg.initial_conc, 60.0);
    const ChemicalField chem_slab = run_vbf_step(
        domain, vcfg, mcfg, "slab", 0.0, mcfg.initial_conc, 60.0);

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

  const Domain domain = make_unit_domain();
  VBFConfig vcfg_static = make_dynamic_vbf_config(25.0e-6, false);
  vcfg_static.mucin_liberation = 5.0e-5;
  vcfg_static.carbon_sink_vmax = 5.5e-5;
  vcfg_static.carbon_sink_km = 1.0e-4;

  const MucinConfig mcfg_off = make_dynamic_mucin_config(false);
  const MucinConfig mcfg_on = make_dynamic_mucin_config(true);

  VbfFluxTotals totals_baseline;
  const ChemicalField chem_baseline = run_vbf_step(
      domain, vcfg_static, mcfg_off, "replicated", 1.0e-3, 1.0e-2, 60.0,
      &totals_baseline);

  // Dynamic toggle disabled while mucin field is enabled.
  VbfFluxTotals totals_a;
  const ChemicalField chem_a = run_vbf_step(
      domain, vcfg_static, mcfg_on, "replicated", 1.0e-3, 1.0e-2, 60.0,
      &totals_a);

  // Dynamic toggle enabled while mucin field is disabled.
  VBFConfig vcfg_dyn_flag = vcfg_static;
  vcfg_dyn_flag.use_dynamic_mucin = true;
  VbfFluxTotals totals_b;
  const ChemicalField chem_b = run_vbf_step(
      domain, vcfg_dyn_flag, mcfg_off, "replicated", 1.0e-3, 1.0e-2, 60.0,
      &totals_b);

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
