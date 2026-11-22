/* -----------------------------------------------------------------------
   GutIBM – Transient QSSA vs dynamic advection-diffusion-decay PDE bounds
   Quantifies receptor-occupancy and kill-probability error across
   bio_dt = 60 s (~300 s burst) in 60 um and 200 um geometries (Issue #472).
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "domain.h"
#include "greens_function.h"
#include "input_parser.h"
#include "plasmid.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

using namespace gutibm;

namespace {

// 16-point Gauss-Legendre nodes and weights on [-1, 1].
constexpr std::array<Real, 16> kGlNodes = {
    -0.98940093499164993260, -0.94457502307323257608,
    -0.86563120238783174388, -0.75540440835500303390,
    -0.61787624440264374845, -0.45801677765722738634,
    -0.28160355077925891323, -0.09501250983763744019,
     0.09501250983763744019,  0.28160355077925891323,
     0.45801677765722738634,  0.61787624440264374845,
     0.75540440835500303390,  0.86563120238783174388,
     0.94457502307323257608,  0.98940093499164993260};

constexpr std::array<Real, 16> kGlWeights = {
    0.02715245941175409485, 0.06225352393864789286,
    0.09515851168249278481, 0.12462897125553387205,
    0.14959598881657673208, 0.16915651939500253819,
    0.18260341504492358887, 0.18945061045506849628,
    0.18945061045506849628, 0.18260341504492358887,
    0.16915651939500253819, 0.14959598881657673208,
    0.12462897125553387205, 0.09515851168249278481,
    0.06225352393864789286, 0.02715245941175409485};

// Exact 3D heat kernel with sealed Neumann reflections at z=0 and z=H and
// uniform wall-parallel advection U = (ux, 0, 0).
Real bounded_heat_kernel(const Vec3& src, const Vec3& dst, Real s,
                         Real d_eff, Real ux, Real decay_rate, Real h) {
  if (s <= 0.0) return 0.0;
  const Real four_d_s = 4.0 * d_eff * s;
  const Real prefactor =
      std::exp(-decay_rate * s) / std::pow(std::numbers::pi * four_d_s, 1.5);
  const Real dx = (dst[0] - src[0]) - ux * s;
  const Real dy = dst[1] - src[1];
  const Real r2_xy = dx * dx + dy * dy;
  const Real xy_factor = std::exp(-r2_xy / four_d_s);
  if (xy_factor <= 0.0) return 0.0;

  // Determine image shells needed for < 1e-14 tail truncation.
  const Real max_dz = std::sqrt(36.0 * four_d_s);
  const int n_shells = std::clamp(
      static_cast<int>(std::ceil((max_dz + h) / (2.0 * h))) + 1, 2, 64);

  Real z_sum = 0.0;
  for (int n = -n_shells; n <= n_shells; ++n) {
    const Real z_plus = src[2] + 2.0 * static_cast<Real>(n) * h;
    const Real z_minus = -src[2] + 2.0 * static_cast<Real>(n) * h;
    const Real dz_p = dst[2] - z_plus;
    const Real dz_m = dst[2] - z_minus;
    z_sum += std::exp(-(dz_p * dz_p) / four_d_s)
        + std::exp(-(dz_m * dz_m) / four_d_s);
  }
  return prefactor * xy_factor * z_sum;
}

// Exact transient PDE concentration at time t for exponentially decaying burst
// Q(t') = Q0 * exp(-t' / release_tau) (or constant Q0 when release_tau <= 0).
Real transient_pde_concentration(const Vec3& src, const Vec3& dst, Real t,
                                 Real q0, Real release_tau, Real d_eff,
                                 Real ux, Real decay_rate, Real h) {
  if (t <= 0.0) return 0.0;
  // Change of variable s = u^2, ds = 2 u du removes small-s stiffness.
  // Geometric panels resolve both the early diffusive arrival peak and long
  // steady-state tails simultaneously.
  const Real u_max = std::sqrt(t);
  const Real dx = dst[0] - src[0];
  const Real dy = dst[1] - src[1];
  const Real dz = dst[2] - src[2];
  const Real r = std::sqrt(dx * dx + dy * dy + dz * dz);
  const Real u_char = std::max(r / std::sqrt(4.0 * std::max(d_eff, 1.0e-30)), 1.0e-3);
  const Real u_first = std::min(0.1 * u_char, 0.05 * u_max);
  constexpr int kGeomPanels = 24;
  std::array<Real, kGeomPanels + 2> edges{};
  edges[0] = 0.0;
  edges[1] = u_first;
  const Real ratio = std::pow(u_max / u_first, 1.0 / static_cast<Real>(kGeomPanels));
  for (int p = 1; p <= kGeomPanels; ++p) {
    edges[static_cast<size_t>(p + 1)] =
        (p == kGeomPanels) ? u_max : edges[static_cast<size_t>(p)] * ratio;
  }

  Real integral = 0.0;
  for (size_t p = 0; p + 1 < edges.size(); ++p) {
    const Real a = edges[p];
    const Real b = edges[p + 1];
    const Real half = 0.5 * (b - a);
    const Real mid = 0.5 * (a + b);
    for (size_t k = 0; k < kGlNodes.size(); ++k) {
      const Real u = mid + half * kGlNodes[k];
      const Real s = u * u;
      if (s <= 0.0) continue;
      const Real release_weight =
          release_tau > 0.0 ? std::exp(-(t - s) / release_tau) : 1.0;
      const Real g = bounded_heat_kernel(src, dst, s, d_eff, ux, decay_rate, h);
      integral += half * kGlWeights[k] * (2.0 * u) * q0 * release_weight * g;
    }
  }
  return integral;
}

Real receptor_occupancy(Real conc, Real kd_apparent) {
  if (conc <= 0.0) return 0.0;
  return conc / (conc + kd_apparent);
}

struct BurstErrorSummary {
  Real max_occ_err_60s = 0.0;
  Real max_occ_err_300s = 0.0;
  Real pkill_qssa_60s = 0.0;
  Real pkill_pde_60s = 0.0;
  Real pkill_qssa_300s = 0.0;
  Real pkill_pde_300s = 0.0;
  Real peak_c_qssa = 0.0;
  Real peak_c_pde = 0.0;
};

BurstErrorSummary evaluate_burst_error(
    const GreensFunction& gf, const GreensFunctionParams& params,
    const Vec3& src, const Vec3& dst, Real ux, Real h,
    Real release_tau, Real kd_apparent, Real k_kill) {
  constexpr Real bio_dt = 60.0;
  constexpr int n_steps = 5;  // 300 s total burst window
  constexpr int substeps_per_bio_dt = 12;  // 5 s trapezoidal quadrature for PDE
  constexpr Real sub_dt = bio_dt / static_cast<Real>(substeps_per_bio_dt);

  const Real c_gf_base = gf.concentration_bounded(src, dst, params);
  const Real d_eff = params.diff_coeff / params.retardation;
  BurstErrorSummary out;
  out.peak_c_qssa = c_gf_base;

  Real hazard_qssa_60 = 0.0;
  Real hazard_pde_60 = 0.0;
  Real hazard_qssa_300 = 0.0;
  Real hazard_pde_300 = 0.0;

  for (int step = 0; step < n_steps; ++step) {
    const Real t_step_start = static_cast<Real>(step) * bio_dt;
    // GutIBM FixBacteriocin evaluates QSSA at the step with exp(-age / tau).
    const Real c_qssa = c_gf_base * std::exp(-t_step_start / release_tau);
    const Real occ_qssa = receptor_occupancy(c_qssa, kd_apparent);
    const Real step_hazard_qssa = k_kill * occ_qssa * bio_dt;
    if (step == 0) hazard_qssa_60 += step_hazard_qssa;
    hazard_qssa_300 += step_hazard_qssa;

    // Integrate dynamic PDE occupancy over [t_step_start, t_step_start + bio_dt].
    Real step_hazard_pde = 0.0;
    for (int sub = 0; sub < substeps_per_bio_dt; ++sub) {
      const Real t0 = t_step_start + static_cast<Real>(sub) * sub_dt;
      const Real t1 = t0 + sub_dt;
      const Real c0 = transient_pde_concentration(
          src, dst, t0, params.source_rate, release_tau,
          d_eff, ux, params.decay_rate, h);
      const Real c1 = transient_pde_concentration(
          src, dst, t1, params.source_rate, release_tau,
          d_eff, ux, params.decay_rate, h);
      out.peak_c_pde = std::max({out.peak_c_pde, c0, c1});
      const Real occ0 = receptor_occupancy(c0, kd_apparent);
      const Real occ1 = receptor_occupancy(c1, kd_apparent);
      step_hazard_pde += 0.5 * (occ0 + occ1) * k_kill * sub_dt;

      const Real occ_err = std::abs(occ_qssa - occ1);
      if (step == 0) {
        out.max_occ_err_60s = std::max(out.max_occ_err_60s, occ_err);
      }
      out.max_occ_err_300s = std::max(out.max_occ_err_300s, occ_err);
    }
    if (step == 0) hazard_pde_60 += step_hazard_pde;
    hazard_pde_300 += step_hazard_pde;
  }

  out.pkill_qssa_60s = 1.0 - std::exp(-hazard_qssa_60);
  out.pkill_pde_60s = 1.0 - std::exp(-hazard_pde_60);
  out.pkill_qssa_300s = 1.0 - std::exp(-hazard_qssa_300);
  out.pkill_pde_300s = 1.0 - std::exp(-hazard_pde_300);
  return out;
}

void test_steady_state_convergence() {
  // Sanity check: constant-rate transient PDE at large t must converge to
  // GutIBM's GreensFunction::concentration_bounded.
  DomainConfig dom_cfg;
  dom_cfg.lo = {0.0, 0.0, 0.0};
  dom_cfg.hi = {200e-6, 200e-6, 100e-6};
  dom_cfg.grid_dx = 5e-6;
  Domain dom;
  dom.init(dom_cfg);

  AdvectionConfig adv_cfg;
  adv_cfg.radial_turnover = 1.0e12;
  adv_cfg.mucus_thickness = 100e-6;
  adv_cfg.distal_length = 1.0e-3;
  adv_cfg.distal_transit_time = 1000.0;  // uniform wall-parallel ux = 1 um/s
  adv_cfg.profile_alpha = 0.0;
  adv_cfg.peristaltic_enabled = false;
  AdvectionField adv;
  adv.init(adv_cfg, dom);

  GreensFunction gf;
  gf.init(dom, adv);

  const PlasmidEntry* colE2_entry = PlasmidLibrary::find("ColE2");
  assert(colE2_entry != nullptr);
  const BICluster colE2 = colE2_entry->cluster;
  GreensFunctionParams params;
  params.diff_coeff = colE2.diff_coeff;
  params.source_rate = (colE2.burst_size / AVOGADRO) / 300.0;
  params.retardation = colE2.retardation;
  params.decay_rate = std::numbers::ln2 / colE2.protease_half_life;
  const Real d_eff = params.diff_coeff / params.retardation;

  const Vec3 src = {100e-6, 100e-6, 50e-6};
  const Vec3 dst = {115e-6, 100e-6, 50e-6};
  const Real c_gf = gf.concentration_bounded(src, dst, params);
  const Real c_pde_ss = transient_pde_concentration(
      src, dst, 2.0e4, params.source_rate, 0.0, d_eff,
      adv.velocity(src)[0], params.decay_rate, 100e-6);

  const Real rel_diff = std::abs(c_pde_ss - c_gf) / std::max(c_gf, 1.0e-30);
  assert(std::isfinite(c_gf) && c_gf > 0.0);
  assert(rel_diff < 1.0e-3);
  std::cout << "  PASS: steady-state PDE limit matches GreensFunction"
            << " (rel_diff=" << std::scientific << std::setprecision(2)
            << rel_diff << ")\n" << std::fixed;
}

struct GeometrySpec {
  const char* label;
  Real lx;
  Real ly;
  Real lz;
  Real z_src;
  std::vector<Real> distances_um;
};

void run_geometry_audit(const GeometrySpec& geom) {
  DomainConfig dom_cfg;
  dom_cfg.lo = {0.0, 0.0, 0.0};
  dom_cfg.hi = {geom.lx, geom.ly, geom.lz};
  dom_cfg.grid_dx = 5e-6;
  Domain dom;
  dom.init(dom_cfg);

  // Uniform wall-parallel flow matching colonic mucosal distal velocity
  // (~0.83 um/s).
  AdvectionConfig adv_cfg;
  adv_cfg.radial_turnover = 1.0e12;
  adv_cfg.mucus_thickness = geom.lz;
  adv_cfg.distal_length = 1.0e-3;
  adv_cfg.distal_transit_time = 1200.0;  // ux = 1e-3 / 1200 = 0.833 um/s
  adv_cfg.profile_alpha = 0.0;
  adv_cfg.peristaltic_enabled = false;
  AdvectionField adv;
  adv.init(adv_cfg, dom);

  GreensFunction gf;
  gf.init(dom, adv);

  const SimulationConfig def = InputParser::default_config();
  const Real release_tau = def.fixes.bacteriocin.burst_release_tau;
  const Real kd_uncompeted = def.fixes.receptor.kd_colicinE_btuB;
  const Real kd_shipped =
      kd_uncompeted * (1.0 + 1.0e-3 / def.fixes.receptor.kd_b12_btuB);

  const std::array<const char*, 2> plasmids = {"ColE1", "ColE2"};
  const Vec3 src = {0.5 * geom.lx, 0.5 * geom.ly, geom.z_src};
  const Real ux = adv.velocity(src)[0];

  std::cout << "\n  Geometry: " << geom.label
            << " (" << (geom.lx * 1e6) << "x" << (geom.ly * 1e6)
            << "x" << (geom.lz * 1e6) << " um, ux=" << (ux * 1e6)
            << " um/s, release_tau=" << release_tau << " s)\n";

  // Track ColE1 vs ColE2 transient errors at the shared mid-field distance
  // (20 um) to assert the physical diffusivity ordering (ColE2 < ColE1).
  Real colE1_pkill_err_20um = 0.0;
  Real colE2_pkill_err_20um = 0.0;

  for (const char* plasmid_name : plasmids) {
    const PlasmidEntry* entry = PlasmidLibrary::find(plasmid_name);
    assert(entry != nullptr);
    const BICluster bi = entry->cluster;
    GreensFunctionParams params;
    params.diff_coeff = bi.diff_coeff;
    params.source_rate = (bi.burst_size / AVOGADRO) / release_tau;
    params.retardation = bi.retardation;
    params.decay_rate = std::numbers::ln2 / bi.protease_half_life;
    const Real k_kill = def.fixes.receptor.kill_rate_colicin;

    Real worst_occ_60 = 0.0;
    Real worst_occ_300 = 0.0;
    Real worst_pkill_60 = 0.0;
    Real worst_pkill_300 = 0.0;
    Real worst_pkill_300_shipped = 0.0;

    for (const Real dist_um : geom.distances_um) {
      const Real dist = dist_um * 1.0e-6;
      // Downstream target (comet-tail direction, highest QSSA concentration).
      const Vec3 dst = {src[0] + dist, src[1], src[2]};
      const BurstErrorSummary uncomp = evaluate_burst_error(
          gf, params, src, dst, ux, geom.lz, release_tau, kd_uncompeted, k_kill);
      const BurstErrorSummary shipped = evaluate_burst_error(
          gf, params, src, dst, ux, geom.lz, release_tau, kd_shipped, k_kill);

      const Real d_pkill_60 = std::abs(uncomp.pkill_qssa_60s - uncomp.pkill_pde_60s);
      const Real d_pkill_300 = std::abs(uncomp.pkill_qssa_300s - uncomp.pkill_pde_300s);
      const Real d_pkill_300_ship =
          std::abs(shipped.pkill_qssa_300s - shipped.pkill_pde_300s);

      assert(std::isfinite(uncomp.max_occ_err_60s));
      assert(std::isfinite(uncomp.max_occ_err_300s));
      assert(uncomp.max_occ_err_60s >= 0.0 && uncomp.max_occ_err_60s <= 1.0);
      assert(uncomp.max_occ_err_300s >= 0.0 && uncomp.max_occ_err_300s <= 1.0);
      assert(d_pkill_60 >= 0.0 && d_pkill_60 <= 1.0);
      assert(d_pkill_300 >= 0.0 && d_pkill_300 <= 1.0);

      // Causality check: on the first 60 s step, instantaneous QSSA cannot
      // underestimate the burst-origin causal PDE kill probability.
      assert(uncomp.pkill_qssa_60s + 1.0e-9 >= uncomp.pkill_pde_60s);

      worst_occ_60 = std::max(worst_occ_60, uncomp.max_occ_err_60s);
      worst_occ_300 = std::max(worst_occ_300, uncomp.max_occ_err_300s);
      worst_pkill_60 = std::max(worst_pkill_60, d_pkill_60);
      worst_pkill_300 = std::max(worst_pkill_300, d_pkill_300);
      worst_pkill_300_shipped = std::max(worst_pkill_300_shipped, d_pkill_300_ship);

      if (std::abs(dist_um - 20.0) < 1.0e-6) {
        if (std::string(plasmid_name) == "ColE1") {
          colE1_pkill_err_20um = d_pkill_300;
        } else {
          colE2_pkill_err_20um = d_pkill_300;
        }
      }

      std::cout << "    " << plasmid_name
                << " r=" << std::setw(4) << std::setprecision(1) << dist_um << " um"
                << " | occ_err(60s)=" << std::setprecision(4) << uncomp.max_occ_err_60s
                << " occ_err(300s)=" << uncomp.max_occ_err_300s
                << " | Pkill_QSSA(60s)=" << uncomp.pkill_qssa_60s
                << " Pkill_PDE(60s)=" << uncomp.pkill_pde_60s
                << " | Pkill_QSSA(300s)=" << uncomp.pkill_qssa_300s
                << " Pkill_PDE(300s)=" << uncomp.pkill_pde_300s
                << " | dPkill_shipped(300s)=" << d_pkill_300_ship << "\n";
    }

    // Also audit transverse (+y) and upstream (-x) directions at 10 um & 20 um
    // to quantify comet-tail directional dependence.
    for (const Real dist_um : {10.0, 20.0}) {
      const Real dist = dist_um * 1.0e-6;
      const std::array<std::pair<const char*, Vec3>, 2> off_axis = {{
          {"transverse", {src[0], src[1] + dist, src[2]}},
          {"upstream  ", {src[0] - dist, src[1], src[2]}},
      }};
      for (const auto& [dir_label, dst] : off_axis) {
        const BurstErrorSummary uncomp = evaluate_burst_error(
            gf, params, src, dst, ux, geom.lz, release_tau, kd_uncompeted, k_kill);
        assert(std::isfinite(uncomp.max_occ_err_60s));
        assert(uncomp.pkill_qssa_60s + 1.0e-9 >= uncomp.pkill_pde_60s);
        std::cout << "    " << plasmid_name << " " << dir_label
                  << " r=" << std::setw(4) << std::setprecision(1) << dist_um << " um"
                  << " | occ_err(60s)=" << std::setprecision(4) << uncomp.max_occ_err_60s
                  << " occ_err(300s)=" << uncomp.max_occ_err_300s
                  << " | Pkill_QSSA(60s)=" << uncomp.pkill_qssa_60s
                  << " Pkill_PDE(60s)=" << uncomp.pkill_pde_60s
                  << " | Pkill_QSSA(300s)=" << uncomp.pkill_qssa_300s
                  << " Pkill_PDE(300s)=" << uncomp.pkill_pde_300s << "\n";
      }
    }

    // Physical sanity bounds:
    // - Weakly retarded ColE2 (D_eff ~ 1.7e-11 m^2/s) equilibrates on ~10-25 s
    //   scales, so its cumulative 300 s kill-probability error remains modest.
    if (std::string(plasmid_name) == "ColE2") {
      assert(worst_pkill_300 < 0.25);
    }
    // - Strongly retarded ColE1 (D_eff ~ 7.0e-13 m^2/s, diffusion time ~100-500 s)
    //   exhibits substantial transient lag at r >= 10 um under uncompeted Kd,
    //   while shipped B12-competed Kd suppresses occupancy and shrinks dPkill.
    if (std::string(plasmid_name) == "ColE1") {
      assert(worst_occ_60 > 0.05);
      assert(worst_pkill_300_shipped < worst_pkill_300);
    }
  }

  // Diffusivity ordering: higher-diffusivity ColE2 must have smaller transient
  // QSSA error than retarded ColE1 at r = 20 um.
  assert(colE2_pkill_err_20um < colE1_pkill_err_20um);
}

}  // namespace

int main() {
  std::cout << "=== QSSA Transient vs Dynamic PDE Bounds ===\n";
  test_steady_state_convergence();

  const GeometrySpec geom_60 = {
      "60um_compact", 60e-6, 60e-6, 50e-6, 25e-6, {5.0, 10.0, 20.0, 25.0}};
  const GeometrySpec geom_200 = {
      "200um_mucus", 200e-6, 200e-6, 100e-6, 50e-6, {5.0, 10.0, 20.0, 30.0, 50.0}};

  run_geometry_audit(geom_60);
  run_geometry_audit(geom_200);

  std::cout << "\nAll QSSA transient vs PDE bound checks passed.\n";
  return 0;
}
