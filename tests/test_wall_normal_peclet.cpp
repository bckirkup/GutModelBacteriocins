/* -----------------------------------------------------------------------
   GutIBM – Wall-normal Peclet guard, auto-dispatch, and high-Pe_z
   Sturm-Liouville / FVM verification benchmark (Issue #478)
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "chemical_field.h"
#include "domain.h"
#include "error.h"
#include "greens_function.h"
#include "greens_function_gpu.h"
#include "input_parser.h"
#include "neumann_image_series.h"
#include "path_utils.h"
#include "plasmid.h"
#include "qssa_solver.h"
#include "robin_correction_table.h"
#include "simulation.h"
#include "species_names.h"

#ifdef GUTIBM_HDF5
#include <hdf5.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numbers>
#include <string>
#include <string_view>
#include <vector>

using namespace gutibm;

namespace {

constexpr Real kHeight = 100.0e-6;
constexpr Real kSourceRate = 1.0e-18;
constexpr Real kDecayRate = 5.0e-5;
constexpr Real kRadialTurnover = 5400.0;

void require(bool condition, const std::string& message) {
  if (condition) return;
  std::cerr << "FAIL: " << message << "\n";
  std::exit(1);
}

struct SlabSystem {
  Domain domain;
  AdvectionField adv;
  GreensFunction gf;
};

SlabSystem make_uniform_flow_system(Real flow_z, Real flow_x = 0.0,
                                    Real grid_dx = 5.0e-6) {
  SlabSystem sys;
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {1.0e-3, 1.0e-3, kHeight};
  dcfg.grid_dx = grid_dx;
  sys.domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.crypts_enabled = false;
  acfg.peristaltic_enabled = false;
  acfg.mucus_thickness = kHeight;
  acfg.profile_alpha = 0.0;
  acfg.distal_length = flow_x == 0.0 ? 0.0 : 1.0e-3;
  acfg.distal_transit_time = flow_x == 0.0 ? 1.0 : acfg.distal_length / flow_x;
  acfg.radial_turnover = flow_z == 0.0
      ? std::numeric_limits<Real>::infinity()
      : kHeight / flow_z;
  sys.adv.init(acfg, sys.domain);
  sys.gf.init(sys.domain, sys.adv);
  return sys;
}

SlabSystem make_profile_flow_system(Real grid_dx = 5.0e-6) {
  SlabSystem sys;
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {1.0e-3, 1.0e-3, kHeight};
  dcfg.grid_dx = grid_dx;
  sys.domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.crypts_enabled = false;
  acfg.peristaltic_enabled = false;
  acfg.mucus_thickness = kHeight;
  acfg.profile_alpha = 1.5;
  acfg.distal_length = 0.05;
  acfg.distal_transit_time = 43200.0;
  acfg.radial_turnover = kRadialTurnover;
  sys.adv.init(acfg, sys.domain);
  sys.gf.init(sys.domain, sys.adv);
  return sys;
}

GreensFunctionParams make_toxin_params(Real d_free, Real retardation,
                                       bool drift_correction) {
  GreensFunctionParams p;
  p.diff_coeff = d_free;
  p.retardation = retardation;
  p.source_rate = kSourceRate;
  p.decay_rate = kDecayRate;
  p.lumen_transfer_length = robin::kZeroTransferLength;
  p.image_series_relative_tolerance = 1.0e-12;
  p.drift_correction = drift_correction;
  return p;
}

// Converged N-mode Sturm-Liouville eigenfunction expansion for the physical
// sealed slab (-D dC/dz + U_z C = 0 at z = 0 and z = H).
Real sturm_liouville_500_reference(Real dx, Real dy, Real z_s, Real z_t,
                                   Real d_eff, Real decay, Real ux, Real uz,
                                   int mode_count = 500) {
  const Real rho = std::hypot(dx, dy);
  const Real a = uz / (2.0 * d_eff);
  const Real kappa0_sq = decay / d_eff + (ux * ux) / (4.0 * d_eff * d_eff);
  const Real kappa_sq = kappa0_sq + a * a;
  const Real kappa0 = std::sqrt(std::max(kappa0_sq, 0.0));
  const Real lateral_rho = std::max(rho, 1.0e-12);
  const Real lat_expo = (ux * dx) / (2.0 * d_eff);
  const Real uz_expo = a * (z_t - z_s);

  Real sum = 0.0;
  if (kappa0 > 0.0) {
    const Real aH = a * kHeight;
    Real z_factor = 0.0;
    if (std::abs(aH) < 1.0e-12) {
      z_factor = std::exp(2.0 * a * z_t) / kHeight;
    } else if (a > 0.0) {
      z_factor = (2.0 * a / (-std::expm1(-2.0 * aH)))
          * std::exp(2.0 * a * (z_t - kHeight));
    } else {
      z_factor = (-2.0 * a / (-std::expm1(2.0 * aH)))
          * std::exp(2.0 * a * z_t);
    }
    const Real arg0 = kappa0 * lateral_rho;
    if (lat_expo - arg0 > -700.0) {
      sum += z_factor * robin::drift_detail::bessel_k0e(arg0)
          * std::exp(lat_expo - arg0);
    }
  }
  for (int n = 1; n <= mode_count; ++n) {
    const Real beta = static_cast<Real>(n) * std::numbers::pi / kHeight;
    const Real kn = std::sqrt(kappa_sq + beta * beta);
    const Real arg = kn * lateral_rho;
    const Real ex = lat_expo + uz_expo - arg;
    if (ex < -700.0) break;
    const Real phi_s = std::cos(beta * z_s) + (a / beta) * std::sin(beta * z_s);
    const Real phi_t = std::cos(beta * z_t) + (a / beta) * std::sin(beta * z_t);
    const Real norm_n = 0.5 * kHeight * (1.0 + (a * a) / (beta * beta));
    sum += phi_s * phi_t * robin::drift_detail::bessel_k0e(arg)
        * std::exp(ex) / norm_n;
  }
  return kSourceRate / (4.0 * std::numbers::pi * d_eff) * (2.0 * sum);
}

// Bernoulli function B(x) = x / (exp(x) - 1) for Scharfetter-Gummel FVM.
Real bernoulli(Real x) {
  if (std::abs(x) < 1.0e-7) {
    return 1.0 - 0.5 * x + x * x / 12.0;
  }
  return x / std::expm1(x);
}

struct CylindricalFvmSolution {
  int nr = 0;
  int nz = 0;
  std::vector<Real> r_center;
  std::vector<Real> z_center;
  std::vector<Real> conc;  // [ir * nz + iz]
  Real max_residual = 0.0;
};

// Independent fine-grid steady-state axisymmetric cylindrical (r, z) FVM PDE
// solver for:
//   -D * (1/r) d/dr(r dC/dr) + d/dz(-D dC/dz + U_z C) + lambda C = Q delta
// with exact Scharfetter-Gummel zero-total-flux walls J_z(0) = J_z(H) = 0.
CylindricalFvmSolution solve_axisymmetric_fvm(Real z_source, Real d_eff,
                                              Real decay, Real flow_z) {
  constexpr int kNr = 96;
  constexpr int kNz = 64;
  constexpr Real kDr0 = 1.0e-6;
  constexpr Real kRMax = 2500.0e-6;
  const Real dz = kHeight / static_cast<Real>(kNz);

  // Geometrically stretched radial faces from 0 to kRMax.
  const Real growth = std::pow(kRMax / kDr0, 1.0 / static_cast<Real>(kNr - 1));
  std::vector<Real> r_face(static_cast<size_t>(kNr + 1), 0.0);
  for (int i = 1; i <= kNr; ++i) {
    r_face[static_cast<size_t>(i)] =
        kDr0 * (std::pow(growth, static_cast<Real>(i)) - 1.0) / (growth - 1.0);
  }
  const Real scale_r = kRMax / r_face.back();
  for (Real& rf : r_face) rf *= scale_r;

  CylindricalFvmSolution sol;
  sol.nr = kNr;
  sol.nz = kNz;
  sol.r_center.resize(kNr);
  sol.z_center.resize(kNz);
  for (int i = 0; i < kNr; ++i) {
    sol.r_center[static_cast<size_t>(i)] =
        0.5 * (r_face[static_cast<size_t>(i)]
               + r_face[static_cast<size_t>(i + 1)]);
  }
  for (int j = 0; j < kNz; ++j) {
    sol.z_center[static_cast<size_t>(j)] = (static_cast<Real>(j) + 0.5) * dz;
  }

  // Exact cell volumes and face conductances.
  std::vector<Real> area_z(kNr);
  std::vector<Real> vol(kNr);
  std::vector<Real> g_r_lo(kNr, 0.0);
  std::vector<Real> g_r_hi(kNr, 0.0);
  for (int i = 0; i < kNr; ++i) {
    const Real r_lo = r_face[static_cast<size_t>(i)];
    const Real r_hi = r_face[static_cast<size_t>(i + 1)];
    area_z[static_cast<size_t>(i)] =
        std::numbers::pi * (r_hi * r_hi - r_lo * r_lo);
    vol[static_cast<size_t>(i)] = area_z[static_cast<size_t>(i)] * dz;
    if (i > 0) {
      const Real dr_lo = sol.r_center[static_cast<size_t>(i)]
          - sol.r_center[static_cast<size_t>(i - 1)];
      g_r_lo[static_cast<size_t>(i)] =
          2.0 * std::numbers::pi * r_lo * dz * d_eff / dr_lo;
    }
    const Real dr_hi = (i + 1 < kNr)
        ? (sol.r_center[static_cast<size_t>(i + 1)]
           - sol.r_center[static_cast<size_t>(i)])
        : (r_hi - sol.r_center[static_cast<size_t>(i)]);
    g_r_hi[static_cast<size_t>(i)] =
        2.0 * std::numbers::pi * r_hi * dz * d_eff / dr_hi;
  }

  // Scharfetter-Gummel vertical flux coefficients per unit area.
  const Real cell_pe = flow_z * dz / d_eff;
  const Real sg_up = (d_eff / dz) * bernoulli(-cell_pe);
  const Real sg_dn = (d_eff / dz) * bernoulli(cell_pe);

  // Place point source at cell (0, j_s) via linear z-interpolation if off-center.
  const Real z_idx_f = z_source / dz - 0.5;
  const int j_lo = std::clamp(static_cast<int>(std::floor(z_idx_f)), 0, kNz - 1);
  const int j_hi = std::min(j_lo + 1, kNz - 1);
  const Real w_hi = std::clamp(z_idx_f - static_cast<Real>(j_lo), 0.0, 1.0);
  const Real w_lo = 1.0 - w_hi;

  // Assemble strictly diagonally dominant Scharfetter-Gummel FVM matrix in
  // banded storage (ordering idx = i * kNz + j, half-bandwidth bw = kNz)
  // and solve directly to machine precision (~1e-15 residual).
  const int n_unknowns = kNr * kNz;
  const int bw = kNz;
  const int band_stride = 2 * bw + 1;
  std::vector<Real> band(static_cast<size_t>(n_unknowns) * band_stride, 0.0);
  std::vector<Real> rhs(static_cast<size_t>(n_unknowns), 0.0);
  auto mat = [&](int row, int col) -> Real& {
    return band[static_cast<size_t>(row) * band_stride
                + static_cast<size_t>(col - row + bw)];
  };

  for (int i = 0; i < kNr; ++i) {
    const Real az = area_z[static_cast<size_t>(i)];
    const Real v = vol[static_cast<size_t>(i)];
    const Real glo = g_r_lo[static_cast<size_t>(i)];
    const Real ghi = g_r_hi[static_cast<size_t>(i)];
    for (int j = 0; j < kNz; ++j) {
      const int row = i * kNz + j;
      Real diag = decay * v + glo + ghi;
      if (i > 0) mat(row, (i - 1) * kNz + j) = -glo;
      if (i + 1 < kNr) mat(row, (i + 1) * kNz + j) = -ghi;
      if (j > 0) {
        diag += az * sg_dn;
        mat(row, i * kNz + (j - 1)) = -az * sg_up;
      }
      if (j + 1 < kNz) {
        diag += az * sg_up;
        mat(row, i * kNz + (j + 1)) = -az * sg_dn;
      }
      mat(row, row) = diag;
      if (i == 0) {
        if (j == j_lo) rhs[static_cast<size_t>(row)] += kSourceRate * w_lo;
        if (j == j_hi && j_hi != j_lo) {
          rhs[static_cast<size_t>(row)] += kSourceRate * w_hi;
        }
      }
    }
  }

  // Keep a copy of RHS for independent residual verification.
  const std::vector<Real> rhs_orig = rhs;
  std::vector<Real> factor = band;
  auto lu = [&](int row, int col) -> Real& {
    return factor[static_cast<size_t>(row) * band_stride
                  + static_cast<size_t>(col - row + bw)];
  };

  for (int k = 0; k < n_unknowns; ++k) {
    const Real pivot = lu(k, k);
    const int i_end = std::min(n_unknowns, k + bw + 1);
    const int j_end = std::min(n_unknowns, k + bw + 1);
    for (int r = k + 1; r < i_end; ++r) {
      const Real mult = lu(r, k) / pivot;
      if (mult == 0.0) continue;
      lu(r, k) = 0.0;
      for (int c = k + 1; c < j_end; ++c) {
        lu(r, c) -= mult * lu(k, c);
      }
      rhs[static_cast<size_t>(r)] -= mult * rhs[static_cast<size_t>(k)];
    }
  }

  sol.conc.assign(static_cast<size_t>(n_unknowns), 0.0);
  for (int r = n_unknowns - 1; r >= 0; --r) {
    Real sum = rhs[static_cast<size_t>(r)];
    const int c_end = std::min(n_unknowns, r + bw + 1);
    for (int c = r + 1; c < c_end; ++c) {
      sum -= lu(r, c) * sol.conc[static_cast<size_t>(c)];
    }
    sol.conc[static_cast<size_t>(r)] = sum / lu(r, r);
  }

  Real max_res = 0.0;
  for (int r = 0; r < n_unknowns; ++r) {
    Real lhs = 0.0;
    const int c_begin = std::max(0, r - bw);
    const int c_end = std::min(n_unknowns, r + bw + 1);
    for (int c = c_begin; c < c_end; ++c) {
      lhs += mat(r, c) * sol.conc[static_cast<size_t>(c)];
    }
    max_res = std::max(max_res, std::abs(lhs - rhs_orig[static_cast<size_t>(r)]));
  }
  sol.max_residual = max_res / kSourceRate;
  return sol;
}

Real relative_l2_error(const std::vector<Real>& actual,
                       const std::vector<Real>& reference) {
  Real num_sq = 0.0;
  Real den_sq = 0.0;
  for (size_t i = 0; i < actual.size(); ++i) {
    const Real diff = actual[i] - reference[i];
    num_sq += diff * diff;
    den_sq += reference[i] * reference[i];
  }
  return std::sqrt(num_sq / std::max(den_sq, 1.0e-300));
}

// 1. Benchmark ColE1 (R=50.2) and ColE2 (R=1.2, 2.04) across near-epithelium,
// mid-depth, and near-lumen point sources against BOTH a 500-mode Sturm-Liouville
// expansion AND an independent fine-grid axisymmetric FVM PDE solver.
void test_cole1_and_cole2_against_sl_and_fvm() {
  struct ToxinCase {
    const char* name;
    Real d_free;
    Real retardation;
    Real flow_z;
  };
  const Real uz_mid = (kHeight / kRadialTurnover) * std::pow(0.5, 1.5);
  const Real uz_lumen = kHeight / kRadialTurnover;
  const std::array<ToxinCase, 4> cases = {{
      {"ColE1_mid_Pe0.822", 4.0e-11, 50.2, uz_mid},
      {"ColE1_lumen_Pe2.325", 4.0e-11, 50.2, uz_lumen},
      {"ColE2_R2.04_lumen", 3.5e-11, 2.04, uz_lumen},
      {"ColE2_R1.2_lumen", 3.5e-11, 1.2, uz_lumen},
  }};
  const std::array<Real, 3> source_heights = {0.1 * kHeight, 0.5 * kHeight, 0.9 * kHeight};
  const std::array<Real, 5> sl_radii = {5.0e-6, 12.0e-6, 25.0e-6, 60.0e-6, 120.0e-6};
  const std::array<Real, 5> sl_z_fracs = {0.0, 0.25, 0.50, 0.75, 1.0};

  Real worst_sl_interior_l2 = 0.0;
  Real worst_sl_wall_l2 = 0.0;
  Real worst_fvm_interior_l2 = 0.0;
  Real worst_fvm_wall_l2 = 0.0;
  Real worst_wall_flux_norm = 0.0;

  for (const auto& tc : cases) {
    const Real d_eff = tc.d_free / tc.retardation;
    const Real pe_z = neumann::wall_normal_peclet(tc.flow_z, kHeight, d_eff);
    auto sys = make_uniform_flow_system(tc.flow_z, 0.0);
    const auto params = make_toxin_params(tc.d_free, tc.retardation, true);
    const auto legacy_params = make_toxin_params(tc.d_free, tc.retardation, false);

    Real case_sl_int = 0.0;
    Real case_sl_wall = 0.0;
    Real case_leg_wall = 0.0;
    Real case_fvm_int = 0.0;
    Real case_fvm_wall = 0.0;
    Real case_flux = 0.0;

    for (const Real z_s : source_heights) {
      const Vec3 src = {500.0e-6, 500.0e-6, z_s};

      // (a) Compare against 500-mode Sturm-Liouville expansion (interior + walls).
      std::vector<Real> gf_int;
      std::vector<Real> sl_int;
      std::vector<Real> gf_wall;
      std::vector<Real> sl_wall;
      std::vector<Real> leg_wall;
      for (const Real rho : sl_radii) {
        for (const Real zf : sl_z_fracs) {
          const Real z_t = zf * kHeight;
          const Vec3 tgt = {src[0] + rho, src[1], z_t};
          const Real actual = sys.gf.concentration_bounded(src, tgt, params);
          const Real ref = sturm_liouville_500_reference(
              rho, 0.0, z_s, z_t, d_eff, kDecayRate, 0.0, tc.flow_z, 500);
          require(std::isfinite(actual) && actual >= 0.0,
                  std::string("non-finite or negative drift concentration for ")
                      + tc.name);
          if (zf == 0.0 || zf == 1.0) {
            gf_wall.push_back(actual);
            sl_wall.push_back(ref);
            leg_wall.push_back(
                sys.gf.concentration_bounded(src, tgt, legacy_params));
          } else {
            gf_int.push_back(actual);
            sl_int.push_back(ref);
          }
        }
      }
      const Real sl_int_l2 = relative_l2_error(gf_int, sl_int);
      const Real sl_wall_l2 = relative_l2_error(gf_wall, sl_wall);
      const Real leg_wall_l2 = relative_l2_error(leg_wall, sl_wall);
      case_sl_int = std::max(case_sl_int, sl_int_l2);
      case_sl_wall = std::max(case_sl_wall, sl_wall_l2);
      case_leg_wall = std::max(case_leg_wall, leg_wall_l2);
      worst_sl_interior_l2 = std::max(worst_sl_interior_l2, sl_int_l2);
      worst_sl_wall_l2 = std::max(worst_sl_wall_l2, sl_wall_l2);
      require(sl_int_l2 <= 1.0e-2 && sl_wall_l2 <= 1.0e-2,
              std::string("Sturm-Liouville L2 exceeded 1% for ") + tc.name
                  + " (int=" + std::to_string(sl_int_l2)
                  + ", wall=" + std::to_string(sl_wall_l2) + ")");

      // Verify legacy image series actually degrades on high-Pe_z ColE1.
      if (pe_z > 0.5) {
        require(leg_wall_l2 > 0.10,
                std::string("legacy ColE1 wall L2 unexpectedly small: ")
                    + std::to_string(leg_wall_l2));
      }

      // (b) Compare against independent fine-grid axisymmetric FVM PDE solver.
      const CylindricalFvmSolution fvm =
          solve_axisymmetric_fvm(z_s, d_eff, kDecayRate, tc.flow_z);
      require(fvm.max_residual < 1.0e-9,
              std::string("FVM solver failed to converge for ") + tc.name
                  + " (res=" + std::to_string(fvm.max_residual) + ")");

      std::vector<Real> fvm_gf_int;
      std::vector<Real> fvm_ref_int;
      std::vector<Real> fvm_gf_wall;
      std::vector<Real> fvm_ref_wall;
      for (int i = 0; i < fvm.nr; ++i) {
        const Real rho = fvm.r_center[static_cast<size_t>(i)];
        // Exclude r < 12 um (point-source discretization zone) and far outer BC.
        if (rho < 12.0e-6 || rho > 200.0e-6) continue;
        for (int j = 0; j < fvm.nz; ++j) {
          const Real z_t = fvm.z_center[static_cast<size_t>(j)];
          const Real ref = fvm.conc[static_cast<size_t>(i * fvm.nz + j)];
          const Real actual = sys.gf.concentration_bounded(
              src, {src[0] + rho, src[1], z_t}, params);
          if (j == 0 || j == fvm.nz - 1) {
            fvm_gf_wall.push_back(actual);
            fvm_ref_wall.push_back(ref);
          } else {
            fvm_gf_int.push_back(actual);
            fvm_ref_int.push_back(ref);
          }
        }
      }
      const Real fvm_int_l2 = relative_l2_error(fvm_gf_int, fvm_ref_int);
      const Real fvm_wall_l2 = relative_l2_error(fvm_gf_wall, fvm_ref_wall);
      case_fvm_int = std::max(case_fvm_int, fvm_int_l2);
      case_fvm_wall = std::max(case_fvm_wall, fvm_wall_l2);
      worst_fvm_interior_l2 = std::max(worst_fvm_interior_l2, fvm_int_l2);
      worst_fvm_wall_l2 = std::max(worst_fvm_wall_l2, fvm_wall_l2);
      require(fvm_int_l2 <= 1.0e-2 && fvm_wall_l2 <= 1.0e-2,
              std::string("FVM L2 exceeded 1% for ") + tc.name
                  + " (int=" + std::to_string(fvm_int_l2)
                  + ", wall=" + std::to_string(fvm_wall_l2) + ")");

      // (c) Discrete zero-total-flux conservation at z = z_lo and z = z_hi:
      // |J_z(wall)| <= 1e-8 * Q / A with A = H^2 (characteristic slab area).
      // Use a 5-point 4th-order one-sided stencil so FD truncation is < 1e-11.
      const Real flux_ref = kSourceRate / (kHeight * kHeight);
      const Real h_fd = 5.0e-9;
      for (const Real rho : {0.0, 10.0e-6, 40.0e-6}) {
        for (const int wall : {0, 1}) {
          const Real zw = wall == 0 ? 0.0 : kHeight;
          std::array<Real, 5> vals{};
          for (int k = 0; k < 5; ++k) {
            const Real z = wall == 0 ? zw + k * h_fd : zw - k * h_fd;
            vals[static_cast<size_t>(k)] = sys.gf.concentration_bounded(
                src, {src[0] + rho, src[1], z}, params);
          }
          const Real dcdz = wall == 0
              ? (-25.0 * vals[0] + 48.0 * vals[1] - 36.0 * vals[2]
                 + 16.0 * vals[3] - 3.0 * vals[4]) / (12.0 * h_fd)
              : (25.0 * vals[0] - 48.0 * vals[1] + 36.0 * vals[2]
                 - 16.0 * vals[3] + 3.0 * vals[4]) / (12.0 * h_fd);
          const Real jz = -d_eff * dcdz + tc.flow_z * vals[0];
          const Real norm_jz = std::abs(jz) / flux_ref;
          case_flux = std::max(case_flux, norm_jz);
          worst_wall_flux_norm = std::max(worst_wall_flux_norm, norm_jz);
          if (!(norm_jz <= 1.0e-8)) {
            std::cerr << "FAIL: wall flux conservation violated for "
                      << tc.name << " (z_s=" << z_s << ", rho=" << rho
                      << ", wall=" << wall << ", |J_z|/(Q/A)=" << norm_jz
                      << ")\n";
            std::exit(1);
          }
        }
      }
    }
    require(sys.gf.negative_field_count() == 0,
            std::string("negative drift concentration detected for ") + tc.name);
    std::cout << "  " << tc.name << " (Pe_z=" << pe_z
              << "): legacy_wall_L2=" << case_leg_wall
              << ", SL int/wall L2=" << case_sl_int << "/" << case_sl_wall
              << ", FVM int/wall L2=" << case_fvm_int << "/" << case_fvm_wall
              << ", max |J_z|/(Q/A)=" << case_flux << "\n";
  }
  std::cout << "  test_cole1_and_cole2_against_sl_and_fvm: PASSED\n";
}

// 2. Verify 3D lateral flow (distal U_x != 0 + wall-normal U_z != 0) upstream,
// transverse, and downstream against 500-mode Sturm-Liouville and wall flux.
void test_3d_anisotropic_distal_and_wall_normal_flow() {
  const Real uz_mid = (kHeight / kRadialTurnover) * std::pow(0.5, 1.5);
  const Real ux_mid = (0.05 / 43200.0) * std::pow(0.5, 1.5);
  auto sys = make_uniform_flow_system(uz_mid, ux_mid);
  const auto params = make_toxin_params(4.0e-11, 50.2, true);
  const Real d_eff = 4.0e-11 / 50.2;
  const Vec3 src = {500.0e-6, 500.0e-6, 0.5 * kHeight};

  const std::array<std::pair<Real, Real>, 5> offsets = {{
      {15.0e-6, 0.0}, {-15.0e-6, 0.0}, {0.0, 15.0e-6},
      {40.0e-6, 10.0e-6}, {-5.0e-6, 20.0e-6},
  }};
  const Real h_fd = 1.0e-8;
  Real worst_bc_res = 0.0;
  for (const auto& [dx, dy] : offsets) {
    const Real c_down = sys.gf.concentration_bounded(
        src, {src[0] + std::abs(dx) + 10.0e-6, src[1] + dy, 0.5 * kHeight},
        params);
    const Real c_up = sys.gf.concentration_bounded(
        src, {src[0] - std::abs(dx) - 10.0e-6, src[1] + dy, 0.5 * kHeight},
        params);
    require(c_down > c_up && c_up >= 0.0,
            "downstream concentration must exceed upstream under U_x > 0");

    for (const int wall : {0, 1}) {
      const Real zw = wall == 0 ? 0.0 : kHeight;
      std::array<Real, 4> vals{};
      for (int k = 0; k < 4; ++k) {
        const Real z = wall == 0 ? zw + k * h_fd : zw - k * h_fd;
        vals[static_cast<size_t>(k)] = sys.gf.concentration_bounded(
            src, {src[0] + dx, src[1] + dy, z}, params);
      }
      if (vals[0] < 1.0e-250) continue;
      const Real dcdz = wall == 0
          ? (-11.0 * vals[0] + 18.0 * vals[1] - 9.0 * vals[2] + 2.0 * vals[3])
              / (6.0 * h_fd)
          : (11.0 * vals[0] - 18.0 * vals[1] + 9.0 * vals[2] - 2.0 * vals[3])
              / (6.0 * h_fd);
      const Real jz = -d_eff * dcdz + uz_mid * vals[0];
      const Real scale = std::abs(d_eff * dcdz) + std::abs(uz_mid * vals[0])
          + d_eff * vals[0] / kHeight;
      worst_bc_res = std::max(worst_bc_res, std::abs(jz) / scale);
    }
  }
  require(worst_bc_res <= 1.0e-6,
          "3D anisotropic wall-flux residual exceeded 1e-6: "
              + std::to_string(worst_bc_res));
  std::cout << "  test_3d_anisotropic_distal_and_wall_normal_flow: PASSED "
               "(max relative wall residual=" << worst_bc_res << ")\n";
}

// 3. Verify runtime Peclet guard behavior: Pe_z <= 0.05 vs Pe_z > 0.05 under
// "auto", "warn", "error", "allow", graded thresholds, and zero-flow control.
void test_runtime_peclet_guard_and_auto_dispatch() {
  auto sys = make_profile_flow_system();

  // Mid-depth ColE1: Pe_z = 0.822 > 0.05; mid-depth ColE2 (R=1.2): Pe_z = 0.022 <= 0.05.
  const Vec3 mid_pos = {500.0e-6, 500.0e-6, 0.5 * kHeight};
  const Vec3 lumen_pos = {500.0e-6, 500.0e-6, 0.9 * kHeight};
  const Vec3 mid_flow = sys.adv.velocity(mid_pos);
  const Real pe_cole1_mid = neumann::wall_normal_peclet(
      mid_flow[2], kHeight, 4.0e-11 / 50.2);
  const Real pe_cole2_mid = neumann::wall_normal_peclet(
      mid_flow[2], kHeight, 3.5e-11 / 1.2);
  require(pe_cole1_mid > 0.8 && pe_cole1_mid < 0.85,
          "unexpected ColE1 mid-depth Pe_z");
  require(pe_cole2_mid > 0.01 && pe_cole2_mid < 0.03,
          "unexpected ColE2 (R=1.2) mid-depth Pe_z");

  // (a) QSSASolver::init policy matrix for ColE1 (Pe_z = 0.822) vs threshold.
  ChemicalSpec spec;
  spec.name = species::BACTERIOCIN_BTUB;
  spec.diff_coeff = 4.0e-11;
  spec.retardation = 50.2;
  spec.decay_rate = kDecayRate;
  const std::vector<ChemicalSpec> cole1_chem = {spec};
  const RuntimeDriftEnvelopeBasis cole1_basis{"ColE1", 4.0e-11 / 50.2};

  // "error" fails closed when Pe_z > threshold, and passes when threshold > Pe_z.
  {
    QSSAConfig cfg;
    cfg.low_screening_policy = "allow";
    cfg.drift_envelope_policy = "error";
    cfg.peclet_z_threshold = 0.05;
    bool threw = false;
    try {
      QSSASolver qssa;
      qssa.init(cfg, sys.domain, sys.adv, false, &cole1_chem, &cole1_basis);
    } catch (const Error& err) {
      const std::string_view msg(err.what());
      threw = msg.find("ColE1") != std::string_view::npos
          && msg.find("Pe_z=") != std::string_view::npos;
    }
    require(threw, "drift_envelope_policy=error did not fail closed naming ColE1 and Pe_z");

    cfg.peclet_z_threshold = 2.0;
    QSSASolver qssa_high_thresh;
    qssa_high_thresh.init(cfg, sys.domain, sys.adv, false, &cole1_chem,
                          &cole1_basis);
  }

  // "auto" and "warn" and "allow" all initialize without throwing.
  for (const char* policy : {"auto", "warn", "allow"}) {
    QSSAConfig cfg;
    cfg.low_screening_policy = "allow";
    cfg.drift_envelope_policy = policy;
    cfg.peclet_z_threshold = 0.05;
    QSSASolver qssa;
    qssa.init(cfg, sys.domain, sys.adv, false, &cole1_chem, &cole1_basis);
  }

  // (b) Per-source auto-dispatch: graded sweep over threshold [0.01, 0.05, 1.0]
  // with mixed sources: ColE1@0.5H (Pe=0.822), ColE2_R1.2@0.5H (Pe=0.022),
  // ColE2_R2.04@0.9H (Pe=0.092).
  const std::vector<Vec3> mixed_sources = {mid_pos, mid_pos, lumen_pos};
  std::vector<GreensFunctionParams> mixed_params = {
      make_toxin_params(4.0e-11, 50.2, false),
      make_toxin_params(3.5e-11, 1.2, false),
      make_toxin_params(3.5e-11, 2.04, false),
  };
  for (auto& p : mixed_params) {
    p.drift_auto_dispatch = true;
  }

  const std::array<std::pair<Real, uint64_t>, 3> threshold_expected = {{
      {0.01, 3ULL},   // all 3 exceed 0.01
      {0.05, 2ULL},   // ColE1@0.5H (0.822) and ColE2@0.9H (0.092) exceed 0.05
      {1.00, 0ULL},   // none exceed 1.0
  }};
  for (const auto& [thresh, expected_dispatched] : threshold_expected) {
    sys.gf.reset_low_screening_diagnostics();
    for (auto& p : mixed_params) p.peclet_z_threshold = thresh;
    std::vector<Real> grid;
    sys.gf.superpose_to_grid(mixed_sources, mixed_params, grid, 30.0e-6);
    require(sys.gf.peclet_auto_dispatched_sources() == expected_dispatched,
            "threshold=" + std::to_string(thresh) + " expected "
                + std::to_string(expected_dispatched)
                + " auto-dispatched sources, got "
                + std::to_string(sys.gf.peclet_auto_dispatched_sources()));
    if (expected_dispatched > 0) {
      require(sys.gf.drift_modal_evaluations() > 0,
              "auto-dispatched sources did not record drift_modal_evaluations");
    } else {
      require(sys.gf.drift_modal_evaluations() == 0,
              "zero auto-dispatched sources should have zero drift_modal_evaluations");
    }
  }

  // (c) Negative control: U_z = 0 never auto-dispatches even with threshold = 0.0.
  auto zero_uz_sys = make_uniform_flow_system(0.0, 1.0e-6);
  for (auto& p : mixed_params) {
    p.peclet_z_threshold = 0.0;
    p.drift_auto_dispatch = true;
  }
  std::vector<Real> zero_grid;
  zero_uz_sys.gf.superpose_to_grid(mixed_sources, mixed_params, zero_grid, 30.0e-6);
  require(zero_uz_sys.gf.peclet_auto_dispatched_sources() == 0,
          "U_z=0 negative control unexpectedly auto-dispatched");
  require(zero_uz_sys.gf.drift_modal_evaluations() == 0,
          "U_z=0 negative control unexpectedly ran modal drift evaluations");

  // (d) Verify robin_host_fallback_sources returns empty for sealed drift-corrected
  // sources (native GPU eligibility, zero forced host fallback).
  std::vector<GreensFunctionParams> sealed_drift_params = {
      make_toxin_params(4.0e-11, 50.2, true),
      make_toxin_params(3.5e-11, 2.04, true),
  };
  const std::vector<Vec3> sealed_sources = {
      {500.0e-6, 500.0e-6, 2.0e-6},
      {500.0e-6, 500.0e-6, 98.0e-6},
  };
  const auto fallback = robin_host_fallback_sources(
      sys.domain, sealed_sources, sealed_drift_params);
  require(fallback.empty(),
          "sealed drift-corrected sources must not fall back to host on GPU");

  std::cout << "  test_runtime_peclet_guard_and_auto_dispatch: PASSED\n";
}

// 4. End-to-end Simulation + HDF5 /run_provenance/ verification under
// qssa.drift_envelope_policy=auto and finite lumen_transfer_length Robin BC.
void test_end_to_end_auto_dispatch_and_provenance() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.time.total_time = 600.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 600.0;
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {80.0e-6, 80.0e-6, 100.0e-6};
  cfg.domain.grid_dx = 10.0e-6;
  cfg.advection.mucus_thickness = 100.0e-6;
  cfg.advection.radial_turnover = kRadialTurnover;
  cfg.advection.distal_length = 0.05;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff = 40.0e-6;
  cfg.qssa.drift_envelope_policy = "auto";
  cfg.qssa.peclet_z_threshold = 0.05;
  cfg.qssa.drift_correction = false;
  cfg.fixes.bacteriocin.sos_lysis_prob = 0.0;
  cfg.fixes.bacteriocin.sos_basal_rate = 1.0;  // force rapid ColE1 SOS induction
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain producer;
  producer.type = 1;
  producer.count = 6;
  producer.mu_max = 5.0e-4;
  producer.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(producer);
  SimulationConfig::InitialStrain bystander;
  bystander.type = 2;
  bystander.count = 6;
  bystander.mu_max = 5.0e-4;
  bystander.receptor_expression["BtuB"] = 0.0;  // resistant bystander keeps run alive
  cfg.initial_strains.push_back(bystander);
  cfg.initial_population.placement = "z_slab";
  cfg.initial_population.z_min = 45.0e-6;
  cfg.initial_population.z_max = 55.0e-6;
  cfg.seed = 123;

#ifdef GUTIBM_HDF5
  const std::string h5_path = resolve_test_h5_path(
      "GUTIBM_TEST_WALL_NORMAL_PECLET_H5", "wall_normal_peclet_provenance");
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = h5_path;
#else
  cfg.hdf5.enabled = false;
#endif

  Simulation sim;
  sim.init(cfg);
  sim.run();

  require(sim.neumann_peclet_auto_dispatched_sources() > 0,
          "ColE1 auto-dispatch simulation recorded 0 peclet_auto_dispatched_sources");
  require(sim.neumann_drift_modal_evaluations() > 0,
          "ColE1 auto-dispatch simulation recorded 0 drift_modal_evaluations");
  require(sim.neumann_negative_field_count() == 0,
          "ColE1 auto-dispatch simulation recorded negative toxin concentration");

#ifdef GUTIBM_HDF5
  hid_t file = H5Fopen(h5_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  require(file >= 0, "failed to open provenance HDF5 file");
  auto read_u64 = [file](const char* path) {
    unsigned long long value = 0;
    hid_t dset = H5Dopen2(file, path, H5P_DEFAULT);
    require(dset >= 0, std::string("missing HDF5 dataset ") + path);
    const herr_t st = H5Dread(dset, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, &value);
    H5Dclose(dset);
    require(st >= 0, std::string("failed reading HDF5 dataset ") + path);
    return static_cast<uint64_t>(value);
  };
  const uint64_t h5_modal = read_u64("run_provenance/drift_modal_evaluations");
  const uint64_t h5_dispatched =
      read_u64("run_provenance/peclet_auto_dispatched_sources");
  H5Fclose(file);
  std::remove(h5_path.c_str());
  require(h5_modal == sim.neumann_drift_modal_evaluations() && h5_modal > 0,
          "HDF5 drift_modal_evaluations mismatch");
  require(h5_dispatched == sim.neumann_peclet_auto_dispatched_sources()
              && h5_dispatched > 0,
          "HDF5 peclet_auto_dispatched_sources mismatch");
#endif

  // Also check finite lumen_transfer_length compatibility under auto-dispatch.
  auto robin_sys = make_uniform_flow_system(
      (kHeight / kRadialTurnover) * std::pow(0.5, 1.5), 0.0);
  auto robin_params = make_toxin_params(4.0e-11, 50.2, false);
  robin_params.drift_auto_dispatch = true;
  robin_params.lumen_transfer_length = 100.0e-6;
  const Vec3 src = {500.0e-6, 500.0e-6, 0.5 * kHeight};
  const Real c_robin = robin_sys.gf.concentration_bounded(
      src, {520.0e-6, 500.0e-6, 0.9 * kHeight}, robin_params);
  auto sealed_params = robin_params;
  sealed_params.lumen_transfer_length = robin::kZeroTransferLength;
  const Real c_sealed = robin_sys.gf.concentration_bounded(
      src, {520.0e-6, 500.0e-6, 0.9 * kHeight}, sealed_params);
  require(std::isfinite(c_robin) && c_robin > 0.0 && c_robin < c_sealed,
          "finite lumen_transfer_length under drift auto-dispatch did not reduce near-lumen field");

  std::cout << "  test_end_to_end_auto_dispatch_and_provenance: PASSED "
               "(dispatched=" << sim.neumann_peclet_auto_dispatched_sources()
            << ", modal_evals=" << sim.neumann_drift_modal_evaluations()
            << ")\n";
}

// 5. Config parser validation for qssa.drift_envelope_policy=auto and
// qssa.peclet_z_threshold.
void test_config_parser_peclet_keys() {
  SimulationConfig cfg = InputParser::default_config();
  require(InputParser::apply_flat_key(cfg, "qssa.drift_envelope_policy", "auto"),
          "failed to parse qssa.drift_envelope_policy=auto");
  require(cfg.qssa.drift_envelope_policy == "auto",
          "qssa.drift_envelope_policy did not store 'auto'");
  require(InputParser::apply_flat_key(cfg, "qssa.peclet_z_threshold", "0.025"),
          "failed to parse qssa.peclet_z_threshold");
  require(std::abs(cfg.qssa.peclet_z_threshold - 0.025) < 1.0e-15,
          "qssa.peclet_z_threshold value mismatch");

  bool rejected_neg = false;
  try {
    InputParser::apply_flat_key(cfg, "qssa.peclet_z_threshold", "-0.1");
  } catch (const ConfigError&) {
    rejected_neg = true;
  }
  require(rejected_neg, "negative qssa.peclet_z_threshold was not rejected");

  bool rejected_policy = false;
  try {
    InputParser::apply_flat_key(cfg, "qssa.drift_envelope_policy", "invalid");
  } catch (const ConfigError&) {
    rejected_policy = true;
  }
  require(rejected_policy, "invalid qssa.drift_envelope_policy was not rejected");

  std::cout << "  test_config_parser_peclet_keys: PASSED\n";
}

}  // namespace

int main() {
  std::cout << "=== Wall-Normal Peclet Guard & High-Pe_z Kernel Tests (#478) ===\n";
  test_config_parser_peclet_keys();
  test_cole1_and_cole2_against_sl_and_fvm();
  test_3d_anisotropic_distal_and_wall_normal_flow();
  test_runtime_peclet_guard_and_auto_dispatch();
  test_end_to_end_auto_dispatch_and_provenance();
  std::cout << "All wall-normal Peclet tests passed.\n";
  return 0;
}
