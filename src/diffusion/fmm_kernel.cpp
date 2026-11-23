/* -----------------------------------------------------------------------
   GutIBM – FMM kernel Taylor coefficient implementation
   ----------------------------------------------------------------------- */

#include "fmm_kernel.h"
#include <array>
#include <cmath>
#include <limits>

namespace gutibm {

namespace fmm_detail {

int coeff_index(int ox, int oy, int oz, int order);
Real binomial(int n, int k);

// Central finite-difference step for a k-th order derivative, scaled to the
// problem's characteristic length. The total error of a 2nd-order central
// stencil for the k-th derivative is ~ (h/L)^2 (truncation) + eps*(L/h)^k
// (round-off), minimized at h ≈ L * eps^{1/(k+2)}. A single tiny step (the
// former 1.0e-9) makes the order-3 stencil divide by ~2h^3 = 2e-27, so
// round-off is amplified ~1e11 — catastrophic cancellation. (Spec 0 §3.)
Real fd_step(int deriv_order, Real length_scale) {
  const Real eps = std::numeric_limits<Real>::epsilon();
  Real rel = std::cbrt(eps);  // order 1: eps^(1/3)
  if (deriv_order == 2) {
    rel = std::pow(eps, 0.25);  // eps^(1/4)
  } else if (deriv_order >= 3) {
    rel = std::pow(eps, 0.2);  // eps^(1/5)
  }
  const Real L = (length_scale > 0.0) ? length_scale : 1.0;
  return L * rel;
}

Real separation_length(const Vec3& a, const Vec3& b) {
  const Real dx = a[0] - b[0];
  const Real dy = a[1] - b[1];
  const Real dz = a[2] - b[2];
  const Real r = std::sqrt(dx * dx + dy * dy + dz * dz);
  // Floor keeps the step positive for the degenerate (coincident) case, which
  // the well-separated FMM never actually evaluates.
  constexpr Real kMinLength = 1.0e-9;
  return (r > kMinLength) ? r : kMinLength;
}

template<typename Fn>
void for_each_multi_index(int order, Fn&& fn) {
  int idx = 0;
  for (int total = 0; total <= order; ++total) {
    for (int iz = 0; iz <= total; ++iz) {
      for (int iy = 0; iy <= total - iz; ++iy) {
        const int ix = total - iy - iz;
        fn(ix, iy, iz, idx++);
      }
    }
  }
}

Real dr_monomial(Real base, int ix, int iy, int iz, const Vec3& dr) {
  Real term = base;
  for (int d = 0; d < ix; ++d) term *= dr[0];
  for (int d = 0; d < iy; ++d) term *= dr[1];
  for (int d = 0; d < iz; ++d) term *= dr[2];
  return term;
}

Real delta_power(int exp_x, int exp_y, int exp_z, const Vec3& delta) {
  Real pow = 1.0;
  for (int d = 0; d < exp_x; ++d) pow *= delta[0];
  for (int d = 0; d < exp_y; ++d) pow *= delta[1];
  for (int d = 0; d < exp_z; ++d) pow *= delta[2];
  return pow;
}

Real shift_term_contribution(const std::vector<Real>& coeffs,
                             int order,
                             int aix, int aiy, int aiz,
                             int bix, int biy, int biz,
                             const Vec3& delta) {
  if (bix > aix || biy > aiy || biz > aiz) return 0.0;
  const int bidx = coeff_index(bix, biy, biz, order);
  const Real coeff = binomial(aix, bix) * binomial(aiy, biy) * binomial(aiz, biz);
  const Real pow = delta_power(aix - bix, aiy - biy, aiz - biz, delta);
  return coeff * coeffs[bidx] * pow;
}

template<typename Fn>
void for_b_multi_index(int atotal, Fn&& fn) {
  for (int btotal = 0; btotal <= atotal; ++btotal) {
    for (int biz = 0; biz <= btotal; ++biz) {
      for (int biy = 0; biy <= btotal - biz; ++biy) {
        const int bix = btotal - biy - biz;
        fn(bix, biy, biz);
      }
    }
  }
}

Real accumulate_shift_sum(const std::vector<Real>& coeffs,
                          int order,
                          int atotal,
                          int aix, int aiy, int aiz,
                          const Vec3& delta) {
  Real sum = 0.0;
  for_b_multi_index(atotal, [&sum, &coeffs, order, aix, aiy, aiz, delta](int bix, int biy, int biz) {
    sum += shift_term_contribution(coeffs, order, aix, aiy, aiz,
                                   bix, biy, biz, delta);
  });
  return sum;
}

void shift_expansion(std::vector<Real>& out,
                     const std::vector<Real>& in,
                     int order,
                     const Vec3& delta) {
  for_each_multi_index(order, [&out, order, &in, delta](int aix, int aiy, int aiz, int aidx) {
    out[aidx] = accumulate_shift_sum(in, order, aix + aiy + aiz,
                                     aix, aiy, aiz, delta);
  });
}

void shift_local_coeffs(std::vector<Real>& out,
                        const std::vector<Real>& in,
                        int order,
                        const Vec3& delta) {
  // L2L polynomial re-expansion:
  //   sum_beta L_beta ((x - c_child) + delta)^beta
  //     = sum_alpha [ sum_{beta >= alpha} C(beta, alpha) L_beta delta^{beta - alpha} ] (x - c_child)^alpha
  // where delta = c_child - c_parent.
  for_each_multi_index(order, [&out, order, &in, delta](int aix, int aiy, int aiz, int aidx) {
    Real sum = 0.0;
    for_each_multi_index(order, [&sum, &in, aix, aiy, aiz, delta](int bix, int biy, int biz, int bidx) {
      if (bix < aix || biy < aiy || biz < aiz) return;
      const Real coeff = binomial(bix, aix) * binomial(biy, aiy) * binomial(biz, aiz);
      const Real pow = delta_power(bix - aix, biy - aiy, biz - aiz, delta);
      sum += coeff * in[bidx] * pow;
    });
    out[aidx] = sum;
  });
}

template <typename FieldFn>
void fill_first_derivatives(std::vector<Real>& values,
                            int order,
                            const Vec3& center,
                            Real h,
                            FieldFn&& field) {
  for (int d = 0; d < 3; ++d) {
    Vec3 plus = center;
    Vec3 minus = center;
    plus[d] += h;
    minus[d] -= h;
    const Real deriv = (field(plus) - field(minus)) / (2.0 * h);
    const int idx = coeff_index(d == 0 ? 1 : 0, d == 1 ? 1 : 0, d == 2 ? 1 : 0,
                                order);
    values[idx] = deriv;
  }
}

template <typename FieldFn>
Real second_derivative(int d1, int d2,
                       const Vec3& center,
                       Real center_val,
                       Real h,
                       FieldFn&& field) {
  if (d1 == d2) {
    Vec3 plus = center;
    Vec3 minus = center;
    plus[d1] += h;
    minus[d1] -= h;
    return (field(plus) - 2.0 * center_val + field(minus)) / (h * h);
  }

  Vec3 pp = center;
  Vec3 pm = center;
  Vec3 mp = center;
  Vec3 mm = center;
  pp[d1] += h; pp[d2] += h;
  pm[d1] += h; pm[d2] -= h;
  mp[d1] -= h; mp[d2] += h;
  mm[d1] -= h; mm[d2] -= h;
  return (field(pp) - field(pm) - field(mp) + field(mm)) / (4.0 * h * h);
}

template <typename FieldFn>
void fill_second_derivatives(std::vector<Real>& values,
                             int order,
                             const Vec3& center,
                             Real center_val,
                             Real h,
                             const FieldFn& field) {
  for (int d1 = 0; d1 < 3; ++d1) {
    for (int d2 = d1; d2 < 3; ++d2) {
      const Real val = second_derivative(d1, d2, center, center_val, h, field);
      const int ox = (d1 == 0 ? 1 : 0) + (d2 == 0 ? 1 : 0);
      const int oy = (d1 == 1 ? 1 : 0) + (d2 == 1 ? 1 : 0);
      const int oz = (d1 == 2 ? 1 : 0) + (d2 == 2 ? 1 : 0);
      const int idx = coeff_index(ox, oy, oz, order);
      // Taylor coefficient 1/(ox! oy! oz!): 0.5 for diagonal (d1 == d2),
      // 1.0 for upper-triangular mixed partial (d1 < d2).
      const Real scale = (d1 == d2) ? 0.5 : 1.0;
      values[idx] = val * scale;
    }
  }
}

template <typename FieldFn>
void fill_third_derivatives(std::vector<Real>& values,
                            int order,
                            const Vec3& center,
                            Real h,
                            FieldFn&& field,
                            Real scale = 1.0 / 6.0) {
  for (int d = 0; d < 3; ++d) {
    Vec3 plus2 = center;
    Vec3 minus2 = center;
    plus2[d] += 2.0 * h;
    minus2[d] -= 2.0 * h;
    Vec3 plus1 = center;
    Vec3 minus1 = center;
    plus1[d] += h;
    minus1[d] -= h;
    const Real val = (field(plus2) - 2.0 * field(plus1)
                      + 2.0 * field(minus1) - field(minus2))
                     / (2.0 * h * h * h);
    const int ox = (d == 0 ? 3 : 0);
    const int oy = (d == 1 ? 3 : 0);
    const int oz = (d == 2 ? 3 : 0);
    const int idx = coeff_index(ox, oy, oz, order);
    values[idx] = val * scale;
  }
}

int num_coefficients(int order) {
  int p = order + 3;
  return (p * (p - 1) * (p - 2)) / 6;
}

int coeff_index(int ox, int oy, int oz, int order) {
  int result = 0;
  bool found = false;
  for_each_multi_index(order, [&result, ox, oy, oz, &found](int ix, int iy, int iz, int idx) {
    if (found) return;
    if (ix == ox && iy == oy && iz == oz) {
      result = idx;
      found = true;
    }
  });
  return result;
}

void multi_index(int idx, int order, int& ox, int& oy, int& oz) {
  int cursor = 0;
  bool found = false;
  for_each_multi_index(order, [&ox, &oy, &oz, &cursor, &found, idx](int ix, int iy, int iz, int /*coeff_idx*/) {
    if (found) return;
    if (cursor == idx) {
      ox = ix;
      oy = iy;
      oz = iz;
      found = true;
    }
    cursor++;
  });
  if (!found) {
    ox = 0;
    oy = 0;
    oz = 0;
  }
}

static Real factorial(int n) {
  Real f = 1.0;
  for (int i = 2; i <= n; ++i) f *= static_cast<Real>(i);
  return f;
}

Real binomial(int n, int k) {
  if (k < 0 || k > n) return 0.0;
  return factorial(n) / (factorial(k) * factorial(n - k));
}

void add_particle(std::vector<Real>& moments,
                  int order,
                  Real charge,
                  const Vec3& position,
                  const Vec3& center) {
  if (const int n = num_coefficients(order); static_cast<int>(moments.size()) < n) {
    moments.assign(n, 0.0);
  }

  const Vec3 dr = {position[0] - center[0],
                   position[1] - center[1],
                   position[2] - center[2]};

  for_each_multi_index(order, [&moments, charge, dr](int ix, int iy, int iz, int idx) {
    moments[idx] += dr_monomial(charge, ix, iy, iz, dr);
  });
}

std::vector<Real> shift_moments(const std::vector<Real>& child,
                                int order,
                                const Vec3& child_center,
                                const Vec3& parent_center) {
  const int n = num_coefficients(order);
  std::vector parent(n, 0.0);

  const Vec3 delta = {child_center[0] - parent_center[0],
                      child_center[1] - parent_center[1],
                      child_center[2] - parent_center[2]};
  shift_expansion(parent, child, order, delta);
  return parent;
}

void add_shifted_moments(std::vector<Real>& parent,
                         const std::vector<Real>& child,
                         int order,
                         const Vec3& child_center,
                         const Vec3& parent_center) {
  std::vector<Real> shifted = shift_moments(child, order, child_center, parent_center);
  for (size_t i = 0; i < shifted.size(); ++i) parent[i] += shifted[i];
}

}  // namespace fmm_detail

Real fd_axis_derivative(const std::function<Real(const Vec3&)>& field,
                        const Vec3& center, int axis, int order,
                        Real length_scale) {
  const Real h = fmm_detail::fd_step(order, length_scale);
  const auto at = [&field, &center, axis](Real offset) {
    Vec3 p = center;
    p[axis] += offset;
    return field(p);
  };
  if (order == 1) {
    return (at(h) - at(-h)) / (2.0 * h);
  }
  if (order == 2) {
    return (at(h) - 2.0 * field(center) + at(-h)) / (h * h);
  }
  // order 3 (five-point antisymmetric stencil)
  return (at(2.0 * h) - 2.0 * at(h) + 2.0 * at(-h) - at(-2.0 * h)) / (2.0 * h * h * h);
}

KernelTaylorCoeffs kernel_taylor_at_source(
    const GreensFunction& gf,
    const Vec3& source,
    const Vec3& target,
    const GreensFunctionParams& params,
    int max_order) {
  KernelTaylorCoeffs out;
  out.max_order = max_order;
  out.values.assign(fmm_detail::num_coefficients(max_order), 0.0);

  GreensFunctionParams unit = params;
  unit.source_rate = 1.0;

  const Real L = fmm_detail::separation_length(source, target);
  const auto eval = [&gf, &unit, target](const Vec3& src) {
    return gf.concentration_bounded(src, target, unit);
  };

  out.values[0] = eval(source);
  if (max_order < 1) return out;

  fmm_detail::fill_first_derivatives(out.values, max_order, source, fmm_detail::fd_step(1, L), eval);
  if (max_order < 2) return out;

  fmm_detail::fill_second_derivatives(
      out.values, max_order, source, out.values[0], fmm_detail::fd_step(2, L), eval);
  if (max_order < 3) return out;

  fmm_detail::fill_third_derivatives(out.values, max_order, source, fmm_detail::fd_step(3, L), eval);
  return out;
}

Real evaluate_multipole(const std::vector<Real>& moments,
                        int order,
                        const KernelTaylorCoeffs& kernel) {
  Real sum = 0.0;
  const int n = fmm_detail::num_coefficients(order);
  for (int i = 0; i < n; ++i) sum += moments[i] * kernel.values[i];
  return sum;
}

Real evaluate_multipole_at_source(const GreensFunction& gf,
                                  const std::vector<Real>& moments,
                                  int order,
                                  const Vec3& source,
                                  const Vec3& target,
                                  const GreensFunctionParams& params) {
  GreensFunctionParams unit = params;
  unit.source_rate = 1.0;

  const Real L = fmm_detail::separation_length(source, target);
  const auto eval = [&gf, &unit, target](const Vec3& src) {
    return gf.concentration_bounded(src, target, unit);
  };

  const Real g0 = eval(source);
  Real sum = moments[0] * g0;
  if (order < 1) return sum;

  // Order 1: skip identically zero centroid dipole moments.
  const Real h1 = fmm_detail::fd_step(1, L);
  for (int d = 0; d < 3; ++d) {
    const Real m1 = moments[1 + d];
    if (std::abs(m1) <= 0.0) continue;
    Vec3 plus = source;
    Vec3 minus = source;
    plus[d] += h1;
    minus[d] -= h1;
    sum += m1 * ((eval(plus) - eval(minus)) / (2.0 * h1));
  }
  if (order < 2) return sum;

  // Order 2: skip zero quadrupole moments and reuse g0.
  const Real h2 = fmm_detail::fd_step(2, L);
  const std::array<int, 3> diag_idx = {4, 6, 9};
  for (int d = 0; d < 3; ++d) {
    const Real m2 = moments[static_cast<size_t>(diag_idx[static_cast<size_t>(d)])];
    if (std::abs(m2) <= 0.0) continue;
    sum += m2 * 0.5 * fmm_detail::second_derivative(d, d, source, g0, h2, eval);
  }
  const std::array<int, 3> md1 = {0, 0, 1};
  const std::array<int, 3> md2 = {1, 2, 2};
  const std::array<int, 3> midx = {5, 7, 8};
  for (int m = 0; m < 3; ++m) {
    const Real m2 = moments[static_cast<size_t>(midx[static_cast<size_t>(m)])];
    if (std::abs(m2) <= 0.0) continue;
    sum += m2 * fmm_detail::second_derivative(
        md1[static_cast<size_t>(m)], md2[static_cast<size_t>(m)], source, g0, h2, eval);
  }
  if (order < 3) return sum;

  // Order 3: pure single-axis derivatives with 1/3! = 1/6 Taylor factor.
  const Real h3 = fmm_detail::fd_step(3, L);
  const std::array<int, 3> third_idx = {10, 13, 19};
  for (int d = 0; d < 3; ++d) {
    const Real m3 = moments[static_cast<size_t>(third_idx[static_cast<size_t>(d)])];
    if (std::abs(m3) <= 0.0) continue;
    Vec3 p2 = source;
    Vec3 m2 = source;
    Vec3 p1 = source;
    Vec3 m1 = source;
    p2[d] += 2.0 * h3;
    m2[d] -= 2.0 * h3;
    p1[d] += h3;
    m1[d] -= h3;
    const Real deriv = (eval(p2) - 2.0 * eval(p1) + 2.0 * eval(m1) - eval(m2))
                     / (12.0 * h3 * h3 * h3);
    sum += m3 * deriv;
  }
  return sum;
}

std::vector<Real> multipole_to_local(
    const std::vector<Real>& moments,
    int order,
    const Vec3& source_center,
    const Vec3& target_center,
    const GreensFunction& gf,
    const GreensFunctionParams& avg_params) {
  const int n = fmm_detail::num_coefficients(order);
  std::vector local(n, 0.0);

  const Real L = fmm_detail::separation_length(source_center, target_center);
  const auto multipole_field = [&gf, &moments, order, source_center, avg_params](const Vec3& eval_target) {
    return evaluate_multipole_at_source(
        gf, moments, order, source_center, eval_target, avg_params);
  };

  local[0] = multipole_field(target_center);
  if (order < 1) return local;

  fmm_detail::fill_first_derivatives(local, order, target_center, fmm_detail::fd_step(1, L), multipole_field);
  if (order < 2) return local;

  fmm_detail::fill_second_derivatives(
      local, order, target_center, local[0], fmm_detail::fd_step(2, L), multipole_field);
  if (order < 3) return local;

  fmm_detail::fill_third_derivatives(local, order, target_center, fmm_detail::fd_step(3, L), multipole_field);
  return local;
}

std::vector<Real> shift_local_expansion(const std::vector<Real>& parent_local,
                                        int order,
                                        const Vec3& parent_center,
                                        const Vec3& child_center) {
  const int n = fmm_detail::num_coefficients(order);
  std::vector child(n, 0.0);

  const Vec3 delta = {child_center[0] - parent_center[0],
                      child_center[1] - parent_center[1],
                      child_center[2] - parent_center[2]};
  fmm_detail::shift_local_coeffs(child, parent_local, order, delta);
  return child;
}

Real evaluate_local(const std::vector<Real>& local,
                    int order,
                    const Vec3& center,
                    const Vec3& target) {
  const Vec3 dr = {target[0] - center[0],
                   target[1] - center[1],
                   target[2] - center[2]};

  Real sum = 0.0;
  fmm_detail::for_each_multi_index(order, [&sum, &local, dr](int ix, int iy, int iz, int idx) {
    sum += fmm_detail::dr_monomial(local[idx], ix, iy, iz, dr);
  });
  return sum;
}

}  // namespace gutibm
