/* -----------------------------------------------------------------------
   GutIBM – FMM kernel Taylor coefficient implementation
   ----------------------------------------------------------------------- */

#include "fmm_kernel.h"
#include <cmath>

namespace gutibm {

namespace fmm_detail {

int coeff_index(int ox, int oy, int oz, int order);
Real binomial(int n, int k);

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
  for_b_multi_index(atotal, [&](int bix, int biy, int biz) {
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
                       Real h,
                       FieldFn&& field) {
  if (d1 == d2) {
    Vec3 plus = center;
    Vec3 minus = center;
    plus[d1] += h;
    minus[d1] -= h;
    return (field(plus) - 2.0 * field(center) + field(minus)) / (h * h);
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
                             Real h,
                             FieldFn&& field,
                             Real scale = 1.0) {
  for (int d1 = 0; d1 < 3; ++d1) {
    for (int d2 = d1; d2 < 3; ++d2) {
      const Real val = second_derivative(d1, d2, center, h, field);
      const int ox = (d1 == 0 ? 1 : 0) + (d2 == 0 ? 1 : 0);
      const int oy = (d1 == 1 ? 1 : 0) + (d2 == 1 ? 1 : 0);
      const int oz = (d1 == 2 ? 1 : 0) + (d2 == 2 ? 1 : 0);
      const int idx = coeff_index(ox, oy, oz, order);
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
                            Real scale = 1.0) {
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
  for_each_multi_index(order, [&result, ox, oy, oz, &found, order](int ix, int iy, int iz, int idx) {
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
  for_each_multi_index(order, [&ox, &oy, &oz, &cursor, &found, order, idx](int ix, int iy, int iz, int /*coeff_idx*/) {
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

  for_each_multi_index(order, [&moments, order, charge, dr](int ix, int iy, int iz, int idx) {
    moments[idx] += dr_monomial(charge, ix, iy, iz, dr);
  });
}

std::vector<Real> shift_moments(const std::vector<Real>& child,
                                int order,
                                const Vec3& child_center,
                                const Vec3& parent_center) {
  const int n = num_coefficients(order);
  std::vector<Real> parent(n, 0.0);

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

  const Real h = 1.0e-9;
  const auto eval = [&gf, &unit, target](const Vec3& src) {
    return gf.concentration_bounded(src, target, unit);
  };

  out.values[0] = eval(source);
  if (max_order < 1) return out;

  fmm_detail::fill_first_derivatives(out.values, max_order, source, h, eval);
  if (max_order < 2) return out;

  fmm_detail::fill_second_derivatives(out.values, max_order, source, h, eval);
  if (max_order < 3) return out;

  fmm_detail::fill_third_derivatives(out.values, max_order, source, h, eval);
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

std::vector<Real> multipole_to_local(
    const std::vector<Real>& moments,
    int order,
    const Vec3& source_center,
    const Vec3& target_center,
    const GreensFunction& gf,
    const GreensFunctionParams& avg_params) {
  const int n = fmm_detail::num_coefficients(order);
  std::vector<Real> local(n, 0.0);

  const Real h = 1.0e-9;
  const auto multipole_field = [&gf, &moments, order, source_center, avg_params](const Vec3& eval_target) {
    KernelTaylorCoeffs k = kernel_taylor_at_source(
        gf, source_center, eval_target, avg_params, order);
    return evaluate_multipole(moments, order, k);
  };

  local[0] = multipole_field(target_center);
  if (order < 1) return local;

  fmm_detail::fill_first_derivatives(local, order, target_center, h, multipole_field);
  if (order < 2) return local;

  fmm_detail::fill_second_derivatives(local, order, target_center, h, multipole_field, 0.5);
  if (order < 3) return local;

  fmm_detail::fill_third_derivatives(local, order, target_center, h, multipole_field, 1.0 / 6.0);
  return local;
}

std::vector<Real> shift_local_expansion(const std::vector<Real>& parent_local,
                                        int order,
                                        const Vec3& parent_center,
                                        const Vec3& child_center) {
  const int n = fmm_detail::num_coefficients(order);
  std::vector<Real> child(n, 0.0);

  const Vec3 delta = {child_center[0] - parent_center[0],
                      child_center[1] - parent_center[1],
                      child_center[2] - parent_center[2]};
  fmm_detail::shift_expansion(child, parent_local, order, delta);
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
  fmm_detail::for_each_multi_index(order, [&sum, order, &local, dr](int ix, int iy, int iz, int idx) {
    sum += fmm_detail::dr_monomial(local[idx], ix, iy, iz, dr);
  });
  return sum;
}

}  // namespace gutibm
