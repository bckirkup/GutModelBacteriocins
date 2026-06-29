/* -----------------------------------------------------------------------
   GutIBM – FMM kernel Taylor coefficient implementation
   ----------------------------------------------------------------------- */

#include "fmm_kernel.h"
#include <cmath>

namespace gutibm {

namespace fmm_detail {

int num_coefficients(int order) {
  int p = order + 3;
  return (p * (p - 1) * (p - 2)) / 6;
}

int coeff_index(int ox, int oy, int oz, int order) {
  int idx = 0;
  for (int total = 0; total <= order; ++total) {
    for (int iz = 0; iz <= total; ++iz) {
      for (int iy = 0; iy <= total - iz; ++iy) {
        if (int ix = total - iy - iz; ix == ox && iy == oy && iz == oz) return idx;
        idx++;
      }
    }
  }
  return idx;
}

void multi_index(int idx, int order, int& ox, int& oy, int& oz) {
  int cursor = 0;
  for (int total = 0; total <= order; ++total) {
    for (int iz = 0; iz <= total; ++iz) {
      for (int iy = 0; iy <= total - iz; ++iy) {
        int ix = total - iy - iz;
        if (cursor == idx) {
          ox = ix;
          oy = iy;
          oz = iz;
          return;
        }
        cursor++;
      }
    }
  }
  ox = oy = oz = 0;
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

  Vec3 dr = {position[0] - center[0],
             position[1] - center[1],
             position[2] - center[2]};

  for (int total = 0; total <= order; ++total) {
    for (int iz = 0; iz <= total; ++iz) {
      for (int iy = 0; iy <= total - iz; ++iy) {
        int ix = total - iy - iz;
        int idx = coeff_index(ix, iy, iz, order);
        Real term = charge;
        for (int d = 0; d < ix; ++d) term *= dr[0];
        for (int d = 0; d < iy; ++d) term *= dr[1];
        for (int d = 0; d < iz; ++d) term *= dr[2];
        moments[idx] += term;
      }
    }
  }
}

std::vector<Real> shift_moments(const std::vector<Real>& child,
                                int order,
                                const Vec3& child_center,
                                const Vec3& parent_center) {
  const int n = num_coefficients(order);
  std::vector<Real> parent(n, 0.0);

  Vec3 delta = {child_center[0] - parent_center[0],
                child_center[1] - parent_center[1],
                child_center[2] - parent_center[2]};

  for (int atotal = 0; atotal <= order; ++atotal) {
    for (int aiz = 0; aiz <= atotal; ++aiz) {
      for (int aiy = 0; aiy <= atotal - aiz; ++aiy) {
        int aix = atotal - aiy - aiz;
        int aidx = coeff_index(aix, aiy, aiz, order);
        Real sum = 0.0;

        for (int btotal = 0; btotal <= atotal; ++btotal) {
          for (int biz = 0; biz <= btotal; ++biz) {
            for (int biy = 0; biy <= btotal - biz; ++biy) {
              int bix = btotal - biy - biz;
              if (bix > aix || biy > aiy || biz > aiz) continue;
              int bidx = coeff_index(bix, biy, biz, order);
              Real coeff = binomial(aix, bix) * binomial(aiy, biy) * binomial(aiz, biz);
              Real pow = 1.0;
              for (int d = 0; d < aix - bix; ++d) pow *= delta[0];
              for (int d = 0; d < aiy - biy; ++d) pow *= delta[1];
              for (int d = 0; d < aiz - biz; ++d) pow *= delta[2];
              sum += coeff * child[bidx] * pow;
            }
          }
        }
        parent[aidx] = sum;
      }
    }
  }
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

  auto eval = [&](const Vec3& src) {
    return gf.concentration_bounded(src, target, unit);
  };

  out.values[0] = eval(source);

  if (max_order < 1) return out;

  for (int d = 0; d < 3; ++d) {
    Vec3 sp = source;
    Vec3 sm = source;
    sp[d] += h;
    sm[d] -= h;
    Real dG = (eval(sp) - eval(sm)) / (2.0 * h);
    int idx = fmm_detail::coeff_index(d == 0 ? 1 : 0, d == 1 ? 1 : 0, d == 2 ? 1 : 0,
                                      max_order);
    out.values[idx] = dG;
  }

  if (max_order < 2) return out;

  for (int d1 = 0; d1 < 3; ++d1) {
    for (int d2 = d1; d2 < 3; ++d2) {
      Real val;
      if (d1 == d2) {
        Vec3 sp = source;
        Vec3 sm = source;
        sp[d1] += h;
        sm[d1] -= h;
        val = (eval(sp) - 2.0 * eval(source) + eval(sm)) / (h * h);
      } else {
        Vec3 pp = source;
        Vec3 pm = source;
        Vec3 mp = source;
        Vec3 mm = source;
        pp[d1] += h; pp[d2] += h;
        pm[d1] += h; pm[d2] -= h;
        mp[d1] -= h; mp[d2] += h;
        mm[d1] -= h; mm[d2] -= h;
        val = (eval(pp) - eval(pm) - eval(mp) + eval(mm)) / (4.0 * h * h);
      }
      int ox = (d1 == 0 ? 1 : 0) + (d2 == 0 ? 1 : 0);
      int oy = (d1 == 1 ? 1 : 0) + (d2 == 1 ? 1 : 0);
      int oz = (d1 == 2 ? 1 : 0) + (d2 == 2 ? 1 : 0);
      int idx = fmm_detail::coeff_index(ox, oy, oz, max_order);
      out.values[idx] = val;
    }
  }

  if (max_order < 3) return out;

  for (int d = 0; d < 3; ++d) {
    Vec3 sp = source;
    Vec3 sm = source;
    sp[d] += 2.0 * h;
    sm[d] -= 2.0 * h;
    Vec3 sp1 = source;
    Vec3 sm1 = source;
    sp1[d] += h;
    sm1[d] -= h;
    Real val = (eval(sp) - 2.0 * eval(sp1) + 2.0 * eval(sm1) - eval(sm))
               / (2.0 * h * h * h);
    int ox = (d == 0 ? 3 : 0);
    int oy = (d == 1 ? 3 : 0);
    int oz = (d == 2 ? 3 : 0);
    int idx = fmm_detail::coeff_index(ox, oy, oz, max_order);
    out.values[idx] = val;
  }

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

  auto multipole_field = [&](const Vec3& eval_target) {
    KernelTaylorCoeffs k = kernel_taylor_at_source(
        gf, source_center, eval_target, avg_params, order);
    return evaluate_multipole(moments, order, k);
  };

  local[0] = multipole_field(target_center);

  if (order < 1) return local;

  for (int d = 0; d < 3; ++d) {
    Vec3 tp = target_center;
    Vec3 tm = target_center;
    tp[d] += h;
    tm[d] -= h;
    Real deriv = (multipole_field(tp) - multipole_field(tm)) / (2.0 * h);
    int idx = fmm_detail::coeff_index(d == 0 ? 1 : 0, d == 1 ? 1 : 0, d == 2 ? 1 : 0,
                                      order);
    local[idx] = deriv;
  }

  if (order < 2) return local;

  for (int d1 = 0; d1 < 3; ++d1) {
    for (int d2 = d1; d2 < 3; ++d2) {
      Real val;
      if (d1 == d2) {
        Vec3 tp = target_center;
        Vec3 tm = target_center;
        tp[d1] += h;
        tm[d1] -= h;
        val = (multipole_field(tp) - 2.0 * multipole_field(target_center)
               + multipole_field(tm)) / (h * h);
      } else {
        Vec3 pp = target_center;
        Vec3 pm = target_center;
        Vec3 mp = target_center;
        Vec3 mm = target_center;
        pp[d1] += h; pp[d2] += h;
        pm[d1] += h; pm[d2] -= h;
        mp[d1] -= h; mp[d2] += h;
        mm[d1] -= h; mm[d2] -= h;
        val = (multipole_field(pp) - multipole_field(pm)
               - multipole_field(mp) + multipole_field(mm)) / (4.0 * h * h);
      }
      int ox = (d1 == 0 ? 1 : 0) + (d2 == 0 ? 1 : 0);
      int oy = (d1 == 1 ? 1 : 0) + (d2 == 1 ? 1 : 0);
      int oz = (d1 == 2 ? 1 : 0) + (d2 == 2 ? 1 : 0);
      int idx = fmm_detail::coeff_index(ox, oy, oz, order);
      local[idx] = val / 2.0;
    }
  }

  if (order < 3) return local;

  for (int d = 0; d < 3; ++d) {
    Vec3 tp = target_center;
    Vec3 tm = target_center;
    tp[d] += 2.0 * h;
    tm[d] -= 2.0 * h;
    Vec3 tp1 = target_center;
    Vec3 tm1 = target_center;
    tp1[d] += h;
    tm1[d] -= h;
    Real val = (multipole_field(tp) - 2.0 * multipole_field(tp1)
                + 2.0 * multipole_field(tm1) - multipole_field(tm))
               / (2.0 * h * h * h);
    int ox = (d == 0 ? 3 : 0);
    int oy = (d == 1 ? 3 : 0);
    int oz = (d == 2 ? 3 : 0);
    int idx = fmm_detail::coeff_index(ox, oy, oz, order);
    local[idx] = val / 6.0;
  }

  return local;
}

std::vector<Real> shift_local_expansion(const std::vector<Real>& parent_local,
                                        int order,
                                        const Vec3& parent_center,
                                        const Vec3& child_center) {
  const int n = fmm_detail::num_coefficients(order);
  std::vector<Real> child(n, 0.0);

  Vec3 delta = {child_center[0] - parent_center[0],
                child_center[1] - parent_center[1],
                child_center[2] - parent_center[2]};

  for (int atotal = 0; atotal <= order; ++atotal) {
    for (int aiz = 0; aiz <= atotal; ++aiz) {
      for (int aiy = 0; aiy <= atotal - aiz; ++aiy) {
        int aix = atotal - aiy - aiz;
        int aidx = fmm_detail::coeff_index(aix, aiy, aiz, order);
        Real sum = 0.0;

        for (int btotal = 0; btotal <= atotal; ++btotal) {
          for (int biz = 0; biz <= btotal; ++biz) {
            for (int biy = 0; biy <= btotal - biz; ++biy) {
              int bix = btotal - biy - biz;
              if (bix > aix || biy > aiy || biz > aiz) continue;
              int bidx = fmm_detail::coeff_index(bix, biy, biz, order);
              Real coeff = fmm_detail::binomial(aix, bix)
                         * fmm_detail::binomial(aiy, biy)
                         * fmm_detail::binomial(aiz, biz);
              Real pow = 1.0;
              for (int d = 0; d < aix - bix; ++d) pow *= delta[0];
              for (int d = 0; d < aiy - biy; ++d) pow *= delta[1];
              for (int d = 0; d < aiz - biz; ++d) pow *= delta[2];
              sum += coeff * parent_local[bidx] * pow;
            }
          }
        }
        child[aidx] = sum;
      }
    }
  }
  return child;
}

Real evaluate_local(const std::vector<Real>& local,
                    int order,
                    const Vec3& center,
                    const Vec3& target) {
  Vec3 dr = {target[0] - center[0],
             target[1] - center[1],
             target[2] - center[2]};

  Real sum = 0.0;
  for (int total = 0; total <= order; ++total) {
    for (int iz = 0; iz <= total; ++iz) {
      for (int iy = 0; iy <= total - iz; ++iy) {
        int ix = total - iy - iz;
        int idx = fmm_detail::coeff_index(ix, iy, iz, order);
        Real term = local[idx];
        for (int d = 0; d < ix; ++d) term *= dr[0];
        for (int d = 0; d < iy; ++d) term *= dr[1];
        for (int d = 0; d < iz; ++d) term *= dr[2];
        sum += term;
      }
    }
  }
  return sum;
}

}  // namespace gutibm
