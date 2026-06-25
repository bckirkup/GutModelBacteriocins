/* -----------------------------------------------------------------------
   GutIBM – FMM kernel Taylor coefficients for the QSSA Green's function

   Provides concentration_bounded() and its source-position partial
   derivatives up to order 3, used by the kernel-independent FMM.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FMM_KERNEL_H
#define GUTIBM_FMM_KERNEL_H

#include "types.h"
#include "greens_function.h"
#include <vector>

namespace gutibm {

struct KernelTaylorCoeffs {
  // Flat storage for multi-index (ix, iy, iz) with ix+iy+iz <= max_order.
  // Indexing matches FMM::num_coefficients / fmm_detail::coeff_index.
  std::vector<Real> values;
  int max_order = 0;
};

// Evaluate G and source-position Taylor coefficients at (source, target).
// coeffs.values[k] = (1/(ix! iy! iz!)) * d^{ix+iy+iz} G / d sx^ix d sy^iy d sz^iz
KernelTaylorCoeffs kernel_taylor_at_source(
    const GreensFunction& gf,
    const Vec3& source,
    const Vec3& target,
    const GreensFunctionParams& params,
    int max_order);

// Evaluate multipole expansion: sum_alpha M_alpha * kernel_taylor_alpha.
Real evaluate_multipole(const std::vector<Real>& moments,
                        int order,
                        const KernelTaylorCoeffs& kernel);

// Convert remote multipole at source_center into local expansion at target_center.
std::vector<Real> multipole_to_local(const std::vector<Real>& moments,
                                     int order,
                                     const Vec3& source_center,
                                     const Vec3& target_center,
                                     const GreensFunction& gf,
                                     const GreensFunctionParams& avg_params);

// Shift local expansion from parent_center to child_center.
std::vector<Real> shift_local_expansion(const std::vector<Real>& parent_local,
                                        int order,
                                        const Vec3& parent_center,
                                        const Vec3& child_center);

// Evaluate local expansion at a point.
Real evaluate_local(const std::vector<Real>& local,
                    int order,
                    const Vec3& center,
                    const Vec3& target);

namespace fmm_detail {

int num_coefficients(int order);
int coeff_index(int ox, int oy, int oz, int order);
void multi_index(int idx, int order, int& ox, int& oy, int& oz);

// P2M: accumulate particle charge at position into moments about center.
void add_particle(std::vector<Real>& moments,
                  int order,
                  Real charge,
                  const Vec3& position,
                  const Vec3& center);

// M2M: shift child moments to parent center.
std::vector<Real> shift_moments(const std::vector<Real>& child,
                                int order,
                                const Vec3& child_center,
                                const Vec3& parent_center);

// Add shifted child moments into parent.
void add_shifted_moments(std::vector<Real>& parent,
                         const std::vector<Real>& child,
                         int order,
                         const Vec3& child_center,
                         const Vec3& parent_center);

}  // namespace fmm_detail

}  // namespace gutibm

#endif  // GUTIBM_FMM_KERNEL_H
