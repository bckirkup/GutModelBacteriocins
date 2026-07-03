/* -----------------------------------------------------------------------
   GutIBM – Analytical Green's function kernels for QSSA diffusion
   
   Instead of explicit FTCS grid-based PDE solvers, we compute
   steady-state concentration fields from point sources using
   analytical solutions to the advection-diffusion equation.
   
   For a point source of strength Q (mol/s) in a uniform flow U
   with effective diffusivity D, the steady-state concentration is:
   
     C(r) = (Q / 4*pi*D*r) * exp(-U*(r - x_downstream) / (2*D))
   
   where r is distance from source and x_downstream is the
   downstream projection.  This gives the characteristic
   "comet-tail" pattern without any grid timestep constraint.
   
   Method of Images enforces no-flux boundaries at z=0 (epithelium)
   and z=h (lumen), preventing artificial mass loss.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_GREENS_FUNCTION_H
#define GUTIBM_GREENS_FUNCTION_H

#include "types.h"
#include <vector>

namespace gutibm {

class Domain;
class AdvectionField;

struct GreensFunctionParams {
  Real diff_coeff;       // effective D (m^2/s) = D_free / retardation
  Real source_rate;      // Q (mol/s)
  Real pI;               // isoelectric point (determines retardation)
  Real retardation;      // mucin retardation factor

  // NOTE: bacteriocin pI classification lives in a single source of truth,
  // `classify_by_pI()` in src/genome/plasmid.h (pI > 8.5 → CORE, pI < 7.0 →
  // HALO, else NEUTRAL). The Green's function / QSSA code consumes the
  // pre-computed `BICluster.bclass` and must never re-classify from pI here.
};

class GreensFunction {
 public:
  GreensFunction() = default;

  void init(const Domain& domain, const AdvectionField& adv);

  // Steady-state concentration at `target` from a point source at `source`
  // using the advection-diffusion Green's function.
  Real concentration(const Vec3& source, const Vec3& target,
                     const GreensFunctionParams& params) const;

  // Same but with Method of Images for bounded z-domain
  Real concentration_bounded(const Vec3& source, const Vec3& target,
                              const GreensFunctionParams& params) const;

  // Superpose contributions from multiple sources onto grid cells.
  // Uses spatial hashing cutoff to limit O(N*M) to O(N*k).
  void superpose_to_grid(const std::vector<Vec3>& sources,
                          const std::vector<GreensFunctionParams>& params,
                          std::vector<Real>& grid_conc,
                          Real cutoff_radius) const;

  // Peclet number at position: Pe = U*L/D
  Real peclet(const Vec3& pos, Real D_eff, Real length_scale) const;

 private:
  void require_init() const;

  // Single image contribution
  Real single_kernel(const Vec3& src, const Vec3& tgt,
                      Real D_eff, Real Q,
                      const Vec3& flow_vel) const;

  const Domain* domain_    = nullptr;
  const AdvectionField* adv_ = nullptr;

  // Method of Images parameters
  Real z_lo_ = 0.0;
  Real z_hi_ = 100.0e-6;
  static constexpr int N_IMAGES = 3;  // number of image pairs
};

}  // namespace gutibm

#endif  // GUTIBM_GREENS_FUNCTION_H
