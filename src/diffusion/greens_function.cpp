/* -----------------------------------------------------------------------
   GutIBM – Analytical Green's function implementation
   
   Steady-state solution to:  ∇·(D∇C) - U·∇C + Q*delta(x-x_s) = 0
   
   In 3D with uniform advection U = (Ux, 0, Uz):
     C(r) = (Q / 4π D_eff |r|) * exp( (U · (r - r_s)) / (2 D_eff) )
          * exp( -|U|·|r - r_s| / (2 D_eff) )
   
   This is the infinite-domain fundamental solution.
   Bounded domains use Method of Images to enforce ∂C/∂n=0 at walls.
   ----------------------------------------------------------------------- */

#include "greens_function.h"
#include "domain.h"
#include "advection.h"
#include <cmath>
#include <algorithm>

namespace gutibm {

void GreensFunction::init(const Domain& domain, const AdvectionField& adv) {
  domain_ = &domain;
  adv_    = &adv;
  z_lo_   = domain.lo()[2];
  z_hi_   = domain.hi()[2];
}

Real GreensFunction::single_kernel(const Vec3& src, const Vec3& tgt,
                                    Real D_eff, Real Q,
                                    const Vec3& flow_vel) const {
  Vec3 delta = domain_->min_image_delta(src, tgt);
  Real r = std::sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);

  if (r < 1.0e-9) return 0.0;  // self-interaction cutoff

  Real U_mag = std::sqrt(flow_vel[0]*flow_vel[0] +
                          flow_vel[1]*flow_vel[1] +
                          flow_vel[2]*flow_vel[2]);

  // Dot product U · delta (downstream projection)
  Real U_dot_r = flow_vel[0]*delta[0] + flow_vel[1]*delta[1] + flow_vel[2]*delta[2];

  // Advection-diffusion Green's function (3D steady-state)
  Real prefactor = Q / (4.0 * PI * D_eff * r);
  Real exponent  = (U_dot_r - U_mag * r) / (2.0 * D_eff);

  // Clamp exponent to avoid overflow
  exponent = std::max(exponent, -500.0);

  return prefactor * std::exp(exponent);
}

Real GreensFunction::concentration(const Vec3& source, const Vec3& target,
                                    const GreensFunctionParams& params) const {
  Real D_eff = params.diff_coeff / params.retardation;
  Vec3 flow  = adv_->velocity(source);
  return single_kernel(source, target, D_eff, params.source_rate, flow);
}

Real GreensFunction::concentration_bounded(const Vec3& source, const Vec3& target,
                                            const GreensFunctionParams& params) const {
  Real D_eff = params.diff_coeff / params.retardation;
  Vec3 flow  = adv_->velocity(source);
  Real Q     = params.source_rate;

  Real total = single_kernel(source, target, D_eff, Q, flow);

  // Method of Images: reflect source across z=z_lo and z=z_hi
  for (int n = 1; n <= N_IMAGES; ++n) {
    // Image below z_lo
    Vec3 img_lo = source;
    img_lo[2] = 2.0 * z_lo_ - source[2] - 2.0 * n * (z_hi_ - z_lo_);
    total += single_kernel(img_lo, target, D_eff, Q, flow);

    // Image above z_hi
    Vec3 img_hi = source;
    img_hi[2] = 2.0 * z_hi_ - source[2] + 2.0 * n * (z_hi_ - z_lo_);
    total += single_kernel(img_hi, target, D_eff, Q, flow);

    // Image reflected below then above
    Vec3 img_lo2 = source;
    img_lo2[2] = 2.0 * z_lo_ - source[2] + 2.0 * n * (z_hi_ - z_lo_);
    total += single_kernel(img_lo2, target, D_eff, Q, flow);

    // Image reflected above then below
    Vec3 img_hi2 = source;
    img_hi2[2] = 2.0 * z_hi_ - source[2] - 2.0 * n * (z_hi_ - z_lo_);
    total += single_kernel(img_hi2, target, D_eff, Q, flow);
  }

  return std::max(total, 0.0);
}

void GreensFunction::superpose_to_grid(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    std::vector<Real>& grid_conc,
    Real cutoff_radius) const {

  Int nx = domain_->nx(), ny = domain_->ny(), nz = domain_->nz();
  Int ncells = domain_->ncells();

  grid_conc.assign(ncells, 0.0);

  // For each source, only contribute to grid cells within cutoff
  for (size_t s = 0; s < sources.size(); ++s) {
    const Vec3& src = sources[s];
    const GreensFunctionParams& p = params[s];

    Int src_ix, src_iy, src_iz;
    domain_->pos_to_grid(src, src_ix, src_iy, src_iz);

    Int span = static_cast<Int>(std::ceil(cutoff_radius / domain_->dx()));

    for (Int dz = -span; dz <= span; ++dz) {
      Int iz = src_iz + dz;
      if (iz < 0 || iz >= nz) continue;
      for (Int dy = -span; dy <= span; ++dy) {
        Int iy = src_iy + dy;
        // Handle PBC in y
        if (domain_->config().periodic[1]) {
          iy = ((iy % ny) + ny) % ny;
        } else if (iy < 0 || iy >= ny) continue;

        for (Int dx = -span; dx <= span; ++dx) {
          Int ix = src_ix + dx;
          if (domain_->config().periodic[0]) {
            ix = ((ix % nx) + nx) % nx;
          } else if (ix < 0 || ix >= nx) continue;

          Vec3 tgt = domain_->cell_center(ix, iy, iz);
          Real c = concentration_bounded(src, tgt, p);
          Int idx = domain_->cell_index(ix, iy, iz);
          grid_conc[idx] += c;
        }
      }
    }
  }
}

Real GreensFunction::peclet(const Vec3& pos, Real D_eff, Real length_scale) const {
  Vec3 v = adv_->velocity(pos);
  Real U = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  return (U * length_scale) / D_eff;
}

}  // namespace gutibm
