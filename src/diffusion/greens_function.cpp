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
#include "greens_function_gpu.h"
#include "dispatch.h"
#include <cmath>
#include <algorithm>
#include "error.h"
#ifdef GUTIBM_OPENMP
#include <omp.h>
#endif

namespace gutibm {

namespace {

bool in_periodic_grid(Int& idx, Int count, bool periodic) {
  if (periodic) {
    idx = ((idx % count) + count) % count;
    return true;
  }
  return idx >= 0 && idx < count;
}

struct SuperposeGridContext {
  Int nx;
  Int ny;
  Int nz;
  Int span;
  bool periodic_x;
  bool periodic_y;
};

void accumulate_cutoff_cell(const Domain& domain,
                            const GreensFunction& gf,
                            const Vec3& src,
                            const GreensFunctionParams& p,
                            Int ix, Int iy, Int iz,
                            std::vector<Real>& grid_conc) {
  const Vec3 tgt = domain.cell_center(ix, iy, iz);
  const Real c = gf.concentration_bounded(src, tgt, p);
  const auto idx = domain.cell_index(ix, iy, iz);
  grid_conc[idx] += c;
}

void accumulate_cutoff_row(const Domain& domain,
                           const GreensFunction& gf,
                           const Vec3& src,
                           const GreensFunctionParams& p,
                           Int src_ix, Int iy, Int iz,
                           const SuperposeGridContext& grid,
                           std::vector<Real>& grid_conc) {
  for (Int dx = -grid.span; dx <= grid.span; ++dx) {
    Int ix = src_ix + dx;
    if (!in_periodic_grid(ix, grid.nx, grid.periodic_x)) continue;
    accumulate_cutoff_cell(domain, gf, src, p, ix, iy, iz, grid_conc);
  }
}

void accumulate_source_cutoff(const Domain& domain,
                              const GreensFunction& gf,
                              const Vec3& src,
                              const GreensFunctionParams& p,
                              const SuperposeGridContext& grid,
                              std::vector<Real>& grid_conc) {
  Int src_ix = 0;
  Int src_iy = 0;
  Int src_iz = 0;
  domain.pos_to_grid(src, src_ix, src_iy, src_iz);

  for (Int dz = -grid.span; dz <= grid.span; ++dz) {
    Int iz = src_iz + dz;
    if (iz < 0 || iz >= grid.nz) continue;

    for (Int dy = -grid.span; dy <= grid.span; ++dy) {
      Int iy = src_iy + dy;
      if (!in_periodic_grid(iy, grid.ny, grid.periodic_y)) continue;
      accumulate_cutoff_row(domain, gf, src, p, src_ix, iy, iz, grid, grid_conc);
    }
  }
}

SuperposeGridContext make_superpose_grid(const Domain& domain, Real cutoff_radius) {
  return {
    domain.nx(),
    domain.ny(),
    domain.nz(),
    static_cast<Int>(std::ceil(cutoff_radius / domain.dx())),
    domain.config().periodic[0],
    domain.config().periodic[1],
  };
}

#ifdef GUTIBM_CUDA
bool try_gpu_superpose(const Domain& domain,
                       const AdvectionField& adv,
                       const std::vector<Vec3>& sources,
                       const std::vector<GreensFunctionParams>& params,
                       std::vector<Real>& grid_conc,
                       Real cutoff_radius) {
  if (!gpu_runtime_enabled()) return false;
  return gpu_superpose_to_grid(domain, adv, sources, params, grid_conc, cutoff_radius);
}
#endif

struct SuperposeSourcesContext {
  const Domain& domain;
  const GreensFunction& gf;
  Real cutoff_radius;
  SuperposeGridContext grid;
};

#ifndef GUTIBM_OPENMP
void superpose_sources_serial(const std::vector<Vec3>& sources,
                              const std::vector<GreensFunctionParams>& params,
                              const SuperposeSourcesContext& ctx,
                              std::vector<Real>& grid_conc) {
  for (size_t s = 0; s < sources.size(); ++s) {
    accumulate_source_cutoff(ctx.domain, ctx.gf, sources[s], params[s],
                             ctx.grid, grid_conc);
  }
}
#endif

#ifdef GUTIBM_OPENMP
struct SuperposeOpenmpContext {
  SuperposeSourcesContext sources;
  Int ncells;
};

void superpose_sources_openmp(const std::vector<Vec3>& sources,
                              const std::vector<GreensFunctionParams>& params,
                              const SuperposeOpenmpContext& ctx,
                              std::vector<Real>& grid_conc) {
  #pragma omp parallel
  {
    std::vector<Real> local_conc(ctx.ncells, 0.0);
    #pragma omp for schedule(dynamic)
    for (size_t s = 0; s < sources.size(); ++s) {
      accumulate_source_cutoff(ctx.sources.domain, ctx.sources.gf,
                               sources[s], params[s], ctx.sources.grid, local_conc);
    }
    #pragma omp critical
    {
      for (Int c = 0; c < ctx.ncells; ++c) {
        grid_conc[c] += local_conc[c];
      }
    }
  }
}
#endif

}  // namespace

void GreensFunction::init(const Domain& domain, const AdvectionField& adv) {
  domain_ = &domain;
  adv_    = &adv;
  z_lo_   = domain.lo()[2];
  z_hi_   = domain.hi()[2];
}

void GreensFunction::require_init() const {
  if (!domain_ || !adv_) {
    throw SimulationError(
        "GreensFunction::init() must be called before concentration queries");
  }
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
  require_init();
  Real D_eff = params.diff_coeff / params.retardation;
  Vec3 flow  = adv_->velocity(source);
  return single_kernel(source, target, D_eff, params.source_rate, flow);
}

Real GreensFunction::concentration_bounded(const Vec3& source, const Vec3& target,
                                            const GreensFunctionParams& params) const {
  require_init();
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
  require_init();

  const Int ncells = domain_->ncells();
  grid_conc.assign(ncells, 0.0);

#ifdef GUTIBM_CUDA
  if (adv_ && domain_ && try_gpu_superpose(*domain_, *adv_, sources, params,
                                           grid_conc, cutoff_radius)) {
    return;
  }
#endif

  const SuperposeGridContext grid = make_superpose_grid(*domain_, cutoff_radius);
  const SuperposeSourcesContext serial_ctx{*domain_, *this, cutoff_radius, grid};

#ifdef GUTIBM_OPENMP
  const SuperposeOpenmpContext omp_ctx{serial_ctx, ncells};
  superpose_sources_openmp(sources, params, omp_ctx, grid_conc);
#else
  superpose_sources_serial(sources, params, serial_ctx, grid_conc);
#endif
}

Real GreensFunction::peclet(const Vec3& pos, Real D_eff, Real length_scale) const {
  require_init();
  Vec3 v = adv_->velocity(pos);
  Real U = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  return (U * length_scale) / D_eff;
}

}  // namespace gutibm
