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
#include "neumann_image_series.h"
#include <cassert>
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
  Int span_x;
  Int span_y;
  Int span_z;
  Int x_begin;
  Int x_end;
  Int storage_nx;
  Int halo_width;
  bool periodic_x;
  bool periodic_y;
};

Int storage_cell_index(const SuperposeGridContext& grid,
                       Int ix, Int iy, Int iz) {
  const Int local_ix = ix - grid.x_begin + grid.halo_width;
  return iz * grid.storage_nx * grid.ny + iy * grid.storage_nx + local_ix;
}

void accumulate_cutoff_cell(const Domain& domain,
                            const GreensFunction& gf,
                            const Vec3& src,
                            const GreensFunctionParams& p,
                            Int ix, Int iy, Int iz,
                            const SuperposeGridContext& grid,
                            std::vector<Real>& grid_conc) {
  const Vec3 tgt = domain.cell_center(ix, iy, iz);
  const Real c = gf.concentration_bounded(src, tgt, p);
  const auto idx = storage_cell_index(grid, ix, iy, iz);
  grid_conc[idx] += c;
}

void accumulate_cutoff_row(const Domain& domain,
                           const GreensFunction& gf,
                           const Vec3& src,
                           const GreensFunctionParams& p,
                           Int src_ix, Int iy, Int iz,
                           const SuperposeGridContext& grid,
                           std::vector<Real>& grid_conc) {
  for (Int dx = -grid.span_x; dx <= grid.span_x; ++dx) {
    if (!is_first_periodic_offset(dx, grid.nx, grid.span_x, grid.periodic_x)) {
      continue;
    }
    Int ix = src_ix + dx;
    if (!in_periodic_grid(ix, grid.nx, grid.periodic_x)) continue;
    if (ix < grid.x_begin || ix >= grid.x_end) continue;
    accumulate_cutoff_cell(domain, gf, src, p, ix, iy, iz, grid, grid_conc);
  }
}

void accumulate_source_cutoff(const Domain& domain,
                              const GreensFunction& gf,
                              const Vec3& src,
                              GreensFunctionParams p,
                              Real strength,
                              const SuperposeGridContext& grid,
                              std::vector<Real>& grid_conc) {
  p.source_rate *= strength;
  Int src_ix = 0;
  Int src_iy = 0;
  Int src_iz = 0;
  domain.pos_to_grid(src, src_ix, src_iy, src_iz);

  for (Int dz = -grid.span_z; dz <= grid.span_z; ++dz) {
    Int iz = src_iz + dz;
    if (iz < 0 || iz >= grid.nz) continue;

    for (Int dy = -grid.span_y; dy <= grid.span_y; ++dy) {
      if (!is_first_periodic_offset(dy, grid.ny, grid.span_y,
                                    grid.periodic_y)) {
        continue;
      }
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
    static_cast<Int>(std::ceil(cutoff_radius / domain.dx_x())),
    static_cast<Int>(std::ceil(cutoff_radius / domain.dx_y())),
    static_cast<Int>(std::ceil(cutoff_radius / domain.dx_z())),
    0,
    domain.nx(),
    domain.nx(),
    0,
    domain.config().periodic[0],
    domain.config().periodic[1],
  };
}

SuperposeGridContext make_superpose_grid(
    const Domain& domain, Real cutoff_radius, Int x_begin, Int x_end,
    Int storage_nx, Int halo_width) {
  auto grid = make_superpose_grid(domain, cutoff_radius);
  grid.x_begin = x_begin;
  grid.x_end = x_end;
  grid.storage_nx = storage_nx;
  grid.halo_width = halo_width;
  return grid;
}

#ifdef GUTIBM_CUDA
bool try_gpu_superpose(const Domain& domain,
                       const AdvectionField& adv,
                       const std::vector<Vec3>& sources,
                       const std::vector<GreensFunctionParams>& params,
                       std::vector<Real>& grid_conc,
                       Real cutoff_radius, uint64_t* cap_hits) {
  if (!gpu_runtime_enabled()) return false;
  return gpu_superpose_to_grid(domain, adv, sources, params, {}, grid_conc,
                               cutoff_radius, cap_hits);
}
#endif

struct SuperposeSourcesContext {
  const Domain& domain;
  const GreensFunction& gf;
  Real cutoff_radius;
  SuperposeGridContext grid;
  const std::vector<Real>* strength_factors;
};

Real source_strength(const SuperposeSourcesContext& ctx, size_t source) {
  if (ctx.strength_factors == nullptr) return 1.0;
  return (*ctx.strength_factors)[source];
}

#ifndef GUTIBM_OPENMP
void superpose_sources_serial(const std::vector<Vec3>& sources,
                              const std::vector<GreensFunctionParams>& params,
                              const SuperposeSourcesContext& ctx,
                              std::vector<Real>& grid_conc) {
  for (size_t s = 0; s < sources.size(); ++s) {
    accumulate_source_cutoff(ctx.domain, ctx.gf, sources[s], params[s],
                             source_strength(ctx, s),
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
    std::vector local_conc(ctx.ncells, 0.0);
    #pragma omp for schedule(dynamic)
    for (size_t s = 0; s < sources.size(); ++s) {
      accumulate_source_cutoff(ctx.sources.domain, ctx.sources.gf,
                               sources[s], params[s],
                               source_strength(ctx.sources, s),
                               ctx.sources.grid, local_conc);
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

void superpose_cpu(const Domain& domain,
                   const GreensFunction& gf,
                   const std::vector<Vec3>& sources,
                   const std::vector<GreensFunctionParams>& params,
                   const std::vector<Real>* strength_factors,
                   std::vector<Real>& grid_conc,
                   Real cutoff_radius) {
  const SuperposeGridContext grid = make_superpose_grid(domain, cutoff_radius);
  const SuperposeSourcesContext ctx{
      domain, gf, cutoff_radius, grid, strength_factors};

#ifdef GUTIBM_OPENMP
  const Int ncells = domain.ncells();
  const SuperposeOpenmpContext omp_ctx{ctx, ncells};
  superpose_sources_openmp(sources, params, omp_ctx, grid_conc);
#else
  superpose_sources_serial(sources, params, ctx, grid_conc);
#endif
}

void superpose_cpu_local(
    const Domain& domain, const GreensFunction& gf,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>* strength_factors,
    std::vector<Real>& grid_conc, Real cutoff_radius,
    Int x_begin, Int x_end, Int storage_nx, Int halo_width) {
  const SuperposeGridContext grid = make_superpose_grid(
      domain, cutoff_radius, x_begin, x_end, storage_nx, halo_width);
  const SuperposeSourcesContext ctx{
      domain, gf, cutoff_radius, grid, strength_factors};

#ifdef GUTIBM_OPENMP
  const Int ncells = static_cast<Int>(grid_conc.size());
  const SuperposeOpenmpContext omp_ctx{ctx, ncells};
  superpose_sources_openmp(sources, params, omp_ctx, grid_conc);
#else
  superpose_sources_serial(sources, params, ctx, grid_conc);
#endif
}

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
                                    Real D_eff, Real Q, Real decay_rate,
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
  const Real screened_speed = decay_rate <= 0.0
      ? U_mag
      : std::sqrt(U_mag * U_mag + 4.0 * D_eff * decay_rate);
  Real exponent  = (U_dot_r - screened_speed * r) / (2.0 * D_eff);

  // Clamp exponent to avoid overflow
  exponent = std::max(exponent, -500.0);

  return prefactor * std::exp(exponent);
}

Real GreensFunction::concentration(const Vec3& source, const Vec3& target,
                                    const GreensFunctionParams& params) const {
  require_init();
  Real D_eff = params.diff_coeff / params.retardation;
  Vec3 flow  = adv_->velocity(source);
  return single_kernel(source, target, D_eff, params.source_rate,
                       params.decay_rate, flow);
}

Real GreensFunction::concentration_bounded(const Vec3& source, const Vec3& target,
                                            const GreensFunctionParams& params) const {
  require_init();
  Real D_eff = params.diff_coeff / params.retardation;
  Vec3 flow  = adv_->velocity(source);
  const Real Q = params.source_rate;
  const auto evaluate_image = [this, &source, &target, D_eff, Q,
                               decay_rate = params.decay_rate, flow](
                                  Real image_z, int reflected) {
    Vec3 image = source;
    image[2] = image_z;
    Vec3 image_flow = flow;
    if (reflected != 0) {
      image_flow[2] = -image_flow[2];
    }
    return single_kernel(image, target, D_eff, Q, decay_rate, image_flow);
  };
  int cap_hit = 0;
  // The image construction is exact only for flow uniform in z. The existing
  // radial_velocity(z) profile varies in z, so this retains its pre-existing
  // uniform-flow approximation.
  const Real total = neumann::sum_image_series(
      source[2], z_lo_, z_hi_, evaluate_image,
      D_eff, params.decay_rate,
      std::sqrt(flow[0] * flow[0] + flow[1] * flow[1] + flow[2] * flow[2]),
      std::abs(flow[2]),
      neumann::kRelativeTolerance, nullptr, &cap_hit);
  if (cap_hit != 0) {
#ifdef GUTIBM_OPENMP
#pragma omp atomic update
#endif
    ++image_series_cap_hits_;
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
                                           grid_conc, cutoff_radius,
                                           &image_series_cap_hits_)) {
    return;
  }
#endif

  superpose_cpu(*domain_, *this, sources, params, nullptr, grid_conc,
                cutoff_radius);
}

void GreensFunction::superpose_to_grid(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    std::vector<Real>& grid_conc,
    Real cutoff_radius) const {
  require_init();
  assert(strength_factors.size() == sources.size());
  grid_conc.assign(domain_->ncells(), 0.0);
  superpose_cpu(*domain_, *this, sources, params, &strength_factors,
                grid_conc, cutoff_radius);
}

void GreensFunction::superpose_to_local_grid(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    std::vector<Real>& grid_conc, Real cutoff_radius, Int x_begin, Int x_end,
    Int storage_nx, Int halo_width) const {
  require_init();
  assert(strength_factors.size() == sources.size());
  superpose_cpu_local(*domain_, *this, sources, params, &strength_factors,
                      grid_conc, cutoff_radius, x_begin, x_end, storage_nx,
                      halo_width);
}

Real GreensFunction::peclet(const Vec3& pos, Real D_eff, Real length_scale) const {
  require_init();
  Vec3 v = adv_->velocity(pos);
  Real U = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  return (U * length_scale) / D_eff;
}

}  // namespace gutibm
