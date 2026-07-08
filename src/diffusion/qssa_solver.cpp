/* -----------------------------------------------------------------------
   GutIBM – QSSA solver implementation
   ----------------------------------------------------------------------- */

#include "qssa_solver.h"
#include "species_names.h"
#include "fmm.h"
#include "domain.h"
#include "advection.h"
#include "chemical_field.h"
#include "agent.h"
#include <cmath>
#include <numbers>
#include <numeric>
#ifdef GUTIBM_OPENMP
#include <omp.h>
#endif

namespace gutibm {

namespace {

constexpr Real k_ln2 = std::numbers::ln2;

bool in_periodic_grid(Int& idx, Int count, bool periodic) {
  if (periodic) {
    idx = ((idx % count) + count) % count;
    return true;
  }
  return idx >= 0 && idx < count;
}

struct NearFieldGridContext {
  const Domain& domain;
  Int nx;
  Int ny;
  Int nz;
  Int span;
  bool periodic_x;
  bool periodic_y;
};

struct FarFieldGridContext {
  const Domain& domain;
  Real fmm_theta;
  Int nx;
  Int ny;
  Int nz;
};

struct MicrocinSourceBuffers {
  std::vector<Vec3>& sources;
  std::vector<GreensFunctionParams>& params;
  std::vector<Real>& strength_factors;
  std::vector<bool>& is_nuclease;
  std::vector<ReceptorType>& targets;
};

void accumulate_near_field_cell(const Domain& domain,
                                const GreensFunction& gf,
                                const Vec3& src,
                                const GreensFunctionParams& p,
                                Int ix, Int iy, Int iz,
                                std::vector<Real>& toxin_conc) {
  const Vec3 tgt = domain.cell_center(ix, iy, iz);
  const Real c = gf.concentration_bounded(src, tgt, p);
  const Int idx = domain.cell_index(ix, iy, iz);
  toxin_conc[idx] += c;
}

void accumulate_near_field_row(const GreensFunction& gf,
                               const Vec3& src,
                               const GreensFunctionParams& p,
                               Int src_ix, Int src_iy, Int iz,
                               const NearFieldGridContext& grid,
                               std::vector<Real>& toxin_conc) {
  for (Int dx = -grid.span; dx <= grid.span; ++dx) {
    Int ix = src_ix + dx;
    if (!in_periodic_grid(ix, grid.nx, grid.periodic_x)) continue;
    accumulate_near_field_cell(grid.domain, gf, src, p, ix, src_iy, iz, toxin_conc);
  }
}

void accumulate_near_field(const Domain& domain,
                           const GreensFunction& gf,
                           const std::vector<Vec3>& sources,
                           const std::vector<GreensFunctionParams>& params,
                           const std::vector<Real>& strength_factors,
                           const NearFieldGridContext& grid,
                           std::vector<Real>& toxin_conc) {
  for (size_t s = 0; s < sources.size(); ++s) {
    const Vec3& src = sources[s];
    GreensFunctionParams p = params[s];
    p.source_rate *= strength_factors[s];

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
        accumulate_near_field_row(gf, src, p, src_ix, iy, iz, grid, toxin_conc);
      }
    }
  }
}

NearFieldGridContext make_near_field_grid(const Domain& domain, Real cutoff_radius) {
  return {
    domain,
    domain.nx(),
    domain.ny(),
    domain.nz(),
    static_cast<Int>(std::ceil(cutoff_radius / domain.dx())),
    domain.config().periodic[0],
    domain.config().periodic[1],
  };
}

GreensFunctionParams weighted_avg_params(
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    const QSSAConfig& cfg,
    std::vector<Real>& strengths) {
  strengths.resize(params.size());
  GreensFunctionParams avg_params{};
  Real total_s = 0.0;
  for (size_t i = 0; i < params.size(); ++i) {
    const Real s = params[i].source_rate * strength_factors[i];
    strengths[i] = s;
    avg_params.diff_coeff  += s * params[i].diff_coeff;
    avg_params.pI          += s * params[i].pI;
    avg_params.retardation += s * params[i].retardation;
    total_s += s;
  }
  if (total_s > 0.0) {
    avg_params.diff_coeff  /= total_s;
    avg_params.pI          /= total_s;
    avg_params.retardation /= total_s;
  } else {
    avg_params.diff_coeff  = cfg.fallback_diff_coeff;
    avg_params.pI          = cfg.fallback_pI;
    avg_params.retardation = cfg.fallback_retardation;
  }
  avg_params.source_rate = 0.0;
  return avg_params;
}

void accumulate_far_field(const FMM& fmm,
                          const GreensFunction& gf,
                          const GreensFunctionParams& avg_params,
                          const FarFieldGridContext& grid,
                          std::vector<Real>& toxin_conc) {
  for (Int iz = 0; iz < grid.nz; ++iz) {
    for (Int iy = 0; iy < grid.ny; ++iy) {
      for (Int ix = 0; ix < grid.nx; ++ix) {
        const Vec3 tgt = grid.domain.cell_center(ix, iy, iz);
        const Int idx = grid.domain.cell_index(ix, iy, iz);
        const Real total = fmm.evaluate_total_field(tgt, grid.fmm_theta, gf, avg_params);
        const Real far = std::max(0.0, total - toxin_conc[idx]);
        toxin_conc[idx] += far;
      }
    }
  }
}

FarFieldGridContext make_far_field_grid(const Domain& domain, Real fmm_theta) {
  return {domain, fmm_theta, domain.nx(), domain.ny(), domain.nz()};
}

void deposit_to_chemical_field(ChemicalField& chem,
                               Int toxin_species_idx,
                               const std::vector<Real>& concentrations) {
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (size_t c = 0; c < concentrations.size(); ++c) {
    chem.conc(toxin_species_idx, static_cast<Int>(c)) = concentrations[c];
  }
}

Real microcin_steady_decay_factor(Real decay_rate,
                                  Real washout,
                                  Real fallback_dilution) {
  const Real k_dilution = std::max(washout, fallback_dilution);
  if (decay_rate <= 0.0) return 1.0;
  return 1.0 / (1.0 + decay_rate / k_dilution);
}

void collect_microcin_sources(const AgentPool& agents,
                              const QSSAConfig& cfg,
                              const ProteaseConfig& protease,
                              const AdvectionField& adv,
                              MicrocinSourceBuffers& out) {
  for (const Agent& a : agents) {
    if (a.state == PhenoState::DEAD || a.state == PhenoState::SOS_INDUCED) continue;

    for (const auto& bi : a.genome.bi_loci) {
      if (bi.release_mode != ReleaseMode::CONTINUOUS) continue;

      GreensFunctionParams gfp;
      gfp.diff_coeff   = bi.diff_coeff;
      gfp.retardation  = bi.retardation;
      gfp.pI           = bi.pI;
      gfp.source_rate  = cfg.microcin_secretion;

      Real factor = 1.0;
      if (protease.enabled && bi.protease_half_life > 0.0) {
        const Real decay_rate = k_ln2 / bi.protease_half_life;
        factor = microcin_steady_decay_factor(
            decay_rate, adv.washout_rate(a.x[2]), protease.dilution_rate);
      }

      out.sources.push_back(a.x);
      out.params.push_back(gfp);
      out.strength_factors.push_back(factor);
      out.is_nuclease.push_back(bi.is_nuclease);
      out.targets.push_back(bi.target);
    }
  }
}

void append_burst_sources(const std::vector<ToxinBurstSource>& bursts,
                          Real current_time,
                          const ProteaseConfig& protease,
                          MicrocinSourceBuffers& out) {
  for (const ToxinBurstSource& burst : bursts) {
    Real factor = 1.0;
    if (protease.enabled && burst.decay_rate > 0.0) {
      const Real age = std::max(0.0, current_time - burst.creation_time);
      factor = std::exp(-burst.decay_rate * age);
      if (factor < 1.0e-12) continue;
    }

    out.sources.push_back(burst.pos);
    out.params.push_back(burst.params);
    out.strength_factors.push_back(factor);
    out.is_nuclease.push_back(burst.is_nuclease);
    out.targets.push_back(burst.target);
  }
}

void filter_sources_by_target(const std::vector<Vec3>& sources,
                              const std::vector<GreensFunctionParams>& params,
                              const std::vector<Real>& strength_factors,
                              const std::vector<ReceptorType>& targets,
                              ReceptorType target,
                              std::vector<Vec3>& out_sources,
                              std::vector<GreensFunctionParams>& out_params,
                              std::vector<Real>& out_strengths) {
  for (size_t i = 0; i < sources.size(); ++i) {
    if (targets[i] != target) continue;
    out_sources.push_back(sources[i]);
    out_params.push_back(params[i]);
    out_strengths.push_back(strength_factors[i]);
  }
}

void zero_species_field(ChemicalField& chem, Int species_idx) {
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Int c = 0; c < chem.ncells(); ++c) {
    chem.conc(species_idx, c) = 0.0;
  }
}

}  // namespace

void QSSASolver::init(const QSSAConfig& cfg, const Domain& domain,
                       const AdvectionField& adv) {
  cfg_    = cfg;
  domain_ = &domain;
  adv_    = &adv;
  gf_.init(domain, adv);
}

void QSSASolver::solve_bacteriocin_field(
    const AgentPool& agents,
    const std::vector<ToxinBurstSource>& bursts,
    Real current_time,
    const ProteaseConfig& protease,
    const AdvectionField& adv,
    ChemicalField& chem,
    Int toxin_species_idx,
    ReceptorType target) const {

  std::vector<Vec3> all_sources;
  std::vector<GreensFunctionParams> all_params;
  std::vector<Real> all_strengths;
  std::vector<bool> is_nuclease;
  std::vector<ReceptorType> all_targets;
  MicrocinSourceBuffers buffers{all_sources, all_params, all_strengths, is_nuclease, all_targets};

  collect_microcin_sources(agents, cfg_, protease, adv, buffers);
  append_burst_sources(bursts, current_time, protease, buffers);

  std::vector<Vec3> sources;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strength_factors;
  filter_sources_by_target(all_sources, all_params, all_strengths, all_targets, target,
                           sources, params, strength_factors);

  if (sources.empty()) {
    zero_species_field(chem, toxin_species_idx);
    return;
  }

  if (cfg_.use_fmm) {
    solve_bacteriocin_field_fmm(sources, params, strength_factors,
                                chem, toxin_species_idx);
    return;
  }

  std::vector<Real> toxin_conc(domain_->ncells(), 0.0);
  const NearFieldGridContext grid = make_near_field_grid(*domain_, cfg_.toxin_cutoff);
  accumulate_near_field(*domain_, gf_, sources, params, strength_factors,
                        grid, toxin_conc);
  deposit_to_chemical_field(chem, toxin_species_idx, toxin_conc);
}

void QSSASolver::solve_all_bacteriocin_fields(
    const AgentPool& agents,
    const std::vector<ToxinBurstSource>& bursts,
    Real current_time,
    const ProteaseConfig& protease,
    const AdvectionField& adv,
    ChemicalField& chem) const {
  for (ReceptorType target : species::BACTERIOCIN_RECEPTOR_TARGETS) {
    const char* name = species::bacteriocin_species_for(target);
    if (name == nullptr) continue;
    Int idx = chem.find(name);
    if (idx < 0) continue;
    solve_bacteriocin_field(agents, bursts, current_time, protease, adv, chem, idx, target);
  }
}

void QSSASolver::solve_bacteriocin_field_fmm(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    ChemicalField& chem,
    Int toxin_species_idx) const {

  std::vector<Real> strengths;
  const GreensFunctionParams avg_params =
      weighted_avg_params(params, strength_factors, cfg_, strengths);

  FMM fmm;
  fmm.build(sources, strengths, *domain_, cfg_.fmm_expansion_order);
  fmm.compute_local_expansions(cfg_.fmm_theta, gf_, avg_params);

  const Int ncells = domain_->ncells();
  const NearFieldGridContext near_grid = make_near_field_grid(*domain_, cfg_.toxin_cutoff);
  const FarFieldGridContext far_grid = make_far_field_grid(*domain_, cfg_.fmm_theta);

  std::vector<Real> toxin_conc(ncells, 0.0);
  accumulate_near_field(*domain_, gf_, sources, params, strength_factors,
                        near_grid, toxin_conc);
  accumulate_far_field(fmm, gf_, avg_params, far_grid, toxin_conc);
  deposit_to_chemical_field(chem, toxin_species_idx, toxin_conc);
}

void QSSASolver::solve_nutrient_depletion(
    const AgentPool& agents,
    ChemicalField& chem,
    const OxygenConfig& oxygen) const {

  Int i_oxygen = chem.find(species::OXYGEN);

  const Real cell_vol = domain_->dx() * domain_->dx() * domain_->dx();

  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (const Agent& a : agents) {
    if (a.state == PhenoState::DEAD) continue;

    Int cell = a.grid_cell;
    if (cell < 0) continue;

    // Spec 6 — per-agent carbon/iron/B12 depletion is handled solely by the
    // metabolism Fix (yield-based uptake in FixMetabolism::grow_agent). The
    // carbon/iron/B12 terms formerly applied here duplicated that uptake
    // (double-counting); they have been removed. Corrinoid (B12) is no longer
    // depleted at all. This function now applies only agent O2 respiration,
    // which has no counterpart in the metabolism Fix.
    if (oxygen.enabled && i_oxygen >= 0 && cell_vol > 0.0) {
      const Real o2_use = oxygen.q_consumption * std::max(a.mu_realized, 0.0);
      #ifdef GUTIBM_OPENMP
      #pragma omp atomic
      #endif
      chem.reac(i_oxygen, cell) -= o2_use / cell_vol;
    }
  }
}

Real QSSASolver::point_concentration(
    const Vec3& target,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors) const {

  Real total = 0.0;
  for (size_t s = 0; s < sources.size(); ++s) {
    GreensFunctionParams p = params[s];
    p.source_rate *= strength_factors[s];
    Real d2 = domain_->min_image_dist_sq(sources[s], target);
    if (Real cutoff2 = cfg_.toxin_cutoff * cfg_.toxin_cutoff; d2 > cutoff2) continue;
    total += gf_.concentration_bounded(sources[s], target, p);
  }
  return total;
}

}  // namespace gutibm
