/* -----------------------------------------------------------------------
   GutIBM – QSSA solver implementation
   ----------------------------------------------------------------------- */

#include "qssa_solver.h"
#include "fmm.h"
#include "domain.h"
#include "advection.h"
#include "chemical_field.h"
#include "agent.h"
#include <cmath>
#include <numeric>
#ifdef GUTIBM_OPENMP
#include <omp.h>
#endif

namespace gutibm {

namespace {

constexpr Real k_ln2 = 0.6931471805599453;

bool in_periodic_grid(Int& idx, Int count, bool periodic) {
  if (periodic) {
    idx = ((idx % count) + count) % count;
    return true;
  }
  return idx >= 0 && idx < count;
}

void accumulate_near_field(const Domain& domain,
                           const GreensFunction& gf,
                           const std::vector<Vec3>& sources,
                           const std::vector<GreensFunctionParams>& params,
                           const std::vector<Real>& strength_factors,
                           Real cutoff_radius,
                           Int nx, Int ny, Int nz,
                           std::vector<Real>& toxin_conc) {
  const auto span = static_cast<Int>(std::ceil(cutoff_radius / domain.dx()));
  const bool periodic_x = domain.config().periodic[0];
  const bool periodic_y = domain.config().periodic[1];

  for (size_t s = 0; s < sources.size(); ++s) {
    const Vec3& src = sources[s];
    GreensFunctionParams p = params[s];
    p.source_rate *= strength_factors[s];

    Int src_ix = 0;
    Int src_iy = 0;
    Int src_iz = 0;
    domain.pos_to_grid(src, src_ix, src_iy, src_iz);

    for (Int dz = -span; dz <= span; ++dz) {
      Int iz = src_iz + dz;
      if (iz < 0 || iz >= nz) continue;

      for (Int dy = -span; dy <= span; ++dy) {
        Int iy = src_iy + dy;
        if (!in_periodic_grid(iy, ny, periodic_y)) continue;

        for (Int dx = -span; dx <= span; ++dx) {
          Int ix = src_ix + dx;
          if (!in_periodic_grid(ix, nx, periodic_x)) continue;

          const Vec3 tgt = domain.cell_center(ix, iy, iz);
          const Real c = gf.concentration_bounded(src, tgt, p);
          const Int idx = domain.cell_index(ix, iy, iz);
          toxin_conc[idx] += c;
        }
      }
    }
  }
}

GreensFunctionParams weighted_avg_params(
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
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
    avg_params.diff_coeff  = 4e-11;
    avg_params.pI          = 7.0;
    avg_params.retardation = 5.0;
  }
  avg_params.source_rate = 0.0;
  return avg_params;
}

void accumulate_far_field(FMM& fmm,
                          const Domain& domain,
                          Real fmm_theta,
                          const GreensFunction& gf,
                          const GreensFunctionParams& avg_params,
                          Int nx, Int ny, Int nz,
                          std::vector<Real>& toxin_conc) {
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        const Vec3 tgt = domain.cell_center(ix, iy, iz);
        const Int idx = domain.cell_index(ix, iy, iz);
        const Real total = fmm.evaluate_total_field(tgt, fmm_theta, gf, avg_params);
        const Real far = std::max(0.0, total - toxin_conc[idx]);
        toxin_conc[idx] += far;
      }
    }
  }
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
                              std::vector<Vec3>& sources,
                              std::vector<GreensFunctionParams>& params,
                              std::vector<Real>& strength_factors,
                              std::vector<bool>& is_nuclease) {
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

      sources.push_back(a.x);
      params.push_back(gfp);
      strength_factors.push_back(factor);
      is_nuclease.push_back(bi.is_nuclease);
    }
  }
}

void append_burst_sources(const std::vector<ToxinBurstSource>& bursts,
                          Real current_time,
                          const ProteaseConfig& protease,
                          std::vector<Vec3>& sources,
                          std::vector<GreensFunctionParams>& params,
                          std::vector<Real>& strength_factors,
                          std::vector<bool>& is_nuclease) {
  for (const ToxinBurstSource& burst : bursts) {
    Real factor = 1.0;
    if (protease.enabled && burst.decay_rate > 0.0) {
      const Real age = std::max(0.0, current_time - burst.creation_time);
      factor = std::exp(-burst.decay_rate * age);
      if (factor < 1.0e-12) continue;
    }

    sources.push_back(burst.pos);
    params.push_back(burst.params);
    strength_factors.push_back(factor);
    is_nuclease.push_back(burst.is_nuclease);
  }
}

void filter_nuclease_sources(const std::vector<Vec3>& sources,
                             const std::vector<GreensFunctionParams>& params,
                             const std::vector<Real>& strength_factors,
                             const std::vector<bool>& is_nuclease,
                             std::vector<Vec3>& out_sources,
                             std::vector<GreensFunctionParams>& out_params,
                             std::vector<Real>& out_strengths) {
  for (size_t i = 0; i < sources.size(); ++i) {
    if (!is_nuclease[i]) continue;
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
    Int toxin_species_idx) const {

  std::vector<Vec3> sources;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strength_factors;
  std::vector<bool> is_nuclease;

  collect_microcin_sources(agents, cfg_, protease, adv,
                           sources, params, strength_factors, is_nuclease);
  append_burst_sources(bursts, current_time, protease,
                       sources, params, strength_factors, is_nuclease);

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
  const Int nx = domain_->nx();
  const Int ny = domain_->ny();
  const Int nz = domain_->nz();
  accumulate_near_field(*domain_, gf_, sources, params, strength_factors,
                        cfg_.toxin_cutoff, nx, ny, nz, toxin_conc);
  deposit_to_chemical_field(chem, toxin_species_idx, toxin_conc);
}

void QSSASolver::solve_bacteriocin_field_fmm(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    ChemicalField& chem,
    Int toxin_species_idx) const {

  std::vector<Real> strengths;
  const GreensFunctionParams avg_params =
      weighted_avg_params(params, strength_factors, strengths);

  FMM fmm;
  fmm.build(sources, strengths, *domain_, cfg_.fmm_expansion_order);
  fmm.compute_local_expansions(cfg_.fmm_theta, gf_, avg_params);

  const Int nx = domain_->nx();
  const Int ny = domain_->ny();
  const Int nz = domain_->nz();
  const Int ncells = domain_->ncells();

  std::vector<Real> toxin_conc(ncells, 0.0);
  accumulate_near_field(*domain_, gf_, sources, params, strength_factors,
                        cfg_.toxin_cutoff, nx, ny, nz, toxin_conc);
  accumulate_far_field(fmm, *domain_, cfg_.fmm_theta, gf_, avg_params,
                       nx, ny, nz, toxin_conc);
  deposit_to_chemical_field(chem, toxin_species_idx, toxin_conc);
}

void QSSASolver::solve_nuclease_toxin_field(
    const AgentPool& agents,
    const std::vector<ToxinBurstSource>& bursts,
    Real current_time,
    const ProteaseConfig& protease,
    const AdvectionField& adv,
    ChemicalField& chem,
    Int nuclease_species_idx) const {

  std::vector<Vec3> sources;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strength_factors;
  std::vector<bool> is_nuclease;

  collect_microcin_sources(agents, cfg_, protease, adv,
                           sources, params, strength_factors, is_nuclease);
  append_burst_sources(bursts, current_time, protease,
                       sources, params, strength_factors, is_nuclease);

  std::vector<Vec3> nuc_sources;
  std::vector<GreensFunctionParams> nuc_params;
  std::vector<Real> nuc_strengths;
  filter_nuclease_sources(sources, params, strength_factors, is_nuclease,
                          nuc_sources, nuc_params, nuc_strengths);

  if (nuc_sources.empty()) {
    zero_species_field(chem, nuclease_species_idx);
    return;
  }

  if (cfg_.use_fmm) {
    solve_bacteriocin_field_fmm(nuc_sources, nuc_params, nuc_strengths,
                                chem, nuclease_species_idx);
    return;
  }

  std::vector<Real> toxin_conc(domain_->ncells(), 0.0);
  const Int nx = domain_->nx();
  const Int ny = domain_->ny();
  const Int nz = domain_->nz();
  accumulate_near_field(*domain_, gf_, nuc_sources, nuc_params, nuc_strengths,
                        cfg_.toxin_cutoff, nx, ny, nz, toxin_conc);
  deposit_to_chemical_field(chem, nuclease_species_idx, toxin_conc);
}

void QSSASolver::solve_nutrient_depletion(
    const AgentPool& agents,
    ChemicalField& chem,
    const OxygenConfig& oxygen) const {

  Int i_iron   = chem.find("iron");
  Int i_b12    = chem.find("b12");
  Int i_carbon = chem.find("carbon");
  Int i_oxygen = chem.find("oxygen");

  const Real cell_vol = domain_->dx() * domain_->dx() * domain_->dx();

  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (const Agent& a : agents) {
    if (a.state == PhenoState::DEAD) continue;

    Int cell = a.grid_cell;
    if (cell < 0) continue;

    Real consumption = a.mu_realized * a.biomass;

    if (i_iron >= 0) {
      #ifdef GUTIBM_OPENMP
      #pragma omp atomic
      #endif
      chem.reac(i_iron, cell) -= consumption * 1.0e-6;
    }
    if (i_b12 >= 0) {
      #ifdef GUTIBM_OPENMP
      #pragma omp atomic
      #endif
      chem.reac(i_b12, cell) -= consumption * 1.0e-9;
    }
    if (i_carbon >= 0) {
      #ifdef GUTIBM_OPENMP
      #pragma omp atomic
      #endif
      chem.reac(i_carbon, cell) -= consumption * 0.5;
    }
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
