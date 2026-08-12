/* -----------------------------------------------------------------------
   GutIBM – QSSA solver implementation
   ----------------------------------------------------------------------- */

#include "qssa_solver.h"
#include "species_names.h"
#include "fmm.h"
#include "fmm_gpu.h"
#include "domain.h"
#include "advection.h"
#include "chemical_field.h"
#include "chemical_field_gpu.h"
#include "greens_function_gpu.h"
#include "dispatch.h"
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
    // FMM approximates mixed-toxin screening with a source-weighted k;
    // sqrt(k) is nonlinear, so heterogeneous clusters are approximate.
    avg_params.decay_rate  += s * params[i].decay_rate;
    total_s += s;
  }
  if (total_s > 0.0) {
    avg_params.diff_coeff  /= total_s;
    avg_params.pI          /= total_s;
    avg_params.retardation /= total_s;
    avg_params.decay_rate  /= total_s;
  } else {
    avg_params.diff_coeff  = cfg.fallback_diff_coeff;
    avg_params.pI          = cfg.fallback_pI;
    avg_params.retardation = cfg.fallback_retardation;
    avg_params.decay_rate = 0.0;
  }
  avg_params.source_rate = 0.0;
  return avg_params;
}

void accumulate_far_field(const FMM& fmm,
                          const GreensFunction& gf,
                          const GreensFunctionParams& avg_params,
                          const FarFieldGridContext& grid,
                          std::vector<Real>& toxin_conc,
                          Real toxin_cutoff,
                          bool near_field_on_device) {
  for (Int iz = 0; iz < grid.nz; ++iz) {
    for (Int iy = 0; iy < grid.ny; ++iy) {
      for (Int ix = 0; ix < grid.nx; ++ix) {
        const Vec3 tgt = grid.domain.cell_center(ix, iy, iz);
        const Int idx = grid.domain.cell_index(ix, iy, iz);
        Real contribution = 0.0;
        if (near_field_on_device) {
          contribution = fmm.evaluate_far_field(
              tgt, grid.fmm_theta, toxin_cutoff, gf, avg_params);
        } else {
          const Real total = fmm.evaluate_total_field(
              tgt, grid.fmm_theta, gf, avg_params);
          contribution = std::max(0.0, total - toxin_conc[idx]);
        }
        toxin_conc[idx] += contribution;
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

bool try_gpu_near_field(const Domain& domain,
                        const AdvectionField& adv,
                        const std::vector<Vec3>& sources,
                        const std::vector<GreensFunctionParams>& params,
                        const std::vector<Real>& strength_factors,
                        Real cutoff_radius,
                        ChemicalField& chem,
                        Int toxin_species_idx,
                        ChemicalFieldGpu* chem_gpu,
                        bool defer_host_sync) {
  if (chem_gpu == nullptr || !chem_gpu->active()) return false;
  double* d_conc = chem_gpu->conc_device(toxin_species_idx);
  if (d_conc == nullptr) return false;
  if (!gpu_superpose_to_device(domain, adv, sources, params, strength_factors,
                               d_conc, cutoff_radius)) {
    return false;
  }
  if (!defer_host_sync) {
    chem_gpu->sync_species_concentrations_to_host(chem, toxin_species_idx);
  }
  return true;
}

bool accumulate_near_field_gpu_or_cpu(const Domain& domain,
                                      const GreensFunction& gf,
                                      const AdvectionField& adv,
                                      const std::vector<Vec3>& sources,
                                      const std::vector<GreensFunctionParams>& params,
                                      const std::vector<Real>& strength_factors,
                                      Real cutoff_radius,
                                      std::vector<Real>& toxin_conc,
                                      ChemicalField& chem,
                                      Int toxin_species_idx,
                                      ChemicalFieldGpu* chem_gpu,
                                      bool defer_host_sync) {
  if (try_gpu_near_field(domain, adv, sources, params, strength_factors,
                         cutoff_radius, chem, toxin_species_idx, chem_gpu,
                         defer_host_sync)) {
    return true;
  }
  gf.superpose_to_grid(sources, params, strength_factors, toxin_conc,
                       cutoff_radius);
  return false;
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
      const Real protease_decay = (protease.enabled
                                   && bi.protease_half_life > 0.0)
          ? k_ln2 / bi.protease_half_life : 0.0;
      const Real dilution_decay = std::max(
          adv.washout_rate(a.x[2]), protease.dilution_rate);
      gfp.decay_rate = protease_decay + dilution_decay;

      out.sources.push_back(a.x);
      out.params.push_back(gfp);
      out.strength_factors.push_back(1.0);
      out.is_nuclease.push_back(bi.is_nuclease);
      out.targets.push_back(bi.target);
    }
  }
}

void append_burst_sources(const std::vector<ToxinBurstSource>& bursts,
                          Real current_time,
                          MicrocinSourceBuffers& out) {
  for (const ToxinBurstSource& burst : bursts) {
    if (burst.release_tau <= 0.0) continue;
    const Real age = std::max(0.0, current_time - burst.creation_time);
    const Real factor = std::exp(-age / burst.release_tau);
    if (factor < 1.0e-12) continue;

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
    ReceptorType target,
    ChemicalFieldGpu* chem_gpu) const {

  std::vector<Vec3> all_sources;
  std::vector<GreensFunctionParams> all_params;
  std::vector<Real> all_strengths;
  std::vector<bool> is_nuclease;
  std::vector<ReceptorType> all_targets;
  MicrocinSourceBuffers buffers{all_sources, all_params, all_strengths, is_nuclease, all_targets};

  collect_microcin_sources(agents, cfg_, protease, adv, buffers);
  append_burst_sources(bursts, current_time, buffers);

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
    solve_bacteriocin_field_fmm(sources, params, strength_factors, adv,
                                chem, toxin_species_idx, chem_gpu);
    return;
  }

  std::vector toxin_conc(domain_->ncells(), 0.0);
  const bool defer_sync = chem_gpu != nullptr && chem_gpu->active();
  if (!accumulate_near_field_gpu_or_cpu(*domain_, gf_, adv, sources, params,
                                      strength_factors, cfg_.toxin_cutoff,
                                      toxin_conc, chem, toxin_species_idx,
                                      chem_gpu, defer_sync)) {
    deposit_to_chemical_field(chem, toxin_species_idx, toxin_conc);
  }
}

void QSSASolver::solve_all_bacteriocin_fields(
    const AgentPool& agents,
    const std::vector<ToxinBurstSource>& bursts,
    Real current_time,
    const ProteaseConfig& protease,
    const AdvectionField& adv,
    ChemicalField& chem,
    ChemicalFieldGpu* chem_gpu) const {
  std::vector<Int> solved_indices;
  solved_indices.reserve(species::BACTERIOCIN_RECEPTOR_TARGETS.size());

  for (ReceptorType target : species::BACTERIOCIN_RECEPTOR_TARGETS) {
    const char* name = species::bacteriocin_species_for(target);
    if (name == nullptr) continue;
    Int idx = chem.find(name);
    if (idx < 0) continue;
    solve_bacteriocin_field(agents, bursts, current_time, protease, adv, chem,
                              idx, target, chem_gpu);
    if (chem_gpu != nullptr && chem_gpu->active()) {
      solved_indices.push_back(idx);
    }
  }

  if (chem_gpu != nullptr && chem_gpu->active()) {
    for (Int spec : solved_indices) {
      chem_gpu->sync_species_concentrations_to_host(chem, spec);
    }
  }
}

void QSSASolver::solve_bacteriocin_field_fmm(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    const AdvectionField& adv,
    ChemicalField& chem,
    Int toxin_species_idx,
    ChemicalFieldGpu* chem_gpu) const {

  std::vector<Real> strengths;
  const GreensFunctionParams avg_params =
      weighted_avg_params(params, strength_factors, cfg_, strengths);

  FMM fmm;
  fmm.build(sources, strengths, *domain_, cfg_.fmm_expansion_order);
  fmm.compute_local_expansions(cfg_.fmm_theta, gf_, avg_params);

  const Int ncells = domain_->ncells();
  const FarFieldGridContext far_grid = make_far_field_grid(*domain_, cfg_.fmm_theta);

  std::vector toxin_conc(ncells, 0.0);
  const bool defer_sync = chem_gpu != nullptr && chem_gpu->active();
  if (const bool near_on_gpu = accumulate_near_field_gpu_or_cpu(
          *domain_, gf_, adv, sources, params, strength_factors,
          cfg_.toxin_cutoff, toxin_conc, chem, toxin_species_idx, chem_gpu,
          defer_sync);
      near_on_gpu) {
    chem_gpu->sync_species_concentrations_to_host(chem, toxin_species_idx);
    for (Int c = 0; c < ncells; ++c) {
      toxin_conc[static_cast<size_t>(c)] = chem.conc(toxin_species_idx, c);
    }
  }
  if (fmm.locals_ready()
      && gpu_accumulate_far_field_local(
          fmm, *domain_, cfg_.fmm_expansion_order, toxin_conc, toxin_conc)) {
    // GPU far-field deposit complete.
  } else {
    accumulate_far_field(fmm, gf_, avg_params, far_grid, toxin_conc,
                         cfg_.toxin_cutoff, false);
  }
  deposit_to_chemical_field(chem, toxin_species_idx, toxin_conc);
  if (chem_gpu != nullptr && chem_gpu->active()) {
    chem_gpu->sync_species_concentrations_to_device(chem, toxin_species_idx);
  }
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
      // Pirt respiration: growth-associated + basal maintenance. The
      // maintenance term is applied per living cell regardless of growth, so
      // the O2 field tracks agent density (a non-growing cell still respires).
      const Real o2_use = oxygen.q_consumption * std::max(a.mu_realized, 0.0)
                        + oxygen.q_maintenance;
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
