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
                           Real cutoff_radius,
                           Int nx, Int ny, Int nz,
                           std::vector<Real>& toxin_conc) {
  const auto span = static_cast<Int>(std::ceil(cutoff_radius / domain.dx()));
  const bool periodic_x = domain.config().periodic[0];
  const bool periodic_y = domain.config().periodic[1];

  auto param_it = params.begin();
  for (const Vec3& src : sources) {
    const GreensFunctionParams& p = *param_it++;

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
    std::vector<Real>& strengths) {
  strengths.resize(params.size());
  GreensFunctionParams avg_params{};
  Real total_s = 0.0;
  size_t i = 0;
  for (const GreensFunctionParams& p : params) {
    const Real s = p.source_rate;
    strengths[i++] = s;
    avg_params.diff_coeff  += s * p.diff_coeff;
    avg_params.pI          += s * p.pI;
    avg_params.retardation += s * p.retardation;
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

void collect_toxin_sources(const AgentPool& agents,
                           const QSSAConfig& cfg,
                           std::vector<Vec3>& sources,
                           std::vector<GreensFunctionParams>& params) {
  for (const Agent& a : agents) {
    if (a.state == PhenoState::DEAD) continue;

    for (const auto& bi : a.genome.bi_loci) {
      GreensFunctionParams gfp;
      gfp.diff_coeff   = bi.diff_coeff;
      gfp.retardation  = bi.retardation;
      gfp.pI           = bi.pI;
      gfp.source_rate = (a.state == PhenoState::SOS_INDUCED)
                        ? cfg.colicin_release_rate
                        : cfg.microcin_secretion;

      sources.push_back(a.x);
      params.push_back(gfp);
    }
  }
}

}  // namespace

void QSSASolver::init(const QSSAConfig& cfg, const Domain& domain,
                       const AdvectionField& adv) {
  cfg_    = cfg;
  domain_ = &domain;
  gf_.init(domain, adv);
}

void QSSASolver::solve_bacteriocin_field(
    const AgentPool& agents,
    ChemicalField& chem,
    Int toxin_species_idx) const {

  std::vector<Vec3> sources;
  std::vector<GreensFunctionParams> params;
  collect_toxin_sources(agents, cfg_, sources, params);
  if (sources.empty()) return;

  if (cfg_.use_fmm) {
    solve_bacteriocin_field_fmm(sources, params, chem, toxin_species_idx);
    return;
  }

  std::vector<Real> toxin_conc;
  gf_.superpose_to_grid(sources, params, toxin_conc, cfg_.toxin_cutoff);
  deposit_to_chemical_field(chem, toxin_species_idx, toxin_conc);
}

void QSSASolver::solve_bacteriocin_field_fmm(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    ChemicalField& chem,
    Int toxin_species_idx) const {

  std::vector<Real> strengths;
  const GreensFunctionParams avg_params = weighted_avg_params(params, strengths);

  FMM fmm;
  fmm.build(sources, strengths, *domain_, cfg_.fmm_expansion_order);
  fmm.compute_local_expansions(cfg_.fmm_theta, gf_, avg_params);

  const Int nx = domain_->nx();
  const Int ny = domain_->ny();
  const Int nz = domain_->nz();
  const Int ncells = domain_->ncells();

  std::vector<Real> toxin_conc(ncells, 0.0);
  accumulate_near_field(*domain_, gf_, sources, params, cfg_.toxin_cutoff,
                      nx, ny, nz, toxin_conc);
  accumulate_far_field(fmm, *domain_, cfg_.fmm_theta, gf_, avg_params,
                       nx, ny, nz, toxin_conc);
  deposit_to_chemical_field(chem, toxin_species_idx, toxin_conc);
}

void QSSASolver::solve_nutrient_depletion(
    const AgentPool& agents,
    ChemicalField& chem) const {

  // Each living agent consumes nutrients from its grid cell
  // This is handled in fix_metabolism; here we just compute the
  // steady-state depletion halo around dense colonies
  Int i_iron   = chem.find("iron");
  Int i_b12    = chem.find("b12");
  Int i_carbon = chem.find("carbon");

  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (const Agent& a : agents) {
    if (a.state == PhenoState::DEAD) continue;

    Int cell = a.grid_cell;
    if (cell < 0) continue;

    // Localized consumption based on biomass and growth rate
    Real consumption = a.mu_realized * a.biomass;

    if (i_iron >= 0) {
      #ifdef GUTIBM_OPENMP
      #pragma omp atomic
      #endif
      chem.reac(i_iron, cell) -= consumption * 1.0e-6;  // iron stoichiometry
    }
    if (i_b12 >= 0) {
      #ifdef GUTIBM_OPENMP
      #pragma omp atomic
      #endif
      chem.reac(i_b12, cell) -= consumption * 1.0e-9;   // B12 stoichiometry
    }
    if (i_carbon >= 0) {
      #ifdef GUTIBM_OPENMP
      #pragma omp atomic
      #endif
      chem.reac(i_carbon, cell) -= consumption * 0.5;    // carbon stoichiometry
    }
  }
}

Real QSSASolver::point_concentration(
    const Vec3& target,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params) const {

  Real total = 0.0;
  auto param_it = params.begin();
  for (const Vec3& src : sources) {
    const GreensFunctionParams& p = *param_it++;
    Real d2 = domain_->min_image_dist_sq(src, target);
    if (Real cutoff2 = cfg_.toxin_cutoff * cfg_.toxin_cutoff; d2 > cutoff2) continue;

    total += gf_.concentration_bounded(src, target, p);
  }
  return total;
}

}  // namespace gutibm
