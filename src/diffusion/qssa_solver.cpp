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

  // Collect all active toxin sources (lysed cells and microcin producers)
  std::vector<Vec3> sources;
  std::vector<GreensFunctionParams> params;

  for (Int i = 0; i < agents.size(); ++i) {
    const Agent& a = agents[i];
    if (a.state == PhenoState::DEAD) continue;

    for (const auto& bi : a.genome.bi_loci) {
      GreensFunctionParams gfp;
      gfp.diff_coeff   = bi.diff_coeff;
      gfp.retardation  = bi.retardation;
      gfp.pI           = bi.pI;

      // SOS-lysed cells release a burst; living producers secrete continuously
      if (a.state == PhenoState::SOS_INDUCED) {
        gfp.source_rate = cfg_.colicin_release_rate;
      } else {
        gfp.source_rate = cfg_.microcin_secretion;
      }

      sources.push_back(a.x);
      params.push_back(gfp);
    }
  }

  if (sources.empty()) return;

  if (cfg_.use_fmm) {
    solve_bacteriocin_field_fmm(sources, params, chem, toxin_species_idx);
    return;
  }

  // Superpose onto grid
  std::vector<Real> toxin_conc;
  gf_.superpose_to_grid(sources, params, toxin_conc, cfg_.toxin_cutoff);

  // Deposit onto chemical field
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Int c = 0; c < chem.ncells(); ++c) {
    chem.conc(toxin_species_idx, c) = toxin_conc[c];
  }
}

void QSSASolver::solve_bacteriocin_field_fmm(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    ChemicalField& chem,
    Int toxin_species_idx) const {

  // Build source strengths for the FMM tree (use source_rate as strength)
  std::vector<Real> strengths(sources.size());
  for (size_t i = 0; i < sources.size(); ++i)
    strengths[i] = params[i].source_rate;

  // Source-rate-weighted average kernel parameters for far-field expansion
  GreensFunctionParams avg_params{};
  Real total_s = 0.0;
  for (size_t i = 0; i < params.size(); ++i) {
    Real s = strengths[i];
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

  FMM fmm;
  fmm.build(sources, strengths, *domain_, cfg_.fmm_expansion_order);
  fmm.compute_local_expansions(cfg_.fmm_theta, gf_, avg_params);

  const Int nx = domain_->nx();
  const Int ny = domain_->ny();
  const Int nz = domain_->nz();
  const Int ncells = domain_->ncells();

  // Near-field: exact evaluation within cutoff via spatial hash approach
  std::vector<Real> toxin_conc(ncells, 0.0);

  for (size_t s = 0; s < sources.size(); ++s) {
    const Vec3& src = sources[s];
    const GreensFunctionParams& p = params[s];

    Int src_ix = 0;
    Int src_iy = 0;
    Int src_iz = 0;
    domain_->pos_to_grid(src, src_ix, src_iy, src_iz);

    auto span = static_cast<Int>(std::ceil(cfg_.toxin_cutoff / domain_->dx()));

    for (Int dz = -span; dz <= span; ++dz) {
      Int iz = src_iz + dz;
      if (iz < 0 || iz >= nz) continue;
      for (Int dy = -span; dy <= span; ++dy) {
        Int iy = src_iy + dy;
        if (domain_->config().periodic[1]) {
          iy = ((iy % ny) + ny) % ny;
        } else if (iy < 0 || iy >= ny) continue;

        for (Int dx = -span; dx <= span; ++dx) {
          Int ix = src_ix + dx;
          if (domain_->config().periodic[0]) {
            ix = ((ix % nx) + nx) % nx;
          } else if (ix < 0 || ix >= nx) continue;

          Vec3 tgt = domain_->cell_center(ix, iy, iz);
          Real c = gf_.concentration_bounded(src, tgt, p);
          Int idx = domain_->cell_index(ix, iy, iz);
          toxin_conc[idx] += c;
        }
      }
    }
  }

  // Far-field: FMM multipole for each grid cell.  When local expansions are
  // precomputed, use total - near to avoid double-counting near sources.
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        Vec3 tgt = domain_->cell_center(ix, iy, iz);
        Int idx = domain_->cell_index(ix, iy, iz);
        Real total = fmm.evaluate_total_field(tgt, cfg_.fmm_theta, gf_, avg_params);
        Real far = std::max(0.0, total - toxin_conc[idx]);
        toxin_conc[idx] += far;
      }
    }
  }

  for (Int c = 0; c < ncells; ++c) {
    chem.conc(toxin_species_idx, c) = toxin_conc[c];
  }
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
  for (Int i = 0; i < agents.size(); ++i) {
    const Agent& a = agents[i];
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
  for (size_t s = 0; s < sources.size(); ++s) {
    Real d2 = domain_->min_image_dist_sq(sources[s], target);
    if (Real cutoff2 = cfg_.toxin_cutoff * cfg_.toxin_cutoff; d2 > cutoff2) continue;

    total += gf_.concentration_bounded(sources[s], target, params[s]);
  }
  return total;
}

}  // namespace gutibm
