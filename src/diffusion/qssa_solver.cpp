/* -----------------------------------------------------------------------
   GutIBM – QSSA solver implementation
   ----------------------------------------------------------------------- */

#include "qssa_solver.h"
#include "domain.h"
#include "advection.h"
#include "chemical_field.h"
#include "agent.h"
#include <cmath>

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

  // Superpose onto grid
  std::vector<Real> toxin_conc;
  gf_.superpose_to_grid(sources, params, toxin_conc, cfg_.toxin_cutoff);

  // Deposit onto chemical field
  for (Int c = 0; c < chem.ncells(); ++c) {
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

  for (Int i = 0; i < agents.size(); ++i) {
    const Agent& a = agents[i];
    if (a.state == PhenoState::DEAD) continue;

    Int cell = a.grid_cell;
    if (cell < 0) continue;

    // Localized consumption based on biomass and growth rate
    Real consumption = a.mu_realized * a.biomass;

    if (i_iron >= 0) {
      chem.reac(i_iron, cell) -= consumption * 1.0e-6;  // iron stoichiometry
    }
    if (i_b12 >= 0) {
      chem.reac(i_b12, cell) -= consumption * 1.0e-9;   // B12 stoichiometry
    }
    if (i_carbon >= 0) {
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
    Real cutoff2 = cfg_.toxin_cutoff * cfg_.toxin_cutoff;
    if (d2 > cutoff2) continue;

    total += gf_.concentration_bounded(sources[s], target, params[s]);
  }
  return total;
}

}  // namespace gutibm
