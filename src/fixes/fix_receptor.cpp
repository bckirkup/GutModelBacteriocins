/* -----------------------------------------------------------------------
   GutIBM – Receptor competitive binding implementation
   ----------------------------------------------------------------------- */

#include "fix_receptor.h"
#include "simulation.h"
#include <cmath>
#include <algorithm>

namespace gutibm {

FixReceptor::FixReceptor(Simulation& sim, const ReceptorConfig& cfg)
    : Fix("receptor", sim), cfg_(cfg) {}

void FixReceptor::compute(Real dt) {
  auto& agents = sim_.agents();
  auto& rng    = sim_.rng();

  for (Int i = 0; i < agents.size(); ++i) {
    Agent& a = agents[i];
    if (a.state == PhenoState::DEAD || a.state == PhenoState::SOS_INDUCED)
      continue;

    Real p_kill = compute_kill_prob(a, dt);

    if (p_kill > 0.0 && rng.bernoulli(p_kill)) {
      // Check immunity
      bool immune = false;
      // Agent is immune to its own toxins (cognate immunity)
      // For external toxins, check if agent carries matching immunity
      // (simplified: if agent has BI cluster targeting same receptor,
      //  the immunity protein protects)
      // Full immunity check is handled here
      if (!immune) {
        a.state = PhenoState::DEAD;
      }
    }
  }
}

Real FixReceptor::compute_kill_prob(const Agent& agent, Real dt) const {
  auto& chem = sim_.chemical_field();
  Int cell = agent.grid_cell;
  if (cell < 0) return 0.0;

  // Get local toxin concentration
  Int i_toxin = chem.find("bacteriocin");
  if (i_toxin < 0) return 0.0;

  Real tox_conc = chem.conc(i_toxin, cell);
  if (tox_conc <= 0.0) return 0.0;

  Real total_kill = 0.0;

  // BtuB-mediated killing (colicin E family)
  {
    Real expr = agent.receptor_expr[static_cast<int>(ReceptorType::BtuB)];
    Int i_b12 = chem.find("b12");
    Real ligand = (i_b12 >= 0) ? chem.conc(i_b12, cell) : 0.0;

    Real occ = toxin_occupancy(tox_conc, ligand,
                                cfg_.kd_colicinE_btuB,
                                cfg_.kd_b12_btuB,
                                expr);
    // Check immunity — best-matching BI cluster determines cross-protection.
    // immunity_binding_affinity < 1 indicates reduced cross-neutralization
    // (e.g. against immunity-escape super-killer toxins, VADI §57).
    Real eff = 1.0;
    for (const auto& bi : agent.genome.bi_loci) {
      if (bi.target == ReceptorType::BtuB) {
        Real candidate = cfg_.immunity_factor * bi.immunity_binding_affinity;
        if (candidate < eff) eff = candidate;
      }
    }
    total_kill += cfg_.kill_rate_colicin * occ * eff * dt;
  }

  // FepA-mediated killing (colicin B/D)
  {
    Real expr = agent.receptor_expr[static_cast<int>(ReceptorType::FepA)];
    Int i_iron = chem.find("iron");
    Real ligand = (i_iron >= 0) ? chem.conc(i_iron, cell) : 0.0;

    Real occ = toxin_occupancy(tox_conc, ligand,
                                cfg_.kd_colicinB_fepA,
                                cfg_.kd_enterobactin,
                                expr);
    Real eff = 1.0;
    for (const auto& bi : agent.genome.bi_loci) {
      if (bi.target == ReceptorType::FepA) {
        Real candidate = cfg_.immunity_factor * bi.immunity_binding_affinity;
        if (candidate < eff) eff = candidate;
      }
    }
    total_kill += cfg_.kill_rate_colicin * occ * eff * dt;
  }

  // CirA-mediated killing (colicin Ia, microcin V)
  {
    Real expr = agent.receptor_expr[static_cast<int>(ReceptorType::CirA)];
    // CirA transports linearized enterobactin — use iron field as proxy
    Int i_iron = chem.find("iron");
    Real ligand = (i_iron >= 0) ? chem.conc(i_iron, cell) * 0.1 : 0.0;

    Real occ = toxin_occupancy(tox_conc, ligand,
                                cfg_.kd_colicinIa_cirA,
                                cfg_.kd_lin_enterobactin,
                                expr);
    Real eff = 1.0;
    for (const auto& bi : agent.genome.bi_loci) {
      if (bi.target == ReceptorType::CirA) {
        Real candidate = cfg_.immunity_factor * bi.immunity_binding_affinity;
        if (candidate < eff) eff = candidate;
      }
    }
    total_kill += cfg_.kill_rate_microcin * occ * eff * dt;
  }

  return std::min(1.0 - std::exp(-total_kill), 1.0);
}

Real FixReceptor::toxin_occupancy(Real tox_conc, Real ligand_conc,
                                   Real kd_tox, Real kd_ligand,
                                   Real receptor_expr) const {
  // Competitive binding: Michaelis-Menten with competitive inhibition
  // Apparent Kd = Kd_tox * (1 + [ligand]/Kd_ligand)
  Real competitive_factor = 1.0 + ligand_conc / kd_ligand;
  Real apparent_kd = kd_tox * competitive_factor;

  // Occupancy = receptor_expr * [Tox] / (Kd_app + [Tox])
  return receptor_expr * tox_conc / (apparent_kd + tox_conc);
}

}  // namespace gutibm
