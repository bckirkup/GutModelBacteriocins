/* -----------------------------------------------------------------------
   GutIBM – Receptor competitive binding implementation
   ----------------------------------------------------------------------- */

#include "fix_receptor.h"
#include "species_names.h"
#include "simulation.h"
#include "receptor_gpu.h"
#include <cmath>
#include <algorithm>
#include <vector>
#ifdef GUTIBM_OPENMP
#include <omp.h>
#endif

namespace gutibm {

FixReceptor::FixReceptor(Simulation& sim, const ReceptorConfig& cfg)
    : Fix("receptor", sim), cfg_(cfg) {}

void FixReceptor::compute(Real dt) {
  using enum PhenoState;
  auto& agents = sim_.agents();
  auto& rng    = sim_.rng();
  Int n = agents.size();

  std::vector kill_probs(n, 0.0);

  bool gpu_ok = false;
  if (sim_.gpu_active()) {
    auto& ag = sim_.agents_gpu();
    ag.sync_from_host(agents);
    std::vector<double> gpu_probs;
    if (gpu_compute_receptor_kill_probs_host_packed(
            ag, agents, sim_.chem_gpu(), sim_.chemical_field(), cfg_, dt,
            gpu_probs)) {
      for (Int i = 0; i < n; ++i) {
        kill_probs[i] = gpu_probs[static_cast<size_t>(i)];
      }
      gpu_ok = true;
    }
  }

  if (!gpu_ok) {
    #ifdef GUTIBM_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (Int i = 0; i < n; ++i) {
      const Agent& a = agents[i];
      if (a.state == DEAD || a.state == SOS_INDUCED)
        continue;
      kill_probs[i] = compute_kill_prob(a, dt);
    }
  }

  // Apply kills serially (RNG is not thread-safe)
  for (Int i = 0; i < n; ++i) {
    Agent& a = agents[i];
    if (a.state == DEAD || a.state == SOS_INDUCED)
      continue;

    if (kill_probs[i] > 0.0 && rng.bernoulli(kill_probs[i])) {
      const KillAssessment assessment = assess_kill(a, dt);
      a.state = DEAD;
      sim_.step_events().colicin_kills++;
      if (sim_.provenance_enabled()) {
        KillProvenanceEvent event;
        event.victim_id = a.identity.tag;
        event.position = a.x;
        event.strain = a.identity.type;
        event.cause = ProvenanceCause::COLICIN;
        event.toxin_concentration = assessment.concentration;
        event.toxin_occupancy = assessment.occupancy;
        event.toxin_hazard = assessment.hazard;
        sim_.record_kill_provenance(event);
      }
    }
  }
}

Real FixReceptor::local_toxin_conc(const ChemicalField& chem, Int cell,
                                    const char* species_name) const {
  Int idx = chem.find(species_name);
  if (idx < 0 || cell < 0) return 0.0;
  return chem.conc(idx, cell);
}

Real FixReceptor::compute_kill_prob(const Agent& agent, Real dt) const {
  return compute_kill_prob(agent, dt, nullptr);
}

FixReceptor::KillAssessment FixReceptor::assess_kill(const Agent& agent,
                                                      Real dt) const {
  KillAssessment assessment;
  assessment.probability = compute_kill_prob(agent, dt, &assessment);
  return assessment;
}

Real FixReceptor::compute_kill_prob(const Agent& agent, Real dt,
                                    KillAssessment* diagnostics) const {
  const auto& chem = sim_.chemical_field();
  Int cell = agent.grid_cell;
  if (cell < 0) return 0.0;

  Real tox_btuB = local_toxin_conc(chem, cell, species::BACTERIOCIN_BTUB);
  Real tox_fepA = local_toxin_conc(chem, cell, species::BACTERIOCIN_FEPA);
  Real tox_cirA = local_toxin_conc(chem, cell, species::BACTERIOCIN_CIRA);
  Real tox_fhuA = local_toxin_conc(chem, cell, species::BACTERIOCIN_FHUA);
  if (diagnostics != nullptr) {
    diagnostics->concentration = {tox_btuB, tox_fepA, tox_cirA, tox_fhuA};
  }
  if (tox_btuB <= 0.0 && tox_fepA <= 0.0 && tox_cirA <= 0.0 && tox_fhuA <= 0.0) {
    return 0.0;
  }

  Real total_kill = 0.0;

  // BtuB-mediated killing (colicin E family)
  if (tox_btuB > 0.0) {
    int ri = to_underlying(ReceptorType::BtuB);
    Real expr = agent.receptor_expr[ri];
    Int i_b12 = chem.find(species::B12);
    auto ligand = (i_b12 >= 0) ? chem.conc(i_b12, cell) : 0.0;

    Real occ = toxin_occupancy(tox_btuB, ligand,
                                cfg_.kd_colicinE_btuB,
                                cfg_.kd_b12_btuB,
                                expr,
                                agent.genome.toxin_affinity[ri],
                                agent.genome.ligand_affinity[ri]);
    Real eff = 1.0;
    for (const auto& bi : agent.genome.bi_loci) {
      if (bi.target == ReceptorType::BtuB) {
        Real candidate = cfg_.immunity_factor * bi.immunity_binding_affinity;
        if (candidate < eff) eff = candidate;
      }
    }
    Real hazard = cfg_.kill_rate_colicin * occ * eff * dt;
    if (diagnostics != nullptr) {
      diagnostics->occupancy[0] = occ;
      diagnostics->hazard[0] = hazard;
    }
    total_kill += hazard;
  }

  // FepA-mediated killing (colicin B/D)
  if (tox_fepA > 0.0) {
    int ri = to_underlying(ReceptorType::FepA);
    Real expr = agent.receptor_expr[ri];
    Int i_ferric_enterobactin = chem.find(species::FERRIC_ENTEROBACTIN);
    Real ligand = (i_ferric_enterobactin >= 0)
        ? chem.conc(i_ferric_enterobactin, cell) : 0.0;

    Real occ = toxin_occupancy(tox_fepA, ligand,
                                cfg_.kd_colicinB_fepA,
                                cfg_.kd_enterobactin,
                                expr,
                                agent.genome.toxin_affinity[ri],
                                agent.genome.ligand_affinity[ri]);
    Real eff = 1.0;
    for (const auto& bi : agent.genome.bi_loci) {
      if (bi.target == ReceptorType::FepA) {
        Real candidate = cfg_.immunity_factor * bi.immunity_binding_affinity;
        if (candidate < eff) eff = candidate;
      }
    }
    Real hazard = cfg_.kill_rate_colicin * occ * eff * dt;
    if (diagnostics != nullptr) {
      diagnostics->occupancy[1] = occ;
      diagnostics->hazard[1] = hazard;
    }
    total_kill += hazard;
  }

  // CirA-mediated killing (colicin Ia, microcin V)
  if (tox_cirA > 0.0) {
    int ri = to_underlying(ReceptorType::CirA);
    Real expr = agent.receptor_expr[ri];
    Int i_ferric_enterobactin = chem.find(species::FERRIC_ENTEROBACTIN);
    Real ligand = (i_ferric_enterobactin >= 0)
        ? chem.conc(i_ferric_enterobactin, cell)
          * cfg_.cirA_linearized_fraction
        : 0.0;

    Real occ = toxin_occupancy(tox_cirA, ligand,
                                cfg_.kd_colicinIa_cirA,
                                cfg_.kd_lin_enterobactin,
                                expr,
                                agent.genome.toxin_affinity[ri],
                                agent.genome.ligand_affinity[ri]);
    Real eff = 1.0;
    for (const auto& bi : agent.genome.bi_loci) {
      if (bi.target == ReceptorType::CirA) {
        Real candidate = cfg_.immunity_factor * bi.immunity_binding_affinity;
        if (candidate < eff) eff = candidate;
      }
    }
    Real hazard = cfg_.kill_rate_microcin * occ * eff * dt;
    if (diagnostics != nullptr) {
      diagnostics->occupancy[2] = occ;
      diagnostics->hazard[2] = hazard;
    }
    total_kill += hazard;
  }

  // FhuA-mediated killing (colicin M)
  if (tox_fhuA > 0.0) {
    int ri = to_underlying(ReceptorType::FhuA);
    Real expr = agent.receptor_expr[ri];
    Int i_ferrichrome = chem.find(species::FERRICHROME);
    Real ligand = (i_ferrichrome >= 0)
        ? chem.conc(i_ferrichrome, cell) : 0.0;

    Real occ = toxin_occupancy(tox_fhuA, ligand,
                                cfg_.kd_colicinM_fhuA,
                                cfg_.kd_ferrichrome,
                                expr,
                                agent.genome.toxin_affinity[ri],
                                agent.genome.ligand_affinity[ri]);
    Real eff = 1.0;
    for (const auto& bi : agent.genome.bi_loci) {
      if (bi.target == ReceptorType::FhuA) {
        Real candidate = cfg_.immunity_factor * bi.immunity_binding_affinity;
        if (candidate < eff) eff = candidate;
      }
    }
    Real hazard = cfg_.kill_rate_colicin * occ * eff * dt;
    if (diagnostics != nullptr) {
      diagnostics->occupancy[3] = occ;
      diagnostics->hazard[3] = hazard;
    }
    total_kill += hazard;
  }

  return std::min(1.0 - std::exp(-total_kill), 1.0);
}

Real FixReceptor::toxin_occupancy(Real tox_conc, Real ligand_conc,
                                   Real kd_tox, Real kd_ligand,
                                   Real receptor_expr,
                                   Real toxin_aff,
                                   Real ligand_aff) const {
  // Partial resistance scales the effective Kd values:
  //   - Low toxin_aff  → higher effective Kd_tox (toxin binds worse)
  //   - Low ligand_aff → higher effective Kd_ligand (ligand binds worse)
  Real eff_kd_tox    = kd_tox    / std::max(toxin_aff,  1.0e-6);
  Real eff_kd_ligand = kd_ligand / std::max(ligand_aff, 1.0e-6);

  // Competitive binding: Michaelis-Menten with competitive inhibition
  Real competitive_factor = 1.0 + ligand_conc / eff_kd_ligand;
  Real apparent_kd = eff_kd_tox * competitive_factor;

  // Occupancy = receptor_expr * [Tox] / (Kd_app + [Tox])
  return receptor_expr * tox_conc / (apparent_kd + tox_conc);
}

}  // namespace gutibm
