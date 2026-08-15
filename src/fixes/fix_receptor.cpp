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
      kill_probs[i] = compute_kill_prob(a, i, dt);
    }
  }

  // Apply kills serially (RNG is not thread-safe)
  for (Int i = 0; i < n; ++i) {
    Agent& a = agents[i];
    if (a.state == DEAD || a.state == SOS_INDUCED)
      continue;

    if (kill_probs[i] > 0.0 && rng.bernoulli(kill_probs[i])) {
      const KillAssessment assessment = assess_kill(a, i, dt);
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

Real FixReceptor::local_toxin_conc(const Agent& agent, Int agent_index,
                                    const char* species_name) const {
  const auto& chem = sim_.chemical_field();
  Int idx = chem.find(species_name);
  if (idx < 0 || agent.grid_cell < 0) return 0.0;
  if (sim_.qssa().agent_sampling()) {
    return sim_.qssa().sampled_toxin_conc(agent_index, idx);
  }
  return chem.conc_global(idx, agent.grid_cell);
}

Real FixReceptor::compute_kill_prob(const Agent& agent, Int agent_index,
                                    Real dt) const {
  return compute_kill_prob(agent, agent_index, dt, nullptr);
}

FixReceptor::KillAssessment FixReceptor::assess_kill(
    const Agent& agent, Int agent_index, Real dt) const {
  KillAssessment assessment;
  assessment.probability = compute_kill_prob(
      agent, agent_index, dt, &assessment);
  return assessment;
}

Real FixReceptor::compute_kill_prob(const Agent& agent, Int agent_index, Real dt,
                                    KillAssessment* diagnostics) const {
  const Int cell = agent.grid_cell;
  if (cell < 0) return 0.0;

  const std::array<ReceptorDescriptor, 4> descriptors = {{
      {ReceptorType::BtuB, species::BACTERIOCIN_BTUB, species::B12,
       &ReceptorConfig::kd_colicinE_btuB, &ReceptorConfig::kd_b12_btuB,
       &ReceptorConfig::kill_rate_colicin, 1.0, 0},
      {ReceptorType::FepA, species::BACTERIOCIN_FEPA,
       species::FERRIC_ENTEROBACTIN, &ReceptorConfig::kd_colicinB_fepA,
       &ReceptorConfig::kd_enterobactin, &ReceptorConfig::kill_rate_colicin,
       1.0, 1},
      {ReceptorType::CirA, species::BACTERIOCIN_CIRA,
       species::FERRIC_ENTEROBACTIN, &ReceptorConfig::kd_colicinIa_cirA,
       &ReceptorConfig::kd_lin_enterobactin, &ReceptorConfig::kill_rate_microcin,
       cfg_.cirA_linearized_fraction, 2},
      {ReceptorType::FhuA, species::BACTERIOCIN_FHUA, species::FERRICHROME,
       &ReceptorConfig::kd_colicinM_fhuA, &ReceptorConfig::kd_ferrichrome,
       &ReceptorConfig::kill_rate_colicin, 1.0, 3},
  }};

  Real total_kill = 0.0;
  bool any_toxin = false;
  for (const auto& descriptor : descriptors) {
    const Real toxin_concentration =
        local_toxin_conc(agent, agent_index, descriptor.toxin_species);
    if (diagnostics != nullptr) {
      diagnostics->concentration[descriptor.diagnostic_index] =
          toxin_concentration;
    }
    if (toxin_concentration <= 0.0) continue;
    any_toxin = true;
    total_kill += compute_receptor_hazard(
        agent, dt, descriptor, toxin_concentration, diagnostics);
  }
  if (!any_toxin) {
    return 0.0;
  }

  return std::min(1.0 - std::exp(-total_kill), 1.0);
}

Real FixReceptor::compute_receptor_hazard(
    const Agent& agent, Real dt, const ReceptorDescriptor& descriptor,
    Real toxin_concentration, KillAssessment* diagnostics) const {
  const Int receptor_index = to_underlying(descriptor.receptor);
  const Int ligand_index =
      sim_.chemical_field().find(descriptor.ligand_species);
  const Real ligand_concentration = ligand_index >= 0
      ? sim_.chemical_field().conc_global(ligand_index, agent.grid_cell)
          * descriptor.ligand_scale
      : 0.0;
  const Real occupancy = toxin_occupancy(
      toxin_concentration, ligand_concentration, cfg_.*descriptor.kd_toxin,
      cfg_.*descriptor.kd_ligand, agent.receptor_expr[receptor_index],
      agent.genome.toxin_affinity[receptor_index],
      agent.genome.ligand_affinity[receptor_index]);
  Real immunity_factor = 1.0;
  for (const auto& bi : agent.genome.bi_loci) {
    if (bi.target == descriptor.receptor) {
      const Real candidate =
          cfg_.immunity_factor * bi.immunity_binding_affinity;
      immunity_factor = std::min(immunity_factor, candidate);
    }
  }
  const Real hazard =
      cfg_.*descriptor.kill_rate * occupancy * immunity_factor * dt;
  if (diagnostics != nullptr) {
    diagnostics->occupancy[descriptor.diagnostic_index] = occupancy;
    diagnostics->hazard[descriptor.diagnostic_index] = hazard;
  }
  return hazard;
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
