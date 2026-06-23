/* -----------------------------------------------------------------------
   GutIBM – Mutation fix implementation
   ----------------------------------------------------------------------- */

#include "fix_mutation.h"
#include "plasmid.h"
#include "simulation.h"
#include <algorithm>

namespace gutibm {

FixMutation::FixMutation(Simulation& sim, const MutationConfig& cfg)
    : Fix("mutation", sim), cfg_(cfg) {}

void FixMutation::compute(Real dt) {
  // Mutations occur during division events.
  // This fix is called after metabolism's division pass.
  // We check recently-divided agents (age ~0).
  auto& agents = sim_.agents();

  for (Int i = 0; i < agents.size(); ++i) {
    Agent& a = agents[i];
    if (a.state == PhenoState::DEAD) continue;
    if (a.age > dt) continue;  // only newly divided cells

    mutate_on_division(a);
  }
}

void FixMutation::mutate_on_division(Agent& agent) {
  auto& rng = sim_.rng();

  // BI locus duplication
  if (rng.bernoulli(cfg_.bi_duplication_rate)) {
    duplicate_bi_locus(agent);
  }

  // BI locus recombination
  if (rng.bernoulli(cfg_.bi_recombination_rate)) {
    recombine_bi_locus(agent);
  }

  // Receptor downregulation mutation
  if (rng.bernoulli(cfg_.receptor_mutation_rate)) {
    mutate_receptor(agent);
  }

  // Partial resistance: extracellular loop missense mutation
  if (rng.bernoulli(cfg_.partial_resistance_rate)) {
    partial_resistance_mutation(agent);
  }

  // Super-killer emergence
  if (rng.bernoulli(cfg_.super_killer_rate)) {
    generate_super_killer(agent);
  }

  // Compensatory chromosomal mutation that ameliorates plasmid cost
  if (rng.bernoulli(cfg_.compensatory_rate)) {
    compensatory_mutation(agent);
  }
}

void FixMutation::duplicate_bi_locus(Agent& agent) {
  if (agent.genome.bi_loci.empty()) return;
  if (static_cast<Int>(agent.genome.bi_loci.size()) >= cfg_.max_bi_loci) return;

  auto& rng = sim_.rng();
  Int idx = rng.randint(0, static_cast<Int>(agent.genome.bi_loci.size()) - 1);
  agent.genome.bi_loci.push_back(agent.genome.bi_loci[idx]);
  agent.genome.mutations++;

  sim_.lineage_tracker().record_mutation(agent.tag, "bi_duplication",
                                          agent.genome.lineage_id);
}

void FixMutation::recombine_bi_locus(Agent& agent) {
  if (agent.genome.bi_loci.size() < 2) return;

  auto& rng = sim_.rng();
  Int i1 = rng.randint(0, static_cast<Int>(agent.genome.bi_loci.size()) - 1);
  Int i2 = rng.randint(0, static_cast<Int>(agent.genome.bi_loci.size()) - 1);
  if (i1 == i2) return;

  // Swap immunity between two BI clusters
  std::swap(agent.genome.bi_loci[i1].immunity_id,
            agent.genome.bi_loci[i2].immunity_id);
  agent.genome.mutations++;

  sim_.lineage_tracker().record_mutation(agent.tag, "bi_recombination",
                                          agent.genome.lineage_id);
}

void FixMutation::mutate_receptor(Agent& agent) {
  auto& rng = sim_.rng();
  Int receptor = rng.randint(0, NUM_RECEPTORS - 1);

  agent.receptor_expr[receptor] =
      std::max(0.0, agent.receptor_expr[receptor] - cfg_.receptor_reduction);
  agent.genome.receptor_expression[receptor] = agent.receptor_expr[receptor];
  agent.genome.mutations++;

  // If expression is very low, mark as resistant
  if (agent.receptor_expr[receptor] < 0.2) {
    agent.state = PhenoState::RESISTANT;
  }

  sim_.lineage_tracker().record_mutation(agent.tag, "receptor_downreg",
                                          agent.genome.lineage_id);
}

void FixMutation::partial_resistance_mutation(Agent& agent) {
  auto& rng = sim_.rng();
  Int receptor = rng.randint(0, NUM_RECEPTORS - 1);

  // Reduce toxin binding affinity by 10–100x (Kd multiplier 0.01–0.1)
  agent.genome.toxin_affinity[receptor] = rng.uniform(0.01, 0.1);

  // Ligand affinity stays within 2x of wild-type (0.5–1.0)
  agent.genome.ligand_affinity[receptor] = rng.uniform(0.5, 1.0);

  // Receptor expression is NOT changed (distinct from full downregulation)
  agent.genome.mutations++;

  sim_.lineage_tracker().record_mutation(agent.tag, "partial_resistance",
                                          agent.genome.lineage_id);
}

void FixMutation::generate_super_killer(Agent& agent) {
  if (agent.genome.bi_loci.empty()) return;

  auto& rng = sim_.rng();
  Int idx = rng.randint(0, static_cast<Int>(agent.genome.bi_loci.size()) - 1);

  BICluster novel = create_novel_toxin(agent.genome.bi_loci[idx]);

  if (static_cast<Int>(agent.genome.bi_loci.size()) < cfg_.max_bi_loci) {
    agent.genome.bi_loci.push_back(novel);
  }
  agent.genome.mutations++;

  const char* label = (novel.immunity_binding_affinity < 1.0)
                          ? "super_killer_escape"
                          : "super_killer";
  sim_.lineage_tracker().record_mutation(agent.tag, label,
                                          agent.genome.lineage_id);
}

void FixMutation::compensatory_mutation(Agent& agent) {
  if (agent.genome.bi_loci.empty()) return;

  // Compensatory chromosomal mutations rapidly ameliorate the metabolic
  // drain of plasmid maintenance over time (VADI §79).
  // This reduces the effective per-locus cost in fix_metabolism.
  agent.genome.plasmid_cost_amelioration += cfg_.compensatory_reduction;
  // Cap: cannot eliminate more than 75% of original cost
  agent.genome.plasmid_cost_amelioration =
      std::min(agent.genome.plasmid_cost_amelioration, 0.015);
  agent.genome.mutations++;

  sim_.lineage_tracker().record_mutation(agent.tag, "compensatory",
                                          agent.genome.lineage_id);
}

BICluster FixMutation::create_novel_toxin(const BICluster& parent) {
  auto& rng = sim_.rng();
  BICluster novel = parent;

  // New toxin ID (unique)
  novel.toxin_id = static_cast<uint16_t>(rng.randint(1000, 65000));

  // Slightly altered properties
  novel.pI += rng.gaussian(0.0, 0.5);
  novel.pI = std::clamp(novel.pI, 3.0, 12.0);

  // Reclassify based on new pI
  novel.bclass = classify_by_pI(novel.pI);

  // Altered target receptor (may hijack a different TBDT)
  if (rng.bernoulli(0.3)) {
    novel.target = static_cast<ReceptorType>(rng.randint(0, NUM_RECEPTORS - 1));
  }

  // Immunity-escape mutation: reduce cognate immunity binding affinity
  // so that existing immunity proteins provide less cross-protection
  // against this novel toxin variant (VADI §57)
  if (rng.bernoulli(cfg_.immunity_escape_prob)) {
    novel.immunity_binding_affinity =
        rng.uniform(cfg_.escape_affinity_lo, cfg_.escape_affinity_hi);
  }

  return novel;
}

}  // namespace gutibm
