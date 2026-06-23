/* -----------------------------------------------------------------------
   GutIBM – Stochastic mutation at BI loci and receptor genes
   
   During cell division, agents have probability of:
   1. BI locus duplication → expanded toxin range
   2. BI locus recombination → novel toxin/immunity combos
   3. Receptor downregulation mutation → resistance + metabolic penalty
   4. Super-killer emergence → toxin that bypasses cognate immunity
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_MUTATION_H
#define GUTIBM_FIX_MUTATION_H

#include "fix.h"
#include "agent.h"

namespace gutibm {

struct MutationConfig {
  // Per-division mutation rates
  Real bi_duplication_rate     = 1.0e-5;   // BI locus tandem duplication
  Real bi_recombination_rate   = 5.0e-6;   // BI locus recombination
  Real receptor_mutation_rate  = 1.0e-7;   // receptor downregulation
  Real super_killer_rate       = 1.0e-8;   // novel toxin variant

  // Receptor mutation effect
  Real receptor_reduction      = 0.1;      // expression drops by 0.1 per mutation

  // Max BI clusters per genome
  Int  max_bi_loci             = 8;

  // Immunity-escape: probability that a super-killer also reduces cognate
  // immunity binding affinity (VADI §57 affinity-neutralization matrix)
  Real immunity_escape_prob    = 0.5;     // fraction of super-killers with escape
  Real escape_affinity_lo      = 0.01;    // lower bound of reduced affinity
  Real escape_affinity_hi      = 0.3;     // upper bound (3–100× reduced protection)

  // Compensatory mutation: ameliorates plasmid metabolic cost over time
  Real compensatory_rate       = 1.0e-6;  // per division
  Real compensatory_reduction  = 0.005;   // per-locus cost reduced by this
};

class FixMutation : public Fix {
 public:
  FixMutation(Simulation& sim, const MutationConfig& cfg);

  void compute(Real dt) override;

 private:
  void mutate_on_division(Agent& agent);
  void duplicate_bi_locus(Agent& agent);
  void recombine_bi_locus(Agent& agent);
  void mutate_receptor(Agent& agent);
  void generate_super_killer(Agent& agent);
  void compensatory_mutation(Agent& agent);

  BICluster create_novel_toxin(const BICluster& parent);

  MutationConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_MUTATION_H
