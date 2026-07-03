/* -----------------------------------------------------------------------
   GutIBM – Conjugation fix implementation
   ----------------------------------------------------------------------- */

#include "fix_conjugation.h"
#include "simulation.h"
#include <cmath>

namespace gutibm {

FixConjugation::FixConjugation(Simulation& sim, const ConjugationConfig& cfg)
    : Fix("conjugation", sim), cfg_(cfg) {}

void FixConjugation::compute(Real dt) {
  auto& agents = sim_.agents();
  const auto& hash   = sim_.domain().spatial_hash();
  const auto& adv    = sim_.advection();
  auto& rng    = sim_.rng();

  for (Int i = 0; i < agents.size(); ++i) {
    Agent& donor = agents[i];
    if (donor.state == PhenoState::DEAD) continue;
    if (!donor.genome.has_conjugative_plasmid) continue;
    if (donor.genome.bi_loci.empty()) continue;

    // Determine effective pilus reach for this conjugation attempt
    Real effective_radius = cfg_.pili_length;
    if (cfg_.pili_heterogeneity) {
      effective_radius = rng.uniform(cfg_.pili_length_min, cfg_.pili_length_max);
    }

    // Find nearby potential recipients within pili reach
    auto neighbors = hash.query_radius(donor.x, effective_radius);

    for (Int j : neighbors) {
      if (j == i) continue;
      Agent& recipient = agents[j];
      if (recipient.state == PhenoState::DEAD) continue;

      // Check distance against sampled pilus length
      if (Real d2 = sim_.domain().min_image_dist_sq(donor.x, recipient.x);
          d2 > effective_radius * effective_radius) continue;

      // Local shear rate
      Real shear = adv.shear_rate(donor.x);

      attempt_transfer(donor, recipient, shear, dt);
    }
  }
}

void FixConjugation::attempt_transfer(Agent& donor, Agent& recipient,
                                        Real shear, Real dt) {
  auto& rng = sim_.rng();

  // MPS probability: inversely proportional to shear
  Real p_mps = std::exp(-shear / cfg_.shear_critical);
  if (Real p_transfer = cfg_.base_transfer_rate * p_mps * dt;
      !rng.bernoulli(p_transfer)) return;

  // Transfer a random BI cluster from donor to recipient
  if (donor.genome.bi_loci.empty()) return;

  Int idx = rng.randint(0, static_cast<Int>(donor.genome.bi_loci.size()) - 1);
  const BICluster& cluster = donor.genome.bi_loci[idx];

  // Check if recipient already has this cluster
  for (const auto& bi : recipient.genome.bi_loci) {
    if (bi.toxin_id == cluster.toxin_id && bi.immunity_id == cluster.immunity_id) {
      return;  // already present
    }
  }

  // Transfer
  recipient.genome.bi_loci.push_back(cluster);
  recipient.genome.has_conjugative_plasmid = true;

  // Record HGT event
  sim_.lineage_tracker().record_hgt(donor.identity.tag, recipient.identity.tag,
                                     cluster.toxin_id);
}

}  // namespace gutibm
