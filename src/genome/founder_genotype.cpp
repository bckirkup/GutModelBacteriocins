#include "founder_genotype.h"

#include "agent.h"
#include "error.h"
#include "plasmid.h"
#include "random.h"
#include "types.h"

#include <cmath>
#include <cstdint>
#include <iostream>

namespace gutibm {

namespace {

// Deterministic [0, 1) draw from (seed, tag, salt) when no RNG is supplied
// (e.g. tabulated Layer-2/Layer-3 founder construction).
Real deterministic_unit_draw(uint64_t seed, TagID tag, uint64_t salt) {
  uint64_t z = seed ^ (static_cast<uint64_t>(tag) * 0x9e3779b97f4a7c15ULL) ^ salt;
  z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
  z = z ^ (z >> 31U);
  constexpr Real kScale = 1.0 / 9007199254740992.0;  // 2^-53
  return static_cast<Real>(z >> 11U) * kScale;
}

bool draw_prophage_bernoulli(Real prob, RNG* rng, uint64_t seed,
                             TagID tag, uint64_t salt) {
  if (prob <= 0.0) return false;
  if (prob >= 1.0) return true;
  if (rng != nullptr) {
    return rng->bernoulli(prob);
  }
  return deterministic_unit_draw(seed, tag, salt) < prob;
}

void apply_founder_prophage_state(Agent& agent,
                                  const SimulationConfig::InitialStrain& strain,
                                  const SimulationConfig& cfg,
                                  RNG* rng) {
  const bool explicit_override = strain.lysogen_prevalence.has_value()
      || strain.active_inducible_fraction.has_value();
  if (!cfg.prophage.enabled && !explicit_override) {
    return;
  }

  const Real prevalence =
      strain.lysogen_prevalence.value_or(cfg.prophage.lysogen_prevalence_init);
  const Real active_frac = strain.active_inducible_fraction.value_or(
      cfg.prophage.active_inducible_fraction);
  if (!std::isfinite(prevalence) || prevalence < 0.0 || prevalence > 1.0) {
    throw ConfigError("founder lysogen_prevalence must be finite and in [0, 1]");
  }
  if (!std::isfinite(active_frac) || active_frac < 0.0 || active_frac > 1.0) {
    throw ConfigError(
        "founder active_inducible_fraction must be finite and in [0, 1]");
  }

  const bool has_immunity = draw_prophage_bernoulli(
      prevalence, rng, cfg.seed, agent.identity.tag, 0x50524f5048414745ULL);
  agent.flags.has_prophage_immunity = has_immunity;
  agent.flags.has_inducible_prophage = has_immunity
      && draw_prophage_bernoulli(active_frac, rng, cfg.seed,
                                 agent.identity.tag, 0x494e44554349424cULL);
  agent.flags.prophage_induced = false;
  agent.timers.prophage_lysis_timer = -1.0;
}

}  // namespace

BICluster configured_plasmid_cluster(const SimulationConfig& cfg,
                                     const PlasmidEntry& entry) {
  BICluster cluster = entry.cluster;
  cluster.retardation =
      retardation_from_pI(cluster.pI, cfg.fixes.bacteriocin.mucin_charge);
  if (cluster.release_mode != ReleaseMode::CONTINUOUS
      && (cfg.fixes.bacteriocin.D_free_colicin_explicit
          || std::abs(cfg.fixes.bacteriocin.D_free_colicin
                      - kDefaultColicinFreeDiffCoeff) > 0.0)) {
    cluster.diff_coeff = cfg.fixes.bacteriocin.D_free_colicin;
  }
  if (const auto it = cfg.plasmid_overrides.find(entry.name);
      it != cfg.plasmid_overrides.end()) {
    const auto& values = it->second;
    if (values.retardation.has_value()) {
      cluster.retardation = *values.retardation;
    }
    if (values.diff_coeff.has_value()) {
      cluster.diff_coeff = *values.diff_coeff;
    }
    if (values.burst_size.has_value()) {
      cluster.burst_size = *values.burst_size;
    }
  }
  return cluster;
}

void assign_plasmids(Agent& agent, const std::vector<std::string>& plasmids,
                     const SimulationConfig& cfg, int rank) {
  for (const auto& pname : plasmids) {
    const PlasmidEntry* entry = PlasmidLibrary::find(pname);
    if (!entry) {
      if (rank == 0) {
        std::cerr << "Warning: unknown plasmid '" << pname
                  << "' — agent spawned without BI locus\n";
      }
      continue;
    }
    BICluster cluster = configured_plasmid_cluster(cfg, *entry);
    agent.genome.bi_loci.push_back(cluster);
    if (entry->conjugative) {
      agent.genome.has_conjugative_plasmid = true;
    }
  }
}

void apply_founder_genotype(Agent& agent,
                            const SimulationConfig::InitialStrain& strain,
                            const SimulationConfig& cfg,
                            RNG* rng) {
  assign_plasmids(agent, strain.plasmids, cfg, /*rank=*/0);
  for (const auto& [name, expression] : strain.receptor_expression) {
    const auto receptor = receptor_type_from_name(name);
    if (!receptor.has_value()) {
      throw ConfigError("unknown receptor name: " + name);
    }
    const Int index = to_underlying(*receptor);
    agent.receptor_expr_base[index] = expression;
    agent.receptor_expr[index] = expression;
    agent.genome.receptor_expression[index] = expression;
    if (receptor_expression_is_resistant(expression)) {
      agent.state = PhenoState::RESISTANT;
    }
  }
  agent.genome.cdi_type = strain.cdi_type;
  agent.genome.cdi_immunity = strain.cdi_immunity;
  apply_founder_prophage_state(agent, strain, cfg, rng);
}

}  // namespace gutibm
