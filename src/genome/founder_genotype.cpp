#include "founder_genotype.h"

#include "agent.h"
#include "error.h"
#include "plasmid.h"
#include "types.h"

#include <iostream>

namespace gutibm {

BICluster configured_plasmid_cluster(const SimulationConfig& cfg,
                                     const PlasmidEntry& entry) {
  BICluster cluster = entry.cluster;
  cluster.retardation =
      retardation_from_pI(cluster.pI, cfg.fixes.bacteriocin.mucin_charge);
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
                            const SimulationConfig& cfg) {
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
}

}  // namespace gutibm
