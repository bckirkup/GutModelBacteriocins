/* ------------------------------------------------------------------
 * founder_genotype.h — single-source founder-genotype application.
 *
 * Every founder-construction site — live Simulation founders, Layer-2
 * patch founders and initial pool cells, and Layer-3 lumen exemplars —
 * must run the same declaration so a configured phenotype (declared
 * receptor expression, plasmid clusters, CDI type/immunity, conjugative
 * plasmid flag) survives every supported route: tabulated construction,
 * deposits/transfers, and checkpoint/restart.
 * ------------------------------------------------------------------ */

#ifndef GUTIBM_FOUNDER_GENOTYPE_H
#define GUTIBM_FOUNDER_GENOTYPE_H

#include "input_parser.h"
#include <string>
#include <vector>

namespace gutibm {

class Agent;
struct BICluster;
struct PlasmidEntry;

// The cluster a plasmid declaration actually carries: the library entry
// with pI-derived retardation (mucin-charge context), then any
// plasmid_overrides.<name>.* fields on top.
[[nodiscard]] BICluster configured_plasmid_cluster(
    const SimulationConfig& cfg, const PlasmidEntry& entry);

// Attach the declared plasmids to the agent's genome. Unknown names warn
// once (rank 0) and are skipped; a conjugative library entry sets the
// conjugative-plasmid flag.
void assign_plasmids(Agent& agent, const std::vector<std::string>& plasmids,
                     const SimulationConfig& cfg, int rank);

// Apply one initial_strains[] declaration to a freshly created agent:
// declared receptor expression (including the PhenoState::RESISTANT
// marker), plasmid BI clusters via configured_plasmid_cluster, and CDI
// type/immunity. Throws ConfigError on an unknown receptor name.
void apply_founder_genotype(Agent& agent,
                            const SimulationConfig::InitialStrain& strain,
                            const SimulationConfig& cfg);

}  // namespace gutibm

#endif  // GUTIBM_FOUNDER_GENOTYPE_H
