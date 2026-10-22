/* ------------------------------------------------------------------
 * mechanism_capability.h — the tabulated layers' capability limit.
 *
 * Layer-2 patches and the Layer-3 chain do not run a live Simulation:
 * no fix plugins execute, no chemistry/diffusion runs, no motility or
 * contact killing happens. A run config that sets mechanism-dynamics
 * parameters on those layers is silently sweeping unused parameters —
 * the SPEC13 selection campaign's failure mode. The check below names
 * every configured parameter that has no consumer at the tabulated
 * layers so it is reported instead of accepted.
 *
 * Founder genotype state (initial_strains plasmids, receptor_expression,
 * CDI fields) is NOT listed: it is carried state that propagates through
 * transfers and restart — a phenotype declaration, not a mechanism knob.
 * ------------------------------------------------------------------ */

#ifndef GUTIBM_MECHANISM_CAPABILITY_H
#define GUTIBM_MECHANISM_CAPABILITY_H

#include "input_parser.h"
#include <string>
#include <vector>

namespace gutibm {

// Dotted config keys whose configured value differs from
// InputParser::default_config() and has no consumer on a tabulated
// layer. Empty when the config stays inside the tabulated surface.
[[nodiscard]] std::vector<std::string> unexecuted_mechanism_params(
    const SimulationConfig& cfg);

}  // namespace gutibm

#endif  // GUTIBM_MECHANISM_CAPABILITY_H
