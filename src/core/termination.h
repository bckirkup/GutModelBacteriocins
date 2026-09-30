/* -----------------------------------------------------------------------
   GutIBM – Run termination causes
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_TERMINATION_H
#define GUTIBM_TERMINATION_H

#include "types.h"

#include <string_view>

namespace gutibm {

enum class TerminationCause : Int {
  HorizonReached = 0,
  DysbiosisGuard = 1,
  PopulationStop = 2,
  StopRequested = 3,
  ClosureViolation = 4,
  IncompleteUnknown = 5,
  // Spec 13 Layer-2 run-level causes (per-patch status stays distinct,
  // see PatchStatus in layer2_config.h).
  PatchInvalidLimit = 6,
  AuditDivergence = 7,
};

constexpr std::string_view termination_cause_name(TerminationCause cause) {
  using enum TerminationCause;
  switch (cause) {
    case HorizonReached: return "horizon_reached";
    case DysbiosisGuard: return "dysbiosis_guard";
    case PopulationStop: return "population_stop";
    case StopRequested: return "stop_requested";
    case ClosureViolation: return "closure_violation";
    case IncompleteUnknown: return "incomplete_unknown";
    case PatchInvalidLimit: return "patch_invalid_limit";
    case AuditDivergence: return "audit_divergence";
  }
  return "incomplete_unknown";
}

}  // namespace gutibm

#endif  // GUTIBM_TERMINATION_H
