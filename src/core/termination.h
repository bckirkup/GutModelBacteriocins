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
};

constexpr std::string_view termination_cause_name(TerminationCause cause) {
  switch (cause) {
    case TerminationCause::HorizonReached: return "horizon_reached";
    case TerminationCause::DysbiosisGuard: return "dysbiosis_guard";
    case TerminationCause::PopulationStop: return "population_stop";
    case TerminationCause::StopRequested: return "stop_requested";
    case TerminationCause::ClosureViolation: return "closure_violation";
    case TerminationCause::IncompleteUnknown: return "incomplete_unknown";
  }
  return "incomplete_unknown";
}

}  // namespace gutibm

#endif  // GUTIBM_TERMINATION_H
