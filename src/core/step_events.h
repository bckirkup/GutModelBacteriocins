/* -----------------------------------------------------------------------
   GutIBM – Per-step event counters for HDF5 summary layer (Spec 4)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_STEP_EVENTS_H
#define GUTIBM_STEP_EVENTS_H

#include "types.h"

namespace gutibm {

struct StepEvents {
  Int sos_inductions = 0;
  Int phage_inductions = 0;
  Int colicin_kills = 0;
  Int cdi_kills = 0;
  Int washout_deaths = 0;
  Int boundary_deaths = 0;
  Int starvation_deaths = 0;
  Int divisions = 0;
  Int conjugation_transfers = 0;
  Int mutations = 0;

  void reset() { *this = StepEvents{}; }
};

}  // namespace gutibm

#endif  // GUTIBM_STEP_EVENTS_H
