/* -----------------------------------------------------------------------
   GutIBM – Per-step event counters for HDF5 summary layer (Spec 4)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_STEP_EVENTS_H
#define GUTIBM_STEP_EVENTS_H

#include "types.h"
#include <array>

namespace gutibm {

enum class ProvenanceCause : Int {
  COLICIN = 0,
  CDI = 1,
  WASHOUT = 2,
  BOUNDARY = 3,
  STARVATION = 4,
};

struct KillProvenanceEvent {
  TagID victim_id = 0;
  Vec3 position{};
  Int strain = 0;
  ProvenanceCause cause = ProvenanceCause::COLICIN;
  TagID cdi_attacker_id = 0;
  bool cdi_attacker_known = false;
  std::array<Real, 4> toxin_concentration{};
  std::array<Real, 4> toxin_occupancy{};
  std::array<Real, 4> toxin_hazard{};
};

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
