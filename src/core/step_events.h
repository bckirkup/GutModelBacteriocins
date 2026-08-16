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
  LYSIS = 5,
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
  Int lysis_deaths = 0;
  Int divisions = 0;
  Int conjugation_transfers = 0;
  Int mutations = 0;
  Int immigrations = 0;

  void add(const StepEvents& other) {
    sos_inductions += other.sos_inductions;
    phage_inductions += other.phage_inductions;
    colicin_kills += other.colicin_kills;
    cdi_kills += other.cdi_kills;
    washout_deaths += other.washout_deaths;
    boundary_deaths += other.boundary_deaths;
    starvation_deaths += other.starvation_deaths;
    lysis_deaths += other.lysis_deaths;
    divisions += other.divisions;
    conjugation_transfers += other.conjugation_transfers;
    mutations += other.mutations;
    immigrations += other.immigrations;
  }

  void reset() { *this = StepEvents{}; }
};

struct MechanicsStats {
  Int displacement_clamps = 0;

  void add(const MechanicsStats& other) {
    displacement_clamps += other.displacement_clamps;
  }

  void reset() { *this = MechanicsStats{}; }
};

}  // namespace gutibm

#endif  // GUTIBM_STEP_EVENTS_H
