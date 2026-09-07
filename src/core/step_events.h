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
  LYSIS = 5,
  // Value 4 was STARVATION and is retired.
};

struct KillProvenanceEvent {
  TagID victim_id = 0;
  // Simulation clock at the kill, stamped by Simulation::record_kill_provenance.
  // Same clock as ToxinBurstSource::creation_time, so a LYSIS event's
  // event_time_s is the release-window origin of the burst it spawned.
  Int event_step = 0;
  Real event_time_s = 0.0;
  Vec3 position{};
  Int strain = 0;
  ProvenanceCause cause = ProvenanceCause::COLICIN;
  TagID cdi_attacker_id = 0;
  bool cdi_attacker_known = false;
  std::array<Real, 4> toxin_concentration{};
  std::array<Real, 4> toxin_occupancy{};
  std::array<Real, 4> toxin_hazard{};
};

// Matches summary n_by_type width in hdf5_writer.cpp.
static constexpr Int MAX_AGENT_TYPES = 8;

struct StepEvents {
  Int sos_inductions = 0;
  Int phage_inductions = 0;
  Int mortality_colicin = 0;
  Int mortality_cdi = 0;
  Int outflow_washout = 0;
  Int outflow_boundary = 0;
  Int mortality_lysis = 0;
  Int divisions = 0;
  // Per-identity-type mother divisions; divisions == sum(divisions_by_type)
  // when all increments go through record_division().
  std::array<Int, MAX_AGENT_TYPES> divisions_by_type{};
  Int conjugation_transfers = 0;
  Int mutations = 0;
  Int immigrations = 0;
  Real sos_basal_rate = 0.0;
  Real sos_post_division_rate = 0.0;
  Real sos_nuclease_cross_induction_rate = 0.0;
  Real sos_ros_rate = 0.0;

  void record_division(Int agent_type) {
    ++divisions;
    if (agent_type >= 0 && agent_type < MAX_AGENT_TYPES) {
      ++divisions_by_type[static_cast<size_t>(agent_type)];
    }
  }

  void add(const StepEvents& other) {
    sos_inductions += other.sos_inductions;
    phage_inductions += other.phage_inductions;
    mortality_colicin += other.mortality_colicin;
    mortality_cdi += other.mortality_cdi;
    outflow_washout += other.outflow_washout;
    outflow_boundary += other.outflow_boundary;
    mortality_lysis += other.mortality_lysis;
    divisions += other.divisions;
    for (Int i = 0; i < MAX_AGENT_TYPES; ++i) {
      divisions_by_type[static_cast<size_t>(i)] +=
          other.divisions_by_type[static_cast<size_t>(i)];
    }
    conjugation_transfers += other.conjugation_transfers;
    mutations += other.mutations;
    immigrations += other.immigrations;
    sos_basal_rate += other.sos_basal_rate;
    sos_post_division_rate += other.sos_post_division_rate;
    sos_nuclease_cross_induction_rate += other.sos_nuclease_cross_induction_rate;
    sos_ros_rate += other.sos_ros_rate;
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
