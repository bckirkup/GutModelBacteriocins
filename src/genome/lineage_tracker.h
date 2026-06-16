/* -----------------------------------------------------------------------
   GutIBM – Lineage tracker and allele time-series recorder
   
   Tracks every agent's ancestry, mutations, and HGT events to
   enable genomic validation:
   - Sweep-and-stasis pattern detection
   - 70–80% resident strain retention rate
   - Accessory genome complexity of residents vs. transients
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_LINEAGE_TRACKER_H
#define GUTIBM_LINEAGE_TRACKER_H

#include "types.h"
#include <string>
#include <vector>
#include <unordered_map>

namespace gutibm {

struct LineageEvent {
  enum Type { BIRTH, DEATH, LYSIS, MUTATION, HGT, WASHOUT };
  Type type;
  Real time;
  TagID agent_tag;
  TagID lineage_id;
  Vec3 position;
  std::string detail;
};

struct LineageSnapshot {
  Real time;
  std::unordered_map<TagID, Int> lineage_counts;  // lineage_id → count
  Int total_agents;
  Int num_lineages;
  Real dominant_fraction;  // fraction of pop belonging to dominant lineage
};

class LineageTracker {
 public:
  LineageTracker() = default;

  void init(Real snapshot_interval);

  // Record events
  void record_birth(TagID tag, TagID parent, TagID lineage, const Vec3& pos);
  void record_death(TagID tag, TagID lineage, const Vec3& pos);
  void record_lysis(TagID tag, const Vec3& pos, TagID lineage);
  void record_mutation(TagID tag, const std::string& type, TagID lineage);
  void record_hgt(TagID donor, TagID recipient, uint16_t toxin_id);
  void record_washout(TagID tag, TagID lineage, const Vec3& pos);

  // Take a population snapshot
  void take_snapshot(Real time, const std::vector<std::pair<TagID, TagID>>& agent_lineages);

  // Access recorded data
  const std::vector<LineageEvent>& events() const { return events_; }
  const std::vector<LineageSnapshot>& snapshots() const { return snapshots_; }

  // Compute resident retention rate: fraction of original lineages
  // still present after N timesteps
  Real resident_retention(Real time_window) const;

  // Identify dominant lineage (super-killer candidate)
  TagID dominant_lineage() const;

 private:
  std::vector<LineageEvent> events_;
  std::vector<LineageSnapshot> snapshots_;
  Real snapshot_interval_ = 3600.0;  // default: every hour
};

}  // namespace gutibm

#endif  // GUTIBM_LINEAGE_TRACKER_H
