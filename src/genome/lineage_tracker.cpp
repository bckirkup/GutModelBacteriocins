/* -----------------------------------------------------------------------
   GutIBM – Lineage tracker implementation
   ----------------------------------------------------------------------- */

#include "lineage_tracker.h"
#include <algorithm>

namespace gutibm {

void LineageTracker::init(Real snapshot_interval) {
  snapshot_interval_ = snapshot_interval;
  events_.clear();
  snapshots_.clear();
}

void LineageTracker::record_birth(TagID tag, TagID parent,
                                    TagID lineage, const Vec3& pos) {
  events_.push_back({LineageEvent::BIRTH, 0.0, tag, lineage, pos,
                     "parent=" + std::to_string(parent)});
}

void LineageTracker::record_death(TagID tag, TagID lineage, const Vec3& pos) {
  events_.push_back({LineageEvent::DEATH, 0.0, tag, lineage, pos, ""});
}

void LineageTracker::record_lysis(TagID tag, const Vec3& pos, TagID lineage) {
  events_.push_back({LineageEvent::LYSIS, 0.0, tag, lineage, pos, ""});
}

void LineageTracker::record_mutation(TagID tag, const std::string& type,
                                       TagID lineage) {
  events_.push_back({LineageEvent::MUTATION, 0.0, tag, lineage,
                     {0, 0, 0}, type});
}

void LineageTracker::record_hgt(TagID donor, TagID recipient,
                                  uint16_t toxin_id) {
  events_.push_back({LineageEvent::HGT, 0.0, recipient, 0,
                     {0, 0, 0},
                     "donor=" + std::to_string(donor) +
                     " toxin=" + std::to_string(toxin_id)});
}

void LineageTracker::record_washout(TagID tag, TagID lineage, const Vec3& pos) {
  events_.push_back({LineageEvent::WASHOUT, 0.0, tag, lineage, pos, ""});
}

void LineageTracker::take_snapshot(
    Real time,
    const std::vector<std::pair<TagID, TagID>>& agent_lineages) {

  LineageSnapshot snap;
  snap.time = time;
  snap.total_agents = static_cast<Int>(agent_lineages.size());

  for (const auto& [tag, lin] : agent_lineages) {
    snap.lineage_counts[lin]++;
  }

  snap.num_lineages = static_cast<Int>(snap.lineage_counts.size());

  // Find dominant lineage fraction
  Int max_count = 0;
  for (const auto& [lin, cnt] : snap.lineage_counts) {
    if (cnt > max_count) max_count = cnt;
  }
  snap.dominant_fraction = (snap.total_agents > 0)
      ? static_cast<Real>(max_count) / snap.total_agents
      : 0.0;

  snapshots_.push_back(std::move(snap));
}

Real LineageTracker::resident_retention(Real time_window) const {
  if (snapshots_.size() < 2) return 1.0;

  const auto& first = snapshots_.front();
  const auto& last  = snapshots_.back();

  if (last.time - first.time < time_window) return 1.0;

  // Count how many original lineages are still present
  Int retained = 0;
  for (const auto& [lin, cnt] : first.lineage_counts) {
    if (last.lineage_counts.count(lin) > 0) {
      retained++;
    }
  }

  return (first.num_lineages > 0)
      ? static_cast<Real>(retained) / first.num_lineages
      : 0.0;
}

TagID LineageTracker::dominant_lineage() const {
  if (snapshots_.empty()) return 0;

  const auto& last = snapshots_.back();
  TagID dom = 0;
  Int max_count = 0;
  for (const auto& [lin, cnt] : last.lineage_counts) {
    if (cnt > max_count) {
      max_count = cnt;
      dom = lin;
    }
  }
  return dom;
}

}  // namespace gutibm
