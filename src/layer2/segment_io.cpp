/* -----------------------------------------------------------------------
   GutIBM – Spec 13 Phase 1: Layer-2 segment run loop + I/O.

   CSV timeseries and JSON provenance are plain stream writes. The
   checkpoint path is HDF5 and reuses agent_transfer_serialize for all
   cell payloads (G4); frozen HDF5 field meanings are untouched because
   the segment writes its own file, never the simulation's datasets.

   Phase-1 checkpoint limitation: a segment containing a live audit slot
   refuses to checkpoint — the embedded Simulation's field state is not
   segment-serializable.
   ----------------------------------------------------------------------- */

#include "error.h"
#include "h5_util.h"
#include "input_parser.h"
#include "patch_table.h"
#include "path_utils.h"
#include "segment.h"
#include "simulation.h"

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#include "agent_transfer.h"
#endif

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace gutibm {

#ifdef GUTIBM_HDF5
using namespace h5;
#endif

namespace {

const char* kTimeseriesHeader =
    "step,time_s,"
    "occ_crypt,occ_exposed_proximal,occ_exposed_distal,"
    "n_patches_crypt,n_patches_exposed_proximal,n_patches_exposed_distal,"
    "n_halted_crypt,n_halted_exposed_proximal,n_halted_exposed_distal,"
    "density_occupied_crypt_cfu_ml,density_occupied_proximal_cfu_ml,"
    "density_occupied_distal_cfu_ml,"
    "segment_mean_cfu_ml,segment_cfu_cm2,"
    "n_occupied_total,patch_cells,pool_cells,pool_packets,"
    "births_cum,washout_cum,contraction_cum,reseeds_cum,distal_cum,"
    "edge_cum,external_cum,deaths_cum,sos_lysis_deaths_cum,toxin_kill_deaths_cum,"
    "pulse_washout_cum,pulse_deaths_cum";

const char* kLiveLineageHeader =
    "step,time_s,patch_id,lineage_id,alive_count,emigrants_cum";

std::string json_escape(const std::string& s) {
  std::string out;
  out.reserve(s.size() + 8);
  for (const char ch : s) {
    switch (ch) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\t': out += "\\t"; break;
      default: out += ch; break;
    }
  }
  return out;
}

void write_live_physics_json(std::ofstream& out, const Layer2Config& c) {
  out << "  \"live_physics\": {\"mode\": \"" << json_escape(c.live_physics_mode)
      << "\", \"crypt_flow_for_crypt_patches\": "
      << (c.resolved_live_crypt_flow() ? "true" : "false")
      << ", \"inherit_oxygen\": "
      << (c.resolved_live_inherit_oxygen() ? "true" : "false")
      << ", \"align_patch_table_delivery\": "
      << (c.resolved_live_align_delivery() ? "true" : "false")
      << ", \"crypt_contraction_shelter_scale\": "
      << c.live_crypt_contraction_shelter_scale
      << ", \"layer2_crypt_disruption_prob\": " << c.crypt.disruption_prob
      << ", \"layer2_exposed_proximal_disruption_prob\": "
      << c.exposed_proximal.disruption_prob
      << ", \"effective_live_crypt_disruption_prob\": "
      << c.effective_live_crypt_disruption_prob()
      << ", \"double_counting_risk\": "
      << (c.live_crypt_double_counting_risk() ? "true" : "false") << "},\n";
}

void write_pulses_json(std::ofstream& out,
                       const std::vector<Layer2Pulse>& pulses,
                       const std::vector<PulseEventRecord>& events) {
  out << "  \"pulses\": [";
  for (size_t i = 0; i < pulses.size(); ++i) {
    const auto& p = pulses[i];
    if (i > 0) out << ", ";
    out << "{\"time_s\": " << p.time_s
        << ", \"washout_fraction\": " << p.washout_fraction
        << ", \"mortality_fraction\": " << p.mortality_fraction
        << ", \"patch_types\": [";
    for (size_t j = 0; j < p.patch_types.size(); ++j) {
      if (j > 0) out << ", ";
      out << '\"' << patch_type_name(p.patch_types[j]) << '\"';
    }
    out << "]}";
  }
  out << "],\n";
  out << "  \"pulse_events\": [";
  for (size_t i = 0; i < events.size(); ++i) {
    const auto& ev = events[i];
    if (i > 0) out << ", ";
    out << "{\"pulse_index\": " << ev.pulse_index
        << ", \"step\": " << ev.step
        << ", \"time_s\": " << ev.time_s
        << ", \"washed_out\": " << ev.washed_out
        << ", \"killed\": " << ev.killed
        << ", \"patches_affected\": " << ev.patches_affected << '}';
  }
  out << "],\n";
}

}  // namespace

int MucusSegment::run() {
  const auto& c = cfg_->layer2;
  const Real dt = cfg_->time.bio_dt;
  Int last_ts_step = -1;
  Int last_lin_step = -1;

  std::ofstream ts;
  if (!c.timeseries_file.empty()) {
    validate_output_file_path(c.timeseries_file);
    ts.open(c.timeseries_file, std::ios::trunc);
    if (!ts) {
      throw IOError("layer2 cannot open timeseries file: " +
                    c.timeseries_file);
    }
    ts << kTimeseriesHeader << '\n';
    emit_timeseries_row(ts);  // t=0 row
    last_ts_step = step_count_;
  }

  std::ofstream lin;
  if (!c.audit.lineage_file.empty() && !live_.empty()) {
    validate_output_file_path(c.audit.lineage_file);
    lin.open(c.audit.lineage_file, std::ios::trunc);
    if (!lin) {
      throw IOError("layer2 cannot open audit lineage file: " +
                    c.audit.lineage_file);
    }
    lin << kLiveLineageHeader << '\n';
    emit_lineage_rows(lin);  // t=0 census
    last_lin_step = step_count_;
  }

  while (time_ < cfg_->time.total_time &&
         termination_cause_ == TerminationCause::IncompleteUnknown) {
    step(dt);
    if (c.summary_interval_steps > 0 &&
        step_count_ % c.summary_interval_steps == 0) {
      if (ts.is_open()) {
        emit_timeseries_row(ts);
        last_ts_step = step_count_;
      }
      if (lin.is_open()) {
        emit_lineage_rows(lin);
        last_lin_step = step_count_;
      }
    }
    if (!c.checkpoint_file.empty() && c.checkpoint_interval_steps > 0 &&
        step_count_ % c.checkpoint_interval_steps == 0) {
      write_checkpoint_now();
    }
  }
  if (termination_cause_ == TerminationCause::IncompleteUnknown &&
      time_ >= cfg_->time.total_time) {
    terminate(TerminationCause::HorizonReached, "reached configured horizon");
  }
  if (ts.is_open() && step_count_ != last_ts_step) {
    emit_timeseries_row(ts);
  }
  if (lin.is_open() && step_count_ != last_lin_step) {
    emit_lineage_rows(lin);
  }
  if (c.checkpoint_final && !c.checkpoint_file.empty()) {
    write_checkpoint_now();
  }
  if (!c.provenance_file.empty()) {
    write_provenance(c.provenance_file);
  }
  std::cout << "layer2 segment terminated: "
            << termination_cause_name(termination_cause_) << " ("
            << termination_detail_ << ")\n";
  return termination_cause_ == TerminationCause::ClosureViolation ? 2 : 0;
}

void MucusSegment::emit_timeseries_row(std::ofstream& out) const {
  const SegmentObservables o = observables();
  out << step_count_ << ',' << time_ << ',';
  for (Int t = 0; t < kPatchTypeCount; ++t) {
    out << o.per_type[t].occupancy << ',';
  }
  for (Int t = 0; t < kPatchTypeCount; ++t) {
    out << o.per_type[t].n_patches << ',';
  }
  for (Int t = 0; t < kPatchTypeCount; ++t) {
    out << o.per_type[t].n_halted << ',';
  }
  for (Int t = 0; t < kPatchTypeCount; ++t) {
    out << o.per_type[t].mean_density_occupied << ',';
  }
  out << o.segment_mean_cfu_ml << ',' << o.segment_cfu_cm2 << ','
      << o.n_occupied_total << ',' << o.patch_cells << ',' << o.pool_cells
      << ',' << o.pool_packets << ',' << ledger_.births << ','
      << ledger_.washout_departures << ',' << ledger_.contraction_departures
      << ',' << ledger_.reseeds << ',' << ledger_.distal_losses << ','
      << ledger_.edge_departures << ',' << ledger_.external_arrivals
      << ',' << ledger_.biological_deaths << ',' << ledger_.sos_lysis_deaths
      << ',' << ledger_.toxin_kill_deaths << ',' << ledger_.pulse_washout
      << ',' << ledger_.pulse_deaths << '\n';
}

void MucusSegment::emit_lineage_rows(std::ofstream& out) const {
  // Per-lineage census of each live patch: alive stocks plus cumulative
  // patch->pool exports, keyed by founder lineage_id. This is the
  // lineage-resolved readout the Layer-1 reference comparison uses.
  for (const auto& [id, st] : live_) {
    std::unordered_map<TagID, Int> alive;
    for (const auto& a : st.sim->agents()) {
      if (a.state != PhenoState::DEAD) ++alive[a.genome.lineage_id];
    }
    std::vector<TagID> keys;
    keys.reserve(alive.size() + st.emigrants_by_lineage.size());
    for (const auto& [lin, n] : alive) {
      (void)n;
      keys.push_back(lin);
    }
    for (const auto& [lin, n] : st.emigrants_by_lineage) {
      (void)n;
      if (!alive.contains(lin)) keys.push_back(lin);
    }
    std::sort(keys.begin(), keys.end());
    for (const TagID lin : keys) {
      const Int n_alive = alive.contains(lin) ? alive.at(lin) : 0;
      const Int n_em = st.emigrants_by_lineage.contains(lin)
                           ? st.emigrants_by_lineage.at(lin)
                           : 0;
      out << step_count_ << ',' << time_ << ',' << id << ',' << lin << ','
          << n_alive << ',' << n_em << '\n';
    }
  }
}

void MucusSegment::write_provenance(const std::string& path) const {
  validate_output_file_path(path);
  std::ofstream out(path, std::ios::trunc);
  if (!out) {
    throw IOError("layer2 cannot open provenance file: " + path);
  }
  const auto& c = cfg_->layer2;
  const SegmentObservables o = observables();
  out << "{\n";
  out << "  \"spec\": \"spec13_phase1_layer2\",\n";
  out << "  \"seed\": " << cfg_->seed << ",\n";
  out << "  \"bio_dt_s\": " << cfg_->time.bio_dt << ",\n";
  out << "  \"steps\": " << step_count_ << ",\n";
  out << "  \"time_s\": " << time_ << ",\n";
  out << "  \"n_patches\": " << c.n_patches << ",\n";
  out << "  \"patch_volume_mL\": " << patch_volume_mL() << ",\n";
  out << "  \"contraction_rate_per_min\": " << c.contraction_rate_per_min
      << ",\n";
  out << "  \"agent_loss_fraction\": " << c.agent_loss_fraction << ",\n";
  out << "  \"transit_half_life_s\": " << c.transit_half_life_s << ",\n";
  out << "  \"reattach_prob_per_transit\": " << c.reattach_prob_per_transit
      << ",\n";
  out << "  \"establish_prob_single\": " << c.establish_prob_single << ",\n";
  out << "  \"establish_ratio\": " << c.establish_ratio << ",\n";
  write_live_physics_json(out, c);
  out << "  \"termination_cause\": \""
      << termination_cause_name(termination_cause_) << "\",\n";
  out << "  \"termination_detail\": \""
      << json_escape(termination_detail_) << "\",\n";
  out << "  \"ledger_closed\": " << (ledger_closed() ? "true" : "false")
      << ",\n";
  out << "  \"ledger\": {\"initial\": " << ledger_.initial_cells
      << ", \"births\": " << ledger_.births
      << ", \"washout\": " << ledger_.washout_departures
      << ", \"contraction\": " << ledger_.contraction_departures
      << ", \"edge\": " << ledger_.edge_departures
      << ", \"reseeds\": " << ledger_.reseeds
      << ", \"distal\": " << ledger_.distal_losses
      << ", \"external\": " << ledger_.external_arrivals
      << ", \"deaths\": " << ledger_.biological_deaths
      << ", \"biological_deaths\": " << ledger_.biological_deaths
      << ", \"sos_lysis_deaths\": " << ledger_.sos_lysis_deaths
      << ", \"toxin_kill_deaths\": " << ledger_.toxin_kill_deaths
      << ", \"pulse_washout\": " << ledger_.pulse_washout
      << ", \"pulse_deaths\": " << ledger_.pulse_deaths << "},\n";
  write_pulses_json(out, c.pulses, pulse_events_);
  if (live_.size() == 1) {
    const auto& [id, st] = *live_.begin();
    out << "  \"audit_live\": {\"patch\": " << id
        << ", \"live_seed\": " << st.live_seed << "},\n";
  }
  out << "  \"audit_live_patches\": [";
  {
    size_t i = 0;
    for (const auto& [id, st] : live_) {
      if (i++ > 0) out << ", ";
      out << "{\"patch\": " << id
          << ", \"live_seed\": " << st.live_seed << '}';
    }
  }
  out << "],\n";
  out << "  \"halted_patches\": [";
  for (size_t i = 0; i < halted_patches_.size(); ++i) {
    const auto& h = halted_patches_[i];
    if (i > 0) out << ", ";
    out << "{\"patch_id\": " << h.patch_id
        << ", \"status\": \"" << patch_status_name(h.status) << "\""
        << ", \"halt_step\": " << h.halt_step
        << ", \"halt_time_s\": " << h.halt_time_s << '}';
  }
  out << "],\n";
  out << "  \"occupancy\": {";
  for (Int t = 0; t < kPatchTypeCount; ++t) {
    const auto& to = o.per_type[t];
    if (t > 0) out << ", ";
    out << '\"' << patch_type_name(static_cast<PatchType>(t))
        << "\": {\"patches\": " << to.n_patches
        << ", \"occupied\": " << to.n_occupied
        << ", \"halted\": " << to.n_halted
        << ", \"occ\": " << to.occupancy
        << ", \"density_occupied_cfu_ml\": " << to.mean_density_occupied
        << '}';
  }
  out << "},\n";
  out << "  \"segment_mean_cfu_ml\": " << o.segment_mean_cfu_ml << ",\n";
  out << "  \"segment_cfu_cm2\": " << o.segment_cfu_cm2 << ",\n";
  out << "  \"audit_records\": [";
  for (size_t i = 0; i < audit_records_.size(); ++i) {
    const auto& r = audit_records_[i];
    if (i > 0) out << ", ";
    out << "{\"time\": " << r.time << ", \"patch\": " << r.patch_id
        << ", \"n_live\": " << r.n_live << ", \"n_shadow\": " << r.n_shadow
        << ", \"rel_diff\": " << r.rel_diff << '}';
  }
  out << "]\n}\n";
}

void MucusSegment::write_checkpoint_now() const {
  write_checkpoint(cfg_->layer2.checkpoint_file);
}

void MucusSegment::write_checkpoint(const std::string& path) const {
  write_checkpoint_group(path, "", true);
}

void MucusSegment::write_checkpoint_group(const std::string& path,
                                          const std::string& prefix,
                                          bool truncate) const {
#ifndef GUTIBM_HDF5
  (void)path;
  (void)prefix;
  (void)truncate;
  throw ConfigError("layer2 checkpoint requires a GUTIBM_HDF5 build");
#else
  if (!live_.empty()) {
    throw ConfigError(
        "layer2 checkpoint with a live audit slot is not supported in "
        "Phase 1 (embedded Simulation field state is not segment-"
        "serializable)");
  }
  validate_output_file_path(path);

  hid_t fid = -1;
  if (truncate) {
    fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                    H5P_DEFAULT);
  } else {
    fid = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  }
  if (fid < 0) {
    throw HDF5Error("cannot open layer2 checkpoint for write: " + path);
  }
  struct FileGuard {
    hid_t f;
    ~FileGuard() { H5Fclose(f); }
  } guard{fid};

  write_state(fid, prefix);
#endif  // GUTIBM_HDF5
}

#ifdef GUTIBM_HDF5
void MucusSegment::write_state(hid_t fid,
                               const std::string& prefix) const {
  const auto p = [&prefix](const char* leaf) {
    return prefix + leaf;
  };

  std::vector<std::vector<char>> patch_blobs(patches_.size());
  std::vector<int64_t> patch_counts(patches_.size());
  std::vector<int64_t> patch_ids(patches_.size());
  std::vector<int32_t> patch_types(patches_.size());
  std::vector<int32_t> patch_status(patches_.size());
  std::vector<int64_t> patch_baseline(patches_.size());
  std::vector<int64_t> patch_next_tag(patches_.size());
  std::vector<double> patch_bloom_since(patches_.size());
  std::vector<std::string> patch_rng(patches_.size());
  std::vector<int64_t> patch_offsets(patches_.size() + 1, 0);
  for (size_t i = 0; i < patches_.size(); ++i) {
    const Patch& p = patches_[i];
    agent_transfer_serialize(p.colonists, patch_blobs[i]);
    patch_counts[i] = static_cast<int64_t>(p.colonists.size());
    patch_ids[i] = p.id;
    patch_types[i] = static_cast<int32_t>(p.type);
    patch_status[i] = static_cast<int32_t>(p.status);
    patch_baseline[i] = p.baseline_count;
    patch_next_tag[i] = static_cast<int64_t>(p.next_cell_tag);
    patch_bloom_since[i] = static_cast<double>(p.bloom_excess_since);
    patch_rng[i] = rng_to_string(p.rng);
    patch_offsets[i + 1] =
        patch_offsets[i] + static_cast<int64_t>(patch_blobs[i].size());
  }
  std::vector<char> patch_blob;
  for (const auto& blob : patch_blobs) {
    patch_blob.insert(patch_blob.end(), blob.begin(), blob.end());
  }

  const size_t n_packets = pool_.size();
  std::vector<std::vector<char>> packet_blobs(n_packets);
  std::vector<int64_t> packet_counts(n_packets);
  std::vector<int64_t> packet_ids(n_packets);
  std::vector<int32_t> packet_kinds(n_packets);
  std::vector<double> packet_times(n_packets);
  std::vector<std::string> packet_rng(n_packets);
  std::vector<int64_t> packet_offsets(n_packets + 1, 0);
  std::vector<char> packet_blob;
  {
    size_t i = 0;
    for (const auto& packet : pool_) {
      agent_transfer_serialize(packet.cells, packet_blobs[i]);
      packet_counts[i] = static_cast<int64_t>(packet.cells.size());
      packet_ids[i] = static_cast<int64_t>(packet.id);
      packet_kinds[i] = static_cast<int32_t>(packet.kind);
      packet_times[i] = static_cast<double>(packet.entry_time);
      packet_rng[i] = rng_to_string(packet.rng);
      packet_offsets[i + 1] =
          packet_offsets[i] + static_cast<int64_t>(packet_blobs[i].size());
      packet_blob.insert(packet_blob.end(), packet_blobs[i].begin(),
                         packet_blobs[i].end());
      ++i;
    }
  }

  make_group(fid, p("/ledger"));
  make_group(fid, p("/patches"));
  make_group(fid, p("/pool"));
  make_group(fid, p("/audit"));
  make_group(fid, p("/streams"));
  make_group(fid, p("/mean_window"));
  make_group(fid, p("/pulses"));

  write_scalar_i64(fid, p("/format_version").c_str(), 1);
  write_scalar_f64(fid, p("/time").c_str(), static_cast<double>(time_));
  write_scalar_i64(fid, p("/step_count").c_str(),
                   static_cast<int64_t>(step_count_));
  write_scalar_i64(fid, p("/seed").c_str(),
                   static_cast<int64_t>(cfg_->seed));
  write_scalar_i64(fid, p("/next_tag").c_str(),
                   static_cast<int64_t>(next_tag_));
  write_scalar_i64(fid, p("/next_packet_id").c_str(),
                   static_cast<int64_t>(next_packet_id_));
  write_scalar_f64(fid, p("/audit_next_time").c_str(),
                   static_cast<double>(audit_next_time_));

  write_scalar_i64(fid, p("/ledger/initial_cells").c_str(),
                   ledger_.initial_cells);
  write_scalar_i64(fid, p("/ledger/births").c_str(), ledger_.births);
  write_scalar_i64(fid, p("/ledger/washout_departures").c_str(),
                   ledger_.washout_departures);
  write_scalar_i64(fid, p("/ledger/contraction_departures").c_str(),
                   ledger_.contraction_departures);
  write_scalar_i64(fid, p("/ledger/edge_departures").c_str(),
                   ledger_.edge_departures);
  write_scalar_i64(fid, p("/ledger/reseeds").c_str(), ledger_.reseeds);
  write_scalar_i64(fid, p("/ledger/distal_losses").c_str(),
                   ledger_.distal_losses);
  write_scalar_i64(fid, p("/ledger/external_arrivals").c_str(),
                   ledger_.external_arrivals);
  write_scalar_i64(fid, p("/ledger/biological_deaths").c_str(),
                   ledger_.biological_deaths);
  write_scalar_i64(fid, p("/ledger/sos_lysis_deaths").c_str(),
                   ledger_.sos_lysis_deaths);
  write_scalar_i64(fid, p("/ledger/toxin_kill_deaths").c_str(),
                   ledger_.toxin_kill_deaths);
  write_scalar_i64(fid, p("/ledger/pulse_washout").c_str(),
                   ledger_.pulse_washout);
  write_scalar_i64(fid, p("/ledger/pulse_deaths").c_str(),
                   ledger_.pulse_deaths);

  write_vec(fid, p("/patches/id").c_str(), H5T_NATIVE_INT64, patch_ids);
  write_vec(fid, p("/patches/type").c_str(), H5T_NATIVE_INT32,
            patch_types);
  write_vec(fid, p("/patches/status").c_str(), H5T_NATIVE_INT32,
            patch_status);
  write_vec(fid, p("/patches/baseline_count").c_str(), H5T_NATIVE_INT64,
            patch_baseline);
  write_vec(fid, p("/patches/next_cell_tag").c_str(), H5T_NATIVE_INT64,
            patch_next_tag);
  write_vec(fid, p("/patches/bloom_excess_since").c_str(),
            H5T_NATIVE_DOUBLE, patch_bloom_since);
  write_vec(fid, p("/patches/colonist_counts").c_str(), H5T_NATIVE_INT64,
            patch_counts);
  write_vec(fid, p("/patches/colonist_offsets").c_str(), H5T_NATIVE_INT64,
            patch_offsets);
  write_vec(fid, p("/patches/colonists").c_str(), H5T_STD_I8LE,
            patch_blob);
  write_str_vec(fid, p("/patches/rng_state").c_str(), patch_rng);

  write_vec(fid, p("/pool/packet_id").c_str(), H5T_NATIVE_INT64,
            packet_ids);
  write_vec(fid, p("/pool/kind").c_str(), H5T_NATIVE_INT32, packet_kinds);
  write_vec(fid, p("/pool/entry_time").c_str(), H5T_NATIVE_DOUBLE,
            packet_times);
  write_vec(fid, p("/pool/cell_counts").c_str(), H5T_NATIVE_INT64,
            packet_counts);
  write_vec(fid, p("/pool/blob_offsets").c_str(), H5T_NATIVE_INT64,
            packet_offsets);
  write_vec(fid, p("/pool/cells").c_str(), H5T_STD_I8LE, packet_blob);
  write_str_vec(fid, p("/pool/rng_state").c_str(), packet_rng);

  write_str_vec(fid, p("/streams/state").c_str(),
                {rng_to_string(contraction_stream_),
                 rng_to_string(init_stream_)});

  {
    std::vector<double> audit_time;
    std::vector<int64_t> audit_patch;
    std::vector<int64_t> audit_live;
    std::vector<int64_t> audit_shadow;
    std::vector<double> audit_rel;
    for (const auto& r : audit_records_) {
      audit_time.push_back(static_cast<double>(r.time));
      audit_patch.push_back(r.patch_id);
      audit_live.push_back(r.n_live);
      audit_shadow.push_back(r.n_shadow);
      audit_rel.push_back(static_cast<double>(r.rel_diff));
    }
    write_vec(fid, p("/audit/time").c_str(), H5T_NATIVE_DOUBLE,
              audit_time);
    write_vec(fid, p("/audit/patch").c_str(), H5T_NATIVE_INT64,
              audit_patch);
    write_vec(fid, p("/audit/n_live").c_str(), H5T_NATIVE_INT64,
              audit_live);
    write_vec(fid, p("/audit/n_shadow").c_str(), H5T_NATIVE_INT64,
              audit_shadow);
    write_vec(fid, p("/audit/rel_diff").c_str(), H5T_NATIVE_DOUBLE,
              audit_rel);
  }

  std::vector<double> win_time;
  std::vector<double> win_value;
  for (const auto& [t, v] : mean_window_) {
    win_time.push_back(static_cast<double>(t));
    win_value.push_back(static_cast<double>(v));
  }
  write_vec(fid, p("/mean_window/time").c_str(), H5T_NATIVE_DOUBLE,
            win_time);
  write_vec(fid, p("/mean_window/value").c_str(), H5T_NATIVE_DOUBLE,
            win_value);

  std::vector<int32_t> pulse_fired_i32;
  pulse_fired_i32.reserve(pulse_fired_.size());
  for (const bool fired : pulse_fired_) {
    pulse_fired_i32.push_back(fired ? 1 : 0);
  }
  write_vec(fid, p("/pulses/fired").c_str(), H5T_NATIVE_INT32,
            pulse_fired_i32);
}
#endif  // GUTIBM_HDF5

void MucusSegment::init_from_checkpoint(const SimulationConfig& cfg,
                                        const std::string& h5_file) {
  init_from_checkpoint_group(cfg, h5_file, "");
}

void MucusSegment::init_from_checkpoint_group(
    const SimulationConfig& cfg, const std::string& h5_file,
    const std::string& prefix) {
#ifndef GUTIBM_HDF5
  (void)cfg;
  (void)h5_file;
  (void)prefix;
  throw ConfigError("layer2 checkpoint requires a GUTIBM_HDF5 build");
#else
  validate_input_file_path(h5_file);
  cfg_ = &cfg;
  validate_config();
  table_ = std::make_unique<PatchTable>();

  hid_t fid = H5Fopen(h5_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fid < 0) {
    throw HDF5Error("cannot open layer2 checkpoint: " + h5_file);
  }
  struct FileGuard {
    hid_t f;
    ~FileGuard() { H5Fclose(f); }
  } guard{fid};

  read_state(fid, prefix);
#endif  // GUTIBM_HDF5
}

#ifdef GUTIBM_HDF5
void MucusSegment::read_state(hid_t fid, const std::string& prefix) {
  const auto gp = [&prefix](const char* leaf) {
    return prefix + leaf;
  };

  const int64_t version =
      read_scalar_i64(fid, gp("/format_version").c_str());
  if (version != 1) {
    throw HDF5Error("unsupported layer2 checkpoint format_version " +
                    std::to_string(version));
  }
  if (read_scalar_i64(fid, gp("/seed").c_str()) !=
      static_cast<int64_t>(cfg_->seed)) {
    throw ConfigError(
        "layer2 checkpoint seed does not match configured seed");
  }
  time_ = static_cast<Real>(read_scalar_f64(fid, gp("/time").c_str()));
  step_count_ =
      static_cast<Int>(read_scalar_i64(fid, gp("/step_count").c_str()));
  next_tag_ =
      static_cast<TagID>(read_scalar_i64(fid, gp("/next_tag").c_str()));
  next_packet_id_ = static_cast<uint64_t>(
      read_scalar_i64(fid, gp("/next_packet_id").c_str()));
  audit_next_time_ = static_cast<Real>(
      read_scalar_f64(fid, gp("/audit_next_time").c_str()));

  ledger_.initial_cells =
      read_scalar_i64(fid, gp("/ledger/initial_cells").c_str());
  ledger_.births = read_scalar_i64(fid, gp("/ledger/births").c_str());
  ledger_.washout_departures =
      read_scalar_i64(fid, gp("/ledger/washout_departures").c_str());
  ledger_.contraction_departures =
      read_scalar_i64(fid, gp("/ledger/contraction_departures").c_str());
  ledger_.edge_departures =
      read_scalar_i64_opt(fid, gp("/ledger/edge_departures"), 0);
  ledger_.reseeds = read_scalar_i64(fid, gp("/ledger/reseeds").c_str());
  ledger_.distal_losses =
      read_scalar_i64(fid, gp("/ledger/distal_losses").c_str());
  ledger_.external_arrivals =
      read_scalar_i64_opt(fid, gp("/ledger/external_arrivals"), 0);
  ledger_.biological_deaths =
      read_scalar_i64_opt(fid, gp("/ledger/biological_deaths"), 0);
  ledger_.sos_lysis_deaths =
      read_scalar_i64_opt(fid, gp("/ledger/sos_lysis_deaths"), 0);
  ledger_.toxin_kill_deaths =
      read_scalar_i64_opt(fid, gp("/ledger/toxin_kill_deaths"), 0);
  ledger_.pulse_washout =
      read_scalar_i64_opt(fid, gp("/ledger/pulse_washout"), 0);
  ledger_.pulse_deaths =
      read_scalar_i64_opt(fid, gp("/ledger/pulse_deaths"), 0);

  const auto patch_ids = read_vec<int64_t>(fid, gp("/patches/id").c_str(),
                                           H5T_NATIVE_INT64);
  const auto patch_types = read_vec<int32_t>(
      fid, gp("/patches/type").c_str(), H5T_NATIVE_INT32);
  const auto patch_status = read_vec<int32_t>(
      fid, gp("/patches/status").c_str(), H5T_NATIVE_INT32);
  const auto patch_baseline = read_vec<int64_t>(
      fid, gp("/patches/baseline_count").c_str(), H5T_NATIVE_INT64);
  const auto patch_next_tag = read_vec<int64_t>(
      fid, gp("/patches/next_cell_tag").c_str(), H5T_NATIVE_INT64);
  const auto patch_bloom = read_vec<double>(
      fid, gp("/patches/bloom_excess_since").c_str(), H5T_NATIVE_DOUBLE);
  const auto patch_counts = read_vec<int64_t>(
      fid, gp("/patches/colonist_counts").c_str(), H5T_NATIVE_INT64);
  const auto patch_offsets = read_vec<int64_t>(
      fid, gp("/patches/colonist_offsets").c_str(), H5T_NATIVE_INT64);
  const auto patch_blob = read_vec<char>(
      fid, gp("/patches/colonists").c_str(), H5T_STD_I8LE);
  const auto patch_rng = read_str_vec(fid, gp("/patches/rng_state").c_str());

  const auto live_indices = cfg_->layer2.audit.resolved_patch_indices();
  patches_.clear();
  patches_.resize(patch_ids.size());
  for (size_t i = 0; i < patch_ids.size(); ++i) {
    Patch& p = patches_[i];
    p.id = static_cast<Int>(patch_ids[i]);
    p.type = static_cast<PatchType>(patch_types[i]);
    p.status = static_cast<PatchStatus>(patch_status[i]);
    p.baseline_count = patch_baseline[i];
    p.next_cell_tag = static_cast<TagID>(patch_next_tag[i]);
    p.bloom_excess_since = static_cast<Real>(patch_bloom[i]);
    rng_from_string(patch_rng[i].c_str(), p.rng);
    const auto lo = static_cast<size_t>(patch_offsets[i]);
    const auto hi = static_cast<size_t>(patch_offsets[i + 1]);
    if (hi > lo) {
      const std::vector<char> slice(patch_blob.begin() + lo,
                                    patch_blob.begin() + hi);
      p.colonists = agent_transfer_deserialize(slice);
      if (p.colonists.size() !=
          static_cast<size_t>(patch_counts[i])) {
        throw HDF5Error("layer2 checkpoint patch colonist count mismatch");
      }
    }
    if (p.live || std::find(live_indices.begin(), live_indices.end(), p.id) !=
                      live_indices.end()) {
      throw ConfigError(
          "layer2 checkpoint with a live audit slot is not supported in "
          "Phase 1");
    }
  }
  order_.resize(patches_.size());
  std::iota(order_.begin(), order_.end(), 0);

  const auto packet_ids = read_vec<int64_t>(
      fid, gp("/pool/packet_id").c_str(), H5T_NATIVE_INT64);
  const auto packet_kinds = read_vec<int32_t>(
      fid, gp("/pool/kind").c_str(), H5T_NATIVE_INT32);
  const auto packet_times = read_vec<double>(
      fid, gp("/pool/entry_time").c_str(), H5T_NATIVE_DOUBLE);
  const auto packet_counts = read_vec<int64_t>(
      fid, gp("/pool/cell_counts").c_str(), H5T_NATIVE_INT64);
  const auto packet_offsets = read_vec<int64_t>(
      fid, gp("/pool/blob_offsets").c_str(), H5T_NATIVE_INT64);
  const auto packet_blob = read_vec<char>(fid, gp("/pool/cells").c_str(),
                                        H5T_STD_I8LE);
  const auto packet_rng = read_str_vec(fid, gp("/pool/rng_state").c_str());

  pool_.clear();
  for (size_t i = 0; i < packet_ids.size(); ++i) {
    LuminalPacket packet;
    packet.id = static_cast<uint64_t>(packet_ids[i]);
    packet.kind = static_cast<PacketKind>(packet_kinds[i]);
    packet.entry_time = static_cast<Real>(packet_times[i]);
    rng_from_string(packet_rng[i].c_str(), packet.rng);
    const auto lo = static_cast<size_t>(packet_offsets[i]);
    const auto hi = static_cast<size_t>(packet_offsets[i + 1]);
    if (hi > lo) {
      const std::vector<char> slice(packet_blob.begin() + lo,
                                    packet_blob.begin() + hi);
      packet.cells = agent_transfer_deserialize(slice);
      if (packet.cells.size() !=
          static_cast<size_t>(packet_counts[i])) {
        throw HDF5Error("layer2 checkpoint packet cell count mismatch");
      }
    }
    pool_.push_back(std::move(packet));
  }

  const auto streams = read_str_vec(fid, gp("/streams/state").c_str());
  if (streams.size() == 2) {
    rng_from_string(streams[0].c_str(), contraction_stream_);
    rng_from_string(streams[1].c_str(), init_stream_);
  }

  const auto audit_time = read_vec<double>(
      fid, gp("/audit/time").c_str(), H5T_NATIVE_DOUBLE);
  const auto audit_patch = read_vec<int64_t>(
      fid, gp("/audit/patch").c_str(), H5T_NATIVE_INT64);
  const auto audit_live = read_vec<int64_t>(
      fid, gp("/audit/n_live").c_str(), H5T_NATIVE_INT64);
  const auto audit_shadow = read_vec<int64_t>(
      fid, gp("/audit/n_shadow").c_str(), H5T_NATIVE_INT64);
  const auto audit_rel = read_vec<double>(
      fid, gp("/audit/rel_diff").c_str(), H5T_NATIVE_DOUBLE);
  for (size_t i = 0; i < audit_time.size(); ++i) {
    audit_records_.push_back(
        {static_cast<Real>(audit_time[i]),
         static_cast<Int>(audit_patch[i]),
         static_cast<Int>(audit_live[i]),
         static_cast<Int>(audit_shadow[i]),
         static_cast<Real>(audit_rel[i])});
  }

  const auto win_time = read_vec<double>(
      fid, gp("/mean_window/time").c_str(), H5T_NATIVE_DOUBLE);
  const auto win_value = read_vec<double>(
      fid, gp("/mean_window/value").c_str(), H5T_NATIVE_DOUBLE);
  for (size_t i = 0; i < win_time.size(); ++i) {
    mean_window_.emplace_back(static_cast<Real>(win_time[i]),
                              static_cast<Real>(win_value[i]));
  }

  pulse_fired_.assign(cfg_->layer2.pulses.size(), false);
  pulse_events_.clear();
  if (H5Lexists(fid, gp("/pulses").c_str(), H5P_DEFAULT) > 0 &&
      H5Lexists(fid, gp("/pulses/fired").c_str(), H5P_DEFAULT) > 0) {
    const auto fired_i32 = read_vec<int32_t>(
        fid, gp("/pulses/fired").c_str(), H5T_NATIVE_INT32);
    const size_t n_copy = std::min(fired_i32.size(), pulse_fired_.size());
    for (size_t i = 0; i < n_copy; ++i) {
      pulse_fired_[i] = (fired_i32[i] != 0);
    }
  }

  if (!ledger_closed()) {
    throw HDF5Error(
        "layer2 checkpoint restored a ledger that does not close");
  }
}
#endif  // GUTIBM_HDF5

}  // namespace gutibm
