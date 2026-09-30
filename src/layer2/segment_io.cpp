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
#include <vector>

namespace gutibm {

using namespace h5;

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
    "edge_cum,external_cum";

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

}  // namespace

int MucusSegment::run() {
  const auto& c = cfg_->layer2;
  const Real dt = cfg_->time.bio_dt;

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
  }

  while (time_ < cfg_->time.total_time &&
         termination_cause_ == TerminationCause::IncompleteUnknown) {
    step(dt);
    if (ts.is_open() && c.summary_interval_steps > 0 &&
        step_count_ % c.summary_interval_steps == 0) {
      emit_timeseries_row(ts);
    }
    if (!c.checkpoint_file.empty() && c.checkpoint_interval_steps > 0 &&
        step_count_ % c.checkpoint_interval_steps == 0) {
      write_checkpoint_now();
    }
  }
  if (ts.is_open()) {
    emit_timeseries_row(ts);
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
      << '\n';
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
      << ", \"external\": " << ledger_.external_arrivals << "},\n";
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
    if (p.live || cfg_->layer2.audit.patch_index == p.id) {
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

  if (!ledger_closed()) {
    throw HDF5Error(
        "layer2 checkpoint restored a ledger that does not close");
  }
}
#endif  // GUTIBM_HDF5

}  // namespace gutibm
