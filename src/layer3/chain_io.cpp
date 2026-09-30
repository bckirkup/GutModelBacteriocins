/* -----------------------------------------------------------------------
   GutIBM – Spec 13 Phase 3: chain run loop, timeseries, provenance,
   and checkpoint serialization.

   Chain checkpoint layout:
     /format_version            i64 (= 3; layer-3 chain marker)
     /seed                      i64 (the chain-level seed)
     /time, /step_count
     /ledger/{initial_cells, mucosal_births, luminal_births,
              stool_exports, luminal_deaths, purged}
     /purged_flag               i64
     /regions/<r>/deposit_tag   i64 (chain-minted deposit serial)
     /regions/<r>/segment/...   MucusSegment checkpoint tree
     /regions/<r>/lumen/{carbon_mol_m3, acetate_mol_m3, n_lineages,
                        lineage_ids, counts, exemplar_offsets,
                        exemplars}
   ----------------------------------------------------------------------- */

#include "agent_transfer.h"
#include "chain.h"
#include "error.h"
#include "h5_util.h"
#include "input_parser.h"
#include "path_utils.h"
#include "segment.h"

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

#include <cstring>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace gutibm {

using namespace h5;

namespace {

const char* kTimeseriesHeader =
    "step,time_s,region,occ_frac,patch_cells,pool_cells,"
    "births_cum,washout_cum,contraction_cum,edge_cum,reseeds_cum,"
    "distal_cum,external_cum,mean_cfu_ml,"
    "lumen_cells,lumen_lineages,lumen_mu_per_s,"
    "lumen_carbon_mol_m3,lumen_acetate_mol_m3,"
    "stool_cum,deaths_cum,purged_cum,stocks,expected,ledger_ok,hapc";

std::string json_quote(const std::string& s) {
  std::string out;
  out.reserve(s.size() + 2);
  out.push_back('"');
  for (const char ch : s) {
    switch (ch) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\t': out += "\\t"; break;
      default: out += ch; break;
    }
  }
  out.push_back('"');
  return out;
}

}  // namespace

void ColonicChain::emit_timeseries_row(std::ofstream& out) const {
  const Real stocks_now = stocks();
  const Real expected = ledger_.expected_stocks();
  const bool ok = ledger_closed();
  for (Int r = 0; r < kRegionCount; ++r) {
    const RegionState& rs = regions_[static_cast<size_t>(r)];
    const SegmentObservables obs = rs.segment->observables();
    const Real occ =
        static_cast<Real>(obs.n_occupied_total) /
        std::max(1.0, static_cast<Real>(cfg_->layer3.n_patches_per_region));
    out << step_count_ << ',' << time_ << ',' << region_name(
               static_cast<RegionKind>(r))
        << ',' << occ << ',' << obs.patch_cells << ',' << obs.pool_cells
        << ',' << rs.segment->ledger().births << ','
        << rs.segment->ledger().washout_departures << ','
        << rs.segment->ledger().contraction_departures << ','
        << rs.segment->ledger().edge_departures << ','
        << rs.segment->ledger().reseeds << ','
        << rs.segment->ledger().distal_losses << ','
        << rs.segment->ledger().external_arrivals << ','
        << obs.segment_mean_cfu_ml << ',' << rs.lumen.total_cells() << ','
        << rs.lumen.lineages().size() << ','
        << last_mean_mu_[static_cast<size_t>(r)] << ','
        << rs.lumen.scalars().carbon_mol_m3 << ','
        << rs.lumen.scalars().acetate_mol_m3 << ','
        << ledger_.stool_exports << ',' << ledger_.luminal_deaths << ','
        << ledger_.purged << ',' << stocks_now << ',' << expected << ','
        << (ok ? 1 : 0) << ','
        << (last_hapc_.fired &&
                    last_hapc_.regions[static_cast<size_t>(r)]
                ? 1
                : 0)
        << '\n';
  }
}

int ColonicChain::run() {
  const auto& c = cfg_->layer3;
  const Real dt = cfg_->time.bio_dt;

  std::ofstream ts;
  if (!c.timeseries_file.empty()) {
    validate_output_file_path(c.timeseries_file);
    ts.open(c.timeseries_file);
    if (!ts) {
      throw IOError("cannot open layer3 timeseries file: " +
                    c.timeseries_file);
    }
    ts << kTimeseriesHeader << '\n';
  }

  while (time_ < cfg_->time.total_time &&
         termination_cause_ == TerminationCause::IncompleteUnknown) {
    step(dt);
    if (termination_cause_ != TerminationCause::IncompleteUnknown) break;
    if (time_ >= cfg_->time.total_time) {
      terminate(TerminationCause::HorizonReached, "reached configured horizon");
      break;
    }
    if (ts.is_open() && c.summary_interval_steps > 0 &&
        step_count_ % c.summary_interval_steps == 0) {
      emit_timeseries_row(ts);
    }
    if (c.checkpoint_interval_steps > 0 && !c.checkpoint_file.empty() &&
        step_count_ % c.checkpoint_interval_steps == 0) {
      write_checkpoint(c.checkpoint_file);
    }
  }
  if (ts.is_open()) {
    ts.flush();
  }
  if (c.checkpoint_final && !c.checkpoint_file.empty()) {
    write_checkpoint(c.checkpoint_file);
  }
  if (!c.provenance_file.empty()) {
    write_provenance(c.provenance_file);
  }

  std::cout << std::format(
      "layer3 chain terminated: cause={} steps={} stocks={:.6e} "
      "expected={:.6e} stool_cum={:.6e}\n",
      static_cast<int>(termination_cause_), step_count_, stocks(),
      ledger_.expected_stocks(), ledger_.stool_exports);
  if (!termination_detail_.empty()) {
    std::cout << "  detail: " << termination_detail_ << '\n';
  }
  return static_cast<int>(termination_cause_);
}

void ColonicChain::write_provenance(const std::string& path) const {
  validate_output_file_path(path);
  std::ofstream out(path);
  if (!out) {
    throw IOError("cannot open layer3 provenance file: " + path);
  }
  const auto& c = cfg_->layer3;
  out << "{\n";
  out << "  \"layer3\": {\n";
  out << "    \"enabled\": " << (c.enabled ? "true" : "false") << ",\n";
  out << "    \"uniform_profile\": "
      << (c.uniform_profile ? "true" : "false") << ",\n";
  out << "    \"n_patches_per_region\": " << c.n_patches_per_region
      << ",\n";
  out << "    \"f_edge\": " << c.f_edge << ",\n";
  out << "    \"reattach_p0\": " << c.reattach_p0 << ",\n";
  out << "    \"hapc_rate_per_day\": " << c.hapc_rate_per_day << ",\n";
  out << "    \"hapc_antegrade_fraction\": " << c.hapc_antegrade_fraction
      << ",\n";
  out << "    \"hapc_day_fraction\": " << c.hapc_day_fraction << ",\n";
  out << "    \"alpha_hapc_per_min\": " << c.alpha_hapc_per_min << ",\n";
  out << "    \"stool_g_per_day\": " << c.stool_g_per_day << ",\n";
  out << "    \"flatness_bound\": " << c.flatness_bound << ",\n";
  out << "    \"lumen_seed_cells\": " << c.lumen_seed_cells << ",\n";
  out << "    \"purge_at_s\": " << c.purge_at_s << ",\n";
  out << "    \"purge_fraction\": " << c.purge_fraction << ",\n";
  out << "    \"regions\": [\n";
  for (Int r = 0; r < kRegionCount; ++r) {
    const Layer3RegionConfig& rc = regions_cfg_[static_cast<size_t>(r)];
    out << "      {\"name\": " << json_quote(
               region_name(static_cast<RegionKind>(r)))
        << ", \"transit_h\": " << rc.transit_h
        << ", \"ph\": " << rc.ph
        << ", \"lumen_carbon_mol_m3\": " << rc.lumen_carbon_mol_m3
        << ", \"lumen_volume_L\": " << rc.lumen_volume_L
        << ", \"contraction_k0_per_min\": " << rc.contraction_k0_per_min
        << ", \"outflow_survival\": " << rc.outflow_survival
        << ", \"flora_density_cfu_ml\": " << rc.flora_density_cfu_ml
        << ", \"patch_supply_scale\": " << rc.patch_supply_scale << "}"
        << (r + 1 < kRegionCount ? "," : "") << '\n';
  }
  out << "    ]\n  },\n";
  out << "  \"seed\": " << cfg_->seed << ",\n";
  out << "  \"dt_s\": " << cfg_->time.bio_dt << ",\n";
  out << "  \"steps_completed\": " << step_count_ << ",\n";
  out << "  \"wall_time_s\": " << time_ << ",\n";
  out << "  \"termination\": {\n"
      << "    \"cause\": " << static_cast<int>(termination_cause_)
      << ",\n    \"detail\": " << json_quote(termination_detail_)
      << "\n  },\n";
  out << "  \"ledger\": {\n"
      << "    \"initial_cells\": " << ledger_.initial_cells << ",\n"
      << "    \"mucosal_births\": " << ledger_.mucosal_births << ",\n"
      << "    \"luminal_births\": " << ledger_.luminal_births << ",\n"
      << "    \"stool_exports\": " << ledger_.stool_exports << ",\n"
      << "    \"luminal_deaths\": " << ledger_.luminal_deaths << ",\n"
      << "    \"purged\": " << ledger_.purged << ",\n"
      << "    \"expected_stocks\": " << ledger_.expected_stocks()
      << ",\n    \"stocks\": " << stocks() << ",\n"
      << "    \"closed\": " << (ledger_closed() ? "true" : "false")
      << "\n  },\n";
  out << "  \"observables\": {\n    \"stool_cfu_per_g\": "
      << stool_cfu_per_g() << ",\n    \"regions\": [\n";
  for (Int r = 0; r < kRegionCount; ++r) {
    const RegionState& rs = regions_[static_cast<size_t>(r)];
    const SegmentObservables obs = rs.segment->observables();
    const Real fraction =
        obs.segment_mean_cfu_ml / regions_cfg_[static_cast<size_t>(r)]
                              .flora_density_cfu_ml;
    out << "      {\"name\": "
        << json_quote(region_name(static_cast<RegionKind>(r)))
        << ", \"occupancy_frac\": "
        << static_cast<Real>(obs.n_occupied_total) /
               std::max(1.0, static_cast<Real>(c.n_patches_per_region))
        << ", \"patch_cells\": " << obs.patch_cells
        << ", \"mean_cfu_ml\": " << obs.segment_mean_cfu_ml
        << ", \"enterobacteriaceae_fraction\": " << fraction
        << ", \"lumen_cells\": " << rs.lumen.total_cells()
        << ", \"lumen_lineages\": " << rs.lumen.lineages().size()
        << ", \"lumen_mean_mu_per_s\": "
        << last_mean_mu_[static_cast<size_t>(r)]
        << ", \"lumen_carbon_mol_m3\": "
        << rs.lumen.scalars().carbon_mol_m3
        << ", \"lumen_acetate_mol_m3\": "
        << rs.lumen.scalars().acetate_mol_m3 << "}"
        << (r + 1 < kRegionCount ? "," : "") << '\n';
  }
  out << "    ]\n  }\n}\n";
}

// ── Checkpoint serialization ────────────────────────────────────────────

#ifdef GUTIBM_HDF5

void LumenCompartment::write_state(hid_t fid,
                                   const std::string& prefix) const {
  auto p = [&](const char* leaf) { return prefix + "/" + leaf; };
  make_group(fid, prefix);
  write_scalar_f64(fid, p("carbon_mol_m3").c_str(),
                   scalars_.carbon_mol_m3);
  write_scalar_f64(fid, p("acetate_mol_m3").c_str(),
                   scalars_.acetate_mol_m3);
  write_scalar_i64(fid, p("n_lineages").c_str(),
                   static_cast<int64_t>(lineages_.size()));

  std::vector<int64_t> ids;
  std::vector<double> counts;
  std::vector<int64_t> offsets;
  std::vector<char> blobs;
  for (const auto& [id, lin] : lineages_) {
    ids.push_back(static_cast<int64_t>(id));
    counts.push_back(lin.count);
    offsets.push_back(static_cast<int64_t>(blobs.size()));
    std::vector<char> blob;
    agent_transfer_serialize({lin.exemplar}, blob);
    blobs.insert(blobs.end(), blob.begin(), blob.end());
  }
  write_vec(fid, p("lineage_ids").c_str(), H5T_NATIVE_INT64, ids);
  write_vec(fid, p("counts").c_str(), H5T_NATIVE_DOUBLE, counts);
  write_vec(fid, p("exemplar_offsets").c_str(), H5T_NATIVE_INT64,
            offsets);
  write_vec(fid, p("exemplars").c_str(), H5T_NATIVE_CHAR, blobs);
}

void LumenCompartment::read_state(hid_t fid, const std::string& prefix) {
  auto gp = [&](const char* leaf) { return prefix + "/" + leaf; };
  scalars_.carbon_mol_m3 =
      read_scalar_f64(fid, gp("carbon_mol_m3").c_str());
  scalars_.acetate_mol_m3 =
      read_scalar_f64(fid, gp("acetate_mol_m3").c_str());
  const int64_t n = read_scalar_i64(fid, gp("n_lineages").c_str());
  std::vector<int64_t> ids =
      read_vec<int64_t>(fid, gp("lineage_ids").c_str(), H5T_NATIVE_INT64);
  std::vector<double> counts =
      read_vec<double>(fid, gp("counts").c_str(), H5T_NATIVE_DOUBLE);
  std::vector<int64_t> offsets = read_vec<int64_t>(
      fid, gp("exemplar_offsets").c_str(), H5T_NATIVE_INT64);
  std::vector<char> blobs =
      read_vec<char>(fid, gp("exemplars").c_str(), H5T_NATIVE_CHAR);
  if (static_cast<int64_t>(ids.size()) != n ||
      static_cast<int64_t>(counts.size()) != n ||
      static_cast<int64_t>(offsets.size()) != n) {
    throw HDF5Error("layer3 checkpoint lineage arrays inconsistent");
  }
  lineages_.clear();
  for (int64_t i = 0; i < n; ++i) {
    const size_t start = static_cast<size_t>(offsets[i]);
    const size_t end = (i + 1 < n)
                           ? static_cast<size_t>(offsets[i + 1])
                           : blobs.size();
    LuminalLineage lin;
    std::vector<char> slice(blobs.begin() + start,
                            blobs.begin() + end);
    lin.exemplar = agent_transfer_deserialize(slice).at(0);
    lin.count = counts[i];
    lineages_[static_cast<TagID>(ids[i])] = std::move(lin);
  }
}

void ColonicChain::write_checkpoint(const std::string& path) const {
  validate_output_file_path(path);
  // Segments own their subtrees; the first call creates the file.
  for (Int r = 0; r < kRegionCount; ++r) {
    regions_[static_cast<size_t>(r)].segment->write_checkpoint_group(
        path, "/regions/" + std::to_string(r) + "/segment", r == 0);
  }
  FileGuard fg{H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT)};
  if (fg.f < 0) {
    throw HDF5Error("cannot reopen layer3 checkpoint: " + path);
  }
  write_scalar_i64(fg.f, "/format_version", 3);
  write_scalar_i64(fg.f, "/seed", static_cast<int64_t>(cfg_->seed));
  write_scalar_f64(fg.f, "/time", time_);
  write_scalar_i64(fg.f, "/step_count", step_count_);
  write_scalar_i64(fg.f, "/purged_flag", purged_ ? 1 : 0);
  make_group(fg.f, "/ledger");
  write_scalar_f64(fg.f, "/ledger/initial_cells", ledger_.initial_cells);
  write_scalar_f64(fg.f, "/ledger/mucosal_births", ledger_.mucosal_births);
  write_scalar_f64(fg.f, "/ledger/luminal_births", ledger_.luminal_births);
  write_scalar_f64(fg.f, "/ledger/stool_exports", ledger_.stool_exports);
  write_scalar_f64(fg.f, "/ledger/luminal_deaths",
                   ledger_.luminal_deaths);
  write_scalar_f64(fg.f, "/ledger/purged", ledger_.purged);
  for (Int r = 0; r < kRegionCount; ++r) {
    const RegionState& rs = regions_[static_cast<size_t>(r)];
    const std::string base = "/regions/" + std::to_string(r);
    write_scalar_i64(fg.f, (base + "/deposit_tag").c_str(),
                     static_cast<int64_t>(rs.next_deposit_tag));
    rs.lumen.write_state(fg.f, base + "/lumen");
  }
}

void ColonicChain::init_from_checkpoint(const SimulationConfig& cfg,
                                        const std::string& h5_file) {
  cfg_ = &cfg;
  const Layer3Config& c = cfg.layer3;
  regions_cfg_.assign({c.cecum, c.transverse, c.descending});
  if (c.uniform_profile) {
    regions_cfg_.assign(kRegionCount, c.transverse);
  }
  validate_config();

  physio_.yield_carbon = cfg.fixes.metabolism.yield_carbon;
  physio_.carbon_cost_factor =
      cfg.chem_env.oxygen.anaerobic_carbon_cost_factor;
  physio_.anaerobic_mu_factor = cfg.chem_env.oxygen.anaerobic_mu_factor;
  physio_.ferm_acid_yield = cfg.chem_env.oxygen.ferm_acid_yield;
  physio_.acid_inhibition_max = cfg.fixes.metabolism.acid_inhibition_max;
  physio_.acid_inhibition_ki = cfg.fixes.metabolism.acid_inhibition_Ki;
  physio_.acetate_pka = cfg.fixes.metabolism.acetate_pKa;
  physio_.diet_carbon_mol_m3 =
      regions_cfg_.front().lumen_carbon_mol_m3;

  regions_.clear();
  regions_.resize(kRegionCount);
  region_cfgs_.clear();
  order_.resize(kRegionCount);
  std::iota(order_.begin(), order_.end(), 0);

  FileGuard fg{H5Fopen(h5_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)};
  if (fg.f < 0) {
    throw HDF5Error("cannot open layer3 checkpoint: " + h5_file);
  }
  if (read_scalar_i64(fg.f, "/format_version") != 3) {
    throw HDF5Error(
        "layer3 checkpoint format_version mismatch (expected 3)");
  }
  if (read_scalar_i64(fg.f, "/seed") !=
      static_cast<int64_t>(cfg.seed)) {
    throw HDF5Error("layer3 checkpoint seed does not match config");
  }

  for (Int r = 0; r < kRegionCount; ++r) {
    const auto i = static_cast<size_t>(r);
    RegionState& rs = regions_[i];
    const std::string base = "/regions/" + std::to_string(r);

    SimulationConfig& rcfg = region_cfgs_.emplace_back(cfg);
    rcfg.seed = mix_seed(cfg.seed, static_cast<uint64_t>(0x5EED),
                         static_cast<uint64_t>(r));
    rcfg.layer2.tag_prefix = static_cast<TagID>(r + 1) << 56;
    rcfg.layer2.n_patches = c.n_patches_per_region;
    rcfg.layer2.edge_shed_fraction = c.f_edge;
    rcfg.layer2.crypt.supply_mult *=
        regions_cfg_[i].patch_supply_scale;
    rcfg.layer2.exposed_proximal.supply_mult *=
        regions_cfg_[i].patch_supply_scale;
    rcfg.layer2.exposed_distal.supply_mult *=
        regions_cfg_[i].patch_supply_scale;
    rcfg.layer2.timeseries_file.clear();
    rcfg.layer2.provenance_file.clear();
    rcfg.layer2.checkpoint_file.clear();
    rcfg.layer2.checkpoint_final = false;

    rs.segment = std::make_unique<MucusSegment>();
    rs.segment->init_from_checkpoint_group(rcfg, h5_file,
                                           base + "/segment");
    rs.next_deposit_tag = static_cast<TagID>(
        read_scalar_i64(fg.f, (base + "/deposit_tag").c_str()));
    rs.lumen.read_state(fg.f, base + "/lumen");
    rs.segment->set_pool_exit_sink(
        [this, i](std::vector<Agent>&& cells) {
          regions_[i].lumen.add_cells(cells);
        });
    rs.segment->set_contraction_rate_per_min(
        regions_cfg_[i].contraction_k0_per_min);
  }

  ledger_.initial_cells = read_scalar_f64(fg.f, "/ledger/initial_cells");
  ledger_.mucosal_births =
      read_scalar_f64(fg.f, "/ledger/mucosal_births");
  ledger_.luminal_births =
      read_scalar_f64(fg.f, "/ledger/luminal_births");
  ledger_.stool_exports = read_scalar_f64(fg.f, "/ledger/stool_exports");
  ledger_.luminal_deaths =
      read_scalar_f64(fg.f, "/ledger/luminal_deaths");
  ledger_.purged = read_scalar_f64(fg.f, "/ledger/purged");
  purged_ = read_scalar_i64(fg.f, "/purged_flag") != 0;
  time_ = read_scalar_f64(fg.f, "/time");
  step_count_ = read_scalar_i64(fg.f, "/step_count");
  last_hapc_ = {};
  termination_cause_ = TerminationCause::IncompleteUnknown;
  termination_detail_.clear();
}

#else  // !GUTIBM_HDF5

void ColonicChain::write_checkpoint(const std::string&) const {
  throw HDF5Error("gut_ibm was built without HDF5 support");
}

void ColonicChain::init_from_checkpoint(const SimulationConfig&,
                                        const std::string&) {
  throw HDF5Error("gut_ibm was built without HDF5 support");
}

#endif  // GUTIBM_HDF5

}  // namespace gutibm
