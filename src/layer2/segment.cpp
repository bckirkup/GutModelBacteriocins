/* -----------------------------------------------------------------------
   GutIBM – Spec 13 Phase 1: Layer-2 mucus segment driver.
   See segment.h for the design contract mapping (S1/S2/S5/G2-G6).
   ----------------------------------------------------------------------- */

#include "segment.h"

#include "agent_transfer.h"
#include "error.h"
#include "input_parser.h"
#include "patch_table.h"
#include "plasmid.h"
#include "simulation.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <numeric>

namespace gutibm {

namespace {

constexpr uint64_t kGolden = 0x9e3779b97f4a7c15ULL;

uint64_t splitmix64(uint64_t x) {
  x += kGolden;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

uint64_t tag_hash(const char* tag) {
  uint64_t h = 0xcbf29ce484222325ULL;  // FNV-1a
  for (const char* c = tag; *c != '\0'; ++c) {
    h = (h ^ static_cast<unsigned char>(*c)) * 0x100000001b3ULL;
  }
  return h;
}

uint64_t hash_combine(uint64_t h, uint64_t v) {
  return splitmix64(h ^ v);
}

constexpr Real kLn2 = 0.6931471805599453;

}  // namespace

const char* patch_type_name(PatchType t) {
  switch (t) {
    case PatchType::Crypt: return "crypt";
    case PatchType::ExposedProximal: return "exposed_proximal";
    case PatchType::ExposedDistal: return "exposed_distal";
  }
  return "unknown";
}

const char* patch_status_name(PatchStatus s) {
  switch (s) {
    case PatchStatus::Active: return "active";
    case PatchStatus::BloomHalted: return "bloom_halted";
    case PatchStatus::SpatialHalted: return "spatial_halted";
  }
  return "unknown";
}

const char* packet_kind_name(PacketKind k) {
  switch (k) {
    case PacketKind::Single: return "single";
    case PacketKind::Fragment: return "fragment";
  }
  return "unknown";
}

// ── Seed mixing (G3) ─────────────────────────────────────────────────────

uint64_t MucusSegment::mix_seed(uint64_t seed, uint64_t a, uint64_t b) {
  return splitmix64(splitmix64(seed ^ splitmix64(a)) ^ splitmix64(b));
}

uint64_t MucusSegment::mix_seed(uint64_t seed, uint64_t a) {
  return splitmix64(seed ^ splitmix64(a));
}

uint64_t MucusSegment::mix_seed(uint64_t seed, const char* tag) {
  return splitmix64(seed ^ tag_hash(tag));
}

// ── Lifecycle ────────────────────────────────────────────────────────────

MucusSegment::MucusSegment() = default;
MucusSegment::~MucusSegment() = default;

const Layer2TypeConfig& MucusSegment::type_config(PatchType t) const {
  switch (t) {
    case PatchType::Crypt: return cfg_->layer2.crypt;
    case PatchType::ExposedProximal: return cfg_->layer2.exposed_proximal;
    case PatchType::ExposedDistal: return cfg_->layer2.exposed_distal;
  }
  return cfg_->layer2.exposed_proximal;
}

Real MucusSegment::patch_volume_mL() const {
  const Real v_m3 =
      cfg_->layer2.patch_lateral_m * cfg_->layer2.patch_lateral_m *
      cfg_->layer2.patch_depth_m;
  return v_m3 * 1.0e6;  // m^3 -> mL
}

Int MucusSegment::spatial_capacity_cells() const {
  const Real cell_vol =
      (4.0 / 3.0) * PI * std::pow(cfg_->layer2.cell_radius_m, 3);
  const Real patch_vol = cfg_->layer2.patch_lateral_m *
                         cfg_->layer2.patch_lateral_m *
                         cfg_->layer2.patch_depth_m;
  const Real cap = cfg_->layer2.spatial_packing_fraction * patch_vol / cell_vol;
  return static_cast<Int>(std::floor(cap));
}

void MucusSegment::validate_config() const {
  const Layer2Config& c = cfg_->layer2;
  if (c.n_patches <= 0) {
    throw ConfigError("layer2.n_patches must be > 0");
  }
  const Real fsum = c.crypt.fraction + c.exposed_proximal.fraction +
                    c.exposed_distal.fraction;
  if (fsum <= 0.0) {
    throw ConfigError("layer2 patch-type fractions sum to zero");
  }
  if (c.occupancy_min_agents < 1) {
    throw ConfigError("layer2.occupancy_min_agents must be >= 1");
  }
  if (c.contraction_rate_per_min < 0.0 || c.agent_loss_fraction < 0.0 ||
      c.agent_loss_fraction > 1.0) {
    throw ConfigError("layer2 contraction/loss parameters out of range");
  }
  if (c.establish_prob_single < 0.0 || c.establish_prob_single > 1.0 ||
      c.establish_ratio < 0.0) {
    throw ConfigError("layer2 establishment parameters out of range");
  }
  if (cfg_->initial_strains.empty()) {
    throw ConfigError("layer2 requires at least one initial_strains entry");
  }
  if (c.audit.patch_index >= c.n_patches) {
    throw ConfigError("layer2.audit_patch_index out of range");
  }
}

void MucusSegment::assign_patch_types() {
  const Layer2Config& c = cfg_->layer2;
  const Real fsum = c.crypt.fraction + c.exposed_proximal.fraction +
                    c.exposed_distal.fraction;
  const Real p_crypt = c.crypt.fraction / fsum;
  const Real p_prox = p_crypt + c.exposed_proximal.fraction / fsum;
  for (auto& patch : patches_) {
    const Real u = patch.rng.uniform();
    if (u < p_crypt) {
      patch.type = PatchType::Crypt;
    } else if (u < p_prox) {
      patch.type = PatchType::ExposedProximal;
    } else {
      patch.type = PatchType::ExposedDistal;
    }
  }
}

void MucusSegment::make_founder(Patch& patch, Int strain_type, Real mu_max,
                                const std::vector<std::string>& plasmids,
                                bool conjugative) {
  const Vec3 pos{0.5 * cfg_->layer2.patch_lateral_m,
                 0.5 * cfg_->layer2.patch_lateral_m,
                 0.5 * cfg_->layer2.patch_depth_m};
  Agent a = Agent::create_default(patch.next_cell_tag++, strain_type, pos,
                                  mu_max);
  a.genome.lineage_id = a.identity.tag;  // each founder roots a lineage
  a.genome.parent_id = 0;
  for (const auto& name : plasmids) {
    const PlasmidEntry* entry = PlasmidLibrary::find(name);
    if (entry == nullptr) {
      throw ConfigError("layer2 founder plasmid not found: " + name);
    }
    a.genome.bi_loci.push_back(entry->cluster);
  }
  a.genome.has_conjugative_plasmid = conjugative;
  patch.colonists.push_back(std::move(a));
}

void MucusSegment::seed_initial_colonists() {
  const Layer2Config& c = cfg_->layer2;
  Int total_strain_count = 0;
  for (const auto& s : cfg_->initial_strains) {
    total_strain_count += std::max<Int>(s.count, 0);
  }
  if (total_strain_count <= 0) {
    throw ConfigError("layer2 initial_strains have zero total count");
  }
  for (auto& patch : patches_) {
    const bool seeded = patch.rng.bernoulli(c.initial_patch_fraction);
    if (!seeded) continue;
    // Proportional founder split over the configured strains.
    Int remaining = c.initial_founder_cells;
    for (size_t si = 0; si < cfg_->initial_strains.size(); ++si) {
      const auto& strain = cfg_->initial_strains[si];
      Int n = static_cast<Int>(
          std::llround(static_cast<Real>(c.initial_founder_cells) *
                       static_cast<Real>(strain.count) /
                       static_cast<Real>(total_strain_count)));
      if (si + 1 == cfg_->initial_strains.size()) {
        n = remaining;  // last strain absorbs rounding remainder
      }
      n = std::clamp<Int>(n, 0, remaining);
      remaining -= n;
      for (Int i = 0; i < n; ++i) {
        make_founder(patch, strain.type, strain.mu_max, strain.plasmids,
                     strain.conjugative);
      }
    }
  }
  for (auto& patch : patches_) {
    patch.baseline_count = static_cast<Int>(patch.colonists.size());
  }
}

void MucusSegment::admissibility_check() const {
  const Int cap = spatial_capacity_cells();
  for (const auto& patch : patches_) {
    if (static_cast<Int>(patch.colonists.size()) > cap) {
      throw ConfigError(
          std::format("layer2 startup inadmissible: patch {} seeds {} cells "
                      "above spatial capacity {}",
                      patch.id, patch.colonists.size(), cap));
    }
  }
  const SegmentObservables obs = observables();
  if (cfg_->layer2.segment_dysbiosis_threshold > 0.0 &&
      obs.segment_mean_cfu_ml >= cfg_->layer2.segment_dysbiosis_threshold) {
    throw ConfigError(
        std::format("layer2 startup inadmissible: initial segment mean "
                    "{:.3e} cells/mL exceeds absolute ladder {}",
                    obs.segment_mean_cfu_ml,
                    cfg_->layer2.segment_dysbiosis_threshold));
  }
}

void MucusSegment::init(const SimulationConfig& cfg) {
  cfg_ = &cfg;
  validate_config();
  table_ = std::make_unique<PatchTable>();

  const Layer2Config& c = cfg.layer2;
  contraction_stream_.seed(mix_seed(cfg.seed, "contraction"));
  init_stream_.seed(mix_seed(cfg.seed, "init"));

  patches_.resize(static_cast<size_t>(c.n_patches));
  order_.resize(static_cast<size_t>(c.n_patches));
  std::iota(order_.begin(), order_.end(), 0);
  for (size_t i = 0; i < patches_.size(); ++i) {
    patches_[i].id = static_cast<Int>(i);
    patches_[i].next_cell_tag =
        static_cast<TagID>(static_cast<int64_t>(i) + 1) << 32;
    patches_[i].rng.seed(mix_seed(cfg.seed, static_cast<uint64_t>(i)));
  }
  next_tag_ = (static_cast<TagID>(c.n_patches) + 1) << 32;

  assign_patch_types();
  seed_initial_colonists();

  // Optional pre-seeded singles in the pool.
  for (Int i = 0; i < c.initial_pool_cells; ++i) {
    const auto& strain = cfg.initial_strains.front();
    const Vec3 pos{0, 0, 0};
    Agent a = Agent::create_default(next_tag_++, strain.type, pos,
                                    strain.mu_max);
    a.genome.lineage_id = a.identity.tag;
    const uint64_t key =
        mix_seed(static_cast<uint64_t>(a.identity.tag),
                 tag_hash("pool_origin"));
    push_packet(PacketKind::Single, {std::move(a)}, key);
  }

  init_audit_slot();
  admissibility_check();

  Int stocks = 0;
  for (const auto& patch : patches_) {
    stocks += patch_count(patch);
  }
  for (const auto& packet : pool_) {
    stocks += static_cast<Int>(packet.cells.size());
  }
  ledger_.initial_cells = stocks;

  audit_next_time_ = c.audit.patch_index >= 0 ? c.audit.interval_s : -1;
}

Int MucusSegment::patch_count(const Patch& patch) const {
  if (patch.live) {
    const auto it = live_.find(patch.id);
    return it == live_.end() ? 0
                             : static_cast<Int>(it->second.sim->agents().size());
  }
  return static_cast<Int>(patch.colonists.size());
}

void MucusSegment::push_packet(PacketKind kind, std::vector<Agent> cells,
                               uint64_t stream_key) {
  if (cells.empty()) return;
  LuminalPacket packet;
  packet.id = next_packet_id_++;
  packet.kind = kind;
  packet.cells = std::move(cells);
  packet.entry_time = time_;
  packet.rng.seed(mix_seed(cfg_->seed, stream_key));
  pool_.push_back(std::move(packet));
}

// ── Tabulated patch mechanics ────────────────────────────────────────────

void MucusSegment::grow_tabulated(Patch& patch, Real dt) {
  const Int n = static_cast<Int>(patch.colonists.size());
  if (n == 0) return;
  const auto& tc = type_config(patch.type);
  const Real k = table_->capacity_cells(tc.supply_mult, patch_volume_mL());
  const Real g = table_->growth_rate(tc.supply_mult);
  const Real mean_births =
      g * static_cast<Real>(n) *
      std::max(0.0, 1.0 - static_cast<Real>(n) / k) * dt;
  const Int births = patch.rng.poisson(mean_births);
  for (Int i = 0; i < births; ++i) {
    const Int parent_idx = patch.rng.randint(0, n - 1);
    Agent daughter = patch.colonists[static_cast<size_t>(parent_idx)];
    daughter.identity.tag = patch.next_cell_tag++;
    daughter.genome.parent_id =
        patch.colonists[static_cast<size_t>(parent_idx)].identity.tag;
    daughter.genome.generation++;
    patch.colonists.push_back(std::move(daughter));
  }
  ledger_.births += births;
}

void MucusSegment::contract_tabulated(Patch& patch, bool ledger) {
  const Real f = cfg_->layer2.agent_loss_fraction;
  std::vector<Agent> fragment;
  auto& cells = patch.colonists;
  const size_t before = cells.size();
  size_t write = 0;
  for (size_t i = 0; i < before; ++i) {
    if (patch.rng.bernoulli(f)) {
      fragment.push_back(std::move(cells[i]));
    } else {
      if (write != i) {
        cells[write] = std::move(cells[i]);
      }
      ++write;
    }
  }
  cells.resize(write);
  if (!ledger) return;
  ledger_.contraction_departures += static_cast<Int>(fragment.size());
  const uint64_t key = mix_seed(static_cast<uint64_t>(patch.id),
                                static_cast<uint64_t>(step_count_));
  push_packet(PacketKind::Fragment, std::move(fragment), key);
}

void MucusSegment::washout_tabulated(Patch& patch, Real dt, bool ledger) {
  const Real rate = type_config(patch.type).single_cell_loss_per_s;
  if (rate <= 0.0) return;
  auto& cells = patch.colonists;
  std::vector<std::pair<uint64_t, Agent>> departed;
  size_t write = 0;
  for (size_t i = 0; i < cells.size(); ++i) {
    if (patch.rng.bernoulli(rate * dt)) {
      departed.emplace_back(
          mix_seed(static_cast<uint64_t>(cells[i].identity.tag),
                   static_cast<uint64_t>(patch.id)),
          std::move(cells[i]));
    } else {
      if (write != i) {
        cells[write] = std::move(cells[i]);
      }
      ++write;
    }
  }
  cells.resize(write);
  if (!ledger) return;
  for (auto& [tag, cell] : departed) {
    push_packet(PacketKind::Single, {std::move(cell)}, tag);
    ++ledger_.washout_departures;
  }
}

// ── Live patch (audit gate) ──────────────────────────────────────────────

SimulationConfig MucusSegment::live_config_for(const Patch& patch) const {
  SimulationConfig live = *cfg_;
  const auto& tc = type_config(patch.type);
  live.domain.lo = {0.0, 0.0, 0.0};
  live.domain.hi = {cfg_->layer2.patch_lateral_m, cfg_->layer2.patch_lateral_m,
                    cfg_->layer2.patch_depth_m};
  live.domain.periodic = {true, true, false};
  live.seed = mix_seed(cfg_->seed, static_cast<uint64_t>(patch.id),
                       tag_hash("live"));
  live.carbon_boundary_conc = 5.0e-3 * tc.supply_mult;
  live.advection.mucus_thickness = cfg_->layer2.patch_depth_m;
  live.advection.distal_length = cfg_->layer2.patch_lateral_m;
  live.advection.crypts_enabled = false;  // S2: crypt is patch type only
  live.chem_env.oxygen.enabled = false;   // ladder condition
  live.chem_env.oxygen.k_ROS_respiratory = 0.0;
  live.chem_env.oxygen.k_ROS_funded = 0.0;
  live.dysbiosis_threshold = 0.0;         // guards live at segment level
  live.immigration.enabled = false;
  live.gpu.enabled = false;               // audit slots are CPU-only
  live.hdf5.enabled = false;
  live.restart.enabled = false;
  live.checkpoint = {};
  live.restart = {};
  live.time.total_time = std::numeric_limits<Real>::max();
  return live;
}

void MucusSegment::init_audit_slot() {
  const auto& c = cfg_->layer2;
  if (c.audit.patch_index < 0) return;
  Patch& patch = patches_[static_cast<size_t>(c.audit.patch_index)];

  LivePatchState state;
  state.sim = std::make_unique<Simulation>();
  SimulationConfig live_cfg = live_config_for(patch);
  // Founders mirror the patch's seeded colonist composition.
  live_cfg.initial_strains.clear();
  for (const auto& strain : cfg_->initial_strains) {
    Int count = 0;
    for (const auto& a : patch.colonists) {
      if (a.identity.type == strain.type) {
        ++count;
      }
    }
    SimulationConfig::InitialStrain s = strain;
    s.count = count;
    live_cfg.initial_strains.push_back(s);
  }
  state.sim->init(live_cfg);
  state.sim->set_departure_hook(
      [this, id = patch.id](const Agent& a) {
        live_[id].harvested_singles.push_back(a);
      });
  state.sim->set_late_step_hook([this, id = patch.id]() {
    LivePatchState& st = live_[id];
    if (st.pending_fragment_fraction <= 0.0) return;
    AgentPool& agents = st.sim->agents();
    const Int n = agents.size();
    const Int take = std::min<Int>(
        n, static_cast<Int>(std::ceil(st.pending_fragment_fraction *
                                      static_cast<Real>(n))));
    // Fragment = contiguous downstream chunk: sort indices by x desc.
    std::vector<Int> order(static_cast<size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&agents](Int a, Int b) {
                return agents[a].x[0] > agents[b].x[0];
              });
    for (Int i = 0; i < take; ++i) {
      Agent& a = agents[order[static_cast<size_t>(i)]];
      st.harvested_fragment.push_back(a);
      a.state = PhenoState::DEAD;
    }
    st.pending_fragment_fraction = 0.0;
  });
  // Shadow: tabulated twin starting from the same population.
  state.shadow.id = patch.id;
  state.shadow.type = patch.type;
  state.shadow.rng.seed(mix_seed(cfg_->seed,
                                 static_cast<uint64_t>(patch.id),
                                 tag_hash("shadow")));
  for (const auto& a : state.sim->agents()) {
    state.shadow.colonists.push_back(a);
  }
  state.shadow.baseline_count =
      static_cast<Int>(state.shadow.colonists.size());
  state.shadow_valid = true;
  state.window_start_time = time_;
  state.last_live_birth_count = state.sim->cumulative_events().divisions;

  patch.live = true;
  patch.colonists.clear();  // population lives inside the live sim now
  live_.emplace(patch.id, std::move(state));
}

void MucusSegment::step_live(Patch& patch, Real dt) {
  LivePatchState& st = live_.at(patch.id);
  st.pending_fragment_fraction =
      patch.disrupted ? cfg_->layer2.agent_loss_fraction : 0.0;
  const Int born_before = st.last_live_birth_count;
  st.sim->step(dt);
  const Int divisions = st.sim->cumulative_events().divisions;
  ledger_.births += divisions - born_before;
  st.last_live_birth_count = divisions;

  if (!st.harvested_fragment.empty()) {
    ledger_.contraction_departures +=
        static_cast<Int>(st.harvested_fragment.size());
    const uint64_t key = mix_seed(static_cast<uint64_t>(patch.id),
                                  static_cast<uint64_t>(step_count_));
    push_packet(PacketKind::Fragment, std::move(st.harvested_fragment),
                key);
    st.harvested_fragment.clear();
  }
  for (auto& cell : st.harvested_singles) {
    push_packet(PacketKind::Single, {std::move(cell)},
                mix_seed(static_cast<uint64_t>(cell.identity.tag),
                         static_cast<uint64_t>(patch.id)));
    ++ledger_.washout_departures;
  }
  st.harvested_singles.clear();

  shadow_step(patch, dt);
}

void MucusSegment::shadow_step(Patch& live_patch, Real dt) {
  LivePatchState& st = live_.at(live_patch.id);
  if (!st.shadow_valid) return;
  grow_tabulated(st.shadow, dt);
  if (live_patch.disrupted) {
    contract_tabulated(st.shadow, false);
  }
  washout_tabulated(st.shadow, dt, false);
}

// ── Step driver ──────────────────────────────────────────────────────────

void MucusSegment::step(Real dt) {
  if (termination_cause_ != TerminationCause::IncompleteUnknown) return;

  const Real p_contract = cfg_->layer2.contraction_rate_per_min * dt / 60.0;
  const bool fire = p_contract > 0.0 && contraction_stream_.bernoulli(p_contract);

  for (const Int idx : order_) {
    Patch& patch = patches_[static_cast<size_t>(idx)];
    patch.disrupted =
        fire && patch.status == PatchStatus::Active &&
        patch.rng.bernoulli(type_config(patch.type).disruption_prob);
  }

  for (const Int idx : order_) {
    Patch& patch = patches_[static_cast<size_t>(idx)];
    if (patch.status != PatchStatus::Active) continue;
    if (patch.live) {
      step_live(patch, dt);
      continue;
    }
    grow_tabulated(patch, dt);
    if (patch.disrupted) {
      contract_tabulated(patch, true);
    }
    washout_tabulated(patch, dt, true);
  }

  pool_phase(dt);

  for (const Int idx : order_) {
    guard_update(patches_[static_cast<size_t>(idx)]);
  }
  run_level_guards();
  audit_check();

  time_ += dt;
  ++step_count_;
}

void MucusSegment::pool_phase(Real dt) {
  const auto& c = cfg_->layer2;
  const Real k_decay = kLn2 / std::max(c.transit_half_life_s, 1e-12);
  const Real p_exit = 1.0 - std::exp(-k_decay * dt);
  const Real p_attempt = (1.0 - std::exp(-k_decay * dt)) *
                         c.reattach_prob_per_transit;
  const Real p_frag =
      std::min(1.0, c.establish_prob_single * c.establish_ratio);

  size_t n_packets = pool_.size();
  for (size_t i = 0; i < n_packets; ++i) {
    LuminalPacket packet = std::move(pool_.front());
    pool_.pop_front();
    if (packet.rng.bernoulli(p_exit)) {
      ledger_.distal_losses += static_cast<Int>(packet.cells.size());
      continue;
    }
    if (packet.rng.bernoulli(p_attempt)) {
      const Int target =
          packet.rng.randint(0, static_cast<Int>(patches_.size()) - 1);
      Patch& dest = patches_[static_cast<size_t>(target)];
      const Real p_est = packet.kind == PacketKind::Fragment
                             ? p_frag
                             : c.establish_prob_single;
      if (dest.status == PatchStatus::Active && !dest.live &&
          packet.rng.bernoulli(p_est)) {
        ledger_.reseeds += static_cast<Int>(packet.cells.size());
        for (auto& cell : packet.cells) {
          dest.colonists.push_back(std::move(cell));
        }
        continue;
      }
      if (dest.live && packet.rng.bernoulli(p_est)) {
        ledger_.reseeds += static_cast<Int>(packet.cells.size());
        for (auto& cell : packet.cells) {
          cell.identity.owner_rank = 0;
          cell.grid_cell = -1;
          LivePatchState& st = live_.at(dest.id);
          st.sim->agents().push_back(std::move(cell));
        }
        continue;
      }
    }
    pool_.push_back(std::move(packet));
  }
}

// ── Guards (S4) ──────────────────────────────────────────────────────────

void MucusSegment::guard_update(Patch& patch) {
  if (patch.status != PatchStatus::Active) return;
  const Int n = patch_count(patch);

  if (n > spatial_capacity_cells()) {
    patch.status = PatchStatus::SpatialHalted;
    return;
  }
  const Int bound = static_cast<Int>(
      std::llround(cfg_->layer2.bloom_factor *
                   static_cast<Real>(patch.baseline_count)));
  if (patch.baseline_count > 0 && n > bound) {
    if (patch.bloom_excess_since < 0.0) {
      patch.bloom_excess_since = time_;
    }
    if (time_ - patch.bloom_excess_since >= cfg_->layer2.bloom_sustain_s) {
      patch.status = PatchStatus::BloomHalted;
    }
  } else {
    patch.bloom_excess_since = -1.0;
  }
}

void MucusSegment::run_level_guards() {
  if (termination_cause_ != TerminationCause::IncompleteUnknown) return;

  Int halted = 0;
  for (const auto& patch : patches_) {
    if (patch.status != PatchStatus::Active) ++halted;
  }
  const Real invalid_frac =
      static_cast<Real>(halted) / static_cast<Real>(patches_.size());
  if (invalid_frac > cfg_->layer2.invalid_patch_fraction_stop) {
    terminate(TerminationCause::PatchInvalidLimit,
              std::format("invalid patch fraction {:.3f} exceeds limit {:.3f}",
                          invalid_frac,
                          cfg_->layer2.invalid_patch_fraction_stop));
    return;
  }

  const SegmentObservables stocks = observables();
  if (stocks.patch_cells + stocks.pool_cells == 0) {
    terminate(TerminationCause::PopulationStop,
              "all cells lost from patches and pool");
    return;
  }

  const auto& c = cfg_->layer2;
  if (c.segment_dysbiosis_threshold > 0.0) {
    const SegmentObservables obs = observables();
    mean_window_.emplace_back(time_, obs.segment_mean_cfu_ml);
    while (!mean_window_.empty() &&
           mean_window_.front().first < time_ - c.segment_guard_window_s) {
      mean_window_.pop_front();
    }
    Real sum = 0.0;
    for (const auto& [t, v] : mean_window_) {
      sum += v;
    }
    const Real trailing =
        mean_window_.empty()
            ? 0.0
            : sum / static_cast<Real>(mean_window_.size());
    if (!mean_window_.empty() &&
        time_ - mean_window_.front().first >= 0.9 * c.segment_guard_window_s &&
        trailing >= c.segment_dysbiosis_threshold) {
      terminate(TerminationCause::DysbiosisGuard,
                std::format("segment mean {:.3e} cells/mL sustained above "
                            "ladder {:.3e} over {:.0f}s window",
                            trailing, c.segment_dysbiosis_threshold,
                            c.segment_guard_window_s));
      return;
    }
  }

  if (!ledger_closed()) {
    terminate(TerminationCause::ClosureViolation,
              "population ledger does not close under the two-channel design");
  }
}

void MucusSegment::audit_check() {
  const auto& c = cfg_->layer2;
  if (audit_next_time_ < 0.0 || time_ < audit_next_time_) return;
  audit_next_time_ += c.audit.interval_s;
  for (auto& [id, st] : live_) {
    const Int n_live = static_cast<Int>(st.sim->agents().size());
    const Int n_shadow = static_cast<Int>(st.shadow.colonists.size());
    const Real denom = std::max<Real>(static_cast<Real>(n_shadow), 1.0);
    const Real rel = std::abs(n_live - n_shadow) / denom;
    audit_records_.push_back({time_, id, n_live, n_shadow, rel});
    if (rel > c.audit.tolerance && c.audit.halt_on_exceed) {
      terminate(TerminationCause::AuditDivergence,
                std::format("audit patch {} diverged: live={} shadow={} "
                            "(rel {:.3f} > {:.3f})",
                            id, n_live, n_shadow, rel, c.audit.tolerance));
      return;
    }
    // Reset the comparison window: shadow re-seeded from live state.
    st.shadow.colonists.clear();
    for (const auto& a : st.sim->agents()) {
      st.shadow.colonists.push_back(a);
    }
    st.shadow.baseline_count =
        static_cast<Int>(st.shadow.colonists.size());
  }
}

void MucusSegment::terminate(TerminationCause cause, std::string detail) {
  if (termination_cause_ != TerminationCause::IncompleteUnknown) return;
  termination_cause_ = cause;
  termination_detail_ = std::move(detail);
}

// ── Observables (G5/G6) ──────────────────────────────────────────────────

SegmentObservables MucusSegment::observables() const {
  SegmentObservables obs;
  const Real vol = patch_volume_mL();
  Real frac_sum = 0.0;
  for (Int t = 0; t < kPatchTypeCount; ++t) {
    frac_sum += type_config(static_cast<PatchType>(t)).fraction;
  }

  for (const auto& patch : patches_) {
    TypeObservables& to =
        obs.per_type[to_underlying(patch.type)];
    ++to.n_patches;
    if (patch.status != PatchStatus::Active) {
      ++to.n_halted;
    }
    const Int n = patch_count(patch);
    if (n >= cfg_->layer2.occupancy_min_agents) {
      ++to.n_occupied;
      to.mean_density_occupied += static_cast<Real>(n) / vol;
      ++obs.n_occupied_total;
    }
    obs.patch_cells += n;
  }
  for (const auto& packet : pool_) {
    obs.pool_cells += static_cast<Int>(packet.cells.size());
  }
  obs.pool_packets = static_cast<Int>(pool_.size());

  for (Int t = 0; t < kPatchTypeCount; ++t) {
    TypeObservables& to = obs.per_type[t];
    const auto& tc = type_config(static_cast<PatchType>(t));
    if (to.n_occupied > 0) {
      to.mean_density_occupied /= static_cast<Real>(to.n_occupied);
    }
    to.occupancy =
        to.n_patches > 0
            ? static_cast<Real>(to.n_occupied) / static_cast<Real>(to.n_patches)
            : 0.0;
    const Real f_t = frac_sum > 0.0 ? tc.fraction / frac_sum : 0.0;
    obs.segment_mean_cfu_ml += f_t * to.occupancy * to.mean_density_occupied;
  }
  // G6: cells/cm^3 * layer thickness (cm) = cells/cm^2.
  obs.segment_cfu_cm2 =
      obs.segment_mean_cfu_ml * cfg_->layer2.mucus_thickness_m * 100.0;
  return obs;
}

bool MucusSegment::ledger_closed() const {
  Int stocks = 0;
  for (const auto& patch : patches_) {
    stocks += patch_count(patch);
  }
  for (const auto& packet : pool_) {
    stocks += static_cast<Int>(packet.cells.size());
  }
  return stocks == ledger_.expected_stocks();
}

uint64_t MucusSegment::fingerprint() const {
  uint64_t h = 0x9e3779b97f4a7c15ULL;
  for (const auto& patch : patches_) {
    h = hash_combine(h, static_cast<uint64_t>(patch.id));
    h = hash_combine(h, static_cast<uint64_t>(to_underlying(patch.type)));
    h = hash_combine(h, static_cast<uint64_t>(to_underlying(patch.status)));
    if (patch.live) {
      const auto& st = live_.at(patch.id);
      for (const auto& a : st.sim->agents()) {
        h = hash_combine(h, static_cast<uint64_t>(a.identity.tag));
      }
    } else {
      for (const auto& a : patch.colonists) {
        h = hash_combine(h, static_cast<uint64_t>(a.identity.tag));
        h = hash_combine(h, static_cast<uint64_t>(a.genome.lineage_id));
      }
    }
  }
  // Packet identity for the fingerprint is content (kind + sorted cell
  // tags): the creation-order id label legitimately differs under patch
  // iteration permutations (G3) while the trajectory must not.
  std::vector<uint64_t> packet_keys;
  packet_keys.reserve(pool_.size());
  for (const auto& packet : pool_) {
    std::vector<uint64_t> tags;
    tags.reserve(packet.cells.size());
    for (const auto& a : packet.cells) {
      tags.push_back(static_cast<uint64_t>(a.identity.tag));
    }
    std::sort(tags.begin(), tags.end());
    uint64_t pk = static_cast<uint64_t>(to_underlying(packet.kind));
    for (const uint64_t tag : tags) {
      pk = hash_combine(pk, tag);
    }
    packet_keys.push_back(pk);
  }
  std::sort(packet_keys.begin(), packet_keys.end());
  for (const uint64_t pk : packet_keys) {
    h = hash_combine(h, pk);
  }
  h = hash_combine(h, static_cast<uint64_t>(ledger_.births));
  h = hash_combine(h, static_cast<uint64_t>(ledger_.washout_departures));
  h = hash_combine(h, static_cast<uint64_t>(ledger_.contraction_departures));
  h = hash_combine(h, static_cast<uint64_t>(ledger_.reseeds));
  h = hash_combine(h, static_cast<uint64_t>(ledger_.distal_losses));
  return h;
}

void MucusSegment::permute_patch_order_for_testing(
    const std::vector<Int>& order) {
  if (order.size() != order_.size()) {
    throw ConfigError("patch iteration order has wrong size");
  }
  order_ = order;
}

}  // namespace gutibm
