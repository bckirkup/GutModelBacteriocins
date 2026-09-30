/* -----------------------------------------------------------------------
   GutIBM – Spec 13 Phase 3: three-region colonic chain driver.
   ----------------------------------------------------------------------- */

#include "chain.h"

#include "error.h"
#include "input_parser.h"
#include "plasmid.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <numeric>

namespace gutibm {

namespace {

uint64_t splitmix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

uint64_t hash_combine(uint64_t h, uint64_t v) {
  return splitmix64(h ^ v);
}

uint64_t tag_hash(const char* tag) {
  uint64_t h = 0xcbf29ce484222325ULL;  // FNV-1a
  for (const char* c = tag; *c != '\0'; ++c) {
    h = (h ^ static_cast<unsigned char>(*c)) * 0x100000001b3ULL;
  }
  return h;
}

// Deposit tags mint by the chain for reattached cells: region prefix
// bits 56-63 plus the 0xFF marker in bits 48-55 — disjoint from every
// patch tag space ((patch_id + 1) << 32 < 2^48).
constexpr TagID kDepositMarker = static_cast<TagID>(0xFF) << 48;
// Seed exemplar tags are region-independent (same inoculum genotype in
// every region's lumen): 0xFD marker in bits 48-55, strain index + 1 in
// the low word.
constexpr TagID kSeedMarker = static_cast<TagID>(0xFD) << 48;

}  // namespace

const char* region_name(RegionKind r) {
  switch (r) {
    case RegionKind::CecumAscending: return "cecum";
    case RegionKind::Transverse: return "transverse";
    case RegionKind::DescendingSigmoid: return "descending";
  }
  return "unknown";
}

ColonicChain::ColonicChain() = default;
ColonicChain::~ColonicChain() = default;

uint64_t ColonicChain::mix_seed(uint64_t seed, uint64_t a, uint64_t b) {
  return splitmix64(splitmix64(seed ^ splitmix64(a)) ^ splitmix64(b));
}
uint64_t ColonicChain::mix_seed(uint64_t seed, uint64_t a) {
  return splitmix64(seed ^ splitmix64(a));
}
uint64_t ColonicChain::mix_seed(uint64_t seed, const char* tag) {
  return splitmix64(seed ^ tag_hash(tag));
}

const Layer3RegionConfig& ColonicChain::region_cfg(RegionKind r) const {
  return regions_cfg_[static_cast<size_t>(r)];
}

void ColonicChain::validate_config() const {
  const Layer3Config& c = cfg_->layer3;
  if (c.n_patches_per_region <= 0 || c.n_patches_per_region > 0x10000) {
    throw ConfigError("layer3.n_patches_per_region out of range");
  }
  for (const Layer3RegionConfig& rc : regions_cfg_) {
    if (rc.transit_h <= 0.0 || rc.ph <= 0.0 || rc.lumen_volume_L <= 0.0 ||
        rc.contraction_k0_per_min < 0.0 || rc.outflow_survival < 0.0 ||
        rc.outflow_survival > 1.0 || rc.flora_density_cfu_ml <= 0.0) {
      throw ConfigError("layer3 region profile values out of range");
    }
  }
  if (c.f_edge < 0.0 || c.f_edge > 1.0 || c.reattach_p0 < 0.0 ||
      c.reattach_p0 > 1.0 || c.alpha_hapc_per_min < 0.0 ||
      c.hapc_rate_per_day < 0.0 || c.hapc_antegrade_fraction < 0.0 ||
      c.hapc_antegrade_fraction > 1.0 || c.hapc_day_fraction < 0.0 ||
      c.hapc_day_fraction > 1.0 || c.stool_g_per_day <= 0.0 ||
      c.purge_fraction < 0.0 || c.purge_fraction > 1.0 ||
      c.lumen_seed_cells < 0.0 || c.flatness_bound <= 0.0) {
    throw ConfigError("layer3 global parameters out of range");
  }
  if (cfg_->layer2.audit.patch_index >= 0) {
    throw ConfigError(
        "layer3 does not support the live-patch audit slot (embedded "
        "Simulation state is not chain-serializable)");
  }
  if (cfg_->initial_strains.empty()) {
    throw ConfigError("layer3 requires at least one initial_strains entry");
  }
}

std::vector<Agent> ColonicChain::seed_lumen_exemplars() const {
  // Luminal seed lineages carry the configured founders' genotypes with
  // region-independent seed tags: the same inoculum genotype is one
  // lineage everywhere it appears.
  std::vector<Agent> out;
  Int si = 0;
  for (const auto& strain : cfg_->initial_strains) {
    ++si;
    Agent a = Agent::create_default(kSeedMarker | si, strain.type,
                                    Vec3{0, 0, 0}, strain.mu_max);
    a.genome.lineage_id = a.identity.tag;
    a.genome.parent_id = 0;
    for (const auto& name : strain.plasmids) {
      const PlasmidEntry* entry = PlasmidLibrary::find(name);
      if (entry == nullptr) {
        throw ConfigError("layer3 founder plasmid not found: " + name);
      }
      a.genome.bi_loci.push_back(entry->cluster);
    }
    a.genome.has_conjugative_plasmid = strain.conjugative;
    out.push_back(std::move(a));
  }
  return out;
}

void ColonicChain::init(const SimulationConfig& cfg) {
  cfg_ = &cfg;
  const Layer3Config& c = cfg.layer3;

  // S7: per-region values come from the sourced config rows; the
  // sanity arm collapses every row to transverse.
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
  // Region 0's luminal inflow boundary reuses its own profile carbon —
  // the ileocecal endpoint of the sourced axial profile.
  physio_.diet_carbon_mol_m3 = regions_cfg_.front().lumen_carbon_mol_m3;

  const std::vector<Agent> exemplars = seed_lumen_exemplars();

  regions_.clear();
  regions_.resize(kRegionCount);
  region_cfgs_.clear();
  order_.resize(kRegionCount);
  std::iota(order_.begin(), order_.end(), 0);

  for (Int r = 0; r < kRegionCount; ++r) {
    const auto i = static_cast<size_t>(r);
    RegionState& rs = regions_[i];
    rs.next_deposit_tag = (static_cast<TagID>(r + 1) << 56) |
                          kDepositMarker;

    SimulationConfig& rcfg = region_cfgs_.emplace_back(cfg);
    // Per-region seed so contraction/init/pool streams are distinct;
    // tag prefix so cell tags cannot collide across regions.
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
    rs.segment->init(rcfg);
    rs.segment->set_pool_exit_sink(
        [this, i](std::vector<Agent>&& cells) {
          regions_[i].lumen.add_cells(cells);
        });
    rs.segment->set_contraction_rate_per_min(
        regions_cfg_[i].contraction_k0_per_min);
    rs.lumen.init(regions_cfg_[i], exemplars, c.lumen_seed_cells);
  }

  ledger_ = {};
  for (Int r = 0; r < kRegionCount; ++r) {
    RegionState& rs = regions_[static_cast<size_t>(r)];
    ledger_.initial_cells +=
        static_cast<Real>(rs.segment->ledger().initial_cells) +
        rs.lumen.total_cells();
  }
  purged_ = false;
  last_hapc_ = {};
  time_ = 0.0;
  step_count_ = 0;
  termination_cause_ = TerminationCause::IncompleteUnknown;
  termination_detail_.clear();
}

// ── Step driver ─────────────────────────────────────────────────────────

Real ColonicChain::hapc_rate_per_s(Real time_s) const {
  const Layer3Config& c = cfg_->layer3;
  const Real tod = std::fmod(time_s, 86400.0);
  const Real day_start = c.hapc_day_start_h * 3600.0;
  const bool day = tod >= day_start && tod < day_start + 43200.0;
  const Real per_day =
      day ? c.hapc_rate_per_day * c.hapc_day_fraction
          : c.hapc_rate_per_day * (1.0 - c.hapc_day_fraction);
  return per_day / 43200.0;
}

HapcEvent ColonicChain::draw_hapc(Real dt, Int step_index) const {
  HapcEvent ev;
  const Layer3Config& c = cfg_->layer3;
  RNG draw;
  draw.seed(mix_seed(cfg_->seed, static_cast<uint64_t>(step_index),
                     tag_hash("hapc")));
  if (!draw.bernoulli(hapc_rate_per_s(time_) * dt)) {
    return ev;
  }
  ev.fired = true;
  ev.origin = draw.randint(0, kRegionCount - 1);
  ev.antegrade = draw.bernoulli(c.hapc_antegrade_fraction);
  if (ev.antegrade) {
    for (Int r = ev.origin; r < kRegionCount; ++r) {
      ev.regions[static_cast<size_t>(r)] = true;
    }
  } else {
    for (Int r = 0; r <= ev.origin; ++r) {
      ev.regions[static_cast<size_t>(r)] = true;
    }
  }
  return ev;
}

void ColonicChain::process_reattach(Int r, const LumenStepPlan& plan,
                                    std::map<TagID, Int>& deposited) {
  RegionState& rs = regions_[static_cast<size_t>(r)];
  for (const auto& [id, attempts] : plan.reattach_events) {
    const LuminalLineage& lin = rs.lumen.lineages().at(id);
    Int placed = 0;
    for (Int k = 0; k < attempts; ++k) {
      Agent a = lin.exemplar;
      a.identity.tag = rs.next_deposit_tag++;
      a.genome.parent_id = lin.exemplar.identity.tag;
      // lineage_id unchanged: the reattached cell stays in its lineage.
      const uint64_t key = mix_seed(
          cfg_->seed, static_cast<uint64_t>(step_count_),
          static_cast<uint64_t>(r) * 0x1000003ULL +
              static_cast<uint64_t>(id) +
              static_cast<uint64_t>(k));
      if (rs.segment->try_deposit(std::move(a), key)) {
        ++placed;
      }
    }
    if (placed > 0) {
      deposited[id] = placed;
    }
  }
}

void ColonicChain::step(Real dt) {
  if (termination_cause_ != TerminationCause::IncompleteUnknown) return;
  const Layer3Config& c = cfg_->layer3;

  last_hapc_ = draw_hapc(dt, step_count_);

  // Mucosal phase: every region's segment steps with the HAPC-boosted
  // contraction rate; pool exits are routed into the same region's
  // lumen by the per-region sink.
  for (const Int idx : order_) {
    RegionState& rs = regions_[static_cast<size_t>(idx)];
    rs.segment->set_contraction_rate_per_min(
        regions_cfg_[static_cast<size_t>(idx)].contraction_k0_per_min +
        (last_hapc_.regions[static_cast<size_t>(idx)]
             ? c.alpha_hapc_per_min
             : 0.0));
    rs.segment->step(dt);
    if (rs.segment->termination_cause() !=
        TerminationCause::IncompleteUnknown) {
      terminate(rs.segment->termination_cause(),
                std::format("region {} ({}): {}", idx,
                            region_name(static_cast<RegionKind>(idx)),
                            rs.segment->termination_detail()));
      return;
    }
  }

  // Luminal phase — two-phase update (G3): all plans are computed
  // against the pre-step snapshot before any state is committed, so
  // region processing order cannot affect the result.
  std::vector<LumenStepPlan> plans(kRegionCount);
  std::vector<LuminalScalars> inflow(kRegionCount);
  inflow[0] = {physio_.diet_carbon_mol_m3, 0.0};
  for (Int r = 1; r < kRegionCount; ++r) {
    inflow[static_cast<size_t>(r)] =
        regions_[static_cast<size_t>(r - 1)].lumen.scalars();
  }
  for (Int r = 0; r < kRegionCount; ++r) {
    RegionState& rs = regions_[static_cast<size_t>(r)];
    const SegmentObservables obs = rs.segment->observables();
    const Real occ =
        static_cast<Real>(obs.n_occupied_total) /
        std::max(1.0, static_cast<Real>(c.n_patches_per_region));
    plans[static_cast<size_t>(r)] = rs.lumen.compute(
        dt, regions_cfg_[static_cast<size_t>(r)], c, physio_,
        inflow[static_cast<size_t>(r)], occ, step_count_,
        static_cast<int>(r), cfg_->seed);
  }

  // Apply phase: reattach attempts instantiate genotype-exemplar cells
  // into the region's mucosa; upstream outflow arrives at survival.
  for (const Int idx : order_) {
    RegionState& rs = regions_[static_cast<size_t>(idx)];
    std::map<TagID, Int> deposited;
    process_reattach(idx, plans[static_cast<size_t>(idx)], deposited);
    const LumenStepPlan* upstream =
        idx > 0 ? &plans[static_cast<size_t>(idx - 1)] : nullptr;
    const Real upstream_survival =
        idx > 0 ? regions_cfg_[static_cast<size_t>(idx - 1)]
                      .outflow_survival
                : 0.0;
    rs.lumen.apply(plans[static_cast<size_t>(idx)], upstream,
                   upstream_survival, deposited);
    last_mean_mu_[static_cast<size_t>(idx)] =
        plans[static_cast<size_t>(idx)].mean_mu_per_s;
    ledger_.luminal_births += plans[static_cast<size_t>(idx)].births;
  }

  // Inter-region transfer bookkeeping: outflow cells that do not
  // survive the trip are luminal deaths; region 2's survivors are
  // stool. Arrivals themselves are added inside apply().
  for (Int r = 0; r < kRegionCount; ++r) {
    const Real out = plans[static_cast<size_t>(r)].outflow_cells;
    const Real surv = regions_cfg_[static_cast<size_t>(r)].outflow_survival;
    ledger_.luminal_deaths += out * (1.0 - surv);
    if (r == kRegionCount - 1) {
      ledger_.stool_exports += out * surv;
    }
  }

  // One-time luminal purge (recovery-kinematics probe).
  if (!purged_ && c.purge_at_s >= 0.0 && time_ + dt > c.purge_at_s) {
    purged_ = true;
    for (RegionState& rs : regions_) {
      ledger_.purged += rs.lumen.purge(c.purge_fraction);
    }
  }

  ledger_.mucosal_births = 0.0;
  for (const RegionState& rs : regions_) {
    ledger_.mucosal_births +=
        static_cast<Real>(rs.segment->ledger().births);
  }

  if (!ledger_closed()) {
    terminate(TerminationCause::ClosureViolation,
              std::format("chain ledger diverged: stocks {:.6e} vs "
                          "expected {:.6e}",
                          stocks(), ledger_.expected_stocks()));
    return;
  }

  time_ += dt;
  ++step_count_;
}

void ColonicChain::permute_region_order_for_testing(
    const std::vector<Int>& order) {
  if (order.size() != order_.size()) {
    throw ConfigError("region iteration order has wrong size");
  }
  order_ = order;
}

// ── Observables ──────────────────────────────────────────────────────────

Real ColonicChain::stocks() const {
  Real s = 0.0;
  for (const RegionState& rs : regions_) {
    const SegmentObservables obs = rs.segment->observables();
    s += static_cast<Real>(obs.patch_cells + obs.pool_cells) +
         rs.lumen.total_cells();
  }
  return s;
}

bool ColonicChain::ledger_closed() const {
  const Real expected = ledger_.expected_stocks();
  const Real tol = std::max(0.5, std::abs(expected) * 1e-6);
  return std::abs(stocks() - expected) <= tol;
}

LumenObservables ColonicChain::lumen_observables(Int i) const {
  const RegionState& rs = regions_[static_cast<size_t>(i)];
  LumenObservables o;
  o.total_cells = rs.lumen.total_cells();
  o.n_lineages = static_cast<Int>(rs.lumen.lineages().size());
  o.mean_mu_per_s = last_mean_mu_[static_cast<size_t>(i)];
  o.scalars = rs.lumen.scalars();
  return o;
}

Real ColonicChain::stool_cfu_per_g() const {
  const Layer3Config& c = cfg_->layer3;
  if (time_ <= 0.0) return 0.0;
  const Real days = time_ / 86400.0;
  const Real cells_per_day = ledger_.stool_exports / days;
  return cells_per_day / c.stool_g_per_day;
}

uint64_t ColonicChain::fingerprint() const {
  uint64_t h = 0x9e3779b97f4a7c15ULL;
  for (const RegionState& rs : regions_) {
    h = hash_combine(h, rs.segment->fingerprint());
    h = hash_combine(h, rs.lumen.fingerprint());
    h = hash_combine(h, static_cast<uint64_t>(rs.next_deposit_tag));
  }
  h = hash_combine(h, static_cast<uint64_t>(ledger_.initial_cells * 1e3));
  h = hash_combine(h, static_cast<uint64_t>(ledger_.luminal_births * 1e6));
  h = hash_combine(h, static_cast<uint64_t>(ledger_.stool_exports * 1e6));
  h = hash_combine(h, static_cast<uint64_t>(ledger_.luminal_deaths * 1e6));
  h = hash_combine(h, static_cast<uint64_t>(ledger_.purged * 1e6));
  h = hash_combine(h, purged_ ? 1ULL : 0ULL);
  return h;
}

void ColonicChain::terminate(TerminationCause cause, std::string detail) {
  termination_cause_ = cause;
  termination_detail_ = std::move(detail);
}

}  // namespace gutibm
