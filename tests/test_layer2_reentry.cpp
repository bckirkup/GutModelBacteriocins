#include "input_parser.h"
#include "segment.h"
#include "simulation.h"
#include "test_shard_dispatch.h"

#include <cassert>
#include <cmath>
#include <iostream>

namespace gutibm {

// Stop between production phases: a subsequent live step can export a cell,
// and the subsequent pool phase could immediately route it again.
struct SingleSlotReentryFixture {
  static void advance_live(MucusSegment& seg, Real dt) {
    seg.step_live(seg.patches_.at(0), dt);
  }

  static void export_at_surface(MucusSegment& seg, Real dt, Real depth) {
    AgentPool& cells = seg.live_.at(0).sim->agents();
    assert(cells.size() == 1);
    cells[0].x = {7e-6, 11e-6, std::nextafter(depth, 0.0)};
    advance_live(seg, dt);
    assert(seg.ledger().washout_departures == 1);
    assert(seg.live_.at(0).sim->step_events().outflow_boundary == 1);
  }

  static LuminalPacket& prepare(MucusSegment& seg, PatchStatus status) {
    assert(seg.patches_.size() == 1 && seg.pool_.size() == 1);
    seg.patches_.front().status = status;
    return seg.pool_.front();
  }

  static void route(MucusSegment& seg, Real dt) { seg.pool_phase(dt); }
};
}  // namespace gutibm

namespace {
using namespace gutibm;
using Fixture = SingleSlotReentryFixture;

SimulationConfig single_slot(bool live, bool boundary) {
  auto cfg = InputParser::default_config();
  cfg.seed = 1234;
  cfg.initial_strains = {SimulationConfig::InitialStrain{}};
  cfg.initial_strains.front().count = 1;
  if (live) {
    // An empty allowlist loads all Fixes. Use just an inactive one.
    cfg.enabled_fixes = {"motility"};
    cfg.cell_bio.motility.enabled = false;
  }
  auto& l = cfg.layer2;
  l.enabled = true;
  l.n_patches = 1;
  l.audit.patch_index = live ? 0 : -1;
  l.patch_lateral_m = 40e-6;
  l.patch_depth_m = 30e-6;
  l.initial_founder_cells = 1;
  l.initial_patch_fraction = boundary ? 1.0 : 0.0;
  l.initial_pool_cells = boundary ? 0 : 1;
  l.transit_half_life_s = cfg.time.bio_dt;
  l.reattach_prob_per_transit = 1.0;
  l.establish_prob_single = 1.0;
  return cfg;
}

// Exactly one copy of the tagged cell, in the expected stock.
const Agent& only_cell(const MucusSegment& seg, bool pooled) {
  assert(seg.ledger_closed() && seg.ledger().expected_stocks() == 1);
  const Agent* resident = nullptr;
  seg.for_each_colonist([&](const Agent& a) {
    assert(resident == nullptr);
    resident = &a;
  });
  if (pooled) {
    assert(resident == nullptr && seg.pool().size() == 1);
    assert(seg.pool().front().cells.size() == 1);
    return seg.pool().front().cells.front();
  }
  assert(seg.pool().empty() && resident != nullptr);
  return *resident;
}

void check_identity(const Agent& cell, const Agent& source) {
  assert(cell.identity.tag == source.identity.tag);
  assert(cell.identity.type == source.identity.type);
  assert(cell.genome.lineage_id == source.genome.lineage_id);
  assert(cell.genome.parent_id == source.genome.parent_id);
  assert(cell.state == source.state && cell.state != PhenoState::DEAD);
}

RNG arm_route(LuminalPacket& packet, uint64_t seed) {
  packet.rng.seed(seed);
  RNG draw = packet.rng;
  assert(!draw.bernoulli(0.5));  // dt = transit half-life
  assert(draw.bernoulli(0.5));   // attempts enabled; sole target
  assert(draw.randint(0, 0) == 0);
  const RNG after_target = draw;
  assert(draw.bernoulli(1.0));
  return after_target;
}

bool local_position(const Vec3& x, const Layer2Config& l) {
  return x[0] >= 0 && x[0] < l.patch_lateral_m &&
         x[1] >= 0 && x[1] < l.patch_lateral_m &&
         x[2] >= 0 && x[2] < l.patch_depth_m;
}

void check_interior(bool live, PatchStatus status, uint64_t seed) {
  const auto cfg = single_slot(live, false);
  MucusSegment seg;
  seg.init(cfg);
  auto& packet = Fixture::prepare(seg, status);
  assert(packet.kind == PacketKind::Single && packet.cells.size() == 1);
  auto& cell = packet.cells.front();
  cell.x = {7e-6, 11e-6, 13e-6};
  cell.genome.lineage_id = 431;
  cell.genome.parent_id = 432;
  const Agent source = only_cell(seg, true);
  const auto id = packet.id;
  const RNG after_target = arm_route(packet, seed);

  Fixture::route(seg, cfg.time.bio_dt);
  const bool active = status == PatchStatus::Active;
  const Agent& result = only_cell(seg, !active);
  assert(seg.ledger().reseeds == (active ? 1 : 0));
  assert(seg.ledger().distal_losses == 0);
  check_identity(result, source);
  assert(local_position(result.x, cfg.layer2) && result.x == source.x);
  if (!active) {
    const auto& retained = seg.pool().front();
    assert(retained.id == id && retained.rng.engine() == after_target.engine());
  } else if (live) {
    Fixture::advance_live(seg, cfg.time.bio_dt);
    check_identity(only_cell(seg, false), source);
    assert(seg.ledger().washout_departures == 0);
  }
}

void eligibility() {
  for (const auto seed : {3u, 4u}) {
    for (const bool live : {false, true}) {
      for (const auto status : {PatchStatus::Active, PatchStatus::BloomHalted,
                                PatchStatus::SpatialHalted}) {
        check_interior(live, status, seed);
      }
    }
  }
}

void boundary_witness() {
  const auto cfg = single_slot(true, true);
  for (const auto seed : {3u, 4u}) {
    MucusSegment seg;
    seg.init(cfg);
    Fixture::export_at_surface(seg, cfg.time.bio_dt, cfg.layer2.patch_depth_m);
    const Agent source = only_cell(seg, true);
    assert(source.x[2] >= cfg.layer2.patch_depth_m);
    arm_route(Fixture::prepare(seg, PatchStatus::Active), seed);
    Fixture::route(seg, cfg.time.bio_dt);
    const Agent& inserted = only_cell(seg, false);
    check_identity(inserted, source);
    assert(seg.ledger().reseeds == 1 && seg.ledger().distal_losses == 0);
    // Report the landing predicate without requiring today's invalid z or
    // prescribing a replacement height. Observe it before the next live step.
    std::cout << "packet_seed=" << seed << " target=0 status=active tag="
              << source.identity.tag << " source=(" << source.x[0] << ','
              << source.x[1] << ',' << source.x[2] << ") inserted=("
              << inserted.x[0] << ',' << inserted.x[1] << ',' << inserted.x[2]
              << ") local=" << local_position(inserted.x, cfg.layer2) << '\n';
    const Int before = seg.ledger().washout_departures;
    Fixture::advance_live(seg, cfg.time.bio_dt);
    const bool reexported = !seg.pool().empty();
    check_identity(only_cell(seg, reexported), source);
    assert(seg.ledger().washout_departures - before == (reexported ? 1 : 0));
    std::cout << "next_live_phase=" << (reexported ? "washout" : "resident") << '\n';
  }
}
}  // namespace

int main(int argc, char** argv) {
  const std::array<testshards::NamedTest, 2> tests = {{
      {"eligibility", eligibility}, {"boundary_witness", boundary_witness}}};
  return testshards::shard_main({argv, static_cast<std::size_t>(argc)}, tests,
                                "Single-slot re-entry\n", "Transfers checked.\n");
}
