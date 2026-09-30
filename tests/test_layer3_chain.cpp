/* -----------------------------------------------------------------------
   Spec 13 Phase 3 — Layer-3 colonic-chain unit coverage.

   Required by the Phase-3 validation gate:
     (i)   population-ledger closure across the mucosa/lumen boundary
           and across region boundaries
     (ii)  seed-fixed invariance under region/lineage iteration order
     (iii) checkpoint round-trip of luminal counts + genotypes +
           per-region state
     (iv)  genotype-preserving count->agent instantiation on
           reattachment
   ----------------------------------------------------------------------- */

#include "agent.h"
#include "chain.h"
#include "input_parser.h"
#include "layer3_config.h"
#include "lumen.h"
#include "path_utils.h"
#include "segment.h"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using namespace gutibm;

namespace {

// Deposit tags minted by the chain carry the 0xFF marker in bits
// 48-55; seed exemplar lineage ids carry 0xFD there. See chain.cpp.
constexpr TagID kDepositMarkerMask = static_cast<TagID>(0xFF) << 48;
constexpr TagID kSeedLineageTag =
    (static_cast<TagID>(0xFD) << 48) | 1;  // strain index 0 + 1

bool is_deposit_tag(TagID t) {
  return (t & kDepositMarkerMask) == kDepositMarkerMask;
}

SimulationConfig base_cfg(bool uniform = false) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.hdf5.enabled = false;
  cfg.enabled_fixes = {};
  cfg.time.bio_dt = 60.0;
  cfg.time.total_time = 3600.0;
  cfg.seed = 777;
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 100;
  strain.mu_max = 5.0e-4;
  strain.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(strain);

  auto& l = cfg.layer2;
  l.n_patches = 40;  // overridden per-region by layer3
  l.initial_patch_fraction = 0.5;
  l.initial_founder_cells = 20;
  l.contraction_rate_per_min = 0.5;
  l.agent_loss_fraction = 0.3;
  // Mucus-turnover channel active but slow enough that reseeds
  // replenish it (the run must survive all 200 steps).
  l.crypt = {0.2, 0.10, 0.5 / 3600.0, 2.0};
  l.exposed_proximal = {0.5, 0.5, 1.0 / 3600.0, 1.0};
  l.exposed_distal = {0.3, 0.5, 2.0 / 3600.0, 0.5};
  l.transit_half_life_s = 600.0;
  l.reattach_prob_per_transit = 0.5;
  l.establish_prob_single = 0.2;
  l.establish_ratio = 10.0;
  l.invalid_patch_fraction_stop = 1.0;

  auto& c = cfg.layer3;
  c.enabled = true;
  c.uniform_profile = uniform;
  c.n_patches_per_region = 40;
  c.f_edge = 0.2;
  c.reattach_p0 = 0.5;
  // Fast luminal transit so outflow/stool channels fire every step.
  c.cecum.transit_h = 0.1;
  c.transverse.transit_h = 0.08;
  c.descending.transit_h = 0.12;
  c.cecum.outflow_survival = 0.90;
  c.transverse.outflow_survival = 0.85;
  c.descending.outflow_survival = 0.95;
  c.cecum.contraction_k0_per_min = 0.5;
  c.transverse.contraction_k0_per_min = 0.4;
  c.descending.contraction_k0_per_min = 0.3;
  c.lumen_seed_cells = 3000.0;
  c.hapc_rate_per_day = 20.0;
  c.alpha_hapc_per_min = 0.5;
  c.purge_at_s = 600.0;         // exercise the purge channel once
  c.purge_fraction = 0.9;
  return cfg;
}

// ── (i) Ledger closure across mucosa/lumen + region boundaries ────────
// Chain stocks == initial + seg_births + luminal_births - stool -
// deaths - purged at every step; internal transfers (mucosal
// departures, pool exits, deposits, downstream arrivals) cancel.
void test_ledger_closure() {
  SimulationConfig cfg = base_cfg();
  ColonicChain chain;
  chain.init(cfg);
  const Real dt = cfg.time.bio_dt;
  assert(std::abs(chain.stocks() - chain.ledger().initial_cells) <
         0.5);

  for (int step = 0; step < 200; ++step) {
    chain.step(dt);
    assert(chain.termination_cause() ==
           TerminationCause::IncompleteUnknown);
    assert(chain.ledger_closed());
  }
  // Every channel must have fired for this to be a real test.
  const ChainLedger& led = chain.ledger();
  assert(led.mucosal_births > 0);
  assert(led.luminal_births > 0);
  assert(led.stool_exports > 0);
  assert(led.luminal_deaths > 0);
  assert(led.purged > 0);
  Int external = 0;
  Int departures = 0;
  for (Int r = 0; r < kRegionCount; ++r) {
    const ColonicChain::RegionView rv = chain.region(r);
    external += rv.segment.ledger().external_arrivals;
    departures += rv.segment.ledger().distal_losses +
                  rv.segment.ledger().contraction_departures;
  }
  assert(external > 0);   // lumen -> mucosa reattachment happened
  assert(departures > 0); // mucosa -> lumen transfers happened
  std::cout << "ledger_closure: stocks=" << chain.stocks()
            << " expected=" << led.expected_stocks()
            << " stool=" << led.stool_exports
            << " deaths=" << led.luminal_deaths
            << " purged=" << led.purged
            << " external=" << external << '\n';
}

// ── (ii) Seed-fixed invariance under region/lineage iteration order ───
// Lineage order is key-sorted by construction; region processing order
// is permutable by the test hook. Same seed + permuted order => same
// trajectory fingerprint.
void test_order_invariance() {
  const Real dt = 60.0;
  const int steps = 60;

  auto run_ordered = [&](const std::vector<Int>& order) {
    SimulationConfig cfg = base_cfg();
    ColonicChain chain;
    chain.init(cfg);
    chain.permute_region_order_for_testing(order);
    for (int i = 0; i < steps; ++i) {
      chain.step(dt);
      assert(chain.termination_cause() ==
             TerminationCause::IncompleteUnknown);
    }
    return chain.fingerprint();
  };

  const uint64_t forward = run_ordered({0, 1, 2});
  const uint64_t reversed = run_ordered({2, 1, 0});
  const uint64_t rotated = run_ordered({2, 0, 1});
  const uint64_t repeat = run_ordered({0, 1, 2});
  assert(forward == reversed);
  assert(forward == rotated);
  assert(forward == repeat);
  std::cout << "order_invariance: fingerprint=" << forward << '\n';
}

// ── (iv) Genotype-preserving count->agent reattachment ─────────────────
// Chain-minted deposits carry a deposit marker tag. Every deposit
// landing in mucosa must carry the lineage exemplar's genotype: seed
// lineage id, the founder's BI loci, and parent_id == exemplar tag.
void test_genotype_preserving_reattach() {
  SimulationConfig cfg = base_cfg();
  cfg.layer2.initial_patch_fraction = 0.1;
  cfg.layer3.reattach_p0 = 0.9;
  ColonicChain chain;
  chain.init(cfg);
  const Real dt = cfg.time.bio_dt;
  for (int i = 0; i < 100; ++i) {
    chain.step(dt);
    assert(chain.termination_cause() ==
           TerminationCause::IncompleteUnknown);
  }

  Int deposits = 0;
  Int genotype_ok = 0;
  for (Int r = 0; r < kRegionCount; ++r) {
    chain.region(r).segment.for_each_colonist([&](const Agent& a) {
      if (!is_deposit_tag(a.identity.tag)) return;
      ++deposits;
      if (a.genome.lineage_id == kSeedLineageTag &&
          a.genome.parent_id == kSeedLineageTag &&
          a.genome.bi_loci.size() == 1) {
        ++genotype_ok;
      }
    });
  }
  // With lumen_seed_cells = 3000 and p0 = 0.9 the seed lineage must have
  // deposited intact at least once in 100 steps.
  assert(deposits > 0);
  assert(genotype_ok > 0);
  std::cout << "genotype_reattach: deposits=" << deposits
            << " genotype_ok=" << genotype_ok << '\n';
}

// ── (iii) Checkpoint round-trip ────────────────────────────────────────
#ifdef GUTIBM_HDF5
void test_checkpoint_roundtrip() {
  const std::string path =
      resolve_test_h5_path("GUTIBM_LAYER3_CKPT_H5", "layer3_chain");
  const Real dt = 60.0;

  SimulationConfig cfg = base_cfg();
  ColonicChain chain;
  chain.init(cfg);
  for (int i = 0; i < 50; ++i) {
    chain.step(dt);
  }
  chain.write_checkpoint(path);

  // State comparison before continuing: luminal counts, lineage ids,
  // exemplar genotype sizes, per-region state.
  ColonicChain restored;
  restored.init_from_checkpoint(cfg, path);
  for (Int r = 0; r < kRegionCount; ++r) {
    const auto& a = chain.region(r).lumen.lineages();
    const auto& b = restored.region(r).lumen.lineages();
    assert(a.size() == b.size());
    for (const auto& [id, lin] : a) {
      const auto it = b.find(id);
      assert(it != b.end());
      assert(it->second.count == lin.count);
      assert(it->second.exemplar.genome.lineage_id ==
             lin.exemplar.genome.lineage_id);
      assert(it->second.exemplar.genome.bi_loci.size() ==
             lin.exemplar.genome.bi_loci.size());
    }
    assert(restored.region(r).lumen.scalars().carbon_mol_m3 ==
           chain.region(r).lumen.scalars().carbon_mol_m3);
  }
  assert(restored.step_count() == chain.step_count());
  assert(restored.fingerprint() == chain.fingerprint());

  // Trajectory equivalence after resume.
  for (int i = 0; i < 50; ++i) {
    chain.step(dt);
    restored.step(dt);
    assert(chain.termination_cause() ==
           TerminationCause::IncompleteUnknown);
    assert(restored.termination_cause() ==
           TerminationCause::IncompleteUnknown);
  }
  assert(restored.fingerprint() == chain.fingerprint());
  std::cout << "checkpoint_roundtrip: fingerprint="
            << restored.fingerprint() << '\n';
}
#endif  // GUTIBM_HDF5

}  // namespace

int main() {
  test_ledger_closure();
  test_order_invariance();
  test_genotype_preserving_reattach();
#ifdef GUTIBM_HDF5
  test_checkpoint_roundtrip();
#endif
  std::cout << "test_layer3_chain: all tests passed\n";
  return 0;
}
