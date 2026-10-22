/* -----------------------------------------------------------------------
   Spec 13 Phase 1 — Layer-2 mucus-segment unit coverage.

   Required by the Phase-1 validation gate:
     (i)   population-ledger closure across patches + pool
     (ii)  seed-fixed invariance under patch-order permutation (G3)
     (iii) checkpoint round-trip of pool + occupancy + patch identity
     (iv)  single-debit attribution per departing agent

   Plus bounds/invariant coverage per the observables contract (G5/G6).
   ----------------------------------------------------------------------- */

#include "error.h"
#include "founder_genotype.h"
#include "genotype_test_helpers.h"
#include "input_parser.h"
#include "layer2_config.h"
#include "layer3_config.h"
#include "lumen.h"
#include "mechanism_capability.h"
#include "path_utils.h"
#include "plasmid.h"
#include "segment.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using namespace gutibm;

namespace {

SimulationConfig base_cfg() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.hdf5.enabled = false;
  cfg.enabled_fixes = {};
  cfg.time.bio_dt = 60.0;
  cfg.time.total_time = 3600.0;
  cfg.seed = 1234;
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 100;
  strain.mu_max = 5.0e-4;
  strain.plasmids = {"ColE1"};
  cfg.initial_strains.push_back(strain);

  auto& l = cfg.layer2;
  l.enabled = true;
  l.n_patches = 40;
  l.initial_patch_fraction = 0.5;
  l.initial_founder_cells = 20;
  l.contraction_rate_per_min = 0.5;
  l.agent_loss_fraction = 0.3;
  l.crypt = {0.2, 0.10, 5.6e-6, 2.0};
  l.exposed_proximal = {0.5, 0.5, 2.8e-5, 1.0};
  l.exposed_distal = {0.3, 0.5, 5.6e-5, 0.5};
  l.transit_half_life_s = 3600.0;
  l.reattach_prob_per_transit = 0.5;
  l.establish_prob_single = 0.2;
  l.establish_ratio = 10.0;
  l.invalid_patch_fraction_stop = 1.0;  // never stop mid-test on invalidity
  return cfg;
}

// Drive `n` steps; returns false if the run terminated early.
bool run_steps(MucusSegment& seg, Real dt, int n) {
  for (int i = 0; i < n; ++i) {
    seg.step(dt);
    if (seg.termination_cause() != TerminationCause::IncompleteUnknown) {
      return false;
    }
  }
  return true;
}

// ── (i) Ledger closure across patches + pool ─────────────────────────────
// Every cell is patch resident, in the pool, exported distally, or never
// existed: stocks == initial + births - distal at every step.
void test_ledger_closure() {
  SimulationConfig cfg = base_cfg();
  MucusSegment seg;
  seg.init(cfg);
  const Real dt = cfg.time.bio_dt;

  Int patch_before = 0;
  Int pool_before = 0;
  for (const auto& p : seg.patches()) patch_before += static_cast<Int>(p.colonists.size());
  for (const auto& q : seg.pool()) pool_before += static_cast<Int>(q.cells.size());
  assert(patch_before + pool_before ==
         seg.ledger().initial_cells);

  for (int step = 0; step < 200; ++step) {
    const Int births_before = seg.ledger().births;
    const Int wash_before = seg.ledger().washout_departures;
    const Int con_before = seg.ledger().contraction_departures;
    const Int res_before = seg.ledger().reseeds;
    const Int dist_before = seg.ledger().distal_losses;

    seg.step(dt);
    assert(seg.termination_cause() == TerminationCause::IncompleteUnknown);

    Int patches = 0;
    for (const auto& p : seg.patches()) patches += static_cast<Int>(p.colonists.size());
    Int pool = 0;
    for (const auto& q : seg.pool()) pool += static_cast<Int>(q.cells.size());

    // Closure: total stocks track births minus the only true sink.
    assert(patches + pool ==
           seg.ledger().initial_cells + seg.ledger().births -
               seg.ledger().distal_losses);
    // Per-step patch delta: births enter, both loss channels and the
    // single reseed direction move between patch and pool stocks.
    const Int departures =
        (seg.ledger().washout_departures - wash_before) +
        (seg.ledger().contraction_departures - con_before);
    const Int arrivals = seg.ledger().reseeds - res_before;
    assert(patches - patch_before ==
           (seg.ledger().births - births_before) - departures + arrivals);
    const Int exits = seg.ledger().distal_losses - dist_before;
    assert(pool - pool_before == departures - arrivals - exits);
    patch_before = patches;
    pool_before = pool;
  }
  // All four channels must have been exercised to make this a real test.
  assert(seg.ledger().births > 0);
  assert(seg.ledger().contraction_departures > 0);
  assert(seg.ledger().washout_departures > 0);
  assert(seg.ledger().distal_losses > 0);
  std::cout << "ledger_closure: births=" << seg.ledger().births
            << " washout=" << seg.ledger().washout_departures
            << " contraction=" << seg.ledger().contraction_departures
            << " reseeds=" << seg.ledger().reseeds
            << " distal=" << seg.ledger().distal_losses << '\n';
}

// ── (iv) Single-debit attribution ────────────────────────────────────────
// With both channels firing on the same patch in the same step, the
// number of cells leaving that patch equals contraction + washout
// departures — each departing agent is debited exactly once.
void test_single_debit_attribution() {
  SimulationConfig cfg = base_cfg();
  auto& l = cfg.layer2;
  l.n_patches = 12;
  l.initial_patch_fraction = 1.0;
  l.initial_founder_cells = 50;
  l.contraction_rate_per_min = 60.0;  // fires every step at dt=60 s
  l.agent_loss_fraction = 0.4;
  l.crypt = {1.0, 0.5, 5.6e-4, 2.0};   // all crypt: high single-cell loss
  l.exposed_proximal = {0.0, 0.5, 2.8e-5, 1.0};
  l.exposed_distal = {0.0, 0.5, 5.6e-5, 0.5};
  l.reattach_prob_per_transit = 0.0;   // isolate departures
  l.establish_prob_single = 0.0;

  MucusSegment seg;
  seg.init(cfg);
  const Real dt = cfg.time.bio_dt;

  for (int step = 0; step < 30; ++step) {
    std::vector<Int> before;
    for (const auto& p : seg.patches()) {
      before.push_back(static_cast<Int>(p.colonists.size()));
    }
    const Int wash_before = seg.ledger().washout_departures;
    const Int con_before = seg.ledger().contraction_departures;
    const Int births_before = seg.ledger().births;
    const Int res_before = seg.ledger().reseeds;

    seg.step(dt);

    const Int departures =
        (seg.ledger().washout_departures - wash_before) +
        (seg.ledger().contraction_departures - con_before);
    const Int births = seg.ledger().births - births_before;
    const Int arrivals = seg.ledger().reseeds - res_before;

    // Patch-level single debit: departed cells = ledger increments.
    for (size_t i = 0; i < seg.patches().size(); ++i) {
      const Int n_now =
          static_cast<Int>(seg.patches()[i].colonists.size());
      // Each patch's population can only drop by what left it.
      assert(n_now <= before[i] + births + arrivals);
    }
    // births/arrivals are spread over patches; the aggregate identity is
    // what must hold exactly:  delta(patch_cells) = births + reseeds - losses
    Int patch_now = 0;
    for (const auto& p : seg.patches()) {
      patch_now += static_cast<Int>(p.colonists.size());
    }
    Int patch_then = 0;
    for (const Int v : before) patch_then += v;
    assert(patch_now == patch_then + births + arrivals - departures);
  }
  assert(seg.ledger().contraction_departures > 0);
  assert(seg.ledger().washout_departures > 0);
  std::cout << "single_debit: washout=" << seg.ledger().washout_departures
            << " contraction=" << seg.ledger().contraction_departures
            << '\n';
}

// ── (ii) Seed-fixed invariance under patch-order permutation (G3) ────────
void test_seed_invariance_under_permutation() {
  SimulationConfig cfg = base_cfg();

  MucusSegment seg_a;
  seg_a.init(cfg);
  const bool ok_a = run_steps(seg_a, cfg.time.bio_dt, 120);
  assert(ok_a);

  MucusSegment seg_b;
  seg_b.init(cfg);
  // Reverse permutation: every patch's draws must be keyed to its own
  // stream, so the trajectory is identical under any iteration order.
  std::vector<Int> rev(static_cast<size_t>(cfg.layer2.n_patches));
  std::iota(rev.begin(), rev.end(), 0);
  std::reverse(rev.begin(), rev.end());
  seg_b.permute_patch_order_for_testing(rev);
  const bool ok_b = run_steps(seg_b, cfg.time.bio_dt, 120);
  assert(ok_b);

  assert(seg_a.fingerprint() == seg_b.fingerprint());
  std::cout << "seed_invariance: fingerprint=" << seg_a.fingerprint()
            << '\n';
}

// ── (iii) Checkpoint round-trip: pool + occupancy + patch identity ───────
void test_checkpoint_roundtrip() {
#ifdef GUTIBM_HDF5
  const std::string path =
      resolve_test_h5_path("GUTIBM_LAYER2_CKPT_H5", "layer2_segment");

  SimulationConfig cfg = base_cfg();
  const Real dt = cfg.time.bio_dt;

  MucusSegment seg_a;
  seg_a.init(cfg);
  assert(run_steps(seg_a, dt, 90));
  seg_a.write_checkpoint(path);

  MucusSegment seg_b;
  seg_b.init_from_checkpoint(cfg, path);
  assert(run_steps(seg_b, dt, 90));

  MucusSegment seg_ref;
  seg_ref.init(cfg);
  assert(run_steps(seg_ref, dt, 180));

  assert(seg_b.fingerprint() == seg_ref.fingerprint());
  // Patch identity persisted: same type/status assignment.
  for (size_t i = 0; i < seg_b.patches().size(); ++i) {
    assert(seg_b.patches()[i].type == seg_ref.patches()[i].type);
    assert(seg_b.patches()[i].id == seg_ref.patches()[i].id);
  }
  // Ledger restores identically.
  assert(seg_b.ledger().initial_cells ==
         seg_ref.ledger().initial_cells);
  assert(seg_b.ledger().births == seg_ref.ledger().births);
  assert(seg_b.ledger().washout_departures ==
         seg_ref.ledger().washout_departures);
  assert(seg_b.ledger().contraction_departures ==
         seg_ref.ledger().contraction_departures);
  assert(seg_b.ledger().reseeds == seg_ref.ledger().reseeds);
  assert(seg_b.ledger().distal_losses ==
         seg_ref.ledger().distal_losses);
  std::cout << "checkpoint_roundtrip: fingerprint=" << seg_b.fingerprint()
            << " pool=" << seg_b.pool().size() << '\n';
#else
  std::cout << "checkpoint_roundtrip: skipped (no HDF5 build)\n";
#endif
}

// ── Bounds and observable invariants (G5/G6) ─────────────────────────────
void test_observables_bounds() {
  SimulationConfig cfg = base_cfg();
  MucusSegment seg;
  seg.init(cfg);
  const Real dt = cfg.time.bio_dt;
  for (int i = 0; i < 60; ++i) {
    seg.step(dt);
    const auto o = seg.observables();
    for (Int t = 0; t < kPatchTypeCount; ++t) {
      const auto& to = o.per_type[t];
      assert(to.occupancy >= 0.0 && to.occupancy <= 1.0);
      assert(to.n_occupied <= to.n_patches);
      assert(std::isfinite(to.mean_density_occupied));
    }
    assert(std::isfinite(o.segment_mean_cfu_ml));
    assert(o.segment_mean_cfu_ml >= 0.0);
    // G6 identity: CFU/cm^2 = cells/mL * thickness(cm).
    const Real expected =
        o.segment_mean_cfu_ml * cfg.layer2.mucus_thickness_m * 100.0;
    assert(std::abs(o.segment_cfu_cm2 - expected) < 1e-9 * (expected + 1));
    assert(seg.ledger_closed());
  }
  std::cout << "observables_bounds: mean="
            << seg.observables().segment_mean_cfu_ml << '\n';
}

// ── Config parsing ────────────────────────────────────────────────────────
void test_layer2_config_parsing() {
  const std::string path =
      std::string(GUTIBM_SOURCE_DIR) + "/tests/fixtures/layer2_config.json";
  SimulationConfig cfg = InputParser::parse(path);
  assert(cfg.layer2.enabled);
  assert(cfg.layer2.n_patches == 64);
  assert(std::abs(cfg.layer2.contraction_rate_per_min - 0.7) < 1e-12);
  assert(std::abs(cfg.layer2.crypt.fraction - 0.1) < 1e-12);
  assert(std::abs(cfg.layer2.crypt.disruption_prob - 0.08) < 1e-12);
  assert(std::abs(cfg.layer2.crypt.single_cell_loss_per_s -
                  0.02 / 3600.0) < 1e-15);
  assert(std::abs(cfg.layer2.exposed_proximal.supply_mult - 1.5) < 1e-12);
  assert(std::abs(cfg.layer2.establish_ratio - 25.0) < 1e-12);
  std::cout << "config_parsing: n_patches=" << cfg.layer2.n_patches
            << " ratio=" << cfg.layer2.establish_ratio << '\n';
}

// ── Patch-level guards do not kill the run (S5) ───────────────────────────
void test_bloom_halt_is_patch_local() {
  SimulationConfig cfg = base_cfg();
  auto& l = cfg.layer2;
  l.n_patches = 8;
  l.initial_patch_fraction = 1.0;
  l.initial_founder_cells = 5;
  l.initial_pool_cells = 2000;     // heavy inflow drives a bloom
  l.transit_half_life_s = 600.0;   // frequent transit draws
  l.reattach_prob_per_transit = 1.0;
  l.establish_prob_single = 1.0;
  l.bloom_factor = 1.5;            // bound = 7 cells per seeded patch
  l.bloom_sustain_s = 0.0;
  l.contraction_rate_per_min = 0.0;
  l.invalid_patch_fraction_stop = 1.0;

  MucusSegment seg;
  seg.init(cfg);
  bool any_halted = false;
  for (int i = 0; i < 60 && !any_halted; ++i) {
    seg.step(cfg.time.bio_dt);
    for (const auto& p : seg.patches()) {
      if (p.status != PatchStatus::Active) any_halted = true;
    }
  }
  assert(any_halted);
  assert(seg.termination_cause() ==
         TerminationCause::IncompleteUnknown);  // run survives
  assert(seg.ledger_closed());
  std::cout << "bloom_halt: run survived with halted patches\n";
}

// ── Founder genotype propagation (external audit, step 1) ────────────
// The declared resistant phenotype, configured plasmid cluster and CDI
// fields must reach every construction site — patch founders AND the
// initial pool — not just the live-Simulation path.
void test_founder_genotype_propagation() {
  SimulationConfig cfg = base_cfg();
  cfg.initial_strains.clear();

  SimulationConfig::InitialStrain s_strain;  // sensitive, no plasmids
  s_strain.type = 1;
  s_strain.count = 50;
  s_strain.mu_max = 5.0e-4;

  SimulationConfig::InitialStrain r_strain;  // declared R phenocopy
  r_strain.type = 2;
  r_strain.count = 50;
  r_strain.mu_max = 4.0e-4;
  r_strain.plasmids = {"ColE1"};
  test::declare_resistant_founder(r_strain);
  cfg.initial_strains = {s_strain, r_strain};

  auto& l = cfg.layer2;
  l.n_patches = 8;
  l.initial_patch_fraction = 1.0;  // every patch seeded
  l.initial_founder_cells = 8;     // proportional 4 + 4 split
  l.initial_pool_cells = 4;        // 2 + 2 pool singles

  MucusSegment seg;
  seg.init(cfg);

  // The cluster founders must carry pI-derived retardation (the
  // configured_plasmid_cluster path), not the raw library entry.
  const Real expected_ret = test::declared_cluster_retardation(cfg);
  const Int btub = to_underlying(ReceptorType::BtuB);

  int r_colonists = 0;
  int s_colonists = 0;
  seg.for_each_colonist([&r_colonists, &s_colonists, btub,
                         expected_ret](const Agent& a) {
    if (a.identity.type == 2) {
      ++r_colonists;
      test::expect_declared_r_phenotype(a, btub);
      test::expect_configured_cluster(a, expected_ret);
      assert(a.genome.lineage_id == a.identity.tag);
    } else {
      ++s_colonists;
      assert(a.state != PhenoState::RESISTANT);
      assert(a.genome.cdi_type == 0);
    }
  });
  assert(r_colonists > 0);
  assert(s_colonists > 0);

  int r_pool = 0;
  int s_pool = 0;
  for (const auto& q : seg.pool()) {
    for (const auto& a : q.cells) {
      if (a.identity.type == 2) {
        ++r_pool;
        test::expect_declared_r_phenotype(a, btub);
        test::expect_configured_cluster(a, expected_ret);
        assert(a.genome.lineage_id == a.identity.tag);
      } else {
        ++s_pool;
      }
    }
  }
  assert(r_pool > 0);
  assert(s_pool > 0);
  std::cout << "founder_genotype: r_colonists=" << r_colonists
            << " s_colonists=" << s_colonists << " r_pool=" << r_pool
            << " s_pool=" << s_pool << '\n';
}

// plasmid_overrides change the carried cluster through
// configured_plasmid_cluster — the same semantics the live path applies.
void test_founder_genotype_helper_contract() {
  SimulationConfig cfg = base_cfg();
  auto& strain = cfg.initial_strains.front();
  test::declare_resistant_founder(strain);
  cfg.plasmid_overrides["ColE1"].retardation = 1.2;
  cfg.plasmid_overrides["ColE1"].burst_size = 42.0;

  Agent a = Agent::create_default(1001, strain.type, Vec3{0, 0, 0},
                                  strain.mu_max);
  apply_founder_genotype(a, strain, cfg);
  const Int btub = to_underlying(ReceptorType::BtuB);
  test::expect_declared_r_phenotype(a, btub);
  assert(a.genome.bi_loci.size() == 1);
  assert(a.genome.bi_loci[0].retardation == 1.2);
  assert(a.genome.bi_loci[0].burst_size == 42.0);
  std::cout << "founder_genotype_helper: overrides applied\n";
}

// Declared phenotype survives checkpoint/restart via the agent_transfer
// blobs that already serialize patch + pool cells.
void test_genotype_survives_checkpoint() {
#ifdef GUTIBM_HDF5
  const std::string path =
      resolve_test_h5_path("GUTIBM_LAYER2_CKPT_H5", "layer2_genotype");

  SimulationConfig cfg = base_cfg();
  auto& strain = cfg.initial_strains.front();
  test::declare_resistant_founder(strain);
  cfg.layer2.n_patches = 8;
  cfg.layer2.initial_patch_fraction = 1.0;
  cfg.layer2.initial_founder_cells = 8;
  cfg.layer2.initial_pool_cells = 4;

  MucusSegment seg_a;
  seg_a.init(cfg);
  seg_a.write_checkpoint(path);

  MucusSegment seg_b;
  seg_b.init_from_checkpoint(cfg, path);

  const Int btub = to_underlying(ReceptorType::BtuB);
  int restored = 0;
  seg_b.for_each_colonist([&restored, btub](const Agent& a) {
    ++restored;
    test::expect_declared_r_phenotype(a, btub);
    assert(a.genome.bi_loci.size() == 1);
  });
  assert(restored > 0);
  int pool_restored = 0;
  for (const auto& q : seg_b.pool()) {
    for (const auto& a : q.cells) {
      ++pool_restored;
      test::expect_declared_r_phenotype(a, btub);
    }
  }
  assert(pool_restored > 0);
  std::cout << "genotype_checkpoint: colonists=" << restored
            << " pool=" << pool_restored << '\n';
#else
  std::cout << "genotype_checkpoint: skipped (no HDF5 build)\n";
#endif
}

// ── Capability limit (external audit, step 1) ────────────────────────
// A tabulated-layer config that sets mechanism-dynamics parameters is a
// sweep of unused parameters: validate_config must refuse it and name
// the unsupported keys rather than produce identical outputs.
void test_capability_limit_reports_dead_params() {
  // A clean tabulated config passes, and founder genotype declarations
  // are carried state — they must NOT be flagged as mechanism knobs.
  {
    SimulationConfig cfg = base_cfg();
    cfg.initial_strains[0].plasmids = {"ColE1"};
    test::declare_resistant_founder(cfg.initial_strains[0]);
    assert(unexecuted_mechanism_params(cfg).empty());
    MucusSegment seg;
    seg.init(cfg);
  }

  const int refused = test::expect_dead_knobs_refused(
      {
          {"sos_lysis_prob sweep", "fixes.bacteriocin.sos_lysis_prob",
           test::knob_sos_lysis},
          {"retardation override sweep", "plasmid_overrides",
           test::knob_retardation_override},
          {"fix list", "enabled_fixes", test::knob_fix_list},
          {"receptor affinity", "fixes.receptor.kd_b12_btuB",
           [](SimulationConfig& c) {
             c.fixes.receptor.kd_b12_btuB = 1e-3;
           }},
          {"ROS mortality", "chem_env.oxygen.k_ROS",
           [](SimulationConfig& c) { c.chem_env.oxygen.k_ROS = 1e-6; }},
          {"immigration", "immigration.enabled",
           [](SimulationConfig& c) { c.immigration.enabled = true; }},
          {"motility knob", "cell_bio.motility.enabled",
           [](SimulationConfig& c) {
             c.cell_bio.motility.enabled = false;
           }},
          {"toxin lumping", "qssa.toxin_lumping",
           [](SimulationConfig& c) { c.qssa.toxin_lumping = "lumped"; }},
      },
      base_cfg(),
      [](const SimulationConfig& c) {
        MucusSegment seg;
        seg.init(c);
      });
  std::cout << "capability_limit: " << refused << " dead knobs refused\n";
}

// ── Live audit slot (Edison step 2): mechanism-capable patch ────────────
// One small-domain patch substituted by a real Simulation so each step
// is fast enough for unit coverage.
SimulationConfig live_cfg_for_test() {
  SimulationConfig cfg = base_cfg();
  auto& l = cfg.layer2;
  l.n_patches = 1;
  l.initial_patch_fraction = 1.0;
  l.initial_founder_cells = 20;
  l.patch_lateral_m = 60.0e-6;
  l.patch_depth_m = 60.0e-6;
  l.contraction_rate_per_min = 0.0;
  l.transit_half_life_s = 600.0;
  l.reattach_prob_per_transit = 0.0;
  l.establish_prob_single = 0.0;
  l.audit.patch_index = 0;
  l.audit.halt_on_exceed = false;
  return cfg;
}

// Biological deaths inside the live engine book into
// ledger_.biological_deaths: with basal SOS lysis forcing kills while
// contraction harvests fragments, closure must hold at every step —
// the same cells cannot be counted as both death and departure.
void test_live_patch_death_accounting() {
  SimulationConfig cfg = live_cfg_for_test();
  cfg.layer2.contraction_rate_per_min = 60.0;  // fragment harvest every step
  cfg.fixes.bacteriocin.sos_basal_rate = 0.02; // ~70% SOS induction per step

  MucusSegment seg;
  seg.init(cfg);
  assert(seg.patches()[0].live);

  Int colonists = 0;
  seg.for_each_colonist([&colonists](const Agent&) { ++colonists; });
  assert(colonists == 20);  // founders mirrored into the live engine

  for (int i = 0; i < 15; ++i) {
    seg.step(cfg.time.bio_dt);
    // Every cell lysing or departing early is a legitimate stop; the
    // ledger must still have closed at the final step.
    assert(seg.termination_cause() ==
               TerminationCause::IncompleteUnknown ||
           seg.termination_cause() == TerminationCause::PopulationStop);
    assert(seg.ledger_closed());
    // Explicit conservation identity: cells in == cells out + deaths.
    Int stocks = 0;
    seg.for_each_colonist([&stocks](const Agent&) { ++stocks; });
    for (const auto& q : seg.pool()) {
      stocks += static_cast<Int>(q.cells.size());
    }
    const auto& L = seg.ledger();
    assert(stocks + L.distal_losses + L.biological_deaths ==
           L.initial_cells + L.births + L.external_arrivals);
    if (seg.termination_cause() != TerminationCause::IncompleteUnknown) {
      break;
    }
  }
  assert(seg.ledger().biological_deaths > 0);
  assert(seg.ledger().contraction_departures > 0);
  std::cout << "live_death_accounting: deaths="
            << seg.ledger().biological_deaths
            << " contraction=" << seg.ledger().contraction_departures
            << " births=" << seg.ledger().births << '\n';
}

// Influx: a pre-seeded luminal pool reattaching into the (only, live)
// patch lands inside the embedded Simulation, not the colonist vector.
void test_live_patch_pool_influx() {
  SimulationConfig cfg = live_cfg_for_test();
  auto& l = cfg.layer2;
  l.initial_pool_cells = 60;
  l.transit_half_life_s = 300.0;
  l.reattach_prob_per_transit = 1.0;
  l.establish_prob_single = 1.0;

  MucusSegment seg;
  seg.init(cfg);
  for (int i = 0; i < 60 && seg.ledger().reseeds == 0; ++i) {
    seg.step(cfg.time.bio_dt);
    assert(seg.termination_cause() ==
               TerminationCause::IncompleteUnknown ||
           seg.termination_cause() == TerminationCause::PopulationStop);
    assert(seg.ledger_closed());
  }
  assert(seg.ledger().reseeds > 0);
  // Reseeded cells are visible through the live-aware colonist view.
  Int colonists = 0;
  seg.for_each_colonist([&colonists](const Agent&) { ++colonists; });
  assert(colonists > 20);
  std::cout << "live_pool_influx: reseeds=" << seg.ledger().reseeds
            << " colonists=" << colonists << '\n';
}

// Efflux + handoff: pool distal exits routed through the pool-exit sink
// land in a LumenCompartment keyed by founder lineage_id — the minimal
// segment->chain handoff the Layer-3 driver wires identically.
void test_live_patch_pool_handoff() {
  SimulationConfig cfg = live_cfg_for_test();
  auto& l = cfg.layer2;
  l.n_patches = 2;
  l.crypt = {0.0, 0.0, 0.0, 0.0};
  l.exposed_proximal = {1.0, 0.9, 1e-3, 1.0};  // heavy single-cell washout
  l.exposed_distal = {0.0, 0.9, 1e-3, 1.0};
  l.transit_half_life_s = 300.0;
  l.reattach_prob_per_transit = 0.0;

  Layer3RegionConfig rc;
  LumenCompartment lumen;
  lumen.init(rc, {}, 0.0);
  MucusSegment seg;
  seg.init(cfg);
  seg.set_pool_exit_sink(
      [&lumen](std::vector<Agent>&& cells) { lumen.add_cells(cells); });

  for (int i = 0; i < 60 && seg.ledger().distal_losses == 0; ++i) {
    seg.step(cfg.time.bio_dt);
    assert(seg.termination_cause() ==
               TerminationCause::IncompleteUnknown ||
           seg.termination_cause() == TerminationCause::PopulationStop);
    assert(seg.ledger_closed());
  }
  assert(seg.ledger().distal_losses > 0);
  assert(!lumen.lineages().empty());
  // Every delivered cell resolves to its founder lineage.
  Real total = 0.0;
  for (const auto& [lin, entry] : lumen.lineages()) {
    assert(lin == entry.exemplar.genome.lineage_id);
    total += entry.count;
  }
  assert(std::abs(total - static_cast<Real>(seg.ledger().distal_losses)) <
         1e-9);
  std::cout << "live_pool_handoff: distal=" << seg.ledger().distal_losses
            << " lumen_lineages=" << lumen.lineages().size() << '\n';
}

// Seed invariance for the live slot: two segments with the same config
// produce identical live-patch trajectories (fingerprint covers live
// agent tags + lineage ids now).
void test_live_patch_seed_invariance() {
  SimulationConfig cfg = live_cfg_for_test();
  cfg.layer2.contraction_rate_per_min = 5.0;

  MucusSegment seg_a;
  seg_a.init(cfg);
  MucusSegment seg_b;
  seg_b.init(cfg);
  for (int i = 0; i < 10; ++i) {
    seg_a.step(cfg.time.bio_dt);
    seg_b.step(cfg.time.bio_dt);
    // Identical seeds must terminate identically (population_stop on a
    // fully drained patch is a legitimate identical end).
    assert(seg_a.termination_cause() == seg_b.termination_cause());
  }
  assert(seg_a.fingerprint() == seg_b.fingerprint());
  std::cout << "live_seed_invariance: fingerprint=" << seg_a.fingerprint()
            << '\n';
}

// The audit slot is the sanctioned mechanism-execution hole in the
// capability gate: the same sos_basal_rate that is refused on a fully
// tabulated segment is admitted when a live patch hosts it.
void test_audit_slot_admits_mechanism_config() {
  SimulationConfig live_cfg = live_cfg_for_test();
  live_cfg.fixes.bacteriocin.sos_basal_rate = 0.02;
  {
    MucusSegment seg;
    seg.init(live_cfg);  // audit.patch_index = 0: admitted
  }
  SimulationConfig tab_cfg = live_cfg;
  tab_cfg.layer2.audit.patch_index = -1;
  bool refused = false;
  try {
    MucusSegment seg;
    seg.init(tab_cfg);
  } catch (const ConfigError& e) {
    refused = std::string(e.what()).find("sos_basal_rate") !=
              std::string::npos;
  }
  assert(refused);
  std::cout << "audit_capability: mechanism config admitted only "
               "with live slot\n";
}

}  // namespace

int main() {
  test_ledger_closure();
  test_single_debit_attribution();
  test_seed_invariance_under_permutation();
  test_checkpoint_roundtrip();
  test_observables_bounds();
  test_layer2_config_parsing();
  test_bloom_halt_is_patch_local();
  test_founder_genotype_propagation();
  test_founder_genotype_helper_contract();
  test_genotype_survives_checkpoint();
  test_capability_limit_reports_dead_params();
  test_live_patch_death_accounting();
  test_live_patch_pool_influx();
  test_live_patch_pool_handoff();
  test_live_patch_seed_invariance();
  test_audit_slot_admits_mechanism_config();
  std::cout << "test_layer2_segment: all passed\n";
  return 0;
}
