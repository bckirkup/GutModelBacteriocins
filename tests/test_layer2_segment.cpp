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
#include <cstdio>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
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

  // Nested and flat JSON parsing for audit.patch_indices + validation rejections
  {
    const std::string tmp = secure_temp_file("gutibm_l2_cfg_");
    {
      std::ofstream out(tmp);
      out << R"({"layer2": {"enabled": true, "n_patches": 8, "audit": {"patch_indices": [0, 3, 7], "halt_on_exceed": false, "tolerance": 0.4}}})";
    }
    SimulationConfig parsed = InputParser::parse(tmp);
    assert((parsed.layer2.audit.patch_indices == std::vector<Int>{0, 3, 7}));
    assert(!parsed.layer2.audit.halt_on_exceed);
    assert(std::abs(parsed.layer2.audit.tolerance - 0.4) < 1e-12);

    // Reject specifying both patch_index and patch_indices in JSON
    {
      std::ofstream out(tmp);
      out << R"({"layer2": {"enabled": true, "n_patches": 8, "audit": {"patch_index": 0, "patch_indices": [1]}}})";
    }
    bool rejected_both = false;
    try {
      (void)InputParser::parse(tmp);
    } catch (const ConfigError&) {
      rejected_both = true;
    }
    assert(rejected_both);

    // Reject duplicate entries in patch_indices in JSON
    {
      std::ofstream out(tmp);
      out << R"({"layer2.enabled": true, "layer2.n_patches": 8, "layer2.audit.patch_indices": [2, 2]})";
    }
    bool rejected_dup = false;
    try {
      (void)InputParser::parse(tmp);
    } catch (const ConfigError&) {
      rejected_dup = true;
    }
    assert(rejected_dup);
    std::remove(tmp.c_str());
  }

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

// Minimal RFC-8259 unique-key checker for JSON documents emitted by prov.json.
bool has_unique_json_keys(const std::string& text) {
  size_t pos = 0;
  const auto skip_ws = [&]() {
    while (pos < text.size() &&
           std::isspace(static_cast<unsigned char>(text[pos]))) {
      ++pos;
    }
  };
  const auto parse_str = [&]() -> std::string {
    assert(pos < text.size() && text[pos] == '"');
    ++pos;
    std::string s;
    while (pos < text.size()) {
      const char ch = text[pos++];
      if (ch == '"') return s;
      if (ch == '\\' && pos < text.size()) {
        s.push_back(text[pos++]);
      } else {
        s.push_back(ch);
      }
    }
    assert(false && "unterminated JSON string");
    return s;
  };

  std::function<bool()> parse_val = [&]() -> bool {
    skip_ws();
    assert(pos < text.size());
    const char c = text[pos];
    if (c == '{') {
      ++pos;
      skip_ws();
      if (pos < text.size() && text[pos] == '}') {
        ++pos;
        return true;
      }
      std::set<std::string> seen_keys;
      while (true) {
        skip_ws();
        const std::string key = parse_str();
        if (!seen_keys.insert(key).second) {
          return false;  // duplicate key in same JSON object
        }
        skip_ws();
        assert(pos < text.size() && text[pos] == ':');
        ++pos;
        if (!parse_val()) return false;
        skip_ws();
        assert(pos < text.size());
        if (text[pos] == '}') {
          ++pos;
          return true;
        }
        assert(text[pos] == ',');
        ++pos;
      }
    }
    if (c == '[') {
      ++pos;
      skip_ws();
      if (pos < text.size() && text[pos] == ']') {
        ++pos;
        return true;
      }
      while (true) {
        if (!parse_val()) return false;
        skip_ws();
        assert(pos < text.size());
        if (text[pos] == ']') {
          ++pos;
          return true;
        }
        assert(text[pos] == ',');
        ++pos;
      }
    }
    if (c == '"') {
      (void)parse_str();
      return true;
    }
    while (pos < text.size() && text[pos] != ',' && text[pos] != '}' &&
           text[pos] != ']' &&
           !std::isspace(static_cast<unsigned char>(text[pos]))) {
      ++pos;
    }
    return true;
  };
  return parse_val();
}

// ── Issue #476: terminal deduplication, horizon cause, multi-slot prov.json ─
void test_layer2_io_dedup_horizon_and_provenance() {
  assert(has_unique_json_keys(R"({"a": 1, "b": [{"a": 2}, {"a": 3}]})"));
  assert(!has_unique_json_keys(R"({"audit_live": {"patch": 0}, "audit_live": {"patch": 1}})"));

  const std::string ts_path = secure_temp_file("gutibm_l2_io_ts_");
  const std::string lin_path = secure_temp_file("gutibm_l2_io_lineage_");
  const std::string prov_path = secure_temp_file("gutibm_l2_io_prov_");

  // Case A: total_steps (4) % summary_interval_steps (2) == 0, 2 live patches,
  // plus BloomHalted patch transition.
  {
    SimulationConfig cfg = base_cfg();
    cfg.time.bio_dt = 60.0;
    cfg.time.total_time = 240.0;  // 4 steps
    auto& l = cfg.layer2;
    l.n_patches = 4;
    l.patch_lateral_m = 60.0e-6;
    l.patch_depth_m = 60.0e-6;
    l.initial_patch_fraction = 1.0;
    l.initial_founder_cells = 6;
    l.initial_pool_cells = 300;
    l.transit_half_life_s = 60.0;
    l.reattach_prob_per_transit = 1.0;
    l.establish_prob_single = 1.0;
    l.bloom_factor = 1.2;
    l.bloom_sustain_s = 0.0;
    l.contraction_rate_per_min = 0.0;
    l.invalid_patch_fraction_stop = 1.0;
    l.summary_interval_steps = 2;
    l.audit.patch_indices = {0, 2};
    l.audit.halt_on_exceed = false;
    l.audit.lineage_file = lin_path;
    l.timeseries_file = ts_path;
    l.provenance_file = prov_path;

    MucusSegment seg;
    seg.init(cfg);
    const int rc = seg.run();
    assert(rc == 0);
    assert(seg.termination_cause() == TerminationCause::HorizonReached);
    assert(seg.termination_detail() == "reached configured horizon");
    assert(!seg.halted_patches().empty());

    // Verify ts.csv: steps {0, 2, 4} each emitted exactly once.
    {
      std::ifstream in(ts_path);
      std::string header;
      assert(std::getline(in, header));
      assert(header.find("deaths_cum,sos_lysis_deaths_cum,toxin_kill_deaths_cum") !=
             std::string::npos);
      std::vector<Int> steps;
      std::string line;
      while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string step_tok;
        assert(std::getline(ss, step_tok, ','));
        steps.push_back(std::stoi(step_tok));
      }
      assert((steps == std::vector<Int>{0, 2, 4}));
    }

    // Verify lineage.csv: unique (step, patch_id, lineage_id) rows, steps {0, 2, 4},
    // and both live patches {0, 2} present.
    {
      std::ifstream in(lin_path);
      std::string header;
      assert(std::getline(in, header));
      std::set<std::tuple<Int, Int, TagID>> seen;
      std::set<Int> steps;
      std::set<Int> patches;
      std::string line;
      while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string step_s;
        std::string time_s;
        std::string patch_s;
        std::string lin_s;
        assert(std::getline(ss, step_s, ','));
        assert(std::getline(ss, time_s, ','));
        assert(std::getline(ss, patch_s, ','));
        assert(std::getline(ss, lin_s, ','));
        const Int st = std::stoi(step_s);
        const Int pid = std::stoi(patch_s);
        const TagID lid = std::stoll(lin_s);
        assert(seen.emplace(st, pid, lid).second);
        steps.insert(st);
        patches.insert(pid);
      }
      assert((steps == std::set<Int>{0, 2, 4}));
      assert((patches == std::set<Int>{0, 2}));
    }

    // Verify prov.json: RFC-8259 unique keys, horizon_reached, audit_live_patches,
    // halted_patches, and split death ledger keys.
    {
      std::ifstream in(prov_path);
      std::ostringstream ss;
      ss << in.rdbuf();
      const std::string prov = ss.str();
      assert(has_unique_json_keys(prov));
      assert(prov.find("\"termination_cause\": \"horizon_reached\"") !=
             std::string::npos);
      assert(prov.find("\"incomplete_unknown\"") == std::string::npos);
      assert(prov.find("\"audit_live\":") == std::string::npos);
      assert(prov.find("\"audit_live_patches\": [{\"patch\": 0, \"live_seed\": ") !=
             std::string::npos);
      assert(prov.find("{\"patch\": 2, \"live_seed\": ") != std::string::npos);
      assert(prov.find("\"halted_patches\": [{\"patch_id\": ") !=
             std::string::npos);
      assert(prov.find("\"status\": \"bloom_halted\"") != std::string::npos);
      assert(prov.find("\"sos_lysis_deaths\": ") != std::string::npos);
      assert(prov.find("\"toxin_kill_deaths\": ") != std::string::npos);
    }
  }

  // Case B: total_steps (3) % summary_interval_steps (2) != 0, 1 live patch
  // (backward-compatible "audit_live" + "audit_live_patches"), plus SpatialHalted.
  {
    SimulationConfig cfg = base_cfg();
    cfg.time.bio_dt = 60.0;
    cfg.time.total_time = 180.0;  // 3 steps
    auto& l = cfg.layer2;
    l.n_patches = 2;
    l.patch_lateral_m = 40.0e-6;
    l.patch_depth_m = 30.0e-6;
    l.cell_radius_m = 10.0e-6;  // small spatial capacity (~6 cells)
    l.spatial_packing_fraction = 0.6;
    l.initial_patch_fraction = 1.0;
    l.initial_founder_cells = 4;
    l.initial_pool_cells = 80;
    l.transit_half_life_s = 60.0;
    l.reattach_prob_per_transit = 1.0;
    l.establish_prob_single = 1.0;
    l.bloom_factor = 100.0;
    l.contraction_rate_per_min = 0.0;
    l.invalid_patch_fraction_stop = 1.0;
    l.summary_interval_steps = 2;
    l.audit.patch_indices = {1};
    l.audit.halt_on_exceed = false;
    l.audit.lineage_file = lin_path;
    l.timeseries_file = ts_path;
    l.provenance_file = prov_path;

    MucusSegment seg;
    seg.init(cfg);
    const int rc = seg.run();
    assert(rc == 0);
    assert(seg.termination_cause() == TerminationCause::HorizonReached);

    std::ifstream in(ts_path);
    std::string header;
    assert(std::getline(in, header));
    std::vector<Int> steps;
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty()) continue;
      std::istringstream ss(line);
      std::string step_tok;
      assert(std::getline(ss, step_tok, ','));
      steps.push_back(std::stoi(step_tok));
    }
    assert((steps == std::vector<Int>{0, 2, 3}));

    std::ifstream pin(prov_path);
    std::ostringstream pss;
    pss << pin.rdbuf();
    const std::string prov = pss.str();
    assert(has_unique_json_keys(prov));
    assert(prov.find("\"audit_live\": {\"patch\": 1,") != std::string::npos);
    assert(prov.find("\"audit_live_patches\": [{\"patch\": 1,") !=
           std::string::npos);
    assert(prov.find("\"status\": \"spatial_halted\"") != std::string::npos);
  }

  std::remove(ts_path.c_str());
  std::remove(lin_path.c_str());
  std::remove(prov_path.c_str());
  std::cout << "layer2_io_dedup_horizon_and_provenance: passed\n";
}

// ── Issue #465: single-slot patch_indices: [k] vs patch_index: k parity ──
void test_single_slot_patch_indices_bitwise_parity() {
  // Config validation rejects simultaneous specification, duplicates, and OOB.
  {
    SimulationConfig bad = live_cfg_for_test();
    bad.layer2.n_patches = 4;
    bad.layer2.audit.patch_index = 0;
    bad.layer2.audit.patch_indices = {1};
    bool rejected = false;
    try {
      MucusSegment seg;
      seg.init(bad);
    } catch (const ConfigError&) {
      rejected = true;
    }
    assert(rejected);
  }
  {
    SimulationConfig bad = live_cfg_for_test();
    bad.layer2.n_patches = 4;
    bad.layer2.audit.patch_index = -1;
    bad.layer2.audit.patch_indices = {2, 2};
    bool rejected = false;
    try {
      MucusSegment seg;
      seg.init(bad);
    } catch (const ConfigError&) {
      rejected = true;
    }
    assert(rejected);
  }
  {
    SimulationConfig bad = live_cfg_for_test();
    bad.layer2.n_patches = 4;
    bad.layer2.audit.patch_index = -1;
    bad.layer2.audit.patch_indices = {0, 4};
    bool rejected = false;
    try {
      MucusSegment seg;
      seg.init(bad);
    } catch (const ConfigError&) {
      rejected = true;
    }
    assert(rejected);
  }

  for (const Int k : {0, 2}) {
    SimulationConfig cfg_scalar = live_cfg_for_test();
    cfg_scalar.layer2.n_patches = 4;
    cfg_scalar.layer2.initial_founder_cells = 12;
    cfg_scalar.layer2.initial_pool_cells = 15;
    cfg_scalar.layer2.contraction_rate_per_min = 20.0;
    cfg_scalar.layer2.reattach_prob_per_transit = 0.8;
    cfg_scalar.layer2.establish_prob_single = 0.8;
    cfg_scalar.fixes.bacteriocin.sos_basal_rate = 0.01;
    cfg_scalar.layer2.audit.patch_index = k;
    cfg_scalar.layer2.audit.patch_indices.clear();

    SimulationConfig cfg_vec = cfg_scalar;
    cfg_vec.layer2.audit.patch_index = -1;
    cfg_vec.layer2.audit.patch_indices = {k};

    MucusSegment seg_scalar;
    seg_scalar.init(cfg_scalar);
    MucusSegment seg_vec;
    seg_vec.init(cfg_vec);
    assert(seg_scalar.fingerprint() == seg_vec.fingerprint());

    for (int step = 0; step < 15; ++step) {
      seg_scalar.step(cfg_scalar.time.bio_dt);
      seg_vec.step(cfg_vec.time.bio_dt);
      assert(seg_scalar.fingerprint() == seg_vec.fingerprint());
      assert(seg_scalar.ledger().births == seg_vec.ledger().births);
      assert(seg_scalar.ledger().reseeds == seg_vec.ledger().reseeds);
      assert(seg_scalar.ledger().biological_deaths ==
             seg_vec.ledger().biological_deaths);
      assert(seg_scalar.ledger().sos_lysis_deaths ==
             seg_vec.ledger().sos_lysis_deaths);
      assert(seg_scalar.ledger().toxin_kill_deaths ==
             seg_vec.ledger().toxin_kill_deaths);
    }

    std::vector<Agent> cells_s;
    std::vector<Agent> cells_v;
    seg_scalar.for_each_colonist([&](const Agent& a) { cells_s.push_back(a); });
    seg_vec.for_each_colonist([&](const Agent& a) { cells_v.push_back(a); });
    assert(cells_s.size() == cells_v.size());
    for (size_t i = 0; i < cells_s.size(); ++i) {
      assert(cells_s[i].identity.tag == cells_v[i].identity.tag);
      assert(cells_s[i].genome.lineage_id == cells_v[i].genome.lineage_id);
      assert(cells_s[i].x == cells_v[i].x);
    }
  }
  std::cout << "single_slot_patch_indices_bitwise_parity: passed\n";
}

// ── Issue #465: live try_deposit onto Active vs Halted live patches ──────
void test_live_try_deposit_active_vs_halted() {
  SimulationConfig cfg = live_cfg_for_test();
  cfg.layer2.n_patches = 1;
  cfg.layer2.patch_lateral_m = 40.0e-6;
  cfg.layer2.patch_depth_m = 30.0e-6;
  cfg.domain.grid_dx = 2.0e-6;
  cfg.layer2.initial_founder_cells = 1;
  cfg.layer2.audit.patch_index = -1;
  cfg.layer2.audit.patch_indices = {0};

  MucusSegment seg;
  seg.init(cfg);
  assert(seg.patches()[0].live);
  assert(seg.patches()[0].status == PatchStatus::Active);

  // Upper-boundary arrival (z >= D = 30 um) lands at outermost grid cell
  // center z = L_z - 0.5 * dz = 29 um.
  const TagID deposit_tag = (static_cast<TagID>(0xFF) << 48) | 1;
  Agent boundary_cell =
      Agent::create_default(deposit_tag, 1, Vec3{7.0e-6, 11.0e-6, 30.0e-6},
                            5.0e-4);
  boundary_cell.genome.lineage_id = deposit_tag;
  assert(seg.try_deposit(std::move(boundary_cell), 101));
  assert(seg.ledger().external_arrivals == 1);
  assert(seg.ledger_closed());

  bool found_boundary = false;
  seg.for_each_colonist([&](const Agent& a) {
    if (a.identity.tag == deposit_tag) {
      found_boundary = true;
      assert(std::abs(a.x[0] - 7.0e-6) < 1e-18);
      assert(std::abs(a.x[1] - 11.0e-6) < 1e-18);
      assert(std::abs(a.x[2] - 29.0e-6) < 1e-18);
    }
  });
  assert(found_boundary);

  // Out-of-bounds lateral coordinate is rejected without incrementing ledger.
  Agent bad_cell =
      Agent::create_default(deposit_tag + 1, 1,
                            Vec3{-1.0e-6, 11.0e-6, 15.0e-6}, 5.0e-4);
  assert(!seg.try_deposit(std::move(bad_cell), 102));
  assert(seg.ledger().external_arrivals == 1);
  assert(seg.ledger_closed());

  // Transition patch 0 to BloomHalted and verify try_deposit rejects it.
  SimulationConfig halt_cfg = cfg;
  halt_cfg.layer2.bloom_factor = 0.4;
  halt_cfg.layer2.bloom_sustain_s = 0.0;
  halt_cfg.layer2.invalid_patch_fraction_stop = 1.0;
  MucusSegment halted_seg;
  halted_seg.init(halt_cfg);
  halted_seg.step(halt_cfg.time.bio_dt);
  assert(halted_seg.patches()[0].status == PatchStatus::BloomHalted);
  Agent halted_attempt =
      Agent::create_default(deposit_tag + 2, 1,
                            Vec3{7.0e-6, 11.0e-6, 30.0e-6}, 5.0e-4);
  assert(!halted_seg.try_deposit(std::move(halted_attempt), 103));
  assert(halted_seg.ledger().external_arrivals == 0);
  assert(halted_seg.ledger_closed());
  std::cout << "live_try_deposit_active_vs_halted: passed\n";
}

// ── Issue #465: multi-slot (2, 3, and 10 all-live) C/S/R reseeding & split deaths ─
void run_multi_slot_csr_case(Int n_patches,
                             const std::vector<Int>& live_indices,
                             int n_steps) {
  SimulationConfig cfg = base_cfg();
  cfg.seed = 20261121;
  cfg.b12_initial_conc = 1.0e-5;  // unshielded BtuB so colicin E kills fire
  for (auto& chem : cfg.chemicals) {
    if (chem.name == "b12") {
      chem.initial_conc = 1.0e-5;
      chem.boundary_conc = 1.0e-5;
    }
  }
  cfg.fixes.bacteriocin.sos_basal_rate = 0.005;
  cfg.fixes.bacteriocin.burst_release_tau = 60.0;
  cfg.fixes.receptor.kill_rate_colicin = 0.1;

  cfg.initial_strains.clear();
  // Strain 1: Producer (C) carrying ColE2
  SimulationConfig::InitialStrain c_strain;
  c_strain.type = 1;
  c_strain.count = 35;
  c_strain.mu_max = 8.0e-4;
  c_strain.plasmids = {"ColE2"};
  // Strain 2: Sensitive (S)
  SimulationConfig::InitialStrain s_strain;
  s_strain.type = 2;
  s_strain.count = 35;
  s_strain.mu_max = 8.0e-4;
  // Strain 3: Resistant (R, BtuB = 0.0)
  SimulationConfig::InitialStrain r_strain;
  r_strain.type = 3;
  r_strain.count = 30;
  r_strain.mu_max = 7.5e-4;
  test::declare_resistant_founder(r_strain);
  cfg.initial_strains = {c_strain, s_strain, r_strain};

  auto& l = cfg.layer2;
  l.enabled = true;
  l.n_patches = n_patches;
  l.patch_lateral_m = 40.0e-6;
  l.patch_depth_m = 30.0e-6;
  l.initial_patch_fraction = 1.0;
  l.initial_founder_cells = 20;
  l.initial_pool_cells = 0;
  l.contraction_rate_per_min = 30.0;
  l.agent_loss_fraction = 0.25;
  l.crypt.disruption_prob = 0.6;
  l.exposed_proximal.disruption_prob = 0.6;
  l.exposed_distal.disruption_prob = 0.6;
  l.transit_half_life_s = 120.0;
  l.reattach_prob_per_transit = 1.0;
  l.establish_prob_single = 1.0;
  l.establish_ratio = 10.0;
  l.bloom_factor = 100.0;
  l.invalid_patch_fraction_stop = 1.0;
  l.audit.patch_index = -1;
  l.audit.patch_indices = live_indices;
  l.audit.halt_on_exceed = false;

  MucusSegment seg;
  seg.init(cfg);

  for (const Int idx : live_indices) {
    assert(seg.patches()[static_cast<size_t>(idx)].live);
  }

  Int max_cross_patch_colonists = 0;
  for (int step = 0; step < n_steps; ++step) {
    seg.step(cfg.time.bio_dt);
    assert(seg.termination_cause() == TerminationCause::IncompleteUnknown);
    assert(seg.ledger_closed());

    const SegmentLedger& L = seg.ledger();
    assert(L.biological_deaths == L.sos_lysis_deaths + L.toxin_kill_deaths);

    Int alive_patch = 0;
    seg.for_each_colonist([&](const Agent&) { ++alive_patch; });
    Int alive_pool = 0;
    for (const auto& pkt : seg.pool()) {
      alive_pool += static_cast<Int>(pkt.cells.size());
    }
    // Exact closed cell conservation identity across patches + pool:
    assert(L.initial_cells + L.births + L.external_arrivals ==
           alive_patch + alive_pool + L.distal_losses + L.biological_deaths);
    // Pool balance identity: departures from patches feed pool, reseeds and
    // distal exits drain pool:
    assert(L.washout_departures + L.contraction_departures +
               L.edge_departures ==
           alive_pool + L.reseeds + L.distal_losses);

    std::vector<Int> counts_per_patch;
    counts_per_patch.reserve(seg.patches().size());
    for (const auto& p : seg.patches()) {
      counts_per_patch.push_back(seg.patch_count(p));
    }
    Int cur_patch = 0;
    Int idx_in_patch = 0;
    Int step_cross_patch = 0;
    seg.for_each_colonist([&](const Agent& a) {
      while (cur_patch < static_cast<Int>(counts_per_patch.size()) &&
             idx_in_patch >= counts_per_patch[static_cast<size_t>(cur_patch)]) {
        ++cur_patch;
        idx_in_patch = 0;
      }
      assert(cur_patch < n_patches);
      const Int origin_patch =
          static_cast<Int>((a.genome.lineage_id >> 32) - 1);
      assert(origin_patch >= 0 && origin_patch < n_patches);
      if (origin_patch != cur_patch) {
        ++step_cross_patch;
      }
      ++idx_in_patch;
    });
    max_cross_patch_colonists =
        std::max(max_cross_patch_colonists, step_cross_patch);
  }

  const SegmentLedger& L = seg.ledger();
  assert(L.reseeds > 0);
  assert(L.sos_lysis_deaths > 0);
  assert(L.toxin_kill_deaths > 0);
  assert(L.biological_deaths == L.sos_lysis_deaths + L.toxin_kill_deaths);
  assert(max_cross_patch_colonists > 0);
  std::cout << "multi_slot_csr(n_patches=" << n_patches
            << ", live=" << live_indices.size()
            << "): reseeds=" << L.reseeds
            << " max_cross_patch=" << max_cross_patch_colonists
            << " sos_lysis=" << L.sos_lysis_deaths
            << " toxin_kills=" << L.toxin_kill_deaths
            << " total_deaths=" << L.biological_deaths << '\n';
}

void test_multi_slot_csr_interpatch_reseed_and_split_deaths() {
  run_multi_slot_csr_case(2, {0, 1}, 12);
  run_multi_slot_csr_case(3, {0, 1, 2}, 12);
  run_multi_slot_csr_case(10, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, 8);
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
  test_layer2_io_dedup_horizon_and_provenance();
  test_single_slot_patch_indices_bitwise_parity();
  test_live_try_deposit_active_vs_halted();
  test_multi_slot_csr_interpatch_reseed_and_split_deaths();
  std::cout << "test_layer2_segment: all passed\n";
  return 0;
}
