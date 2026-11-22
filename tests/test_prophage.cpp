/* -----------------------------------------------------------------------
   GutIBM – Spec 14 temperate prophage unit & integration tests (GAP-10)
   ----------------------------------------------------------------------- */

#include "agent_transfer.h"
#include "chain.h"
#include "error.h"
#include "fix_prophage.h"
#include "founder_genotype.h"
#include "hdf5_reader.h"
#include "hdf5_writer.h"
#include "input_parser.h"
#include "mechanism_capability.h"
#include "path_utils.h"
#include "segment.h"
#include "sim_fingerprint.h"
#include "simulation.h"
#include "species_names.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

#ifndef GUTIBM_SOURCE_DIR
#define GUTIBM_SOURCE_DIR "."
#endif

using namespace gutibm;

namespace {

std::string source_path(const std::string& rel) {
  return std::string(GUTIBM_SOURCE_DIR) + "/" + rel;
}

SimulationConfig make_closed_prophage_base(uint64_t seed) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {60.0e-6, 60.0e-6, 40.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.time.bio_dt = 60.0;
  cfg.time.total_time = 300.0;
  cfg.time.output_interval = 300.0;
  cfg.seed = seed;
  cfg.hdf5.enabled = false;
  cfg.dysbiosis_threshold = 0.0;
  cfg.advection.radial_turnover = 1.0e12;
  cfg.advection.distal_transit_time = 1.0e12;
  cfg.advection.peristaltic_enabled = false;
  cfg.advection.crypts_enabled = false;
  cfg.fixes.bacteriocin.sos_basal_rate = 0.0;
  cfg.fixes.bacteriocin.sos_lysis_prob = 0.0;
  cfg.fixes.bacteriocin.sos_cross_induction_rate = 0.0;
  cfg.fixes.mutation.bi_duplication_rate = 0.0;
  cfg.fixes.mutation.bi_recombination_rate = 0.0;
  cfg.fixes.mutation.receptor_mutation_rate = 0.0;
  cfg.fixes.mutation.super_killer_rate = 0.0;
  cfg.fixes.mutation.partial_resistance_rate = 0.0;
  cfg.fixes.conjugation.base_transfer_rate = 0.0;
  cfg.initial_strains.clear();
  return cfg;
}

StepEvents cumulative_total(const Simulation& sim) {
  StepEvents total = sim.cumulative_events();
  total.add(sim.step_events());
  return total;
}

void test_default_off_bitwise_parity() {
  SimulationConfig with_prophage = make_closed_prophage_base(7101);
  with_prophage.prophage.enabled = false;
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 40;
  strain.mu_max = 5.0e-4;
  strain.plasmids = {"ColE1"};
  with_prophage.initial_strains.push_back(strain);
  with_prophage.enabled_fixes = {
      "metabolism", "bacteriocin", "receptor", "prophage", "mechanics"};

  SimulationConfig without_prophage = with_prophage;
  without_prophage.enabled_fixes = {
      "metabolism", "bacteriocin", "receptor", "mechanics"};

  Simulation sim_a;
  sim_a.init(with_prophage);
  assert(sim_a.chemical_field().find(species::VIRION) < 0);

  Simulation sim_b;
  sim_b.init(without_prophage);

  for (int s = 0; s < 5; ++s) {
    sim_a.step(with_prophage.time.bio_dt);
    sim_b.step(without_prophage.time.bio_dt);
  }

  const uint64_t fp_a = test_util::simulation_fingerprint(sim_a);
  const uint64_t fp_b = test_util::simulation_fingerprint(sim_b);
  assert(fp_a == fp_b);

  const StepEvents ev = cumulative_total(sim_a);
  assert(ev.prophage_inductions == 0);
  assert(ev.prophage_lytic_kills == 0);
  assert(ev.prophage_lysogen_conversions == 0);
  assert(std::abs(ev.virions_released) < 1e-15);
  std::cout << "  test_default_off_bitwise_parity: PASSED\n";
}

void test_cryptic_vs_active_lysogen_and_superinfection_immunity() {
  // Part A: Under high DNA damage, cryptic carriers never induce/lyse while
  // active lysogens experience strong induction and virion release.
  {
    SimulationConfig cfg = make_closed_prophage_base(7201);
    cfg.prophage.enabled = true;
    cfg.prophage.k_spont = 0.0;
    cfg.prophage.k_max_damage = 5.0e-2;
    cfg.prophage.K_damage = 1.0e-6;
    cfg.prophage.hill_n = 2.0;
    cfg.prophage.damage_stimulus = 5.0e-5;  // 50 * K_damage
    cfg.prophage.virion_burst_size = 20.0;
    cfg.prophage.k_adsorb = 0.0;  // isolate induction from secondary infection
    cfg.enabled_fixes = {"prophage"};

    SimulationConfig::InitialStrain cryptic;
    cryptic.type = 1;
    cryptic.count = 80;
    cryptic.mu_max = 0.0;
    cryptic.lysogen_prevalence = 1.0;
    cryptic.active_inducible_fraction = 0.0;

    SimulationConfig::InitialStrain active;
    active.type = 2;
    active.count = 80;
    active.mu_max = 0.0;
    active.lysogen_prevalence = 1.0;
    active.active_inducible_fraction = 1.0;

    cfg.initial_strains = {cryptic, active};

    Simulation sim;
    sim.init(cfg);

    Int init_cryptic = 0;
    Int init_active = 0;
    for (const Agent& a : sim.agents()) {
      if (a.identity.type == 1) {
        assert(a.flags.has_prophage_immunity);
        assert(!a.flags.has_inducible_prophage);
        ++init_cryptic;
      } else if (a.identity.type == 2) {
        assert(a.flags.has_prophage_immunity);
        assert(a.flags.has_inducible_prophage);
        ++init_active;
      }
    }
    assert(init_cryptic == 80);
    assert(init_active == 80);

    for (int s = 0; s < 3; ++s) {
      sim.step(cfg.time.bio_dt);
    }

    Int live_cryptic = 0;
    Int live_active = 0;
    for (const Agent& a : sim.agents()) {
      if (a.state == PhenoState::DEAD) continue;
      if (a.identity.type == 1) ++live_cryptic;
      if (a.identity.type == 2) ++live_active;
    }
    const StepEvents ev = cumulative_total(sim);
    assert(live_cryptic == 80);  // Cryptic carriers NEVER lyse under high D
    assert(live_active < 20);    // Active lysogens strongly induced
    assert(ev.prophage_inductions == 80 - live_active);
    assert(ev.mortality_lysis == ev.prophage_inductions);
    assert(ev.prophage_lytic_kills == 0);
  }

  // Part B: Superinfection immunity — seed dense free virions with p_lysogeny=0.
  // Both cryptic and active lysogens adsorb/neutralize virions with ZERO lytic
  // kills and ZERO state change, while susceptible non-lysogens suffer kills.
  {
    SimulationConfig cfg = make_closed_prophage_base(7202);
    cfg.prophage.enabled = true;
    cfg.prophage.k_spont = 0.0;
    cfg.prophage.k_max_damage = 0.0;
    cfg.prophage.damage_stimulus = 0.0;
    cfg.prophage.k_virion_decay = 0.0;
    cfg.prophage.virion_washout_scale = 0.0;
    cfg.prophage.k_adsorb = 1.0e-14;
    cfg.prophage.p_lysogeny = 0.0;
    cfg.prophage.virion_burst_size = 10.0;
    cfg.enabled_fixes = {"prophage"};

    SimulationConfig::InitialStrain cryptic;
    cryptic.type = 1;
    cryptic.count = 60;
    cryptic.mu_max = 0.0;
    cryptic.lysogen_prevalence = 1.0;
    cryptic.active_inducible_fraction = 0.0;

    SimulationConfig::InitialStrain susceptible;
    susceptible.type = 2;
    susceptible.count = 60;
    susceptible.mu_max = 0.0;
    susceptible.lysogen_prevalence = 0.0;
    susceptible.active_inducible_fraction = 0.0;

    cfg.initial_strains = {cryptic, susceptible};

    Simulation sim;
    sim.init(cfg);

    const Int i_virion = sim.chemical_field().find(species::VIRION);
    assert(i_virion >= 0);
    const Real cell_vol = sim.domain().cell_volume();
    // Seed 50 virions per voxel everywhere
    const Real seed_conc = 50.0 / cell_vol;
    for (Int c = 0; c < sim.chemical_field().ncells(); ++c) {
      sim.chemical_field().conc(i_virion, c) = seed_conc;
    }

    for (int s = 0; s < 3; ++s) {
      sim.step(cfg.time.bio_dt);
    }

    Int live_cryptic = 0;
    Int live_susceptible = 0;
    for (const Agent& a : sim.agents()) {
      if (a.state == PhenoState::DEAD) continue;
      if (a.identity.type == 1) {
        assert(a.flags.has_prophage_immunity);
        assert(!a.flags.has_inducible_prophage);
        ++live_cryptic;
      } else if (a.identity.type == 2) {
        ++live_susceptible;
      }
    }
    const StepEvents ev = cumulative_total(sim);
    assert(live_cryptic == 60);       // 100% superinfection immunity
    assert(live_susceptible < 20);    // Susceptibles lytically killed
    assert(ev.prophage_lytic_kills == 60 - live_susceptible);
    assert(ev.virions_adsorbed > 100.0);
  }
  std::cout << "  test_cryptic_vs_active_lysogen_and_superinfection_immunity: PASSED\n";
}

void test_hill_response_monotonicity() {
  ProphageConfig pcfg;
  pcfg.enabled = true;
  pcfg.k_spont = 1.0e-6;
  pcfg.k_max_damage = 1.0e-2;
  pcfg.K_damage = 2.0e-6;
  pcfg.hill_n = 2.0;

  const std::vector<Real> doses = {
      0.0, 0.5 * pcfg.K_damage, pcfg.K_damage, 5.0 * pcfg.K_damage};
  const std::vector<Real> expected_hill = {
      0.0, 0.2, 0.5, 25.0 / 26.0};

  std::vector<Real> analytic_hazards;
  for (size_t i = 0; i < doses.size(); ++i) {
    const Real h = prophage_hill_induction_hazard(doses[i], pcfg);
    const Real expected = pcfg.k_spont + pcfg.k_max_damage * expected_hill[i];
    assert(std::abs(h - expected) < 1e-12 * std::max(1.0, expected));
    analytic_hazards.push_back(h);
  }
  for (size_t i = 1; i < analytic_hazards.size(); ++i) {
    assert(analytic_hazards[i] > analytic_hazards[i - 1]);
  }

  // Graded live simulation sweep over the same four damage levels
  std::vector<Int> inductions;
  for (const Real dose : doses) {
    SimulationConfig cfg = make_closed_prophage_base(7301);
    cfg.prophage = pcfg;
    cfg.prophage.damage_stimulus = dose;
    cfg.prophage.k_adsorb = 0.0;
    cfg.enabled_fixes = {"prophage"};

    SimulationConfig::InitialStrain lysogens;
    lysogens.type = 1;
    lysogens.count = 300;
    lysogens.mu_max = 0.0;
    lysogens.lysogen_prevalence = 1.0;
    lysogens.active_inducible_fraction = 1.0;
    cfg.initial_strains = {lysogens};

    Simulation sim;
    sim.init(cfg);
    sim.step(cfg.time.bio_dt);
    inductions.push_back(cumulative_total(sim).prophage_inductions);
  }

  assert(std::is_sorted(inductions.begin(), inductions.end()));
  assert(inductions.front() < 5);
  assert(inductions.back() > 100);
  assert(inductions[2] > inductions[1]);
  assert(inductions[3] > inductions[2]);
  std::cout << "  test_hill_response_monotonicity: PASSED (inductions="
            << inductions[0] << "," << inductions[1] << ","
            << inductions[2] << "," << inductions[3] << ")\n";
}

void test_cell_and_virion_conservation_ledger() {
  SimulationConfig cfg = make_closed_prophage_base(7401);
  cfg.carbon_boundary_conc = 5.0e-2;
  for (auto& spec : cfg.chemicals) {
    if (spec.name == species::CARBON) {
      spec.initial_conc = 5.0e-2;
      spec.boundary_conc = 5.0e-2;
    }
  }
  cfg.prophage.enabled = true;
  cfg.prophage.k_spont = 2.0e-4;
  cfg.prophage.k_max_damage = 0.0;
  cfg.prophage.virion_burst_size = 15.0;
  cfg.prophage.D_virion = 2.0e-11;
  cfg.prophage.virion_retardation = 1.0;
  cfg.prophage.k_virion_decay = 5.0e-4;
  cfg.prophage.virion_washout_scale = 0.0;
  cfg.prophage.k_adsorb = 2.0e-14;
  cfg.prophage.p_lysogeny = 0.5;
  cfg.enabled_fixes = {"metabolism", "prophage", "mechanics"};

  SimulationConfig::InitialStrain active;
  active.type = 1;
  active.count = 50;
  active.mu_max = 6.0e-4;
  active.lysogen_prevalence = 1.0;
  active.active_inducible_fraction = 1.0;

  SimulationConfig::InitialStrain cryptic;
  cryptic.type = 2;
  cryptic.count = 30;
  cryptic.mu_max = 6.0e-4;
  cryptic.lysogen_prevalence = 1.0;
  cryptic.active_inducible_fraction = 0.0;

  SimulationConfig::InitialStrain susceptible;
  susceptible.type = 3;
  susceptible.count = 80;
  susceptible.mu_max = 6.0e-4;
  susceptible.lysogen_prevalence = 0.0;
  susceptible.active_inducible_fraction = 0.0;

  cfg.initial_strains = {active, cryptic, susceptible};

  Simulation sim;
  sim.init(cfg);

  const Int initial_live = sim.agents().size();
  const Real initial_virions = FixProphage::total_virions_in_field(sim);
  assert(std::abs(initial_virions) < 1e-15);

  for (int s = 0; s < 35; ++s) {
    sim.step(cfg.time.bio_dt);
  }

  const StepEvents ev = cumulative_total(sim);
  const Int final_live = sim.agents().size();
  const Real final_virions = FixProphage::total_virions_in_field(sim);

  // 1. Exact cell count ledger closure (closed patch, no washout/colicin/CDI)
  assert(ev.outflow_washout == 0);
  assert(ev.outflow_boundary == 0);
  assert(ev.mortality_colicin == 0);
  assert(ev.mortality_cdi == 0);
  const Int expected_live =
      initial_live + ev.divisions - ev.mortality_lysis - ev.prophage_lytic_kills;
  assert(final_live == expected_live);

  // 2. Exact virion count conservation ledger closure
  assert(ev.prophage_inductions > 0);
  assert(ev.virions_released > 0.0);
  assert(ev.virions_adsorbed > 0.0);
  assert(ev.virions_decayed > 0.0);
  const Real lhs = final_virions + ev.virions_adsorbed + ev.virions_decayed;
  const Real rhs = initial_virions + ev.virions_released;
  const Real rel_err = std::abs(lhs - rhs) / std::max(1.0, rhs);
  assert(rel_err < 1.0e-9);

  // Zero reaction clipping on virion species
  const Int i_virion = sim.chemical_field().find(species::VIRION);
  assert(i_virion >= 0);
  const Real clip_total =
      sim.chemical_field().flux_accounting().reaction_clip_cumulative[i_virion]
      + sim.chemical_field().flux_accounting().reaction_clip_interval[i_virion];
  assert(clip_total < 1.0e-12);

  // 3. Verify cryptic division inheritance (type 2 daughters stay cryptic)
  Int cryptic_gen_gt0 = 0;
  for (const Agent& a : sim.agents()) {
    if (a.identity.type == 2) {
      assert(a.flags.has_prophage_immunity);
      assert(!a.flags.has_inducible_prophage);
      if (a.genome.generation > 0) ++cryptic_gen_gt0;
    }
  }
  assert(cryptic_gen_gt0 > 0);
  std::cout << "  test_cell_and_virion_conservation_ledger: PASSED (rel_err="
            << rel_err << ", conversions=" << ev.prophage_lysogen_conversions
            << ", lytic_kills=" << ev.prophage_lytic_kills << ")\n";
}

void test_latent_lysis_delay() {
  SimulationConfig cfg = make_closed_prophage_base(7501);
  cfg.prophage.enabled = true;
  cfg.prophage.k_spont = 0.0;
  cfg.prophage.k_max_damage = 1.0;  // guaranteed induction in step 1
  cfg.prophage.K_damage = 1.0e-6;
  cfg.prophage.damage_stimulus = 1.0e-3;
  cfg.prophage.lysis_delay_s = 120.0;  // 2 * bio_dt
  cfg.prophage.virion_burst_size = 12.0;
  cfg.prophage.k_adsorb = 0.0;
  cfg.prophage.k_virion_decay = 0.0;
  cfg.prophage.virion_washout_scale = 0.0;
  cfg.enabled_fixes = {"prophage"};

  SimulationConfig::InitialStrain active;
  active.type = 1;
  active.count = 20;
  active.mu_max = 0.0;
  active.lysogen_prevalence = 1.0;
  active.active_inducible_fraction = 1.0;
  cfg.initial_strains = {active};

  Simulation sim;
  sim.init(cfg);

  // Step 1: all 20 induce, none lyse yet (latent timer = 120 s)
  sim.step(60.0);
  StepEvents ev1 = cumulative_total(sim);
  assert(ev1.prophage_inductions == 20);
  assert(ev1.mortality_lysis == 0);
  assert(std::abs(ev1.virions_released) < 1e-15);
  assert(sim.agents().size() == 20);
  for (const Agent& a : sim.agents()) {
    assert(a.flags.prophage_induced);
  }

  // Step 2: timer decrements 120 -> 60 s (still latent)
  sim.step(60.0);
  StepEvents ev2 = cumulative_total(sim);
  assert(ev2.mortality_lysis == 0);
  assert(sim.agents().size() == 20);

  // Step 3: timer decrements 60 -> 0 s -> all 20 lyse and release bursts
  sim.step(60.0);
  StepEvents ev3 = cumulative_total(sim);
  assert(ev3.prophage_inductions == 20);
  assert(ev3.mortality_lysis == 20);
  assert(std::abs(ev3.virions_released - 240.0) < 1e-9);
  assert(std::abs(FixProphage::total_virions_in_field(sim) - 240.0) < 1e-9);
  assert(sim.agents().size() == 0);
  std::cout << "  test_latent_lysis_delay: PASSED\n";
}

void test_dual_regime_presets_and_realized_daily_rate() {
  const SimulationConfig human = InputParser::parse(
      source_path("examples/prophage_presets/human_commensal_community.json"));
  assert(human.prophage.enabled);
  assert(std::abs(human.prophage.lysogen_prevalence_init - 0.8) < 1e-12);
  assert(std::abs(human.prophage.active_inducible_fraction - 0.375) < 1e-12);
  assert(std::abs(human.prophage.k_spont - 1.2e-7) < 1e-15);
  const Real expected_human_daily =
      human.prophage.lysogen_prevalence_init
      * human.prophage.active_inducible_fraction
      * human.prophage.k_spont * 86400.0;
  assert(expected_human_daily >= 0.001 && expected_human_daily <= 0.01);

  SimulationConfig murine = InputParser::parse(
      source_path("examples/prophage_presets/murine_monoxenic_sos.json"));
  assert(murine.prophage.enabled);
  assert(std::abs(murine.prophage.lysogen_prevalence_init - 1.0) < 1e-12);
  assert(std::abs(murine.prophage.active_inducible_fraction - 1.0) < 1e-12);
  assert(murine.prophage.k_spont >= 1.70e-6
         && murine.prophage.k_spont <= 4.44e-6);

  // Run murine preset briefly to verify realized daily rate emission
  murine.domain.hi = {60.0e-6, 60.0e-6, 40.0e-6};
  murine.domain.grid_dx = 5.0e-6;
  murine.hdf5.enabled = false;
  murine.dysbiosis_threshold = 0.0;
  murine.enabled_fixes = {"prophage"};
  murine.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 400;
  s.mu_max = 0.0;
  murine.initial_strains.push_back(s);

  Simulation sim;
  sim.init(murine);
  for (int step = 0; step < 30; ++step) {
    sim.step(60.0);
  }
  const StepEvents ev = cumulative_total(sim);
  assert(ev.prophage_bacterium_seconds > 0.0);
  const Real recomputed =
      (static_cast<Real>(ev.prophage_inductions)
       / ev.prophage_bacterium_seconds) * 86400.0;
  assert(std::isfinite(ev.prophage_induction_rate_per_bacterium_day));
  assert(std::abs(ev.prophage_induction_rate_per_bacterium_day - recomputed)
         < 1e-12);
  std::cout << "  test_dual_regime_presets_and_realized_daily_rate: PASSED"
            << " (human_expected=" << expected_human_daily
            << "/bact/d, murine_realized="
            << ev.prophage_induction_rate_per_bacterium_day << "/bact/d)\n";
}

struct HenrotRunOutcome {
  Int inductions = 0;
  Int conversions = 0;
  Int lytic_kills = 0;
  Real virions_released = 0.0;
};

HenrotRunOutcome run_henrot_arm(bool mucosal_regime, bool with_colE2,
                                uint64_t seed) {
  SimulationConfig cfg = make_closed_prophage_base(seed);
  cfg.domain.hi = {50.0e-6, 50.0e-6, 30.0e-6};
  cfg.domain.grid_dx = 5.0e-6;
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.time.bio_dt = 60.0;
  cfg.qssa.toxin_cutoff = 100.0e-6;
  cfg.qssa.use_fmm = false;
  // Unshielded BtuB per #464 so ColE2 exposure is active
  cfg.fixes.receptor.btuB_competitor_conc = 0.0;
  cfg.fixes.receptor.kill_rate_colicin = 1.0e-5;  // sub-lethal ColE2 SOS probe
  cfg.fixes.bacteriocin.sos_basal_rate = 5.0e-3;
  cfg.fixes.bacteriocin.sos_lysis_prob = 0.0;
  cfg.fixes.bacteriocin.sos_cross_induction_rate = 0.0;
  cfg.enabled_fixes = {"bacteriocin", "receptor", "prophage"};

  cfg.prophage.enabled = true;
  cfg.prophage.k_spont = 1.0e-6;
  cfg.prophage.k_max_damage = 3.0e-2;
  cfg.prophage.K_damage = 1.0e-8;
  cfg.prophage.hill_n = 2.0;
  cfg.prophage.damage_stimulus = 0.0;
  cfg.prophage.nuclease_damage_weight = 1.0;
  cfg.prophage.virion_burst_size = 40.0;
  cfg.prophage.k_adsorb = 4.0e-14;
  cfg.prophage.p_lysogeny = 0.5;

  if (mucosal_regime) {
    // Mucosal advection/retardation + decay + high cryptic superinfection blockade
    cfg.advection.mucus_thickness = 30.0e-6;
    cfg.advection.radial_turnover = 600.0;
    cfg.advection.distal_transit_time = 7200.0;
    cfg.prophage.D_virion = 4.0e-12;
    cfg.prophage.virion_retardation = 10.0;
    cfg.prophage.k_virion_decay = 2.0e-2;
    cfg.prophage.virion_washout_scale = 1.0;
  } else {
    // Well-mixed in-vitro positive control: high D_virion, zero mucus washout,
    // low decay, no cryptic blockade
    cfg.prophage.D_virion = 5.0e-10;
    cfg.prophage.virion_retardation = 1.0;
    cfg.prophage.k_virion_decay = 1.0e-5;
    cfg.prophage.virion_washout_scale = 0.0;
  }

  SimulationConfig::InitialStrain producer;
  producer.type = 1;
  producer.count = 40;
  producer.mu_max = 0.0;
  producer.plasmids = with_colE2 ? std::vector<std::string>{"ColE2"}
                                 : std::vector<std::string>{};
  producer.lysogen_prevalence = 0.0;

  SimulationConfig::InitialStrain donor_lysogen;
  donor_lysogen.type = 2;
  donor_lysogen.count = 40;
  donor_lysogen.mu_max = 0.0;
  donor_lysogen.lysogen_prevalence = 1.0;
  donor_lysogen.active_inducible_fraction = 1.0;

  SimulationConfig::InitialStrain recipient;
  recipient.type = 3;
  recipient.count = 160;
  recipient.mu_max = 0.0;
  // In vivo: high background cryptic immunity (90% cryptic virion sink);
  // In vitro: 0% cryptic blockade (100% susceptible recipients)
  recipient.lysogen_prevalence = mucosal_regime ? 0.90 : 0.0;
  recipient.active_inducible_fraction = 0.0;

  cfg.initial_strains = {producer, donor_lysogen, recipient};

  Simulation sim;
  sim.init(cfg);
  for (int step = 0; step < 18; ++step) {
    sim.step(cfg.time.bio_dt);
  }
  const StepEvents ev = cumulative_total(sim);
  return {ev.prophage_inductions, ev.prophage_lysogen_conversions,
          ev.prophage_lytic_kills, ev.virions_released};
}

void test_henrot_in_vitro_vs_in_vivo_contrast() {
  const HenrotRunOutcome invitro_colE2 = run_henrot_arm(false, true, 7601);
  const HenrotRunOutcome invitro_null = run_henrot_arm(false, false, 7601);
  const HenrotRunOutcome invivo_colE2 = run_henrot_arm(true, true, 7601);

  // Positive in-vitro control: ColE2 DNase triggers Hill prophage induction
  // and yields measurable lysogen conversion in susceptible recipients.
  assert(invitro_colE2.inductions > invitro_null.inductions + 10);
  assert(invitro_colE2.conversions >= 10);
  assert(invitro_colE2.conversions > invitro_null.conversions + 8);

  // Mucosal in-vivo null/suppressed transfer regime (Henrot et al. 2026):
  // retardation + washout/decay + cryptic superinfection immunity sharply
  // suppresses lysogen conversion relative to the well-mixed in-vitro control.
  assert(invivo_colE2.conversions * 4 < invitro_colE2.conversions);

  std::cout << "  test_henrot_in_vitro_vs_in_vivo_contrast: PASSED"
            << " (in_vitro_ColE2 conv=" << invitro_colE2.conversions
            << ", in_vitro_null conv=" << invitro_null.conversions
            << ", in_vivo_ColE2 conv=" << invivo_colE2.conversions << ")\n";
}

void test_mpi_transfer_and_hdf5_checkpoint_roundtrip() {
  // 1. MPI agent_transfer pack/unpack round-trip
  Agent a = Agent::create_default(91, 2, {12e-6, 18e-6, 8e-6}, 5e-4);
  a.flags.has_prophage_immunity = true;
  a.flags.has_inducible_prophage = false;  // cryptic carrier
  a.flags.prophage_induced = false;
  a.timers.prophage_lysis_timer = -1.0;

  Agent b = Agent::create_default(92, 3, {14e-6, 18e-6, 8e-6}, 5e-4);
  b.flags.has_prophage_immunity = true;
  b.flags.has_inducible_prophage = true;
  b.flags.prophage_induced = true;
  b.timers.prophage_lysis_timer = 180.0;

  std::vector<char> buf;
  agent_transfer_serialize({a, b}, buf);
  const std::vector<Agent> unpacked = agent_transfer_deserialize(buf);
  assert(unpacked.size() == 2);
  assert(unpacked[0].flags.has_prophage_immunity);
  assert(!unpacked[0].flags.has_inducible_prophage);
  assert(!unpacked[0].flags.prophage_induced);
  assert(unpacked[1].flags.has_prophage_immunity);
  assert(unpacked[1].flags.has_inducible_prophage);
  assert(unpacked[1].flags.prophage_induced);
  assert(std::abs(unpacked[1].timers.prophage_lysis_timer - 180.0) < 1e-12);

#ifdef GUTIBM_HDF5
  // 2. HDF5 checkpoint round-trip + backward-compatible legacy default test
  const std::string h5_path = secure_temp_file("gutibm_prophage_ckpt_");
  SimulationConfig cfg = make_closed_prophage_base(7701);
  cfg.prophage.enabled = true;
  cfg.prophage.k_spont = 0.0;
  cfg.prophage.k_max_damage = 1.0;
  cfg.prophage.K_damage = 1.0e-6;
  cfg.prophage.damage_stimulus = 1.0e-3;
  cfg.prophage.lysis_delay_s = 180.0;
  cfg.enabled_fixes = {"prophage"};
  cfg.hdf5.enabled = false;

  SimulationConfig::InitialStrain cryptic;
  cryptic.type = 1;
  cryptic.count = 8;
  cryptic.mu_max = 0.0;
  cryptic.lysogen_prevalence = 1.0;
  cryptic.active_inducible_fraction = 0.0;

  SimulationConfig::InitialStrain active;
  active.type = 2;
  active.count = 8;
  active.mu_max = 0.0;
  active.lysogen_prevalence = 1.0;
  active.active_inducible_fraction = 1.0;

  cfg.initial_strains = {cryptic, active};

  {
    Simulation sim;
    sim.init(cfg);
    sim.step(60.0);
    assert(HDF5Writer::write_closed_restart(
        sim, h5_path, sim.step_count(), sim.time(), cfg.time.bio_dt));
  }

  Simulation restarted;
  restarted.init_from_checkpoint(cfg, h5_path);
  Int restored_cryptic = 0;
  Int restored_active_induced = 0;
  for (const Agent& ag : restarted.agents()) {
    if (ag.identity.type == 1) {
      assert(ag.flags.has_prophage_immunity);
      assert(!ag.flags.has_inducible_prophage);
      assert(!ag.flags.prophage_induced);
      ++restored_cryptic;
    } else if (ag.identity.type == 2) {
      assert(ag.flags.has_prophage_immunity);
      assert(ag.flags.has_inducible_prophage);
      assert(ag.flags.prophage_induced);
      assert(std::abs(ag.timers.prophage_lysis_timer - 180.0) < 1e-12);
      ++restored_active_induced;
    }
  }
  assert(restored_cryptic == 8);
  assert(restored_active_induced == 8);
  assert(restarted.cumulative_events().prophage_inductions == 8);

  // Strip optional prophage datasets to simulate a legacy checkpoint file
  {
    hid_t fid = H5Fopen(h5_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    assert(fid >= 0);
    for (const char* ds : {"agents/step_000001/has_prophage_immunity",
                           "agents/step_000001/has_inducible_prophage",
                           "agents/step_000001/prophage_induced",
                           "agents/step_000001/prophage_lysis_timer"}) {
      if (H5Lexists(fid, ds, H5P_DEFAULT) > 0) {
        assert(H5Ldelete(fid, ds, H5P_DEFAULT) >= 0);
      }
    }
    H5Fclose(fid);
  }
  Simulation legacy_restarted;
  legacy_restarted.init_from_checkpoint(cfg, h5_path);
  for (const Agent& ag : legacy_restarted.agents()) {
    assert(!ag.flags.has_prophage_immunity);
    assert(!ag.flags.has_inducible_prophage);
    assert(!ag.flags.prophage_induced);
    assert(ag.timers.prophage_lysis_timer < 0.0);
  }
  std::remove(h5_path.c_str());
#endif
  std::cout << "  test_mpi_transfer_and_hdf5_checkpoint_roundtrip: PASSED\n";
}

void test_mechanism_capability_refuses_prophage_on_tabulated_layers() {
  SimulationConfig def = InputParser::default_config();
  assert(unexecuted_mechanism_params(def).empty());

  SimulationConfig cfg = def;
  cfg.prophage.enabled = true;
  cfg.prophage.k_spont = 2.0e-6;
  const auto dead = unexecuted_mechanism_params(cfg);
  assert(std::ranges::find(dead, "prophage.enabled") != dead.end());
  assert(std::ranges::find(dead, "prophage.k_spont") != dead.end());

  // Tabulated-only Layer-2 segment must refuse prophage config
  SimulationConfig l2 = def;
  l2.prophage.enabled = true;
  l2.layer2.enabled = true;
  l2.layer2.n_patches = 4;
  l2.time.total_time = 600.0;
  l2.layer2.audit.patch_index = -1;
  l2.layer2.audit.patch_indices.clear();
  bool l2_refused = false;
  try {
    MucusSegment seg;
    seg.init(l2);
  } catch (const ConfigError&) {
    l2_refused = true;
  }
  assert(l2_refused);

  // Tabulated-only Layer-3 chain must also refuse prophage config
  SimulationConfig l3 = def;
  l3.prophage.enabled = true;
  l3.layer3.enabled = true;
  l3.time.total_time = 600.0;
  bool l3_refused = false;
  try {
    ColonicChain chain;
    chain.init(l3);
  } catch (const ConfigError&) {
    l3_refused = true;
  }
  assert(l3_refused);

  std::cout << "  test_mechanism_capability_refuses_prophage_on_tabulated_layers: PASSED\n";
}

void test_cognate_nuclease_immunity_specificity() {
  // ColE2 (DNase targeting BtuB) triggers DNA-damage prophage induction.
  // Carrying ColE1 (pore-former targeting BtuB, is_nuclease=false) must NOT
  // cross-protect against ColE2 DNase induction, whereas carrying ColE2
  // (cognate DNase immunity, is_nuclease=true) suppresses induction.
  auto run_case = [](const std::vector<std::string>& lysogen_plasmids) -> Int {
    SimulationConfig cfg = make_closed_prophage_base(7811);
    cfg.domain.hi = {40.0e-6, 40.0e-6, 30.0e-6};
    cfg.domain.grid_dx = 5.0e-6;
    cfg.domain.hash_cell_size = 10.0e-6;
    cfg.qssa.toxin_cutoff = 100.0e-6;
    cfg.qssa.use_fmm = false;
    cfg.prophage.enabled = true;
    cfg.prophage.k_spont = 0.0;
    cfg.prophage.k_max_damage = 0.05;
    cfg.prophage.K_damage = 1.0e-12;
    cfg.prophage.hill_n = 1.0;
    cfg.prophage.nuclease_damage_weight = 1.0;
    cfg.prophage.ros_damage_weight = 0.0;
    cfg.prophage.damage_stimulus = 0.0;
    cfg.prophage.k_adsorb = 0.0;
    cfg.fixes.bacteriocin.sos_basal_rate = 0.0;
    cfg.fixes.bacteriocin.sos_lysis_prob = 0.0;
    cfg.fixes.bacteriocin.sos_cross_induction_rate = 0.0;
    cfg.fixes.receptor.btuB_competitor_conc = 0.0;
    cfg.fixes.receptor.kill_rate_colicin = 0.0;
    cfg.fixes.receptor.immunity_factor = 0.0;
    cfg.enabled_fixes = {"bacteriocin", "receptor", "prophage"};

    SimulationConfig::InitialStrain donor;
    donor.type = 1;
    donor.count = 4;
    donor.mu_max = 0.0;
    donor.plasmids = {"ColE2"};
    donor.lysogen_prevalence = 0.0;

    SimulationConfig::InitialStrain target_lysogen;
    target_lysogen.type = 2;
    target_lysogen.count = 40;
    target_lysogen.mu_max = 0.0;
    target_lysogen.plasmids = lysogen_plasmids;
    target_lysogen.lysogen_prevalence = 1.0;
    target_lysogen.active_inducible_fraction = 1.0;

    cfg.initial_strains = {donor, target_lysogen};
    InputParser::finalize_config(cfg);

    Simulation sim;
    sim.init(cfg);
    // Prime the ColE2 donor cells to lyse on Step 1 post_step without any
    // plasmid SOS induction on target_lysogen cells.
    for (Agent& a : sim.agents()) {
      if (a.identity.type == 1) {
        a.state = PhenoState::SOS_INDUCED;
        a.timers.sos_timer = 10.0;
      }
    }
    for (int s = 0; s < 6; ++s) {
      sim.step(cfg.time.bio_dt);
    }
    return cumulative_total(sim).prophage_inductions;
  };

  const Int ind_unprotected = run_case({});
  const Int ind_pore_btub = run_case({"ColE1"});
  const Int ind_cognate_dnase = run_case({"ColE2"});

  assert(ind_unprotected >= 25);
  assert(ind_pore_btub == ind_unprotected);
  assert(ind_cognate_dnase == 0);
  std::cout << "  test_cognate_nuclease_immunity_specificity: PASSED"
            << " (none=" << ind_unprotected
            << ", ColE1=" << ind_pore_btub
            << ", ColE2=" << ind_cognate_dnase << ")\n";
}

void test_latent_lysis_arrests_division_and_legacy_transfer_blob() {
  // 1. Induced cells in latent lysis window (prophage_induced=true) must not
  //    divide even when biomass exceeds division threshold.
  {
    SimulationConfig cfg = make_closed_prophage_base(7821);
    cfg.prophage.enabled = true;
    cfg.prophage.k_spont = 0.0;
    cfg.prophage.k_max_damage = 1.0;
    cfg.prophage.K_damage = 1.0e-6;
    cfg.prophage.damage_stimulus = 1.0;
    cfg.prophage.lysis_delay_s = 600.0;
    cfg.fixes.metabolism.iron_uptake_enabled = false;
    cfg.fixes.metabolism.b12_uptake_enabled = false;
    cfg.fixes.metabolism.maintenance_rate = 0.0;
    cfg.enabled_fixes = {"metabolism", "prophage"};
    for (auto& c : cfg.chemicals) {
      if (c.name == species::CARBON) {
        c.initial_conc = 10.0;
        c.boundary_conc = 10.0;
      }
    }

    SimulationConfig::InitialStrain lys;
    lys.type = 1;
    lys.count = 16;
    lys.mu_max = 0.05;
    lys.lysogen_prevalence = 1.0;
    lys.active_inducible_fraction = 1.0;
    cfg.initial_strains = {lys};
    InputParser::finalize_config(cfg);

    Simulation sim;
    sim.init(cfg);
    sim.step(cfg.time.bio_dt);
    const Int divs_step1 = cumulative_total(sim).divisions;
    assert(divs_step1 > 0);
    assert(cumulative_total(sim).prophage_inductions == 16 + divs_step1);
    for (int s = 0; s < 3; ++s) {
      sim.step(cfg.time.bio_dt);
    }
    assert(cumulative_total(sim).divisions == divs_step1);
    assert(cumulative_total(sim).prophage_lytic_kills == 0);
  }

  // 2. Legacy v1 agent_transfer blob (count >= 0 without magic header and
  //    without prophage fields) deserializes cleanly with safe defaults.
  {
    struct LegacyId { int64_t tag; int32_t type; int32_t owner_rank; };
    struct LegacySpatial {
      std::array<double, 3> x;
      std::array<double, 3> v;
      double radius;
      double mass;
      double outer_radius;
    };
    struct LegacyMetab {
      double mu_max;
      double mu_realized;
      double biomass;
      double maintenance;
      double realized_fermentation_fraction;
      double pending_growth_chemistry_biomass;
      double respired_oxygen_rate;
      std::array<double, NUM_RECEPTORS> receptor_expr;
      std::array<double, NUM_RECEPTORS> receptor_expr_base;
      double km_iron;
      double km_b12;
      double km_carbon;
    };
    struct LegacyLife {
      int32_t state;
      double age;
      double sos_timer;
      double death_time;
      int32_t grid_cell;
      int32_t in_crypt;
      int32_t microcin_penalty_applied;
    };
    struct LegacyMot {
      std::array<double, 3> swim_direction;
      double swim_speed;
      double run_timer;
      int32_t is_stopped;
      double stop_timer;
      double prev_carbon;
      double prev_oxygen;
      double prev_ai2;
    };
    struct LegacyGen {
      int64_t lineage_id;
      int64_t parent_id;
      uint32_t generation;
      uint32_t mutations;
      int32_t has_conjugative_plasmid;
      double plasmid_cost_amelioration;
      uint16_t cdi_type;
      uint16_t cdi_immunity;
      std::array<double, NUM_RECEPTORS> genome_receptor_expr;
      std::array<double, NUM_RECEPTORS> toxin_affinity;
      std::array<double, NUM_RECEPTORS> ligand_affinity;
      int32_t num_bi_loci;
    };
    struct LegacyBlobRec {
      LegacyId id;
      LegacySpatial sp;
      LegacyMetab met;
      LegacyLife life;
      LegacyMot mot;
      LegacyGen gen;
    };

    LegacyBlobRec rec{};
    rec.id = {777, 2, 0};
    rec.sp.x = {5e-6, 6e-6, 7e-6};
    rec.sp.radius = 0.5e-6;
    rec.met.mu_max = 4e-4;
    rec.met.biomass = 1e-15;
    rec.life.state = static_cast<int32_t>(to_underlying(PhenoState::NORMAL));
    rec.life.sos_timer = -1.0;
    rec.life.in_crypt = 1;
    rec.gen.lineage_id = 777;
    rec.gen.num_bi_loci = 0;

    const int32_t legacy_count = 1;
    std::vector<char> legacy_buf(sizeof(legacy_count) + sizeof(rec));
    std::memcpy(legacy_buf.data(), &legacy_count, sizeof(legacy_count));
    std::memcpy(legacy_buf.data() + sizeof(legacy_count), &rec, sizeof(rec));

    const std::vector<Agent> out = agent_transfer_deserialize(legacy_buf);
    assert(out.size() == 1);
    assert(out[0].identity.tag == 777);
    assert(out[0].flags.in_crypt);
    assert(!out[0].flags.has_prophage_immunity);
    assert(!out[0].flags.has_inducible_prophage);
    assert(!out[0].flags.prophage_induced);
    assert(std::abs(out[0].timers.prophage_lysis_timer - (-1.0)) < 1e-12);
  }
  std::cout << "  test_latent_lysis_arrests_division_and_legacy_transfer_blob: PASSED\n";
}

}  // namespace

int main() {
  std::cout << "=== Spec 14 Temperate Prophage Tests (GAP-10) ===\n";
  test_default_off_bitwise_parity();
  test_cryptic_vs_active_lysogen_and_superinfection_immunity();
  test_hill_response_monotonicity();
  test_cell_and_virion_conservation_ledger();
  test_latent_lysis_delay();
  test_dual_regime_presets_and_realized_daily_rate();
  test_henrot_in_vitro_vs_in_vivo_contrast();
  test_mpi_transfer_and_hdf5_checkpoint_roundtrip();
  test_mechanism_capability_refuses_prophage_on_tabulated_layers();
  test_cognate_nuclease_immunity_specificity();
  test_latent_lysis_arrests_division_and_legacy_transfer_blob();
  std::cout << "All Spec 14 temperate prophage tests passed.\n";
  return 0;
}
