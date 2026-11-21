/* -----------------------------------------------------------------------
   GutIBM – Unit tests for fix_bacteriocin (Spec 2)
   ----------------------------------------------------------------------- */

#include "fix_bacteriocin.h"
#include "founder_genotype.h"
#include "plasmid.h"
#include "simulation.h"
#include "input_parser.h"
#include "qssa_solver.h"
#include "domain.h"
#include "advection.h"
#include "chemical_field.h"
#include "species_names.h"
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace gutibm;

static Simulation make_empty_sim(uint64_t seed = 42) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.seed = seed;

  Simulation sim;
  sim.init(cfg);
  return sim;
}

static Agent make_agent_at_center(Simulation& sim, Int type) {
  Vec3 center = {
    0.5 * (sim.domain().lo()[0] + sim.domain().hi()[0]),
    0.5 * (sim.domain().lo()[1] + sim.domain().hi()[1]),
    0.5 * (sim.domain().lo()[2] + sim.domain().hi()[2]),
  };
  Agent a = Agent::create_default(sim.agents().next_tag(), type, center, 5e-4);
  Int ix;
  Int iy;
  Int iz;
  sim.domain().pos_to_grid(a.x, ix, iy, iz);
  a.grid_cell = sim.domain().cell_index(ix, iy, iz);
  return a;
}

static ToxinBurstSource burst_from_cluster(
    const BICluster& bi, const Vec3& position) {
  ToxinBurstSource burst;
  burst.pos = position;
  burst.params.diff_coeff = bi.diff_coeff;
  burst.params.retardation = bi.retardation;
  burst.params.pI = bi.pI;
  burst.params.source_rate = 1.0e-18;
  burst.release_tau = 300.0;
  burst.is_nuclease = bi.is_nuclease;
  burst.target = bi.target;
  return burst;
}

void test_pi_diffusion_classes() {
  auto e1 = PlasmidLibrary::colicin_E1();
  auto e2 = PlasmidLibrary::colicin_E2();
  auto col_b = PlasmidLibrary::colicin_B();
  auto mcc = PlasmidLibrary::microcin_V();

  assert(e1.bclass == BacteriocinClass::LETHAL_CORE);
  assert(e1.pI > 8.5);
  assert(e1.release_mode == ReleaseMode::SOS_LYSIS);
  assert(e1.burst_size == 1.0e5);

  assert(e2.is_nuclease);
  assert(e2.release_mode == ReleaseMode::SOS_LYSIS);

  assert(col_b.bclass == BacteriocinClass::LETHAL_HALO);
  assert(col_b.pI < 7.0);
  assert(col_b.release_mode == ReleaseMode::PHAGE_LYSIS);

  assert(mcc.molecular_weight < 10000.0);
  assert(mcc.release_mode == ReleaseMode::CONTINUOUS);

  std::cout << "  test_pi_diffusion_classes: PASSED\n";
}

void test_pi_retardation_sensitivity_and_calibration() {
  const MucinChargeConfig cfg;
  const std::array<Real, 5> pIs = {6.0, 7.0, 8.0, 9.0, 10.0};
  Real previous = retardation_from_pI(pIs.front(), cfg);
  for (size_t i = 1; i < pIs.size(); ++i) {
    const Real current = retardation_from_pI(pIs[i], cfg);
    assert(std::isfinite(current));
    assert(current > previous);
    previous = current;
  }

  for (int i = 0; i <= 18; ++i) {
    const Real pI = 3.0 + static_cast<Real>(i) * 0.5;
    const Real value = retardation_from_pI(pI, cfg);
    assert(std::isfinite(value));
    assert(value >= cfg.r_min);
    assert(value <= cfg.r_min + cfg.amplitude);
  }

  MucinChargeConfig no_charge = cfg;
  no_charge.amplitude = 0.0;
  const Real constant = retardation_from_pI(3.0, no_charge);
  for (const Real pI : pIs) {
    assert(std::abs(retardation_from_pI(pI, no_charge) - constant) < 1e-12);
  }

  MucinChargeConfig alkaline = cfg;
  alkaline.ph = 8.0;
  assert(retardation_from_pI(8.0, alkaline)
         < retardation_from_pI(8.0, cfg));
  assert(retardation_from_pI(8.0, cfg)
         < retardation_from_pI(8.0, MucinChargeConfig{
             .r_min = cfg.r_min,
             .amplitude = cfg.amplitude,
             .dz_half = cfg.dz_half,
             .width = cfg.width,
             .ph = 6.0}));

  struct Calibration {
    Real pI;
    Real library;
    Real diff_coeff;
  };
  const std::array<Calibration, 6> calibration = {{
      {5.0, 1.2, 1.0e-10},
      {5.4, 1.5, 4.0e-11},
      {6.5, 3.0, 3.5e-11},
      {7.2, 5.0, 4.0e-11},
      {9.0, 50.0, 4.0e-11},
      {9.3, 60.0, 5.0e-11},
  }};
  for (const Calibration& row : calibration) {
    const Real fitted = retardation_from_pI(row.pI, cfg);
    // ColE2 uses the toxin-plus-Im2 complex pI, so this is a calibration
    // change-detector with one expected off-curve library entry.
    assert(std::abs(fitted - row.library) / row.library < 0.35);
    const Real d_eff = row.diff_coeff / fitted;
    assert(std::isfinite(d_eff) && d_eff > 0.0);
  }

  std::cout << "  test_pi_retardation_sensitivity_and_calibration: PASSED\n";
}

void test_microcin_mu_penalty() {
  BacteriocinConfig cfg;
  cfg.microcin_mu_penalty = 0.03;

  auto sim = make_empty_sim();
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::microcin_V());
  const Real mu_before = a.mu_max;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);
  fix.compute(60.0);

  const Real mu_after = sim.agents()[0].mu_max;
  assert(std::abs(mu_after - mu_before * (1.0 - cfg.microcin_mu_penalty)) < 1e-12);

  std::cout << "  test_microcin_mu_penalty: PASSED\n";
}

void test_microcin_mu_penalty_persists() {
  BacteriocinConfig cfg;
  cfg.microcin_mu_penalty = 0.03;

  auto sim = make_empty_sim();
  Agent producer = make_agent_at_center(sim, 1);
  producer.genome.bi_loci.push_back(PlasmidLibrary::microcin_V());
  const Real producer_mu_before = producer.mu_max;
  sim.agents().push_back(std::move(producer));
  Agent control = make_agent_at_center(sim, 2);
  const Real control_mu_before = control.mu_max;
  sim.agents().push_back(std::move(control));

  for (int step = 0; step < 10; ++step) {
    sim.step(60.0);
  }

  const Real producer_mu_after = sim.agents()[0].mu_max;
  assert(std::abs(producer_mu_after
                  - producer_mu_before * (1.0 - cfg.microcin_mu_penalty))
          < 1e-12);
  assert(std::abs(producer_mu_after
                  - producer_mu_before
                      * std::pow(1.0 - cfg.microcin_mu_penalty, 10))
          > 1e-8);
  assert(std::abs(sim.agents()[1].mu_max - control_mu_before) < 1e-12);

  std::cout << "  test_microcin_mu_penalty_persists: PASSED\n";
}

void test_sos_induction_requires_bi_loci() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 1.0;

  auto sim = make_empty_sim(1001);
  Agent a = make_agent_at_center(sim, 1);
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::NORMAL);

  std::cout << "  test_sos_induction_requires_bi_loci: PASSED\n";
}

void test_sos_induction_high_basal_rate() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 1.0;

  auto sim = make_empty_sim(2002);
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::SOS_INDUCED);
  assert(sim.agents()[0].timers.sos_timer > 0.0);

  std::cout << "  test_sos_induction_high_basal_rate: PASSED\n";
}

void test_steady_state_qssa_source() {
  // A continuously-secreted microcin (MccV, ReleaseMode::CONTINUOUS) must show up
  // as a steady-state QSSA source at its producer, while a lysis-only colicin
  // (ColE1, ReleaseMode::SOS_LYSIS) contributes NO continuous source without a
  // lysis burst — colicins reach the field via ToxinBurstSource (see the SOS
  // lysis / burst tests), not via collect_microcin_sources.
  const BICluster mcc_v = PlasmidLibrary::microcin_V();
  assert(mcc_v.release_mode == ReleaseMode::CONTINUOUS);

  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {100e-6, 100e-6, 50e-6};
  dcfg.grid_dx = 5e-6;
  Domain domain;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.mucus_thickness = 50e-6;
  AdvectionField adv;
  adv.init(acfg, domain);

  QSSAConfig qcfg;
  qcfg.toxin_cutoff = 80e-6;
  qcfg.microcin_secretion = 1.0e-20;
  QSSASolver qssa;
  qssa.init(qcfg, domain, adv);

  ChemicalSpec toxin;
  toxin.name = species::BACTERIOCIN_CIRA;
  toxin.diff_coeff = 4.0e-11;
  toxin.retardation = 10.0;

  const Vec3 pos = {50e-6, 50e-6, 25e-6};
  Int ix = 0;
  Int iy = 0;
  Int iz = 0;
  domain.pos_to_grid(pos, ix, iy, iz);
  const Int idx = domain.cell_index(ix, iy, iz);

  ProteaseConfig protease;
  protease.enabled = false;
  const std::vector<ToxinBurstSource> no_bursts;

  // Continuous microcin secretor → nonzero steady-state field.
  {
    ChemicalField chem;
    chem.init(domain, {toxin});
    AgentPool agents;
    Agent producer = Agent::create_default(1, 1, pos, 5e-4);
    producer.genome.bi_loci.push_back(mcc_v);
    agents.push_back(std::move(producer));
    qssa.solve_bacteriocin_field(agents, no_bursts, 0.0, protease, adv, chem, 0,
                                 ReceptorType::CirA);
    assert(chem.conc(0, idx) > 0.0);
    std::cout << "  test_steady_state_qssa_source: microcin c=" << chem.conc(0, idx) << "\n";
  }

  // Lysis-only colicin (no burst) → no continuous source.
  {
    const BICluster col_e1 = PlasmidLibrary::colicin_E1();
    assert(col_e1.release_mode == ReleaseMode::SOS_LYSIS);
    ChemicalField chem;
    chem.init(domain, {toxin});
    AgentPool agents;
    Agent producer = Agent::create_default(2, 1, pos, 5e-4);
    producer.genome.bi_loci.push_back(col_e1);
    agents.push_back(std::move(producer));
    qssa.solve_bacteriocin_field(agents, no_bursts, 0.0, protease, adv, chem, 0,
                                 ReceptorType::CirA);
    assert(chem.conc(0, idx) == 0.0);
  }

  std::cout << "  test_steady_state_qssa_source: PASSED\n";
}

void test_sos_induction_on_division() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 0.0;
  cfg.sos_lysis_prob = 1.0;

  auto sim = make_empty_sim(3003);
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  a.flags.just_divided = true;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::SOS_INDUCED);

  std::cout << "  test_sos_induction_on_division: PASSED\n";
}

void test_no_sos_without_division() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 0.0;
  cfg.sos_lysis_prob = 1.0;

  auto sim = make_empty_sim(4004);
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  a.flags.just_divided = false;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::NORMAL);

  std::cout << "  test_no_sos_without_division: PASSED\n";
}

void test_phage_induction() {
  BacteriocinConfig cfg;

  auto sim = make_empty_sim(5005);
  Agent a = make_agent_at_center(sim, 1);
  BICluster col_b = PlasmidLibrary::colicin_B();
  col_b.phage_induction_rate = 100.0;
  a.genome.bi_loci.push_back(col_b);
  a.mu_realized = 5.0e-4;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::SOS_INDUCED);

  std::cout << "  test_phage_induction: PASSED\n";
}

void test_phage_does_not_trigger_sos_path() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 1.0;

  auto sim = make_empty_sim(6006);
  Agent a = make_agent_at_center(sim, 1);
  BICluster col_b = PlasmidLibrary::colicin_B();
  col_b.phage_induction_rate = 0.0;
  a.genome.bi_loci.push_back(col_b);
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::NORMAL);

  std::cout << "  test_phage_does_not_trigger_sos_path: PASSED\n";
}

void test_microcin_no_lysis() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 1.0;

  auto sim = make_empty_sim(7007);
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::microcin_V());
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::NORMAL);

  std::cout << "  test_microcin_no_lysis: PASSED\n";
}

void test_cross_induction() {
  BacteriocinConfig cfg;
  cfg.sos_basal_rate = 0.0;
  cfg.sos_lysis_prob = 0.0;
  cfg.sos_cross_induction_rate = 1.0e3;

  auto sim = make_empty_sim(8008);
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  sim.agents().push_back(std::move(a));
  ToxinBurstSource burst = burst_from_cluster(
      PlasmidLibrary::colicin_E2(), sim.agents()[0].x);
  burst.params.source_rate = 1.0e-12;
  sim.add_toxin_burst(burst);
  sim.qssa().solve_all_bacteriocin_fields(
      sim.agents(), sim.toxin_bursts(), 0.0,
      sim.config().chem_env.protease, sim.advection(),
      sim.chemical_field(), nullptr, false);

  FixBacteriocin fix(sim, cfg);
  assert(fix.nuclease_cross_induction_rate(sim.agents()[0], 0) > 0.0);
  fix.compute(60.0);

  assert(sim.agents()[0].state == PhenoState::SOS_INDUCED);

  std::cout << "  test_cross_induction: PASSED\n";
}

void test_per_colicin_burst_size() {
  BacteriocinConfig cfg;
  cfg.burst_release_tau = 300.0;

  auto sim = make_empty_sim();
  Agent a = make_agent_at_center(sim, 1);
  BICluster e1 = PlasmidLibrary::colicin_E1();
  assert(e1.burst_size == 1.0e5);
  a.genome.bi_loci.push_back(e1);
  a.state = PhenoState::SOS_INDUCED;
  a.timers.sos_timer = 1.0;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.post_step(2.0);

  assert(sim.toxin_bursts().size() == 1);
  const Real expected = (1.0e5 / AVOGADRO) / cfg.burst_release_tau;
  assert(std::abs(sim.toxin_bursts()[0].params.source_rate - expected) < 1e-30);
  assert(std::abs(sim.toxin_bursts()[0].release_tau - cfg.burst_release_tau) < 1e-12);
  assert(sim.toxin_bursts()[0].is_nuclease == false);

  std::cout << "  test_per_colicin_burst_size: PASSED\n";
}

void test_sos_lysis_post_step() {
  BacteriocinConfig cfg;

  auto sim = make_empty_sim();
  Agent a = make_agent_at_center(sim, 1);
  a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
  a.state = PhenoState::SOS_INDUCED;
  a.timers.sos_timer = 30.0;
  sim.agents().push_back(std::move(a));

  FixBacteriocin fix(sim, cfg);
  fix.post_step(60.0);

  assert(sim.agents()[0].state == PhenoState::DEAD);

  std::cout << "  test_sos_lysis_post_step: PASSED\n";
}

static SimulationConfig funded_ros_config() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.oxygen.delivery_uptake_enabled = true;
  cfg.chem_env.oxygen.ros_driver = "funded";
  cfg.fixes.metabolism.uptake_limit = "delivery";
  return cfg;
}

static SimulationConfig ambient_ros_config() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.initial_strains.clear();
  cfg.hdf5.enabled = false;
  cfg.domain.hi = {50e-6, 50e-6, 25e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.oxygen.ros_driver = "ambient";
  return cfg;
}

static Real ambient_ros_rate(Real k_ros) {
  SimulationConfig cfg = ambient_ros_config();
  cfg.chem_env.oxygen.k_ROS = k_ros;
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  Agent agent = make_agent_at_center(sim, 1);
  sim.agents().push_back(std::move(agent));
  Agent& probe = sim.agents()[sim.agents().size() - 1];
  probe.mu_realized = 1.0e-4;
  assert(sim.local_O2(probe) > 0.0);
  return sim.ros_induction_rate(probe);
}

void test_ros_ambient_default_is_disabled() {
  SimulationConfig cfg = ambient_ros_config();
  assert(cfg.chem_env.oxygen.k_ROS == 0.0);
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  Agent agent = make_agent_at_center(sim, 1);
  sim.agents().push_back(std::move(agent));
  Agent& probe = sim.agents()[sim.agents().size() - 1];
  probe.mu_realized = 1.0e-4;
  assert(sim.local_O2(probe) > 0.0);
  // CHANGE DETECTOR: the shipped default intentionally retires ambient ROS mortality.
  assert(sim.ros_induction_rate(probe) == 0.0);

  std::cout << "  test_ros_ambient_default_is_disabled: PASSED\n";
}

void test_ros_ambient_explicit_rate_sensitivity() {
  const std::array<Real, 3> coefficients = {1.0e1, 1.0e2, 2.0e2};
  std::array<Real, 3> rates = {};
  size_t index = 0;
  for (const Real coefficient : coefficients) {
    rates[index++] = ambient_ros_rate(coefficient);
  }

  assert(rates[0] > 0.0);
  assert(rates[1] > 0.0);
  assert(rates[2] > 0.0);
  assert(rates[1] > rates[0] * 9.9);
  assert(rates[2] > rates[1] * 1.9);
  assert(std::abs(rates[1] / rates[0] - 10.0) < 1.0e-12);
  assert(std::abs(rates[2] / rates[1] - 2.0) < 1.0e-12);

  std::cout << "  test_ros_ambient_explicit_rate_sensitivity: rates="
            << rates[0] << "," << rates[1] << "," << rates[2] << "\n";
}

void test_ros_funded_absolute_rate_sensitivity() {
  SimulationConfig cfg = funded_ros_config();
  cfg.chem_env.oxygen.k_ROS_funded = 2.0;
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  Agent agent = make_agent_at_center(sim, 1);
  sim.agents().push_back(std::move(agent));
  Agent& probe = sim.agents()[sim.agents().size() - 1];
  probe.biomass = 1.0e-15;
  const std::array<Real, 3> fluxes = {1.0e-20, 2.0e-20, 4.0e-20};
  std::array<Real, 3> rates = {};
  size_t index = 0;
  for (const Real flux : fluxes) {
    probe.respired_oxygen_rate = flux;
    rates[index++] = sim.ros_induction_rate(probe);
  }
  assert(rates[1] > rates[0] * 1.9);
  assert(rates[2] > rates[1] * 1.9);
  assert(std::abs(rates[2] / rates[0] - 4.0) < 1.0e-12);

  probe.biomass = 1.0e-14;
  probe.respired_oxygen_rate = 2.0e-20;
  const Real larger_biomass_rate = sim.ros_induction_rate(probe);
  probe.biomass = 1.0e-15;
  const Real smaller_biomass_rate = sim.ros_induction_rate(probe);
  assert(std::abs(larger_biomass_rate - smaller_biomass_rate) < 1.0e-30);

  probe.biomass = 0.0;
  assert(sim.ros_induction_rate(probe) > 0.0);
  probe.respired_oxygen_rate = 0.0;
  assert(sim.ros_induction_rate(probe) == 0.0);
  std::cout << "  test_ros_funded_absolute_rate_sensitivity: PASSED\n";
}

void test_ros_funded_specific_rate_preserved() {
  SimulationConfig cfg = funded_ros_config();
  cfg.chem_env.oxygen.k_ROS_respiratory = 2.0;
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  Agent agent = make_agent_at_center(sim, 1);
  sim.agents().push_back(std::move(agent));
  Agent& probe = sim.agents()[sim.agents().size() - 1];
  probe.respired_oxygen_rate = 2.0e-20;
  probe.biomass = 1.0e-15;
  const Real smaller_biomass_rate = sim.ros_induction_rate(probe);
  probe.biomass = 1.0e-14;
  const Real larger_biomass_rate = sim.ros_induction_rate(probe);
  assert(smaller_biomass_rate > larger_biomass_rate * 9.9);
  assert(smaller_biomass_rate < larger_biomass_rate * 10.1);
  std::cout << "  test_ros_funded_specific_rate_preserved: PASSED\n";
}

void test_ros_funded_generation_round_trip() {
  constexpr Real q_o2_per_generation = 1.63e-14;
  constexpr Real generation_time = 8172.0;
  const std::array<Real, 3> probabilities = {0.01, 0.02, 0.05};
  std::vector<Real> hazards;
  hazards.reserve(probabilities.size());
  for (const Real probability : probabilities) {
    SimulationConfig cfg = funded_ros_config();
    cfg.chem_env.oxygen.k_ROS_funded =
        -std::log1p(-probability) / q_o2_per_generation;
    InputParser::finalize_config(cfg);

    Simulation sim;
    sim.init(cfg);
    Agent agent = make_agent_at_center(sim, 1);
    sim.agents().push_back(std::move(agent));
    Agent& probe = sim.agents()[sim.agents().size() - 1];
    probe.respired_oxygen_rate = q_o2_per_generation / generation_time;
    const Real hazard = sim.ros_induction_rate(probe);
    const Real integrated_hazard = hazard * generation_time;
    const Real target_hazard = -std::log1p(-probability);
    assert(std::abs(integrated_hazard - target_hazard) / target_hazard
           < 1.0e-9);
    const Real recovered_probability = 1.0 - std::exp(-integrated_hazard);
    assert(std::abs(recovered_probability - probability) / probability
           < 1.0e-9);
    hazards.push_back(hazard);
  }
  assert(hazards[1] > hazards[0]);
  assert(hazards[2] > hazards[1]);
  std::cout << "  test_ros_funded_generation_round_trip: PASSED\n";
}

void test_sos_division_pair_probability_binomial() {
  constexpr Int num_divisions = 60000;
  constexpr Real bio_dt = 60.0;
  const std::array<Real, 2> target_probs = {0.01, 0.10};

  for (size_t idx = 0; idx < target_probs.size(); ++idx) {
    const Real p = target_probs[idx];
    SimulationConfig cfg = InputParser::default_config();
    cfg.initial_strains.clear();
    cfg.hdf5.enabled = false;
    cfg.domain.hi = {50e-6, 50e-6, 25e-6};
    cfg.domain.grid_dx = 5e-6;
    cfg.time.bio_dt = bio_dt;
    cfg.seed = 9001 + idx * 100;
    cfg.chem_env.oxygen.enabled = false;
    cfg.chem_env.oxygen.k_ROS = 0.0;
    cfg.fixes.bacteriocin.sos_basal_rate = 0.0;
    cfg.fixes.bacteriocin.sos_cross_induction_rate = 0.0;
    cfg.fixes.bacteriocin.sos_lysis_prob = p;
    InputParser::finalize_config(cfg);

    Simulation sim;
    sim.init(cfg);
    const BICluster col_e1 = PlasmidLibrary::colicin_E1();

    for (Int i = 0; i < 2 * num_divisions; ++i) {
      Agent daughter = make_agent_at_center(sim, 1);
      daughter.genome.bi_loci.push_back(col_e1);
      daughter.flags.just_divided = true;
      sim.agents().push_back(std::move(daughter));
    }

    FixBacteriocin fix(sim, sim.config().fixes.bacteriocin);
    fix.compute(bio_dt);

    Int pair_events = 0;
    Int induced_cells = 0;
    for (Int pair = 0; pair < num_divisions; ++pair) {
      const bool d0 =
          sim.agents()[2 * pair].state == PhenoState::SOS_INDUCED;
      const bool d1 =
          sim.agents()[2 * pair + 1].state == PhenoState::SOS_INDUCED;
      if (d0 || d1) ++pair_events;
      if (d0) ++induced_cells;
      if (d1) ++induced_cells;
    }

    const Real empirical_pair_prob =
        static_cast<Real>(pair_events) / static_cast<Real>(num_divisions);
    const Real se =
        std::sqrt(p * (1.0 - p) / static_cast<Real>(num_divisions));
    const Real legacy_unnormalized_pair_prob = 1.0 - std::exp(-2.0 * p);
    assert(std::abs(empirical_pair_prob - p) <= 3.0 * se);
    assert(std::abs(empirical_pair_prob - legacy_unnormalized_pair_prob)
           > 5.0 * se);

    const Real empirical_cells_per_div =
        static_cast<Real>(induced_cells) / static_cast<Real>(num_divisions);
    const Real expected_cells_per_div = 2.0 * (1.0 - std::sqrt(1.0 - p));
    assert(std::abs(empirical_cells_per_div - expected_cells_per_div)
           <= 3.5 * se);
  }

  std::cout << "  test_sos_division_pair_probability_binomial: PASSED\n";
}

void test_sos_division_probability_boundaries() {
  constexpr Real bio_dt = 60.0;
  // p = 0.0 -> zero inductions and zero post-division hazard
  {
    SimulationConfig cfg = InputParser::default_config();
    cfg.initial_strains.clear();
    cfg.hdf5.enabled = false;
    cfg.domain.hi = {50e-6, 50e-6, 25e-6};
    cfg.domain.grid_dx = 5e-6;
    cfg.time.bio_dt = bio_dt;
    cfg.seed = 9401;
    cfg.chem_env.oxygen.enabled = false;
    cfg.fixes.bacteriocin.sos_basal_rate = 0.0;
    cfg.fixes.bacteriocin.sos_cross_induction_rate = 0.0;
    cfg.fixes.bacteriocin.sos_lysis_prob = 0.0;
    InputParser::finalize_config(cfg);

    Simulation sim;
    sim.init(cfg);
    for (Int i = 0; i < 200; ++i) {
      Agent a = make_agent_at_center(sim, 1);
      a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
      a.flags.just_divided = true;
      sim.agents().push_back(std::move(a));
    }
    FixBacteriocin fix(sim, sim.config().fixes.bacteriocin);
    fix.compute(bio_dt);
    assert(sim.step_events().sos_inductions == 0);
    assert(sim.step_events().sos_post_division_rate == 0.0);
  }

  // p = 1.0 -> all post-division cells induce with finite rate (no NaN/Inf)
  {
    SimulationConfig cfg = InputParser::default_config();
    cfg.initial_strains.clear();
    cfg.hdf5.enabled = false;
    cfg.domain.hi = {50e-6, 50e-6, 25e-6};
    cfg.domain.grid_dx = 5e-6;
    cfg.time.bio_dt = bio_dt;
    cfg.seed = 9402;
    cfg.chem_env.oxygen.enabled = false;
    cfg.fixes.bacteriocin.sos_basal_rate = 0.0;
    cfg.fixes.bacteriocin.sos_cross_induction_rate = 0.0;
    cfg.fixes.bacteriocin.sos_lysis_prob = 1.0;
    InputParser::finalize_config(cfg);

    Simulation sim;
    sim.init(cfg);
    for (Int i = 0; i < 200; ++i) {
      Agent a = make_agent_at_center(sim, 1);
      a.genome.bi_loci.push_back(PlasmidLibrary::colicin_E1());
      a.flags.just_divided = true;
      sim.agents().push_back(std::move(a));
    }
    FixBacteriocin fix(sim, sim.config().fixes.bacteriocin);
    fix.compute(bio_dt);
    assert(sim.step_events().sos_inductions == 200);
    assert(std::isfinite(sim.step_events().sos_post_division_rate));
    assert(sim.step_events().sos_post_division_rate > 0.0);
  }

  // Out-of-range values (< 0 or > 1) are rejected by InputParser
  for (const char* bad_val : {"-0.01", "1.01"}) {
    SimulationConfig cfg = InputParser::default_config();
    bool threw = false;
    try {
      InputParser::apply_flat_key(cfg, "sos_lysis_prob", bad_val);
    } catch (const ConfigError&) {
      threw = true;
    }
    assert(threw);
  }

  std::cout << "  test_sos_division_probability_boundaries: PASSED\n";
}

void test_d_free_colicin_runtime_propagation() {
  // 1. Default config preserves library per-plasmid diff_coeff and 4e-11 field spec.
  {
    const SimulationConfig def = InputParser::default_config();
    const PlasmidEntry* e2_entry = PlasmidLibrary::find("ColE2");
    const PlasmidEntry* mcc_entry = PlasmidLibrary::find("MccV");
    assert(e2_entry != nullptr && mcc_entry != nullptr);
    const BICluster e2_def = configured_plasmid_cluster(def, *e2_entry);
    const BICluster mcc_def = configured_plasmid_cluster(def, *mcc_entry);
    assert(std::abs(e2_def.diff_coeff - 3.5e-11) < 1e-20);
    assert(std::abs(mcc_def.diff_coeff - 1.0e-10) < 1e-20);
  }

  // 2. Setting bacteriocin.D_free_colicin propagates to chemical species,
  //    QSSA fallback, configured colicin clusters, and runtime lysis bursts,
  //    while explicit plasmid_overrides.<name>.diff_coeff takes precedence.
  const std::array<Real, 3> d_values = {2.0e-11, 4.0e-11, 7.0e-11};
  const std::array<const char*, 3> d_strings = {"2.0e-11", "4.0e-11", "7.0e-11"};
  std::array<Real, 3> near_source_concs{};

  for (size_t idx = 0; idx < d_values.size(); ++idx) {
    SimulationConfig cfg = InputParser::default_config();
    cfg.initial_strains.clear();
    cfg.hdf5.enabled = false;
    cfg.domain.hi = {50e-6, 50e-6, 25e-6};
    cfg.domain.grid_dx = 5e-6;
    cfg.seed = 9500;
    InputParser::apply_flat_key(
        cfg, "bacteriocin.D_free_colicin", d_strings[idx]);
    cfg.plasmid_overrides["ColE2"].retardation = 2.0;

    SimulationConfig::InitialStrain producer;
    producer.type = 1;
    producer.count = 1;
    producer.mu_max = 5.0e-4;
    producer.plasmids = {"ColE2"};
    cfg.initial_strains.push_back(producer);

    Simulation sim;
    sim.init(cfg);

    for (const char* species_name : {
             species::BACTERIOCIN_BTUB, species::BACTERIOCIN_FEPA,
             species::BACTERIOCIN_CIRA, species::BACTERIOCIN_FHUA}) {
      const Int sp_idx = sim.chemical_field().find(species_name);
      assert(sp_idx >= 0);
      assert(std::abs(sim.chemical_field().spec(sp_idx).diff_coeff
                      - d_values[idx])
             < 1e-20);
    }
    assert(std::abs(sim.config().qssa.fallback_diff_coeff - d_values[idx])
           < 1e-20);

    assert(sim.agents().size() == 1);
    Agent& a = sim.agents()[0];
    assert(a.genome.bi_loci.size() == 1);
    assert(std::abs(a.genome.bi_loci[0].diff_coeff - d_values[idx]) < 1e-20);

    // Force lysis and solve QSSA field to verify runtime field sensitivity to D.
    a.x = {25e-6, 25e-6, 12.5e-6};
    a.state = PhenoState::SOS_INDUCED;
    a.timers.sos_timer = 1.0;
    FixBacteriocin fix(sim, sim.config().fixes.bacteriocin);
    fix.post_step(2.0);
    assert(sim.toxin_bursts().size() == 1);
    assert(std::abs(sim.toxin_bursts()[0].params.diff_coeff - d_values[idx])
           < 1e-20);

    sim.qssa().solve_all_bacteriocin_fields(
        sim.agents(), sim.toxin_bursts(), sim.time(),
        sim.config().chem_env.protease, sim.advection(),
        sim.chemical_field(), nullptr, true);
    const Int btub_idx =
        sim.chemical_field().find(species::BACTERIOCIN_BTUB);
    Int ix = 0;
    Int iy = 0;
    Int iz = 0;
    sim.domain().pos_to_grid({30e-6, 25e-6, 12.5e-6}, ix, iy, iz);
    near_source_concs[idx] =
        sim.chemical_field().conc(btub_idx, sim.domain().cell_index(ix, iy, iz));
    assert(std::isfinite(near_source_concs[idx])
           && near_source_concs[idx] > 0.0);
  }

  // Higher free diffusivity spreads the point source faster, lowering peak
  // concentration near the source (1 / (4 * pi * D_eff * r) scaling).
  assert(near_source_concs[0] > near_source_concs[1]);
  assert(near_source_concs[1] > near_source_concs[2]);

  // 3. Explicit plasmid_overrides.<name>.diff_coeff overrides D_free_colicin
  //    for that plasmid while other colicins inherit D_free_colicin.
  {
    SimulationConfig cfg = InputParser::default_config();
    InputParser::apply_flat_key(cfg, "bacteriocin.D_free_colicin", "7.0e-11");
    InputParser::apply_flat_key(
        cfg, "plasmid_overrides.ColE2.diff_coeff", "9.0e-11");
    InputParser::finalize_config(cfg);

    const PlasmidEntry* e1 = PlasmidLibrary::find("ColE1");
    const PlasmidEntry* e2 = PlasmidLibrary::find("ColE2");
    const PlasmidEntry* mcc = PlasmidLibrary::find("MccV");
    assert(e1 != nullptr && e2 != nullptr && mcc != nullptr);
    assert(std::abs(configured_plasmid_cluster(cfg, *e1).diff_coeff - 7.0e-11)
           < 1e-20);
    assert(std::abs(configured_plasmid_cluster(cfg, *e2).diff_coeff - 9.0e-11)
           < 1e-20);
    assert(std::abs(configured_plasmid_cluster(cfg, *mcc).diff_coeff - 1.0e-10)
           < 1e-20);
  }

  std::cout << "  test_d_free_colicin_runtime_propagation: PASSED\n";
}

int main() {
  std::cout << "=== Bacteriocin Fix Tests ===\n";
  test_pi_diffusion_classes();
  test_pi_retardation_sensitivity_and_calibration();
  test_microcin_mu_penalty();
  test_microcin_mu_penalty_persists();
  test_sos_induction_requires_bi_loci();
  test_sos_induction_high_basal_rate();
  test_steady_state_qssa_source();
  test_sos_induction_on_division();
  test_no_sos_without_division();
  test_phage_induction();
  test_phage_does_not_trigger_sos_path();
  test_microcin_no_lysis();
  test_cross_induction();
  test_per_colicin_burst_size();
  test_sos_lysis_post_step();
  test_ros_ambient_default_is_disabled();
  test_ros_ambient_explicit_rate_sensitivity();
  test_ros_funded_absolute_rate_sensitivity();
  test_ros_funded_specific_rate_preserved();
  test_ros_funded_generation_round_trip();
  test_sos_division_pair_probability_binomial();
  test_sos_division_probability_boundaries();
  test_d_free_colicin_runtime_propagation();
  std::cout << "All bacteriocin fix tests passed.\n";
  return 0;
}
