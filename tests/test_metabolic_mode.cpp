#include "metabolic_mode.h"

#include "hdf5_reader.h"
#include "input_parser.h"
#include "path_utils.h"
#include "simulation.h"
#include "species_names.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

namespace {

using gutibm::Agent;
using gutibm::FixMetabolism;
using gutibm::InputParser;
using gutibm::Real;
using gutibm::Simulation;
using gutibm::SimulationConfig;
namespace species = gutibm::species;

constexpr Real kDt = 60.0;

SimulationConfig simulation_config() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.hi = {100.0e-6, 100.0e-6, 40.0e-6};
  cfg.hdf5.enabled = false;
  cfg.time.total_time = kDt;
  cfg.time.bio_dt = kDt;
  cfg.time.output_interval = kDt;
  cfg.dysbiosis_threshold = 0.0;
  cfg.enabled_fixes = {"metabolism"};
  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain strain;
  strain.type = 1;
  strain.count = 1;
  strain.mu_max = 1.0e-4;
  cfg.initial_strains.push_back(strain);
  return cfg;
}

void set_species(Simulation& sim, const char* name, Real value) {
  const gutibm::Int species = sim.chemical_field().find(name);
  assert(species >= 0);
  for (gutibm::Int cell = 0;
       cell < sim.chemical_field().ncells(); ++cell) {
    sim.chemical_field().conc(species, cell) = value;
  }
}

void set_agent_cells(Simulation& sim) {
  assert(sim.agents().size() >= 2);
  const gutibm::Int first = sim.agents()[0].grid_cell;
  assert(first + 1 < sim.chemical_field().ncells());
  sim.agents()[1].grid_cell = first + 1;
}

struct Trajectory {
  std::vector<Real> carbon;
  std::vector<Real> biomass;
};

struct RespirationProbe {
  Real fraction = 0.0;
  Real growth_demand = 0.0;
  Real funded_growth = 0.0;
};

RespirationProbe run_respiration_probe(
    const std::string& driver, Real oxygen_concentration,
    bool zero_growth = false) {
  SimulationConfig cfg = simulation_config();
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.oxygen.delivery_uptake_enabled = true;
  cfg.chem_env.oxygen.respiration_driver = driver;
  cfg.chem_env.oxygen.metabolic_switch_enabled = true;
  cfg.chem_env.oxygen.tau_metabolic_switch = 1.0;
  cfg.chem_env.oxygen.q_consumption = 1.0e-15;
  cfg.chem_env.oxygen.q_maintenance = 0.0;
  cfg.chem_env.oxygen.vbf_sink = 0.0;
  cfg.chem_env.oxygen.epithelial_conc = oxygen_concentration;
  cfg.fixes.metabolism.uptake_limit = "delivery";
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.fixes.metabolism.division_threshold = 1.0e9;
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  set_species(sim, species::CARBON, zero_growth ? 0.0 : 1.0);
  set_species(sim, species::OXYGEN, oxygen_concentration);
  if (zero_growth) {
    sim.agents()[0].realized_fermentation_fraction = 0.37;
  }
  sim.step(kDt);

  const gutibm::Int oxygen = sim.chemical_field().find(species::OXYGEN);
  const gutibm::Int cell = sim.agents()[0].grid_cell;
  RespirationProbe result;
  result.fraction = sim.agents()[0].realized_fermentation_fraction;
  result.growth_demand = sim.agents()[0].pending_oxygen_growth;
  result.funded_growth = sim.chemical_field().sink_realized_global(
      oxygen, cell);
  return result;
}

void test_funded_respiration_graded_and_bounded() {
  const std::array<Real, 4> concentrations = {
      0.0, 1.0e-5, 1.0e-3, 1.0e-1};
  std::array<Real, 4> fractions{};
  for (size_t i = 0; i < concentrations.size(); ++i) {
    const RespirationProbe probe = run_respiration_probe(
        "funded", concentrations[i]);
    fractions[i] = probe.fraction;
    assert(std::isfinite(probe.fraction));
    assert(probe.fraction >= 0.0 && probe.fraction <= 1.0);
    assert(probe.funded_growth <= probe.growth_demand + 1.0e-30);
  }
  assert(fractions[0] > fractions[1]);
  assert(fractions[1] > fractions[2]);
  assert(fractions[2] > fractions[3]);
  assert(fractions[0] - fractions[3] > 0.5);
  std::cout << "  test_funded_respiration_graded_and_bounded: PASSED\n";
}

void test_funded_respiration_discriminates_from_ambient() {
  const RespirationProbe ambient = run_respiration_probe(
      "ambient", 55.0e-6);
  const RespirationProbe funded = run_respiration_probe(
      "funded", 55.0e-6);
  assert(ambient.fraction < 0.2);
  assert(funded.fraction > ambient.fraction + 0.2);
  std::cout << "  test_funded_respiration_discriminates_from_ambient: PASSED\n";
}

void test_funded_respiration_zero_demand_preserves_state() {
  const RespirationProbe probe = run_respiration_probe(
      "funded", 55.0e-6, true);
  assert(std::abs(probe.fraction - 0.37) < 1.0e-12);
  assert(probe.growth_demand <= 0.0);
  assert(probe.funded_growth <= probe.growth_demand + 1.0e-30);
  std::cout << "  test_funded_respiration_zero_demand_preserves_state: PASSED\n";
}

Trajectory run_compatibility(bool explicit_disabled_keys) {
  SimulationConfig cfg = simulation_config();
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.acetate.enabled = true;
  cfg.time.total_time = 2.0 * kDt;
  cfg.initial_strains[0].count = 2;
  cfg.initial_strains[0].mu_max = 1.0e-4;
  if (explicit_disabled_keys) {
    assert(InputParser::apply_flat_key(
        cfg, "oxygen.metabolic_switch_enabled", "false"));
    assert(InputParser::apply_flat_key(
        cfg, "acid_inhibition_enabled", "false"));
  }
  InputParser::finalize_config(cfg);
  Simulation sim;
  sim.init(cfg);
  sim.run();

  const gutibm::Int carbon = sim.chemical_field().find(species::CARBON);
  Trajectory result;
  result.carbon.reserve(static_cast<size_t>(sim.chemical_field().ncells()));
  for (gutibm::Int cell = 0;
       cell < sim.chemical_field().ncells(); ++cell) {
    result.carbon.push_back(sim.chemical_field().conc(carbon, cell));
  }
  result.biomass.reserve(static_cast<size_t>(sim.agents().size()));
  for (const Agent& agent : sim.agents()) {
    result.biomass.push_back(agent.biomass);
  }
  return result;
}

void test_simulation_backward_compatibility() {
  const Trajectory absent = run_compatibility(false);
  const Trajectory explicit_disabled = run_compatibility(true);
  assert(absent.carbon == explicit_disabled.carbon);
  assert(absent.biomass == explicit_disabled.biomass);
  std::cout << "  test_simulation_backward_compatibility: PASSED\n";
}

void test_simulation_carbon_cost_direction() {
  SimulationConfig cfg = simulation_config();
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.acetate.enabled = false;
  cfg.chem_env.oxygen.metabolic_switch_enabled = true;
  cfg.chem_env.oxygen.tau_metabolic_switch = 1.0;
  cfg.chem_env.oxygen.aerobic_mu_factor = 1.0;
  cfg.chem_env.oxygen.anaerobic_mu_factor = 1.0;
  cfg.initial_strains[0].count = 2;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  set_agent_cells(sim);
  set_species(sim, species::CARBON, 1.0);
  set_species(sim, species::OXYGEN, 0.0);
  const gutibm::Int oxygen = sim.chemical_field().find(species::OXYGEN);
  sim.chemical_field().conc_global(
      oxygen, sim.agents()[0].grid_cell) = 1.0;

  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  fix.compute(kDt);
  const gutibm::Int carbon = sim.chemical_field().find(species::CARBON);
  const Real volume = sim.domain().cell_volume();
  const Real respiring_growth =
      sim.agents()[0].mu_realized * sim.agents()[0].biomass * kDt;
  const Real fermenting_growth =
      sim.agents()[1].mu_realized * sim.agents()[1].biomass * kDt;
  const Real respiring_draw = -sim.chemical_field().reac_global(
      carbon, sim.agents()[0].grid_cell) * volume * kDt;
  const Real fermenting_draw = -sim.chemical_field().reac_global(
      carbon, sim.agents()[1].grid_cell) * volume * kDt;
  assert(sim.agents()[0].realized_fermentation_fraction
         < sim.agents()[1].realized_fermentation_fraction);
  assert(std::abs(respiring_growth - fermenting_growth)
         < respiring_growth * 1.0e-6);
  assert(fermenting_draw / fermenting_growth
         > respiring_draw / respiring_growth);
  std::cout << "  test_simulation_carbon_cost_direction: PASSED\n";
}

void test_mu_crit_overflow_sensitivity() {
  SimulationConfig cfg = simulation_config();
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.oxygen.metabolic_switch_enabled = true;
  cfg.chem_env.oxygen.tau_metabolic_switch = 1.0;
  cfg.chem_env.oxygen.aerobic_mu_factor = 1.0;
  cfg.chem_env.oxygen.anaerobic_mu_factor = 1.0;
  cfg.chem_env.acetate.enabled = true;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.fixes.metabolism.division_threshold = 1.0e9;
  cfg.initial_strains[0].count = 2;
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  set_agent_cells(sim);
  set_species(sim, species::CARBON, 1.0);
  set_species(sim, species::OXYGEN, 1.0e30);
  set_species(sim, species::ACETATE, 0.0);
  sim.agents()[0].mu_max = 5.5e-4;
  sim.agents()[1].mu_max = 1.0e-5;

  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  fix.compute(kDt);

  const gutibm::Int acetate = sim.chemical_field().find(species::ACETATE);
  const Real volume = sim.domain().cell_volume();
  const Real above_acetate = sim.chemical_field().reac_global(
      acetate, sim.agents()[0].grid_cell) * volume * kDt;
  const Real below_acetate = sim.chemical_field().reac_global(
      acetate, sim.agents()[1].grid_cell) * volume * kDt;
  const Real mu_crit = sim.config().chem_env.oxygen.mu_crit;
  int above_threshold = 0;
  for (const Agent& agent : sim.agents()) {
    if (agent.mu_realized > mu_crit) ++above_threshold;
  }
  std::cout << "  mu_crit=" << mu_crit
            << " fraction_above=" << static_cast<Real>(above_threshold)
                / static_cast<Real>(sim.agents().size())
            << " mu_above=" << sim.agents()[0].mu_realized
            << " mu_below=" << sim.agents()[1].mu_realized
            << " ferm_above="
                << sim.agents()[0].realized_fermentation_fraction
            << " ferm_below="
                << sim.agents()[1].realized_fermentation_fraction
            << "\n";
  assert(sim.agents()[0].mu_realized > mu_crit);
  assert(sim.agents()[1].mu_realized < mu_crit);
  assert(sim.agents()[0].realized_fermentation_fraction > 0.1);
  assert(sim.agents()[1].realized_fermentation_fraction < 1.0e-12);
  assert(above_acetate > 0.0);
  assert(std::abs(below_acetate) < 1.0e-30);
  std::cout << "  test_mu_crit_overflow_sensitivity: PASSED\n";
}

Real run_acid_probe(Real acetate, Real ph) {
  SimulationConfig cfg = simulation_config();
  cfg.chem_env.oxygen.enabled = false;
  cfg.chem_env.acetate.enabled = true;
  cfg.fixes.metabolism.acid_inhibition_enabled = true;
  cfg.fixes.bacteriocin.mucin_charge.ph = ph;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  InputParser::finalize_config(cfg);
  Simulation sim;
  sim.init(cfg);
  set_species(sim, species::ACETATE, acetate);
  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  fix.compute(kDt);
  return sim.agents()[0].mu_realized;
}

void test_simulation_acid_inhibition() {
  const Real no_acid = run_acid_probe(0.0, 5.0);
  const Real low_acid = run_acid_probe(50.0, 5.0);
  const Real high_acid = run_acid_probe(100.0, 5.0);
  assert(no_acid > low_acid && low_acid > high_acid);
  assert(no_acid - high_acid > no_acid * 1.0e-3);
  assert(run_acid_probe(100.0, 5.0)
         < run_acid_probe(100.0, 6.0));
  std::cout << "  test_simulation_acid_inhibition: PASSED\n";
}

Real run_maintenance_probe(bool fermentative) {
  SimulationConfig cfg = simulation_config();
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.oxygen.metabolic_switch_enabled = true;
  cfg.chem_env.oxygen.tau_metabolic_switch = 1.0;
  cfg.chem_env.oxygen.anaerobic_maintenance_factor = 15.0;
  cfg.fixes.metabolism.carbon_maintenance_rate = 1.0e-5;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  InputParser::finalize_config(cfg);
  Simulation sim;
  sim.init(cfg);
  set_species(sim, species::CARBON, 1.0);
  set_species(sim, species::OXYGEN, fermentative ? 0.0 : 1.0);
  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  fix.compute(kDt);
  const gutibm::Int carbon = sim.chemical_field().find(species::CARBON);
  return sim.chemical_field().flux_accounting()
      .maintenance_step[static_cast<size_t>(carbon)];
}

void test_simulation_maintenance_coupling() {
  const Real aerobic = run_maintenance_probe(false);
  const Real anaerobic = run_maintenance_probe(true);
  assert(anaerobic > aerobic * 10.0);
  std::cout << "  test_simulation_maintenance_coupling: PASSED\n";
}

void test_simulation_inertia() {
  SimulationConfig cfg = simulation_config();
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.oxygen.metabolic_switch_enabled = true;
  cfg.chem_env.oxygen.tau_metabolic_switch = kDt;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  InputParser::finalize_config(cfg);
  Simulation sim;
  sim.init(cfg);
  set_species(sim, species::OXYGEN, 0.0);
  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  fix.compute(kDt);
  const Real expected_first = 1.0 - std::exp(-1.0);
  assert(std::abs(sim.agents()[0].realized_fermentation_fraction
                  - expected_first) < 1.0e-12);
  set_species(sim, species::OXYGEN, 1.0e30);
  Real previous = sim.agents()[0].realized_fermentation_fraction;
  for (int step = 0; step < 3; ++step) {
    sim.chemical_field().zero_reactions();
    fix.compute(kDt);
    const Real expected = previous * std::exp(-1.0);
    assert(std::abs(sim.agents()[0].realized_fermentation_fraction
                    - expected) < 1.0e-12);
    previous = expected;
  }
  std::cout << "  test_simulation_inertia: PASSED\n";
}

void test_simulation_daughter_inheritance() {
  SimulationConfig cfg = simulation_config();
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.oxygen.metabolic_switch_enabled = true;
  cfg.chem_env.oxygen.tau_metabolic_switch = 1.0;
  cfg.fixes.metabolism.division_threshold = 0.5;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  InputParser::finalize_config(cfg);
  Simulation sim;
  sim.init(cfg);
  sim.agents()[0].realized_fermentation_fraction = 0.73;
  set_species(sim, species::OXYGEN, 0.0);
  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  fix.compute(kDt);
  assert(sim.agents().size() == 2);
  assert(sim.agents()[0].realized_fermentation_fraction
         == sim.agents()[1].realized_fermentation_fraction);
  assert(sim.agents()[1].realized_fermentation_fraction > 0.73);
  std::cout << "  test_simulation_daughter_inheritance: PASSED\n";
}

void test_simulation_bounds_and_acid_stoichiometry() {
  SimulationConfig cfg = simulation_config();
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.acetate.enabled = true;
  cfg.chem_env.oxygen.metabolic_switch_enabled = true;
  cfg.chem_env.oxygen.tau_metabolic_switch = 1.0;
  cfg.chem_env.oxygen.ferm_acid_yield = 1.0;
  cfg.fixes.metabolism.maintenance_rate = 0.0;
  cfg.fixes.metabolism.division_threshold = 1.0e9;
  InputParser::finalize_config(cfg);
  Simulation sim;
  sim.init(cfg);
  set_species(sim, species::CARBON, 1.0);
  set_species(sim, species::OXYGEN, 0.0);
  set_species(sim, species::ACETATE, 0.0);
  FixMetabolism fix(sim, sim.config().fixes.metabolism);
  const gutibm::Int carbon = sim.chemical_field().find(species::CARBON);
  const gutibm::Int acetate = sim.chemical_field().find(species::ACETATE);
  const Real volume = sim.domain().cell_volume();
  for (int step = 0; step < 5; ++step) {
    sim.chemical_field().zero_reactions();
    fix.compute(kDt);
    const Real carbon_consumed = -sim.chemical_field().reac_global(
        carbon, sim.agents()[0].grid_cell) * volume * kDt;
    const Real acid_produced = sim.chemical_field().reac_global(
        acetate, sim.agents()[0].grid_cell) * volume * kDt;
    const Real fraction =
        sim.agents()[0].realized_fermentation_fraction;
    assert(std::isfinite(fraction));
    assert(fraction >= 0.0 && fraction <= 1.0);
    assert(std::isfinite(sim.agents()[0].mu_realized));
    assert(acid_produced >= 0.0);
    assert(acid_produced <= carbon_consumed * (1.0 + 1.0e-12));
  }
  std::cout << "  test_simulation_bounds_and_acid_stoichiometry: PASSED\n";
}

#ifdef GUTIBM_HDF5
bool dataset_exists(hid_t file, const char* path) {
  htri_t exists = 0;
  H5E_BEGIN_TRY {
    exists = H5Lexists(file, path, H5P_DEFAULT);
  }
  H5E_END_TRY;
  return exists > 0;
}

void test_summary_and_checkpoint_round_trip() {
  const std::string filename = gutibm::resolve_test_h5_path(
      "GUTIBM_METABOLIC_MODE_H5", "metabolic_mode");
  SimulationConfig cfg = simulation_config();
  cfg.chem_env.oxygen.enabled = true;
  cfg.chem_env.oxygen.metabolic_switch_enabled = true;
  cfg.hdf5.enabled = true;
  cfg.hdf5.filename = filename;
  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 1;
  cfg.hdf5.schedule.lineage = 1;
  cfg.hdf5.schedule.genome = 1;
  cfg.hdf5.schedule.grid = 1;
  cfg.hdf5.schedule.grid_species = {"all"};
  cfg.time.total_time = kDt;
  cfg.initial_strains[0].count = 2;
  InputParser::finalize_config(cfg);
  Simulation sim;
  sim.init(cfg);
  set_species(sim, species::OXYGEN, 0.0);
  sim.run();

  const auto snapshot = gutibm::HDF5Reader::load_snapshot(
      filename, "step_000001");
  assert(snapshot.agents.realized_fermentation_fraction.size()
         == snapshot.agents.id.size());
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(file >= 0);
  assert(dataset_exists(
      file, "summary/step_000001/mean_realized_fermentation_fraction"));
  H5Fclose(file);

  cfg.hdf5.enabled = false;
  Simulation resumed;
  resumed.init_from_checkpoint(cfg, filename, "step_000001");
  assert(static_cast<size_t>(resumed.agents().size())
         == snapshot.agents.id.size());
  assert(std::abs(resumed.agents()[0].realized_fermentation_fraction
                  - snapshot.agents.realized_fermentation_fraction[0])
         < 1.0e-15);
  std::cout << "  test_summary_and_checkpoint_round_trip: PASSED\n";
}
#endif

}  // namespace

int main() {
  using gutibm::metabolic_mode::acid_inhibition;
  using gutibm::metabolic_mode::fermentation_fraction;
  using gutibm::metabolic_mode::interpolate;
  using gutibm::metabolic_mode::relax;

  std::cout << "=== Metabolic Mode Tests ===\n";

  const double low_oxygen = fermentation_fraction(0.0, 1.0e-4, 3.0e-4);
  const double medium_oxygen = fermentation_fraction(0.5, 1.0e-4, 3.0e-4);
  const double high_oxygen = fermentation_fraction(1.0, 1.0e-4, 3.0e-4);
  assert(low_oxygen >= 0.0 && low_oxygen <= 1.0);
  assert(low_oxygen > medium_oxygen && medium_oxygen > high_oxygen);

  const double overflow = fermentation_fraction(1.0, 6.0e-4, 3.0e-4);
  assert(overflow > high_oxygen);
  assert(interpolate(1.0, 4.1, overflow)
         > interpolate(1.0, 4.1, high_oxygen));
  const double respiring_cost = interpolate(1.0, 4.1, high_oxygen);
  const double fermenting_cost = interpolate(1.0, 4.1, low_oxygen);
  assert(fermenting_cost > respiring_cost);
  assert(fermenting_cost <= 4.1);

  const double relaxed = relax(0.0, 1.0, 3600.0, 3600.0);
  assert(std::abs(relaxed - (1.0 - std::exp(-1.0))) < 1.0e-12);
  const double acid_low = acid_inhibition(20.0, 5.0, 4.76, 50.0, 0.8);
  const double acid_mid = acid_inhibition(80.0, 5.0, 4.76, 50.0, 0.8);
  const double acid_high = acid_inhibition(160.0, 5.0, 4.76, 50.0, 0.8);
  assert(acid_low < acid_mid && acid_mid < acid_high);
  assert(acid_high <= 0.8);
  assert(acid_inhibition(80.0, 5.0, 4.76, 50.0, 0.8)
         > acid_inhibition(80.0, 6.0, 4.76, 50.0, 0.8));
  assert(acid_inhibition(0.0, 6.0, 4.76, 50.0, 0.8) == 0.0);
  std::cout << "  helper_algebra: PASSED\n";

  test_simulation_backward_compatibility();
  test_funded_respiration_graded_and_bounded();
  test_funded_respiration_discriminates_from_ambient();
  test_funded_respiration_zero_demand_preserves_state();
  test_simulation_carbon_cost_direction();
  test_mu_crit_overflow_sensitivity();
  test_simulation_acid_inhibition();
  test_simulation_maintenance_coupling();
  test_simulation_inertia();
  test_simulation_daughter_inheritance();
  test_simulation_bounds_and_acid_stoichiometry();
#ifdef GUTIBM_HDF5
  test_summary_and_checkpoint_round_trip();
#endif

  std::cout << "All metabolic mode tests passed.\n";
  return 0;
}
