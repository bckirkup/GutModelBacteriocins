/* -----------------------------------------------------------------------
   GutIBM – OpenMP parity test
   Verifies that the OpenMP-enabled build produces results consistent
   with the serial build (within floating-point tolerance).
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "sim_fingerprint.h"
#include <array>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>

#ifdef GUTIBM_OPENMP
#include <omp.h>
#endif

using namespace gutibm;

static SimulationConfig make_test_config(int seed) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.time.total_time      = 300.0;
  cfg.time.bio_dt          = 60.0;
  cfg.time.output_interval = 300.0;
  cfg.seed            = seed;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length   = 100e-6;
  cfg.advection.radial_turnover = 5400.0;
  cfg.advection.distal_transit_time = 43200.0;
  cfg.qssa.toxin_cutoff = 50e-6;
  cfg.qssa.nutrient_cutoff = 25e-6;
  cfg.hdf5.enabled = false;
  cfg.cell_bio.fur.enabled = false;
  cfg.cell_bio.cdi.enabled = false;
  cfg.cell_bio.motility.enabled = false;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s1;
  s1.type = 1; s1.count = 20; s1.mu_max = 5e-4;
  s1.plasmids = {"ColE1"}; s1.conjugative = false;
  cfg.initial_strains.push_back(s1);

  SimulationConfig::InitialStrain s2;
  s2.type = 2; s2.count = 10; s2.mu_max = 5e-4;
  s2.plasmids = {}; s2.conjugative = false;
  cfg.initial_strains.push_back(s2);

  return cfg;
}

void test_openmp_compile_flag() {
#ifdef GUTIBM_OPENMP
  int threads = omp_get_max_threads();
  std::cout << "  test_openmp_compile_flag: PASSED (OpenMP enabled, max_threads="
            << threads << ")\n";
#else
  std::cout << "  test_openmp_compile_flag: PASSED (OpenMP disabled, serial mode)\n";
#endif
}

void test_simulation_completes() {
  SimulationConfig cfg = make_test_config(42);
  Simulation sim;
  sim.init(cfg);

  assert(sim.agents().size() == 30);
  sim.run();

  assert(sim.time() > 0.0);
  assert(sim.step_count() > 0);

  Int alive = 0;
  for (const Agent& a : sim.agents()) {
    if (a.state != PhenoState::DEAD) alive++;
  }
  assert(alive > 0);

  std::cout << "  test_simulation_completes: PASSED (alive=" << alive
            << " steps=" << sim.step_count() << ")\n";
}

void test_deterministic_growth() {
  // Run twice with same seed, verify biomass totals match
  std::array<Real, 2> total_biomass = {0.0, 0.0};

  for (int trial : {0, 1}) {
    SimulationConfig cfg = make_test_config(999);

    // Use only non-toxic agents to avoid stochastic kill events
    cfg.initial_strains.clear();
    SimulationConfig::InitialStrain s;
    s.type = 1; s.count = 10; s.mu_max = 5e-4;
    s.plasmids = {}; s.conjugative = false;
    cfg.initial_strains.push_back(s);

    Simulation sim;
    sim.init(cfg);
    sim.run();

    for (const Agent& a : sim.agents()) {
      if (a.state != PhenoState::DEAD)
        total_biomass[trial] += a.biomass;
    }
  }

  Real rel_diff = std::abs(total_biomass[0] - total_biomass[1])
                / std::max(total_biomass[0], 1e-30);

  // With OpenMP, floating-point reduction order may differ slightly
  // Allow 1e-10 relative tolerance
  assert(rel_diff < 1e-10);

  std::cout << "  test_deterministic_growth: PASSED (biomass0="
            << total_biomass[0] << " biomass1=" << total_biomass[1]
            << " rel_diff=" << rel_diff << ")\n";
}

void test_chemical_field_parity() {
  SimulationConfig cfg = make_test_config(123);
  Simulation sim;
  sim.init(cfg);

  // Run a few steps
  for (int step = 0; step < 3; ++step) {
    sim.step(cfg.time.bio_dt);
  }

  const auto& chem = sim.chemical_field();
  Int i_carbon = chem.find("carbon");
  Int i_iron = chem.find("iron");

  if (i_carbon >= 0) {
    Real sum = 0.0;
    for (Real val : chem.conc_data()[static_cast<size_t>(i_carbon)]) {
      assert(!std::isnan(val));
      assert(!std::isinf(val));
      assert(val >= 0.0);
      sum += val;
    }
    assert(sum > 0.0);
    std::cout << "  test_chemical_field_parity: carbon sum=" << sum << "\n";
  }

  if (i_iron >= 0) {
    for (Real val : chem.conc_data()[static_cast<size_t>(i_iron)]) {
      assert(!std::isnan(val));
      assert(!std::isinf(val));
      assert(val >= 0.0);
    }
  }

  std::cout << "  test_chemical_field_parity: PASSED\n";
}

void test_grid_coupling_consistency() {
  SimulationConfig cfg = make_test_config(777);
  Simulation sim;
  sim.init(cfg);
  sim.step(cfg.time.bio_dt);

  // Every alive agent should have a valid grid cell
  for (const Agent& a : sim.agents()) {
    if (a.state == PhenoState::DEAD) continue;
    assert(a.grid_cell >= 0);
    assert(a.grid_cell < sim.chemical_field().ncells());
  }

  std::cout << "  test_grid_coupling_consistency: PASSED\n";
}

void test_cross_build_fingerprint() {
  // Deterministic, non-stochastic scenario for serial vs OpenMP comparison.
  // Emits FINGERPRINT=<hex> for scripts/compare_openmp_parity.sh
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo  = {0, 0, 0};
  cfg.domain.hi  = {80e-6, 80e-6, 40e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.time.total_time      = 180.0;
  cfg.time.bio_dt          = 60.0;
  cfg.time.output_interval = 180.0;
  cfg.seed            = 31415;
  cfg.hdf5.enabled    = false;
  cfg.advection.mucus_thickness = 40e-6;
  cfg.advection.distal_length   = 80e-6;
  cfg.qssa.toxin_cutoff = 40e-6;
  cfg.qssa.nutrient_cutoff = 20e-6;
  cfg.cell_bio.fur.enabled = false;
  cfg.cell_bio.cdi.enabled = false;
  cfg.cell_bio.motility.enabled = false;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = 12;
  s.mu_max = 4e-4;
  s.plasmids = {};
  s.conjugative = false;
  cfg.initial_strains.push_back(s);

  Simulation sim;
  sim.init(cfg);
  sim.run();

  std::string fp = test_util::fingerprint_hex(sim);
  std::cout << "FINGERPRINT=" << fp << "\n";
  assert(!fp.empty());

  std::cout << "  test_cross_build_fingerprint: PASSED (fp=" << fp << ")\n";
}

int main() {
  std::cout << "=== OpenMP Parity Tests ===\n";
  test_openmp_compile_flag();
  test_simulation_completes();
  test_deterministic_growth();
  test_chemical_field_parity();
  test_grid_coupling_consistency();
  test_cross_build_fingerprint();
  std::cout << "All OpenMP parity tests passed.\n";
  return 0;
}
