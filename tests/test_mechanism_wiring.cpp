/* -----------------------------------------------------------------------
   GutIBM – Mechanism-wiring verification tests

   These tests answer the cross-cutting question "are the mechanisms from all
   the specs wired together in a way that makes logical sense?" rather than
   probing any single Fix in isolation. They enforce two invariants drawn from
   the ci-test-design skill:

     1. Species mass balance — every coupled species reaches a *bounded* steady
        state. A species that only has a source (accumulates without bound) or
        only has a sink (drains to zero) signals a missing coupling. This is
        exactly the failure mode Spec 5 §1 (carbon source without sink) and
        §3 (B12 sink without source) describe.

     2. Directional coupling sensitivity — turning a coupling on/up must move
        the downstream quantity in the biologically correct direction (more
        agents ⇒ less O2; larger sink ⇒ less substrate). A coupling that is
        parsed but dead, or wired backwards, is caught here.

   Failures are recorded explicitly (not via assert) so the test still fails
   under NDEBUG/Release builds.
   ----------------------------------------------------------------------- */

#include "simulation.h"
#include "input_parser.h"
#include "vbf.h"
#include "chemical_field.h"
#include "domain.h"
#include "fix_metabolism.h"
#include "chem_environment_config.h"
#include "species_names.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace gutibm;

namespace {

int g_failures = 0;

void expect(bool cond, const std::string& msg) {
  if (!cond) {
    std::cerr << "FAIL: " << msg << "\n";
    ++g_failures;
  }
}

// Integrate one closed-cell step of the VBF reaction terms: zero reactions,
// let the VBF apply its sources/sinks, then apply conc += reac * dt (clamped
// at zero, matching Simulation::module_chemistry). No diffusion/boundaries —
// this isolates the VBF coupling so the mass-balance conclusion is unambiguous.
void step_vbf_closed_cell(VBF& vbf, ChemicalField& chem, const Domain& domain,
                          const OxygenConfig& oxy, const AcetateConfig& ace,
                          const MucinConfig& muc, Real dt) {
  chem.zero_reactions();
  vbf.apply_nutrient_coupling(chem, domain, dt, oxy, ace, muc);
  const Int ns = chem.num_species();
  for (Int s = 0; s < ns; ++s) {
    for (Int c = 0; c < chem.ncells(); ++c) {
      chem.conc(s, c) = std::max(chem.conc(s, c) + chem.reac(s, c) * dt, 0.0);
    }
  }
}

Domain make_domain() {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {20e-6, 20e-6, 20e-6};
  dcfg.grid_dx = 5e-6;
  Domain domain;
  domain.init(dcfg);
  return domain;
}

Real max_conc(const ChemicalField& chem, Int spec) {
  Real m = 0.0;
  for (Int c = 0; c < chem.ncells(); ++c) {
    m = std::max(m, chem.conc(spec, c));
  }
  return m;
}

// ── Spec 5 §1 — carbon source without a sink accumulates without bound; the
// VBF carbon sink must convert that into a bounded steady state. ────────────
void test_carbon_sink_bounds_accumulation() {
  Domain domain = make_domain();

  ChemicalSpec carbon;
  carbon.name = species::CARBON;
  carbon.initial_conc = 0.0;

  OxygenConfig oxy;   // disabled
  AcetateConfig ace;  // disabled
  MucinConfig muc;    // disabled -> static liberation path

  // Uniform carbon liberation (no z-gradient), sink off vs on.
  VBFConfig base;
  base.mucin_liberation = 5.0e-5;   // constant carbon source
  base.mucin_z_gradient_enabled = false;
  base.use_dynamic_mucin = false;

  VBFConfig no_sink = base;
  no_sink.carbon_sink_vmax = 0.0;

  VBFConfig with_sink = base;
  with_sink.carbon_sink_vmax = 1.0e-3;  // >> source, so steady state is small
  with_sink.carbon_sink_km   = 1.0e-3;

  const Real dt = 60.0;
  const int steps = 400;

  ChemicalField chem_a;
  chem_a.init(domain, {carbon});
  VBF vbf_a;
  vbf_a.init(no_sink, domain);

  ChemicalField chem_b;
  chem_b.init(domain, {carbon});
  VBF vbf_b;
  vbf_b.init(with_sink, domain);

  for (int i = 0; i < steps; ++i) {
    step_vbf_closed_cell(vbf_a, chem_a, domain, oxy, ace, muc, dt);
    step_vbf_closed_cell(vbf_b, chem_b, domain, oxy, ace, muc, dt);
  }

  const Int ic = chem_a.find(species::CARBON);
  const Real c_no_sink = max_conc(chem_a, ic);
  const Real c_with_sink = max_conc(chem_b, ic);

  // Analytic unbounded growth: source * dt * steps.
  const Real expected_unbounded = 5.0e-5 * dt * steps;  // = 1.2 mol/m^3
  // Analytic bounded steady state: km * S / (vmax - S).
  const Real steady = 1.0e-3 * 5.0e-5 / (1.0e-3 - 5.0e-5);

  expect(c_no_sink > 0.5 * expected_unbounded,
         "carbon without a sink should accumulate unbounded (Spec 5 §1 gap)");
  expect(c_with_sink < 100.0 * steady && c_with_sink < 1.0e-3,
         "carbon with the VBF sink must reach a small bounded steady state");
  expect(c_with_sink < 0.01 * c_no_sink,
         "enabling the carbon sink must dramatically lower carbon vs no sink");

  std::cout << "  test_carbon_sink_bounds_accumulation: PASSED"
            << " (no_sink=" << c_no_sink << " with_sink=" << c_with_sink
            << " steady~=" << steady << ")\n";
}

// Directional sensitivity: a larger Vmax removes more carbon per step. ──────
void test_carbon_sink_sensitivity() {
  Domain domain = make_domain();

  ChemicalSpec carbon;
  carbon.name = species::CARBON;
  carbon.initial_conc = 1.0e-2;  // fixed, well above Km

  OxygenConfig oxy;
  AcetateConfig ace;
  MucinConfig muc;

  auto sink_reac = [&](Real vmax) {
    VBFConfig cfg;
    cfg.mucin_liberation = 0.0;  // isolate the sink term
    cfg.mucin_z_gradient_enabled = false;
    cfg.carbon_sink_vmax = vmax;
    cfg.carbon_sink_km = 1.0e-3;
    ChemicalField chem;
    chem.init(domain, {carbon});
    VBF vbf;
    vbf.init(cfg, domain);
    chem.zero_reactions();
    vbf.apply_nutrient_coupling(chem, domain, 60.0, oxy, ace, muc);
    return chem.reac(chem.find(species::CARBON), domain.cell_index(0, 0, 0));
  };

  const Real r_off = sink_reac(0.0);
  const Real r_lo = sink_reac(1.0e-4);
  const Real r_hi = sink_reac(1.0e-3);

  expect(r_off == 0.0, "carbon sink disabled (vmax=0) must not touch carbon");
  expect(r_lo < 0.0 && r_hi < 0.0, "carbon sink must be a negative reaction term");
  expect(r_hi < r_lo, "larger carbon sink Vmax must remove more carbon (monotonic)");

  std::cout << "  test_carbon_sink_sensitivity: PASSED"
            << " (off=" << r_off << " lo=" << r_lo << " hi=" << r_hi << ")\n";
}

// Build a small integration config with a chosen number of agents. ──────────
SimulationConfig make_integration_cfg(int n_agents, uint64_t seed) {
  SimulationConfig cfg = InputParser::default_config();
  cfg.domain.lo = {0, 0, 0};
  cfg.domain.hi = {100e-6, 100e-6, 50e-6};
  cfg.domain.grid_dx = 5e-6;
  cfg.domain.hash_cell_size = 10e-6;
  cfg.time.total_time = 600.0;
  cfg.time.bio_dt = 60.0;
  cfg.time.output_interval = 600.0;
  cfg.seed = seed;
  cfg.hdf5.enabled = false;
  cfg.cell_bio.motility.enabled = false;
  cfg.cell_bio.cdi.enabled = false;
  cfg.advection.mucus_thickness = 50e-6;
  cfg.advection.distal_length = 100e-6;
  cfg.qssa.toxin_cutoff = 50e-6;
  cfg.qssa.nutrient_cutoff = 25e-6;

  cfg.initial_strains.clear();
  SimulationConfig::InitialStrain s;
  s.type = 1;
  s.count = n_agents;
  s.mu_max = 5.0e-4;
  s.plasmids = {};
  s.conjugative = false;
  cfg.initial_strains.push_back(s);
  return cfg;
}

Real mean_conc(const ChemicalField& chem, Int spec) {
  Real sum = 0.0;
  for (Int c = 0; c < chem.ncells(); ++c) sum += chem.conc(spec, c);
  return chem.ncells() > 0 ? sum / static_cast<Real>(chem.ncells()) : 0.0;
}

// ── O2 consumption is wired end-to-end (Spec 1): more agents must draw the
// dissolved-oxygen field down further. Proves agent metabolism actually
// couples back into the chemical environment through the full step loop. ────
Real run_o2_and_report_mean(int n_agents, uint64_t seed, bool keep_background) {
  SimulationConfig cfg = make_integration_cfg(n_agents, seed);
  cfg.chem_env.oxygen.enabled = true;
  if (!keep_background) {
    // Isolate the agent<->O2 coupling: the uniform VBF background O2 sink would
    // otherwise dominate the field and mask the agent-count signal.
    cfg.chem_env.oxygen.vbf_sink = 0.0;
  }
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  for (int i = 0; i < 6; ++i) sim.step(60.0);

  const ChemicalField& chem = sim.chemical_field();
  const Int io2 = chem.find(species::OXYGEN);
  return io2 >= 0 ? mean_conc(chem, io2) : -1.0;
}

void test_o2_consumption_wired() {
  const Real epithelial = OxygenConfig{}.epithelial_conc;
  const Real o2_few = run_o2_and_report_mean(5, 111, /*keep_background=*/false);
  const Real o2_many = run_o2_and_report_mean(80, 111, /*keep_background=*/false);

  expect(o2_few >= 0.0 && o2_many >= 0.0, "oxygen species must be registered when enabled");
  expect(o2_few < epithelial,
         "agents must consume O2 (mean field below the epithelial boundary)");
  expect(o2_many < o2_few,
         "more agents must deplete O2 further (directional agent<->chemistry coupling)");

  std::cout << "  test_o2_consumption_wired: PASSED"
            << " (few=" << o2_few << " many=" << o2_many
            << " epithelial=" << epithelial << ")\n";
}

// ── Regression for Edison's observation: with the VBF background O2 sink ON
// (the realistic config), the field must still track agent *density* — not
// just settle at the background-sink equilibrium. This failed for two reasons
// that this test now guards against together:
//   1. the background sink was zero-order, removing more O2 than present and
//      hard-zeroing the interior in one step (masking any agent signal);
//   2. respiration was purely growth-coupled, so non-growing cells drew no O2.
// With the first-order sink + Pirt maintenance term, denser agent populations
// must pull the mean O2 progressively below the agent-free baseline. ─────────
void test_o2_tracks_density_over_background() {
  const Real o2_none = run_o2_and_report_mean(0, 222, /*keep_background=*/true);
  const Real o2_few = run_o2_and_report_mean(10, 222, /*keep_background=*/true);
  const Real o2_many = run_o2_and_report_mean(80, 222, /*keep_background=*/true);

  expect(o2_none > 0.0, "background-only O2 equilibrium must be positive");
  expect(o2_few < o2_none,
         "agents must draw O2 below the background-sink equilibrium");
  expect(o2_many < o2_few,
         "denser agent population must deplete O2 further, even with the "
         "VBF background sink active (density tracking, not just equilibrium)");

  std::cout << "  test_o2_tracks_density_over_background: PASSED"
            << " (none=" << o2_none << " few=" << o2_few
            << " many=" << o2_many << ")\n";
}

// ── Spec 6 §3 — corrinoid (B12) is a constant bioavailable pool, NOT a
// depletable field. The anaerobic majority produces total corrinoids (~1 uM)
// far faster than E. coli consumes them, so neither the VBF nor the agents may
// add or remove B12: the field must stay pinned at its 1 uM initial/boundary
// value across a full multi-step run. (Replaces the Spec 5 §3 VBF-B12-source
// test, which is obsolete now that B12 is neither produced nor consumed.) ────
void test_corrinoid_field_constant() {
  SimulationConfig cfg = make_integration_cfg(60, 4242);
  cfg.chem_env.oxygen.enabled = true;  // exercise the full nutrient path
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);

  const ChemicalField& chem = sim.chemical_field();
  const Int ib12 = chem.find(species::B12);
  expect(ib12 >= 0, "B12/corrinoid species must be registered");
  if (ib12 < 0) return;

  const Real expected = 1.0e-6;  // 1 uM total corrinoid (Spec 6 §3)
  const Real tol = 1.0e-12;      // constant field: no drift permitted

  Real max_dev = 0.0;
  for (int i = 0; i < 10; ++i) {
    sim.step(60.0);
    for (Int c = 0; c < chem.ncells(); ++c) {
      max_dev = std::max(max_dev, std::fabs(chem.conc(ib12, c) - expected));
    }
  }

  expect(max_dev < tol,
         "corrinoid field must stay pinned at 1 uM (no production, no depletion)");

  std::cout << "  test_corrinoid_field_constant: PASSED"
            << " (expected=" << expected << " max_dev=" << max_dev << ")\n";
}

// ── Spec 5 §4 — the dysbiosis safety net must halt a run once density leaves
// the calibrated regime, and must be inert when disabled. ───────────────────
void test_dysbiosis_halt() {
  // Control: threshold disabled -> run proceeds through many steps.
  SimulationConfig ctrl = make_integration_cfg(50, 2024);
  ctrl.dysbiosis_threshold = 0.0;
  Simulation sim_ctrl;
  sim_ctrl.init(ctrl);
  sim_ctrl.run();

  // Halting: threshold well below the initial density -> run stops early.
  SimulationConfig halt = make_integration_cfg(50, 2024);
  // Initial density = 50 cells / (100e-6*100e-6*50e-6 m^3 * 1e9 mL/m^3).
  halt.dysbiosis_threshold = 1.0e4;  // cells/mL, far below actual
  Simulation sim_halt;
  sim_halt.init(halt);
  sim_halt.run();

  expect(sim_ctrl.step_count() > sim_halt.step_count(),
         "dysbiosis threshold must halt the run earlier than the disabled control");
  expect(sim_halt.step_count() >= 1 && sim_halt.step_count() <= 2,
         "dysbiosis halt should trigger within the first step or two");

  std::cout << "  test_dysbiosis_halt: PASSED"
            << " (control_steps=" << sim_ctrl.step_count()
            << " halt_steps=" << sim_halt.step_count() << ")\n";
}

void test_metabolism_uptake_has_rate_units() {
  SimulationConfig cfg = make_integration_cfg(1, 707);
  cfg.time.total_time = 60.0;

  Simulation sim;
  sim.init(cfg);
  ChemicalField& chem = sim.chemical_field();
  chem.zero_reactions();

  const Real dt = 60.0;
  const Real biomass_before = sim.agents()[0].biomass;
  FixMetabolism metabolism(sim, cfg.fixes.metabolism);
  metabolism.compute(dt);

  const Agent& agent = sim.agents()[0];
  const Real biomass_gain = agent.biomass - biomass_before;
  const Real cell_volume = sim.domain().dx() * sim.domain().dx()
      * sim.domain().dx();
  const Int carbon = chem.find(species::CARBON);
  const Real expected_rate = biomass_gain * cfg.fixes.metabolism.yield_carbon
      / (cell_volume * dt);
  const Real actual_rate = -chem.reac(carbon, agent.grid_cell);
  const Real tolerance = std::max(1.0e-15, std::abs(expected_rate) * 1.0e-12);

  expect(std::abs(actual_rate - expected_rate) <= tolerance,
         "metabolism uptake must be a concentration rate, not a per-step amount");
  std::cout << "  test_metabolism_uptake_has_rate_units: PASSED"
            << " (rate=" << actual_rate << ")\n";
}

// ── Overarching mass-balance guard: with oxygen and the carbon sink active
// (and B12 a constant pool), no coupled species may go NaN, negative, or blow
// up over a multi-step run. This is the "everything wired together stays sane"
// net. ──────────────────────────────────────────────────────────────────────
void test_all_species_bounded_steady_state() {
  SimulationConfig cfg = make_integration_cfg(40, 909);
  cfg.chem_env.oxygen.enabled = true;
  cfg.vbf.carbon_sink_vmax = 1.0e-3;
  cfg.vbf.carbon_sink_km = 1.0e-3;
  InputParser::finalize_config(cfg);

  Simulation sim;
  sim.init(cfg);
  for (int i = 0; i < 8; ++i) sim.step(60.0);

  const ChemicalField& chem = sim.chemical_field();
  const std::vector<std::string> tracked = {
      species::CARBON, species::IRON, species::B12, species::OXYGEN};

  for (const std::string& name : tracked) {
    const Int idx = chem.find(name);
    if (idx < 0) continue;
    for (Int c = 0; c < chem.ncells(); ++c) {
      const Real v = chem.conc(idx, c);
      if (std::isnan(v) || std::isinf(v) || v < 0.0 || v > 1.0e3) {
        std::ostringstream m;
        m << "species '" << name << "' left the bounded regime: conc=" << v
          << " at cell " << c;
        expect(false, m.str());
        break;
      }
    }
  }

  // Carbon must stay bounded despite its constant VBF source (sink active).
  const Int ic = chem.find(species::CARBON);
  if (ic >= 0) {
    expect(max_conc(chem, ic) < 1.0,
           "with the sink active carbon must stay bounded across a full run");
  }

  std::cout << "  test_all_species_bounded_steady_state: PASSED\n";
}

}  // namespace

int main() {
  std::cout << "=== Mechanism Wiring Tests ===\n";
  test_carbon_sink_bounds_accumulation();
  test_carbon_sink_sensitivity();
  test_o2_consumption_wired();
  test_o2_tracks_density_over_background();
  test_corrinoid_field_constant();
  test_dysbiosis_halt();
  test_metabolism_uptake_has_rate_units();
  test_all_species_bounded_steady_state();

  if (g_failures == 0) {
    std::cout << "All mechanism wiring tests passed.\n";
    return 0;
  }
  std::cerr << g_failures << " mechanism wiring test(s) FAILED.\n";
  return 1;
}
