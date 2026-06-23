/* -----------------------------------------------------------------------
   GutIBM – Tests for peristaltic mixing (issue #11)
   Verifies oscillatory advection perturbation behavior.
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "domain.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;

static constexpr Real TOL = 1.0e-12;

// Helper: create a minimal domain for testing
static Domain make_test_domain() {
  DomainConfig dcfg;
  dcfg.lo = {0, 0, 0};
  dcfg.hi = {1e-3, 1e-3, 100e-6};
  dcfg.grid_dx = 10e-6;
  dcfg.hash_cell_size = 10e-6;
  dcfg.periodic = {true, true, false};
  Domain d;
  d.init(dcfg);
  return d;
}

void test_peristaltic_disabled_gives_unity() {
  AdvectionConfig cfg;
  cfg.peristaltic_enabled = false;

  Domain dom = make_test_domain();
  AdvectionField field;
  field.init(cfg, dom);

  Vec3 pos = {0.5e-3, 0.5e-3, 50e-6};
  field.set_time(0.0);
  assert(std::abs(field.peristaltic_factor(pos) - 1.0) < TOL);

  field.set_time(10.0);
  assert(std::abs(field.peristaltic_factor(pos) - 1.0) < TOL);

  field.set_time(100.0);
  assert(std::abs(field.peristaltic_factor(pos) - 1.0) < TOL);

  std::cout << "  test_peristaltic_disabled_gives_unity: PASSED\n";
}

void test_peristaltic_disabled_velocity_unchanged() {
  AdvectionConfig cfg;
  cfg.peristaltic_enabled = false;

  Domain dom = make_test_domain();
  AdvectionField field;
  field.init(cfg, dom);

  Vec3 pos = {0.5e-3, 0.5e-3, 75e-6};
  field.set_time(5.0);

  Vec3 v = field.velocity(pos);
  // Should equal base velocities (no modulation)
  Real expected_vx = field.distal_velocity(pos[2]);
  Real expected_vz = field.radial_velocity(pos[2]);

  assert(std::abs(v[0] - expected_vx) < TOL);
  assert(std::abs(v[1]) < TOL);
  assert(std::abs(v[2] - expected_vz) < TOL);

  std::cout << "  test_peristaltic_disabled_velocity_unchanged: PASSED\n";
}

void test_peristaltic_oscillation_uniform() {
  AdvectionConfig cfg;
  cfg.peristaltic_enabled = true;
  cfg.peristaltic_period = 20.0;
  cfg.peristaltic_amplitude = 0.5;
  cfg.peristaltic_wavelength = 0.0;  // uniform (no spatial phase)

  Domain dom = make_test_domain();
  AdvectionField field;
  field.init(cfg, dom);

  Vec3 pos = {0.5e-3, 0.5e-3, 50e-6};

  // At t=0: sin(0) = 0 → factor = 1.0
  field.set_time(0.0);
  assert(std::abs(field.peristaltic_factor(pos) - 1.0) < TOL);

  // At t = T/4 = 5s: sin(pi/2) = 1 → factor = 1.5
  field.set_time(5.0);
  assert(std::abs(field.peristaltic_factor(pos) - 1.5) < TOL);

  // At t = T/2 = 10s: sin(pi) = 0 → factor = 1.0
  field.set_time(10.0);
  assert(std::abs(field.peristaltic_factor(pos) - 1.0) < TOL);

  // At t = 3T/4 = 15s: sin(3pi/2) = -1 → factor = 0.5
  field.set_time(15.0);
  assert(std::abs(field.peristaltic_factor(pos) - 0.5) < TOL);

  // At t = T = 20s: sin(2pi) = 0 → factor = 1.0
  field.set_time(20.0);
  assert(std::abs(field.peristaltic_factor(pos) - 1.0) < TOL);

  std::cout << "  test_peristaltic_oscillation_uniform: PASSED\n";
}

void test_peristaltic_modulates_velocity() {
  AdvectionConfig cfg;
  cfg.peristaltic_enabled = true;
  cfg.peristaltic_period = 20.0;
  cfg.peristaltic_amplitude = 0.5;
  cfg.peristaltic_wavelength = 0.0;

  Domain dom = make_test_domain();
  AdvectionField field;
  field.init(cfg, dom);

  Vec3 pos = {0.5e-3, 0.5e-3, 75e-6};

  Real base_vx = field.distal_velocity(pos[2]);
  Real base_vz = field.radial_velocity(pos[2]);

  // At peak (t=T/4): velocity should be 1.5x base
  field.set_time(5.0);
  Vec3 v_peak = field.velocity(pos);
  assert(std::abs(v_peak[0] - 1.5 * base_vx) < TOL);
  assert(std::abs(v_peak[2] - 1.5 * base_vz) < TOL);

  // At trough (t=3T/4): velocity should be 0.5x base
  field.set_time(15.0);
  Vec3 v_trough = field.velocity(pos);
  assert(std::abs(v_trough[0] - 0.5 * base_vx) < TOL);
  assert(std::abs(v_trough[2] - 0.5 * base_vz) < TOL);

  std::cout << "  test_peristaltic_modulates_velocity: PASSED\n";
}

void test_peristaltic_propagating_wave() {
  AdvectionConfig cfg;
  cfg.peristaltic_enabled = true;
  cfg.peristaltic_period = 20.0;
  cfg.peristaltic_amplitude = 0.5;
  cfg.peristaltic_wavelength = 1.0e-3;  // wavelength = domain length

  Domain dom = make_test_domain();
  AdvectionField field;
  field.init(cfg, dom);

  field.set_time(0.0);

  // At x=0: phase = 0 → sin(0) = 0 → factor = 1.0
  Vec3 pos_origin = {0.0, 0.5e-3, 50e-6};
  assert(std::abs(field.peristaltic_factor(pos_origin) - 1.0) < TOL);

  // At x = wavelength/4 = 0.25e-3: phase = -pi/2 → sin(-pi/2) = -1 → factor = 0.5
  Vec3 pos_quarter = {0.25e-3, 0.5e-3, 50e-6};
  assert(std::abs(field.peristaltic_factor(pos_quarter) - 0.5) < TOL);

  // At x = wavelength/2 = 0.5e-3: phase = -pi → sin(-pi) = 0 → factor = 1.0
  Vec3 pos_half = {0.5e-3, 0.5e-3, 50e-6};
  assert(std::abs(field.peristaltic_factor(pos_half) - 1.0) < TOL);

  // At x = 3*wavelength/4 = 0.75e-3: phase = -3pi/2 → sin(-3pi/2) = 1 → factor = 1.5
  Vec3 pos_three_quarter = {0.75e-3, 0.5e-3, 50e-6};
  assert(std::abs(field.peristaltic_factor(pos_three_quarter) - 1.5) < TOL);

  std::cout << "  test_peristaltic_propagating_wave: PASSED\n";
}

void test_peristaltic_uniform_same_at_different_positions() {
  // With wavelength=0, all positions should get the same factor
  AdvectionConfig cfg;
  cfg.peristaltic_enabled = true;
  cfg.peristaltic_period = 20.0;
  cfg.peristaltic_amplitude = 0.5;
  cfg.peristaltic_wavelength = 0.0;

  Domain dom = make_test_domain();
  AdvectionField field;
  field.init(cfg, dom);

  field.set_time(7.3);

  Vec3 pos1 = {0.0, 0.0, 50e-6};
  Vec3 pos2 = {0.5e-3, 0.5e-3, 50e-6};
  Vec3 pos3 = {1.0e-3, 1.0e-3, 100e-6};

  Real f1 = field.peristaltic_factor(pos1);
  Real f2 = field.peristaltic_factor(pos2);
  Real f3 = field.peristaltic_factor(pos3);

  assert(std::abs(f1 - f2) < TOL);
  assert(std::abs(f1 - f3) < TOL);

  std::cout << "  test_peristaltic_uniform_same_at_different_positions: PASSED\n";
}

int main() {
  std::cout << "=== Advection Peristaltic Tests ===\n";
  test_peristaltic_disabled_gives_unity();
  test_peristaltic_disabled_velocity_unchanged();
  test_peristaltic_oscillation_uniform();
  test_peristaltic_modulates_velocity();
  test_peristaltic_propagating_wave();
  test_peristaltic_uniform_same_at_different_positions();
  std::cout << "All peristaltic tests passed.\n";
  return 0;
}
