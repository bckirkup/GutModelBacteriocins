/* -----------------------------------------------------------------------
   GutIBM – Independent two-wall Neumann image-series tests
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "domain.h"
#include "greens_function.h"
#include "input_parser.h"
#include "neumann_image_series.h"
#include "error.h"
#include "species_names.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>

using namespace gutibm;

namespace {

struct TestSystem {
  Domain domain;
  AdvectionField adv;
  GreensFunction gf;
};

TestSystem make_system(Real grid_dx = 5.0e-6) {
  TestSystem system;
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {1.0e-3, 1.0e-3, 100.0e-6};
  dcfg.grid_dx = grid_dx;
  system.domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 1.0e20;
  acfg.distal_transit_time = 1.0e20;
  acfg.mucus_thickness = 100.0e-6;
  acfg.distal_length = 1.0e-3;
  acfg.taylor_aris_enabled = false;
  system.adv.init(acfg, system.domain);
  system.gf.init(system.domain, system.adv);
  return system;
}

GreensFunctionParams test_params(Real decay_rate) {
  GreensFunctionParams params;
  params.diff_coeff = 4.0e-11;
  params.retardation = 1.0;
  params.source_rate = 1.0e-18;
  params.decay_rate = decay_rate;
  return params;
}

void test_boundary_condition() {
  auto system = make_system();
  const Real height = system.domain.size()[2];
  const Real decay = 1.0e-4;
  const auto params = test_params(decay);
  const std::array<Real, 7> source_fractions =
      {0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95};
  const std::array<Real, 2> lateral_offsets = {20.0e-6, 50.0e-6};
  const Real epsilon = 1.0e-9;
  const Real difference_step = 1.0e-9;
  Real maximum_residual = 0.0;

  for (const Real fraction : source_fractions) {
    const Vec3 source = {500.0e-6, 500.0e-6, fraction * height};
    for (const Real lateral : lateral_offsets) {
      for (const Real wall : {system.domain.lo()[2], system.domain.hi()[2]}) {
        Vec3 target = {source[0] + lateral, source[1], wall + (
            wall == system.domain.lo()[2] ? epsilon : -epsilon)};
        Vec3 lower = target;
        Vec3 upper = target;
        lower[2] -= difference_step;
        upper[2] += difference_step;
        const Real lower_value =
            system.gf.concentration_bounded(source, lower, params);
        const Real upper_value =
            system.gf.concentration_bounded(source, upper, params);
        const Real value =
            system.gf.concentration_bounded(source, target, params);
        const Real derivative =
            (upper_value - lower_value) / (2.0 * difference_step);
        const Real residual = std::abs(derivative * height / value);
        maximum_residual = std::max(maximum_residual, residual);
        if (!(residual < 1.0e-3)) {
          std::cerr << "boundary residual exceeded tolerance: "
                    << residual << "\n";
          std::exit(1);
        }
      }
    }
  }
  std::cout << "  test_boundary_condition: PASSED (max residual="
            << maximum_residual << ")\n";
}

Real mode_expansion(const GreensFunctionParams& params, Real rho,
                    Real z, Real source_z, Real z_lo, Real height) {
  const Real d_eff = params.diff_coeff / params.retardation;
  const Real screening = params.decay_rate / d_eff;
  Real sum = 0.0;
  for (int n = 0; n <= 600; ++n) {
    const Real mode = static_cast<Real>(n) * std::numbers::pi / height;
    const Real kappa = std::sqrt(screening + mode * mode);
    const Real term = (n == 0 ? 1.0 : 2.0)
        * std::cos(mode * (source_z - z_lo))
        * std::cos(mode * (z - z_lo))
        * std::cyl_bessel_k(0, kappa * rho);
    sum += term;
    if (n > 20 && std::abs(term) < 1.0e-15 * std::max(1.0, std::abs(sum))) {
      break;
    }
  }
  return params.source_rate * sum
      / (2.0 * std::numbers::pi * d_eff * height);
}

void test_mode_expansion_agreement() {
  auto system = make_system();
  const Real height = system.domain.size()[2];
  const auto params = test_params(1.0e-4);
  const std::array<Real, 3> rhos = {15.0e-6, 30.0e-6, 60.0e-6};
  const std::array<Real, 3> source_heights = {0.1, 0.35, 0.8};
  const std::array<Real, 3> target_heights = {0.2, 0.55, 0.9};
  Real maximum_relative_error = 0.0;

  for (size_t i = 0; i < rhos.size(); ++i) {
    const Vec3 source = {500.0e-6, 500.0e-6,
                         source_heights[i] * height};
    const Vec3 target = {source[0] + rhos[i], source[1],
                         target_heights[i] * height};
    const Real image_value =
        system.gf.concentration_bounded(source, target, params);
    const Real mode_value = mode_expansion(
        params, rhos[i], target[2], source[2], system.domain.lo()[2], height);
    const Real relative_error =
        std::abs(image_value - mode_value) / std::abs(mode_value);
    maximum_relative_error = std::max(maximum_relative_error, relative_error);
    if (!(relative_error < 1.0e-6)) {
      std::cerr << "mode expansion mismatch: " << relative_error << "\n";
      std::exit(1);
    }
  }
  std::cout << "  test_mode_expansion_agreement: PASSED (max relative error="
            << maximum_relative_error << ")\n";
}

void test_shipped_screening_resolution() {
  auto system = make_system();
  const Real height = system.domain.size()[2];
  const auto params = test_params(1.0e-4);
  const Real source_z = 0.3 * height;
  const Real target_z = 0.7 * height;
  const Real rho = 30.0e-6;
  const Real d_eff = params.diff_coeff / params.retardation;
  const Real screening = std::sqrt(params.decay_rate / d_eff);
  const auto budget = neumann::image_series_budget(
      d_eff, params.decay_rate, 0.0, 0.0, height);
  if (!(screening * height > 0.15 && screening * height < 0.17)
      || budget.max_shells < 60 || budget.max_shells >= neumann::kMaxImageShells) {
    std::cerr << "shipped screening budget was not resolved: kH="
              << screening * height << " shells=" << budget.max_shells << "\n";
    std::exit(1);
  }
  const auto kernel = [=](Real image_z, int) {
    const Real distance = std::hypot(rho, target_z - image_z);
    return params.source_rate
        * std::exp(-screening * distance)
        / (4.0 * std::numbers::pi * d_eff * distance);
  };
  int shell_count = 0;
  int cap_hit = 0;
  const Real image_value = neumann::sum_image_series(
      source_z, 0.0, height, kernel, d_eff, params.decay_rate, 0.0, 0.0,
      neumann::kRelativeTolerance, &shell_count, &cap_hit);
  const Real mode_value = mode_expansion(
      params, rho, target_z, source_z, 0.0, height);
  const Real relative_error =
      std::abs(image_value - mode_value) / std::abs(mode_value);
  if (cap_hit != 0 || shell_count < 50 || !(relative_error < 1.0e-6)) {
    std::cerr << "shipped screening image series under-resolved: shells="
              << shell_count << " cap_hit=" << cap_hit
              << " relative_error=" << relative_error << "\n";
    std::exit(1);
  }
  std::cout << "  test_shipped_screening_resolution: PASSED (kH="
            << screening * height << ", shells=" << shell_count
            << ", relative error=" << relative_error << ")\n";
}

void test_global_mass_balance() {
  auto system = make_system(10.0e-6);
  const Real decay = 1.0e-2;
  const auto params = test_params(decay);
  const Vec3 source = {500.0e-6, 500.0e-6, 50.0e-6};
  Real integral = 0.0;
  for (Int ix = 0; ix < system.domain.nx(); ++ix) {
    for (Int iy = 0; iy < system.domain.ny(); ++iy) {
      for (Int iz = 0; iz < system.domain.nz(); ++iz) {
        integral += system.gf.concentration_bounded(
            source, system.domain.cell_center(ix, iy, iz), params);
      }
    }
  }
  integral *= system.domain.cell_volume();
  const Real relative_error =
      std::abs(decay * integral / params.source_rate - 1.0);
  if (!(relative_error < 0.01)) {
    std::cerr << "mass balance mismatch: " << relative_error << "\n";
    std::exit(1);
  }
  std::cout << "  test_global_mass_balance: PASSED (relative error="
            << relative_error << ")\n";
}

struct ScalarKernel {
  Real target_z;
  Real height;
  Real screening;

  Real operator()(Real image_z, int) const {
    const Real distance = std::abs(image_z - target_z);
    return std::exp(-screening * distance) / (1.0 + distance / height);
  }
};

void test_truncation_and_cap_guard() {
  const ScalarKernel kernel = {37.0e-6, 100.0e-6, 100000.0};
  int shell_count = 0;
  int cap_hit = 0;
  const Real adaptive = neumann::sum_image_series(
      23.0e-6, 0.0, 100.0e-6, kernel, neumann::kRelativeTolerance,
      neumann::kMaxImageShells, &shell_count, &cap_hit);
  if (cap_hit != 0 || shell_count <= 0
      || shell_count >= neumann::kMaxImageShells) {
    std::cerr << "adaptive series did not converge below cap\n";
    std::exit(1);
  }

  Real previous = 0.0;
  for (int m = 0; m <= shell_count; ++m) {
    int ignored_cap = 0;
    const Real partial = neumann::sum_image_series(
        23.0e-6, 0.0, 100.0e-6, kernel, 0.0, m, nullptr,
        &ignored_cap);
    if (m > 0 && !(partial > previous)) {
      std::cerr << "image series was not monotone at shell " << m << "\n";
      std::exit(1);
    }
    previous = partial;
  }
  int ignored_cap = 0;
  const Real next_shell = neumann::sum_image_series(
      23.0e-6, 0.0, 100.0e-6, kernel, 0.0, shell_count + 1, nullptr,
      &ignored_cap) - adaptive;
  if (!(std::abs(next_shell)
        < neumann::kRelativeTolerance * std::abs(adaptive))) {
    std::cerr << "adaptive truncation increment exceeded tolerance\n";
    std::exit(1);
  }

  auto system = make_system();
  const Vec3 source = {500.0e-6, 500.0e-6, 50.0e-6};
  const Vec3 target = {530.0e-6, 500.0e-6, 50.0e-6};
  auto zero_decay = test_params(0.0);
  system.gf.reset_image_series_cap_hits();
  (void)system.gf.concentration_bounded(source, target, zero_decay);
  if (system.gf.image_series_cap_hits() == 0) {
    std::cerr << "zero-decay series did not report a cap hit\n";
    std::exit(1);
  }
  std::cout << "  test_truncation_and_cap_guard: PASSED (adaptive shells="
            << shell_count << ", zero-decay cap hit)\n";
}

void test_parser_screening_guard() {
  auto config = InputParser::default_config();
  for (auto& spec : config.chemicals) {
    if (spec.name == species::BACTERIOCIN_BTUB) {
      spec.decay_rate = 0.0;
      break;
    }
  }
  const char* previous = std::getenv("GUTIBM_STRICT_CONFIG");
  const std::string saved = previous == nullptr ? "" : previous;
  setenv("GUTIBM_STRICT_CONFIG", "1", 1);
  bool threw = false;
  try {
    InputParser::finalize_config(config);
  } catch (const ConfigError&) {
    threw = true;
  }
  if (saved.empty()) {
    unsetenv("GUTIBM_STRICT_CONFIG");
  } else {
    setenv("GUTIBM_STRICT_CONFIG", saved.c_str(), 1);
  }
  if (!threw) {
    std::cerr << "strict low-screening parser guard did not throw\n";
    std::exit(1);
  }
  std::cout << "  test_parser_screening_guard: PASSED\n";
}

}  // namespace

int main() {
  std::cout << "=== Independent Neumann Image-Series Tests ===\n";
  test_boundary_condition();
  test_mode_expansion_agreement();
  test_shipped_screening_resolution();
  test_global_mass_balance();
  test_truncation_and_cap_guard();
  test_parser_screening_guard();
  std::cout << "All independent Neumann image-series tests passed.\n";
  return 0;
}
