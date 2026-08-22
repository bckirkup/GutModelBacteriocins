#include "gpu_kernels.h"
#include "gpu_test_support.h"
#include "metabolic_mode.h"
#include "receptor_utils.h"
#include "vbf_carbon_sink.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#ifdef GUTIBM_CUDA
#include "device_memory.h"
#include <cuda_runtime.h>
#endif

namespace {

#ifdef GUTIBM_CUDA

using gutibm::DeviceBuffer;
using gutibm::gpu::AdvectionParams;
using gutibm::gpu::DomainParams;
using gutibm::gpu::GfSourceParams;
using gutibm::gpu::MechanicsLaunchParams;
using gutibm::gpu::VbfLaunchParams;

constexpr int kNx = 4;
constexpr int kNy = 4;
constexpr int kNz = 4;
constexpr int kCells = kNx * kNy * kNz;
constexpr double kPi = 3.14159265358979323846;
constexpr double kTolerance = 1.0e-10;

template <typename T>
std::vector<T> download(const DeviceBuffer<T>& device, int count) {
  std::vector<T> host(static_cast<size_t>(count));
  device.download(host.data(), static_cast<size_t>(count));
  return host;
}

void synchronize() {
  const cudaError_t status = cudaDeviceSynchronize();
  assert(status == cudaSuccess);
}

bool close(double lhs, double rhs, double tolerance = kTolerance) {
  return std::abs(lhs - rhs) <= tolerance;
}

DomainParams domain_params() {
  DomainParams domain{};
  domain.nx = kNx;
  domain.ny = kNy;
  domain.nz = kNz;
  domain.dx_x = 1.0;
  domain.dx_y = 1.0;
  domain.dx_z = 1.0;
  domain.lo = {0.0, 0.0, 0.0};
  domain.extent = {4.0, 4.0, 4.0};
  domain.periodic = {true, true, false};
  domain.z_lo = 0.0;
  domain.z_hi = 4.0;
  return domain;
}

AdvectionParams zero_advection() {
  AdvectionParams advection{};
  advection.h = 4.0;
  advection.lo_z = 0.0;
  advection.profile_alpha = 1.0;
  return advection;
}

struct PeriodicCoefficients {
  double gamma = 0.0;
  double corner = 0.0;
  double denominator = 1.0;
  std::vector<double> correction;
};

PeriodicCoefficients periodic_coefficients(int size, double alpha) {
  PeriodicCoefficients coefficients;
  const double diagonal_value = 1.0 + 2.0 * alpha;
  coefficients.gamma = -diagonal_value;
  coefficients.corner = -alpha;

  std::vector<double> lower(static_cast<size_t>(size - 1), -alpha);
  std::vector<double> upper(static_cast<size_t>(size - 1), -alpha);
  std::vector<double> diagonal(static_cast<size_t>(size), diagonal_value);
  diagonal.front() -= coefficients.gamma;
  diagonal.back() -=
      coefficients.corner * coefficients.corner / coefficients.gamma;

  coefficients.correction.assign(static_cast<size_t>(size), 0.0);
  coefficients.correction.front() = coefficients.gamma;
  coefficients.correction.back() = coefficients.corner;
  for (int row = 1; row < size; ++row) {
    const double multiplier = lower[static_cast<size_t>(row - 1)]
        / diagonal[static_cast<size_t>(row - 1)];
    diagonal[static_cast<size_t>(row)] -=
        multiplier * upper[static_cast<size_t>(row - 1)];
    coefficients.correction[static_cast<size_t>(row)] -=
        multiplier * coefficients.correction[static_cast<size_t>(row - 1)];
  }
  coefficients.correction.back() /= diagonal.back();
  for (int row = size - 2; row >= 0; --row) {
    coefficients.correction[static_cast<size_t>(row)] =
        (coefficients.correction[static_cast<size_t>(row)] -
         upper[static_cast<size_t>(row)] *
             coefficients.correction[static_cast<size_t>(row + 1)]) /
        diagonal[static_cast<size_t>(row)];
  }
  coefficients.denominator =
      1.0 + coefficients.correction.front() +
      coefficients.corner * coefficients.correction.back() /
          coefficients.gamma;
  return coefficients;
}

void test_field_update() {
  DeviceBuffer<double> concentration(kCells);
  DeviceBuffer<double> reaction(kCells);
  DeviceBuffer<double> clip(1);
  std::vector<double> host_reaction(kCells, 0.0);
  host_reaction[0] = -2.0;
  concentration.upload(std::vector<double>(kCells, 1.0));
  reaction.upload(host_reaction);
  clip.upload(std::vector<double>{0.0});
  gutibm::gpu::launch_field_update_kernel(
      concentration.data(), reaction.data(), kCells, 1, 1.0, clip.data(), 1.0,
      kNx, kNy, kNz, 0, kNx, nullptr);
  synchronize();
  assert(download(concentration, kCells)[0] == 0.0);
  assert(close(download(clip, 1)[0], 1.0));

  gutibm::gpu::launch_field_update_kernel(
      concentration.data(), reaction.data(), kCells, 1, 1.0, clip.data(), 1.0,
      kNx, kNy, kNz, 0, kNx, nullptr);
  synchronize();
  // Change-detector: the second deficit is 2.0, so cumulative clipping is 3.0.
  assert(close(download(clip, 1)[0], 3.0));

  const std::array<double, 3> dts{0.5, 1.0, 2.0};
  double previous = 1.0;
  for (const double dt : dts) {
    concentration.upload(std::vector<double>(kCells, 1.0));
    reaction.upload(std::vector<double>(kCells, 0.5));
    gutibm::gpu::launch_field_update_kernel(
        concentration.data(), reaction.data(), kCells, 1, dt, nullptr, 1.0,
        kNx, kNy, kNz, 0, kNx, nullptr);
    synchronize();
    const double value = download(concentration, kCells)[0];
    assert(value > previous);
    previous = value;
  }
  concentration.upload(std::vector<double>(kCells, 1.0));
  gutibm::gpu::launch_field_update_kernel(
      concentration.data(), reaction.data(), kCells, 1, 0.0, nullptr, 1.0,
      kNx, kNy, kNz, 0, kNx, nullptr);
  synchronize();
  assert(download(concentration, kCells)[0] == 1.0);
}

void test_apply_boundaries() {
  DeviceBuffer<double> concentration(kCells);
  DeviceBuffer<double> boundary(1);
  concentration.upload(std::vector<double>(kCells, 0.25));
  boundary.upload(std::vector<double>{0.75});
  gutibm::gpu::launch_apply_boundaries_kernel(
      concentration.data(), kNx, kNy, kNz, 1, boundary.data(), nullptr);
  synchronize();
  const auto result = download(concentration, kCells);
  for (int iy = 0; iy < kNy; ++iy) {
    for (int ix = 0; ix < kNx; ++ix) {
      const int bottom = iy * kNx + ix;
      const int interior = kNx * kNy + iy * kNx + ix;
      const int top = (kNz - 1) * kNx * kNy + iy * kNx + ix;
      assert(result[bottom] == 0.75);
      assert(result[interior] == 0.25);
      assert(result[top] == 0.25);
    }
  }
}

void test_grid_coupling() {
  constexpr int agents = 5;
  DeviceBuffer<double> x(agents);
  DeviceBuffer<double> y(agents);
  DeviceBuffer<double> z(agents);
  DeviceBuffer<int> cells(agents);
  DeviceBuffer<int> state(agents);
  x.upload(std::vector<double>{0.1, 2.1, 3.9, -1.0, 1.0});
  y.upload(std::vector<double>{1.1, 2.1, 3.9, 5.0, 0.0});
  z.upload(std::vector<double>{0.1, 3.1, 3.9, 4.0, 2.0});
  state.upload(std::vector<int>{0, 0, 0, 0, 3});
  gutibm::gpu::launch_grid_coupling_kernel(
      x.data(), y.data(), z.data(), cells.data(), state.data(), 0.0, 0.0, 0.0,
      1.0, 1.0, 1.0, kNx, kNy, kNz, agents, nullptr);
  synchronize();
  const auto mapped = download(cells, agents);
  assert(mapped[0] == 4);
  assert(mapped[1] == 2 + 2 * kNx + 3 * kNx * kNy);
  assert(mapped[2] == 3 + 3 * kNx + 3 * kNx * kNy);
  assert(mapped[3] == (kNz - 1) * kNx * kNy + (kNy - 1) * kNx);
  assert(mapped[4] == -1);
}

struct MetabolismRun {
  double carbon_reaction = 0.0;
  double carbon_reaction_agent1 = 0.0;
  double uptake = 0.0;
  double fepA_expression = 0.0;
  double noniron_expression = 0.0;
  double acetate_reaction = 0.0;
  double acetate_reaction_agent1 = 0.0;
  double biomass = 0.0;
  double biomass_agent1 = 0.0;
  double mu = 0.0;
  double mu_agent1 = 0.0;
  double fermentation_fraction = 0.0;
  double fermentation_fraction_agent1 = 0.0;
  double demand = 0.0;
  double limited_agents = 0.0;
  double maintenance = 0.0;
  double maintenance_shortfall = 0.0;
  double maintenance_limited_agents = 0.0;
};

MetabolismRun run_metabolism(double seed, double maximum_growth,
                             bool fur_enabled = false,
                             bool acetate_enabled = false,
                             double iron_concentration = 1.0,
                             int uptake_limit_mode = 0,
                             double effective_diffusivity = 0.0,
                             double carbon_maintenance_rate = 0.0,
                             bool metabolic_switch_enabled = false,
                             double oxygen_concentration = 0.0,
                             std::array<double, 2> initial_fermentation = {
                                 0.0, 0.0},
                             bool activate_second_agent = false,
                             double acetate_concentration = 2.0,
                             double acetate_scavenge_rate = 0.1,
                             double maintenance_rate = 0.0) {
  constexpr int agents = 2;
  DeviceBuffer<double> concentration(kCells);
  DeviceBuffer<double> iron(kCells);
  DeviceBuffer<double> b12(kCells);
  DeviceBuffer<double> acetate(kCells);
  DeviceBuffer<double> eut(kCells);
  DeviceBuffer<double> oxygen(kCells);
  DeviceBuffer<double> reaction_carbon(kCells);
  DeviceBuffer<double> reaction_iron(kCells);
  DeviceBuffer<double> reaction_b12(kCells);
  DeviceBuffer<double> reaction_acetate(kCells);
  DeviceBuffer<double> mu(agents);
  DeviceBuffer<double> biomass(agents);
  DeviceBuffer<double> radius(agents);
  DeviceBuffer<double> mass(agents);
  DeviceBuffer<double> age(agents);
  DeviceBuffer<double> fermentation_fraction(agents);
  DeviceBuffer<double> mu_max(agents);
  DeviceBuffer<double> km_b12(agents);
  DeviceBuffer<double> km_carbon(agents);
  DeviceBuffer<double> receptor(gutibm::NUM_RECEPTORS * agents);
  DeviceBuffer<double> receptor_base(gutibm::NUM_RECEPTORS * agents);
  DeviceBuffer<double> ligand(gutibm::NUM_RECEPTORS * agents);
  DeviceBuffer<int> iron_receptor(gutibm::NUM_RECEPTORS);
  DeviceBuffer<int> cells(agents);
  DeviceBuffer<int> state(agents);
  DeviceBuffer<int> loci(agents);
  DeviceBuffer<double> amelioration(agents);
  DeviceBuffer<double> uptake(4);
  DeviceBuffer<double> maintenance_available(kCells);
  DeviceBuffer<double> uptake_limit(5);

  concentration.upload(std::vector<double>(kCells, 1.0));
  iron.upload(std::vector<double>(kCells, iron_concentration));
  b12.upload(std::vector<double>(kCells, 1.0));
  acetate.upload(std::vector<double>(
      kCells, acetate_enabled ? acetate_concentration : 0.0));
  eut.upload(std::vector<double>(kCells, 0.0));
  oxygen.upload(std::vector<double>(kCells, oxygen_concentration));
  reaction_carbon.upload(std::vector<double>(kCells, seed));
  reaction_iron.upload(std::vector<double>(kCells, 0.5));
  reaction_b12.upload(std::vector<double>(kCells, 0.0));
  reaction_acetate.upload(std::vector<double>(kCells, 0.0));
  mu.upload(std::vector<double>(agents, 0.0));
  biomass.upload(std::vector<double>(agents, 1.0));
  radius.upload(std::vector<double>(agents, 1.0));
  mass.upload(std::vector<double>(agents, 1.0));
  age.upload(std::vector<double>(agents, 0.0));
  fermentation_fraction.upload(std::vector<double>{
      initial_fermentation[0], initial_fermentation[1]});
  mu_max.upload(std::vector<double>(agents, maximum_growth));
  km_b12.upload(std::vector<double>(agents, 0.1));
  km_carbon.upload(std::vector<double>(agents, 0.1));
  receptor.upload(std::vector<double>(gutibm::NUM_RECEPTORS * agents, 1.0));
  receptor_base.upload(
      std::vector<double>(gutibm::NUM_RECEPTORS * agents, 2.0));
  ligand.upload(std::vector<double>(gutibm::NUM_RECEPTORS * agents, 1.0));
  std::vector<int> iron_receptor_flags(gutibm::NUM_RECEPTORS, 0);
  for (int receptor = 0; receptor < gutibm::NUM_RECEPTORS; ++receptor) {
    iron_receptor_flags[static_cast<size_t>(receptor)] =
        gutibm::is_iron_receptor(receptor) ? 1 : 0;
  }
  iron_receptor.upload(iron_receptor_flags);
  cells.upload(activate_second_agent ? std::vector<int>{4, 5}
                                     : std::vector<int>{4, -1});
  state.upload(activate_second_agent ? std::vector<int>{0, 0}
                                     : std::vector<int>{0, 3});
  loci.upload(std::vector<int>(agents, 0));
  amelioration.upload(std::vector<double>(agents, 0.0));
  uptake.upload(std::vector<double>(4, 0.0));
  maintenance_available.upload(std::vector<double>(kCells, 1.0));
  uptake_limit.upload(std::vector<double>(5, 0.0));

  gutibm::gpu::launch_metabolism_kernel(
      concentration.data(), iron.data(), b12.data(), acetate.data(), eut.data(),
      oxygen.data(), reaction_carbon.data(), reaction_iron.data(),
      reaction_b12.data(), reaction_acetate.data(),
      mu.data(), biomass.data(), radius.data(), mass.data(), age.data(),
      fermentation_fraction.data(), cells.data(), state.data(), mu_max.data(),
      km_b12.data(),
      km_carbon.data(), receptor.data(), receptor_base.data(), ligand.data(),
      iron_receptor.data(), loci.data(),
      amelioration.data(), agents, agents, agents, gutibm::NUM_RECEPTORS,
      1.0, 1.0, 1.0,
      0.1, 0.1, 0.1, 0.1,
      maintenance_rate, 0.0, 1.0, carbon_maintenance_rate, 0.0, 0.0, 0.2,
      0.1, 0.1, 0.1,
      1, 1, 1, fur_enabled ? 1 : 0, 1.0e-5, 10.0, 5.0,
      acetate_enabled ? 1 : 0, 3.0e-4, 0.25, acetate_scavenge_rate, 2.0,
      metabolic_switch_enabled ? 1 : 0, 0.0, 1.0,
      metabolic_switch_enabled ? 1 : 0, 3.0e-4, 1.0, 0.55, 1.0, 4.1,
      0.0, 1.0, 15.0, 0, 0.8, 50.0, 4.76, 7.0, uptake.data(),
      maintenance_available.data(),
      uptake_limit_mode, effective_diffusivity, effective_diffusivity,
      uptake_limit.data(),
      kNx, kNy, kNx, 0, kNx, 0, nullptr);
  synchronize();
  const auto reaction = download(reaction_carbon, kCells);
  const auto acetate_result = download(reaction_acetate, kCells);
  const auto expression =
      download(receptor, gutibm::NUM_RECEPTORS * agents);
  const auto biomass_result = download(biomass, agents);
  const auto mu_result = download(mu, agents);
  const auto fermentation_result = download(fermentation_fraction, agents);
  const auto uptake_host = download(uptake, 4);
  const auto uptake_limit_host = download(uptake_limit, 5);
  MetabolismRun result;
  result.carbon_reaction = reaction[4];
  result.carbon_reaction_agent1 = reaction[5];
  result.uptake = uptake_host[0];
  result.fepA_expression = expression[1 * agents];
  result.noniron_expression = expression[0];
  result.acetate_reaction = acetate_result[4];
  result.acetate_reaction_agent1 = acetate_result[5];
  result.biomass = biomass_result[0];
  result.biomass_agent1 = biomass_result[1];
  result.mu = mu_result[0];
  result.mu_agent1 = mu_result[1];
  result.fermentation_fraction = fermentation_result[0];
  result.fermentation_fraction_agent1 = fermentation_result[1];
  result.demand = uptake_limit_host[0];
  result.limited_agents = uptake_limit_host[2];
  result.maintenance = uptake_host[2];
  result.maintenance_shortfall = uptake_host[3];
  result.maintenance_limited_agents = uptake_limit_host[4];
  return result;
}

void test_metabolism() {
  const MetabolismRun zero_seed = run_metabolism(0.0, 1.0e-3);
  const MetabolismRun seeded = run_metabolism(0.25, 1.0e-3);
  assert(zero_seed.carbon_reaction < 0.0);
  assert(close(seeded.carbon_reaction - zero_seed.carbon_reaction, 0.25));
  const MetabolismRun low = run_metabolism(0.0, 5.0e-4);
  const MetabolismRun high = run_metabolism(0.0, 2.0e-3);
  assert(low.uptake < zero_seed.uptake);
  assert(zero_seed.uptake < high.uptake);
}

void test_metabolism_negative_mu() {
  const MetabolismRun result = run_metabolism(
      0.0, 1.0e-3, false, false, 1.0, 0, 0.0, 0.0, false, 0.0,
      std::array<double, 2>{0.0, 0.0}, false, 2.0, 0.1, 2.0e-3);
  assert(std::isfinite(result.mu));
  assert(result.mu < 0.0);
  assert(std::isfinite(result.biomass));
  assert(result.biomass < 1.0);
  assert(close(result.carbon_reaction, 0.0));
  assert(close(result.uptake, 0.0));
}

void test_metabolism_uptake_limit() {
  const MetabolismRun unlimited = run_metabolism(0.0, 1.0e-3);
  const MetabolismRun capped =
      run_metabolism(0.0, 1.0e-3, false, false, 1.0, 1, 1.0e-6);
  const MetabolismRun voxel =
      run_metabolism(0.0, 1.0e-3, false, false, 1.0, 2, 1.0e-6);
  assert(unlimited.demand > 0.0);
  assert(close(capped.demand, unlimited.demand));
  assert(unlimited.limited_agents == 0.0);
  assert(capped.limited_agents == 1.0);
  assert(capped.uptake < unlimited.uptake);
  assert(capped.uptake > 0.0);
  assert(capped.mu < unlimited.mu);
  assert(voxel.uptake >= capped.uptake);
  assert(voxel.uptake <= unlimited.uptake);
}

void test_metabolism_maintenance() {
  const MetabolismRun result = run_metabolism(
      0.0, 0.0, false, false, 1.0, 0, 0.0, 1.0e-2);
  assert(result.maintenance > 0.0);
  assert(result.carbon_reaction < 0.0);
  assert(result.maintenance_shortfall == 0.0);
  assert(result.maintenance_limited_agents == 0.0);
}

void test_metabolism_maintenance_uptake_limit() {
  const MetabolismRun result = run_metabolism(
      0.0, 0.0, false, false, 1.0, 1, 0.1, 10.0);
  assert(std::isfinite(result.maintenance));
  assert(std::isfinite(result.maintenance_shortfall));
  assert(std::isfinite(result.maintenance_limited_agents));
  assert(result.maintenance > 0.0);
  assert(result.maintenance_shortfall > 0.0);
  assert(result.maintenance_limited_agents == 1.0);
  assert(result.maintenance + result.maintenance_shortfall > 9.0);
}

void test_metabolism_fur() {
  const MetabolismRun result = run_metabolism(
      0.0, 1.0e-3, true, false, 1.0e-6);
  const double expected_factor = 1.0 + 10.0 * 1.0e-5
      / (1.0e-5 + 1.0e-6);
  assert(result.fepA_expression == 5.0);
  assert(result.noniron_expression == 2.0);
  assert(expected_factor * 2.0 > 5.0);
  const MetabolismRun disabled = run_metabolism(
      0.0, 1.0e-3, false, false, 1.0e-6);
  assert(disabled.fepA_expression == 1.0);
  assert(disabled.noniron_expression == 1.0);
}

void test_metabolism_acetate() {
  const MetabolismRun result = run_metabolism(
      0.0, 1.0e-3, false, true, 1.0);
  const double updated_biomass = std::max(
      result.biomass, 1.0e-20);
  const double expected = result.mu > 3.0e-4
      ? 0.25 * updated_biomass - 0.1 * 2.0 / (2.0 + 2.0)
          * updated_biomass
      : -0.1 * 2.0 / (2.0 + 2.0) * updated_biomass;
  assert(close(result.acetate_reaction, expected));
  assert(result.biomass > 1.0);
}

void test_diffuse_x_periodic() {
  const double alpha = 0.1;
  const PeriodicCoefficients coefficients = periodic_coefficients(kNx, alpha);
  DeviceBuffer<double> correction(kNx);
  correction.upload(coefficients.correction);
  DeviceBuffer<double> field(kCells);
  field.upload(std::vector<double>(kCells, 2.0));
  gutibm::gpu::launch_diffuse_x_periodic(
      field.data(), kNx, kNy, kNz, alpha, coefficients.gamma,
      coefficients.corner, coefficients.denominator, correction.data(), nullptr);
  synchronize();
  const auto uniform = download(field, kCells);
  for (const double value : uniform) assert(close(value, 2.0));
  assert(close(std::accumulate(uniform.begin(), uniform.end(), 0.0),
               2.0 * kCells));

  std::vector<double> mode(kCells, 0.0);
  for (int iz = 0; iz < kNz; ++iz) {
    for (int iy = 0; iy < kNy; ++iy) {
      for (int ix = 0; ix < kNx; ++ix) {
        mode[iz * kNx * kNy + iy * kNx + ix] =
            std::sin(2.0 * kPi * ix / kNx);
      }
    }
  }
  field.upload(mode);
  gutibm::gpu::launch_diffuse_x_periodic(
      field.data(), kNx, kNy, kNz, alpha, coefficients.gamma,
      coefficients.corner, coefficients.denominator, correction.data(), nullptr);
  synchronize();
  const auto result = download(field, kCells);
  const double eigenvalue =
      1.0 + 2.0 * alpha - 2.0 * alpha * std::cos(2.0 * kPi / kNx);
  for (int cell = 0; cell < kCells; ++cell) {
    assert(close(result[cell], mode[cell] / eigenvalue));
  }
}

void test_diffuse_y_periodic() {
  const double alpha = 0.1;
  const PeriodicCoefficients coefficients = periodic_coefficients(kNy, alpha);
  DeviceBuffer<double> correction(kNy);
  correction.upload(coefficients.correction);
  DeviceBuffer<double> field(kCells);
  field.upload(std::vector<double>(kCells, 2.0));
  gutibm::gpu::launch_diffuse_y_periodic(
      field.data(), kNx, kNy, kNz, 0, kNx, alpha, coefficients.gamma,
      coefficients.corner, coefficients.denominator, correction.data(), nullptr);
  synchronize();
  const auto uniform = download(field, kCells);
  for (const double value : uniform) assert(close(value, 2.0));
  assert(close(std::accumulate(uniform.begin(), uniform.end(), 0.0),
               2.0 * kCells));

  std::vector<double> mode(kCells, 0.0);
  for (int iz = 0; iz < kNz; ++iz) {
    for (int iy = 0; iy < kNy; ++iy) {
      for (int ix = 0; ix < kNx; ++ix) {
        mode[iz * kNx * kNy + iy * kNx + ix] =
            std::sin(2.0 * kPi * iy / kNy);
      }
    }
  }
  field.upload(mode);
  const double before = std::accumulate(mode.begin(), mode.end(), 0.0);
  gutibm::gpu::launch_diffuse_y_periodic(
      field.data(), kNx, kNy, kNz, 0, kNx, alpha, coefficients.gamma,
      coefficients.corner, coefficients.denominator, correction.data(), nullptr);
  synchronize();
  const auto result = download(field, kCells);
  const double after = std::accumulate(result.begin(), result.end(), 0.0);
  assert(close(before, after));
}

void test_diffuse_z_bounded() {
  DeviceBuffer<double> field(kCells);
  DeviceBuffer<double> exchange(1);
  field.upload(std::vector<double>(kCells, 2.0));
  exchange.upload(std::vector<double>{0.0});
  gutibm::gpu::launch_diffuse_z_bounded(
      field.data(), kNx, kNy, kNz, 0, kNx, 0.1, 2.0, 0, 0.0, 0.0, 1.0,
      exchange.data(), nullptr);
  synchronize();
  for (const double value : download(field, kCells)) assert(close(value, 2.0));
  assert(close(download(exchange, 1)[0], 0.0));
}

void test_set_epithelial_boundary() {
  constexpr int storage_nx = 6;
  constexpr int owned_x_begin = 1;
  constexpr int owned_x_end = 3;
  constexpr int storage_cells = storage_nx * kNy * kNz;
  constexpr int face_cells = (owned_x_end - owned_x_begin) * kNy;
  const std::array<double, 3> targets{0.0, 0.5, 1.0};
  std::array<double, 3> injected{};
  for (size_t run = 0; run < targets.size(); ++run) {
    DeviceBuffer<double> field(storage_cells);
    DeviceBuffer<double> accounting(1);
    field.upload(std::vector<double>(storage_cells, 0.25));
    accounting.upload(std::vector<double>{0.0});
    gutibm::gpu::launch_set_epithelial_boundary(
        field.data(), storage_nx, kNy, owned_x_begin, owned_x_end, targets[run],
        2.0, accounting.data(), nullptr);
    synchronize();
    const auto result = download(field, storage_cells);
    injected[run] = download(accounting, 1)[0];
    for (int iy = 0; iy < kNy; ++iy) {
      for (int ix = owned_x_begin; ix < owned_x_end; ++ix) {
        const int boundary = iy * storage_nx + ix;
        assert(result[boundary] == targets[run]);
        for (int iz = 1; iz < kNz; ++iz) {
          const int interior =
              iz * kNy * storage_nx + iy * storage_nx + ix;
          assert(result[interior] == 0.25);
        }
      }
    }
  }
  assert(injected[0] < injected[1]);
  assert(injected[1] < injected[2]);
  assert(close(injected[2] / injected[1], 3.0));
  assert(face_cells > 0);
}

void test_set_luminal_neumann() {
  DeviceBuffer<double> field(kCells);
  std::vector<double> initial(kCells, 0.0);
  for (int iy = 0; iy < kNy; ++iy) {
    for (int ix = 0; ix < kNx; ++ix) {
      initial[(kNz - 2) * kNx * kNy + iy * kNx + ix] = 3.0;
      initial[(kNz - 1) * kNx * kNy + iy * kNx + ix] = 7.0;
    }
  }
  field.upload(initial);
  gutibm::gpu::launch_set_luminal_neumann(
      field.data(), kNx, kNy, kNz, 0, kNx, nullptr);
  synchronize();
  const auto result = download(field, kCells);
  for (int iy = 0; iy < kNy; ++iy) {
    for (int ix = 0; ix < kNx; ++ix) {
      assert(result[(kNz - 1) * kNx * kNy + iy * kNx + ix] == 3.0);
      assert(result[kNx + iy * kNx + ix] == 0.0);
    }
  }
}

void test_shift_z_gradient() {
  DeviceBuffer<double> field(kCells);
  field.upload(std::vector<double>(kCells, 0.0));
  gutibm::gpu::launch_shift_z_gradient(
      field.data(), kNx, kNy, kNz, 0, kNx, 1.0, 2.0, 2.0, 5.0, 1.0,
      nullptr);
  synchronize();
  const auto result = download(field, kCells);
  assert(close(result[0], 5.0));
  assert(result[kNx * kNy] > result[2 * kNx * kNy]);
  assert(close(result[3 * kNx * kNy], result[2 * kNx * kNy]));
}

void test_clamp_nonneg() {
  DeviceBuffer<double> field(kCells);
  std::vector<double> initial(kCells, 2.0);
  initial[0] = -1.0;
  field.upload(initial);
  gutibm::gpu::launch_clamp_nonneg(
      field.data(), kNx, kNy, kNz, 0, kNx, nullptr);
  synchronize();
  const auto result = download(field, kCells);
  assert(result[0] == 0.0);
  for (int cell = 1; cell < kCells; ++cell) assert(result[cell] == 2.0);
}

void test_superpose() {
  const DomainParams domain = domain_params();
  const AdvectionParams advection = zero_advection();
  DeviceBuffer<double> x(2);
  DeviceBuffer<double> y(2);
  DeviceBuffer<double> z(2);
  DeviceBuffer<GfSourceParams> params(2);
  DeviceBuffer<double> one(kCells);
  DeviceBuffer<double> two(kCells);
  DeviceBuffer<double> combined(kCells);
  x.upload(std::vector<double>{1.5, 2.5});
  y.upload(std::vector<double>{1.5, 1.5});
  z.upload(std::vector<double>{1.5, 1.5});
  params.upload(std::vector<GfSourceParams>{{1.0, 1.0, 1.0, 0.0},
                                             {1.0, 1.0, 1.0, 0.0}});
  one.upload(std::vector<double>(kCells, 0.0));
  two.upload(std::vector<double>(kCells, 0.0));
  combined.upload(std::vector<double>(kCells, 0.0));
  gutibm::gpu::launch_superpose_kernel(
      x.data(), y.data(), z.data(), params.data(), one.data(), domain, advection,
      1, 1, 1, 1, nullptr);
  gutibm::gpu::launch_superpose_kernel(
      x.data() + 1, y.data() + 1, z.data() + 1, params.data() + 1, two.data(),
      domain, advection, 1, 1, 1, 1, nullptr);
  gutibm::gpu::launch_superpose_kernel(
      x.data(), y.data(), z.data(), params.data(), combined.data(), domain,
      advection, 2, 1, 1, 1, nullptr);
  synchronize();
  const auto one_host = download(one, kCells);
  const auto two_host = download(two, kCells);
  const auto combined_host = download(combined, kCells);
  for (int cell = 0; cell < kCells; ++cell) {
    assert(close(combined_host[cell], one_host[cell] + two_host[cell]));
  }
  const int near_source = kNx * kNy + kNx;
  const int far_source = kNx * kNy + kNx + 3;
  assert(one_host[near_source] > one_host[far_source]);
}

void test_fmm_far_local() {
  DeviceBuffer<double> local(1);
  DeviceBuffer<double> center(3);
  DeviceBuffer<int> cell_leaf(kCells);
  DeviceBuffer<double> output(kCells);
  local.upload(std::vector<double>{1.0});
  center.upload(std::vector<double>{1.5, 1.5, 1.5});
  cell_leaf.upload(std::vector<int>(kCells, 0));
  output.upload(std::vector<double>(kCells, 0.0));
  gutibm::gpu::launch_fmm_far_local_kernel(
      local.data(), center.data(), cell_leaf.data(), nullptr, output.data(),
      kCells, 1, 1, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, kNx, kNy, nullptr);
  synchronize();
  for (const double value : download(output, kCells)) assert(value == 1.0);
  gutibm::gpu::launch_fmm_far_local_kernel(
      local.data(), center.data(), cell_leaf.data(), nullptr, output.data(),
      kCells, 1, 1, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, kNx, kNy, nullptr);
  synchronize();
  for (const double value : download(output, kCells)) assert(value == 2.0);
}

void test_vbf_coupling() {
  constexpr int storage_nx = 6;
  constexpr int storage_cells = storage_nx * kNy * kNz;
  DeviceBuffer<double> reaction_carbon(storage_cells);
  DeviceBuffer<double> reaction_iron(storage_cells);
  DeviceBuffer<double> reaction_oxygen(storage_cells);
  DeviceBuffer<double> reaction_acetate(storage_cells);
  DeviceBuffer<double> reaction_mucin(storage_cells);
  DeviceBuffer<double> concentration(storage_cells);
  reaction_carbon.upload(std::vector<double>(storage_cells, 0.0));
  reaction_iron.upload(std::vector<double>(storage_cells, 0.0));
  reaction_oxygen.upload(std::vector<double>(storage_cells, 0.0));
  reaction_acetate.upload(std::vector<double>(storage_cells, 0.0));
  reaction_mucin.upload(std::vector<double>(storage_cells, 0.0));
  concentration.upload(std::vector<double>(storage_cells, 1.0));
  VbfLaunchParams parameters{};
  parameters.storage_nx = storage_nx;
  parameters.owned_x_begin = 1;
  parameters.owned_x_end = 3;
  parameters.global_nx = kNx;
  parameters.ny = kNy;
  parameters.nz = kNz;
  parameters.dx_x = 1.0;
  parameters.dx_y = 1.0;
  parameters.dx_z = 1.0;
  parameters.mucin_liberation = 0.1;
  parameters.vbf_density = 1.0;
  gutibm::gpu::launch_vbf_coupling_kernel(
      2 * kNy * kNz, parameters, reaction_carbon.data(), concentration.data(),
      reaction_iron.data(), concentration.data(), reaction_oxygen.data(),
      concentration.data(), reaction_acetate.data(), reaction_mucin.data(),
      concentration.data(), nullptr, 1.0, nullptr);
  synchronize();
  const auto result = download(reaction_carbon, storage_cells);
  for (int iz = 0; iz < kNz; ++iz) {
    for (int iy = 0; iy < kNy; ++iy) {
      for (int ix = 1; ix < 3; ++ix) {
        const int cell = iz * storage_nx * kNy + iy * storage_nx + ix;
        assert(result[cell] > 0.0);
      }
    }
  }
}

void test_vbf_implicit_carbon_sink() {
  DeviceBuffer<double> reaction_carbon(kCells);
  DeviceBuffer<double> concentration(kCells);
  reaction_carbon.upload(std::vector<double>(kCells, 0.0));
  concentration.upload(std::vector<double>(kCells, 1.0e-4));

  VbfLaunchParams parameters{};
  parameters.storage_nx = kNx;
  parameters.owned_x_begin = 0;
  parameters.owned_x_end = kNx;
  parameters.global_nx = kNx;
  parameters.ny = kNy;
  parameters.nz = kNz;
  parameters.dx_x = 1.0;
  parameters.dx_y = 1.0;
  parameters.dx_z = 1.0;
  parameters.carbon_sink_vmax = 5.0e-3;
  parameters.carbon_sink_km = 1.0e-4;
  constexpr double dt = 60.0;
  const double expected_sink = gutibm::vbf::implicit_carbon_sink(
      1.0e-4, parameters.carbon_sink_vmax, parameters.carbon_sink_km, dt);
  const double explicit_sink = parameters.carbon_sink_vmax * 1.0e-4
      / (parameters.carbon_sink_km + 1.0e-4);
  assert(explicit_sink > 100.0 * expected_sink);

  gutibm::gpu::launch_vbf_coupling_kernel(
      kCells, parameters, reaction_carbon.data(), concentration.data(),
      nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
      dt, nullptr);
  synchronize();
  for (const double value : download(reaction_carbon, kCells)) {
    assert(std::abs(value + expected_sink) < 1.0e-12);
  }
}

void test_o2_depletion() {
  constexpr int storage_nx = 6;
  constexpr int storage_cells = storage_nx * kNy * kNz;
  DeviceBuffer<double> reaction(storage_cells);
  DeviceBuffer<double> mu(2);
  DeviceBuffer<double> fermentation_fraction(2);
  DeviceBuffer<int> cells(2);
  DeviceBuffer<int> state(2);
  reaction.upload(std::vector<double>(storage_cells, 0.0));
  mu.upload(std::vector<double>{1.0, 1.0});
  fermentation_fraction.upload(std::vector<double>{0.0, 0.0});
  cells.upload(std::vector<int>{1, 1});
  state.upload(std::vector<int>{0, 0});
  gutibm::gpu::launch_o2_depletion_kernel(
      reaction.data(), mu.data(), fermentation_fraction.data(), cells.data(),
      state.data(), 2, 0.5, 0.1, 0, 1.0, kNx, kNy, storage_nx, 1, 3, 1,
      nullptr);
  synchronize();
  const auto result = download(reaction, storage_cells);
  // Change-detector: two agents each consume 0.5*1.0+0.1 = 0.6.
  assert(close(result[1], -1.2));
  assert(result[1] < 0.0);
}

void test_o2_depletion_metabolic_mode() {
  constexpr int storage_nx = 6;
  constexpr int storage_cells = storage_nx * kNy * kNz;
  constexpr double q_consumption = 0.5;
  constexpr double q_maintenance = 0.1;
  constexpr double mu = 1.0;
  DeviceBuffer<double> reaction(storage_cells);
  DeviceBuffer<double> mu_realized(3);
  DeviceBuffer<double> fermentation_fraction(3);
  DeviceBuffer<int> cells(3);
  DeviceBuffer<int> state(3);
  reaction.upload(std::vector<double>(storage_cells, 0.0));
  mu_realized.upload(std::vector<double>(3, mu));
  fermentation_fraction.upload(std::vector<double>{0.0, 0.5, 1.0});
  cells.upload(std::vector<int>{1, 2, 3});
  state.upload(std::vector<int>(3, 0));
  gutibm::gpu::launch_o2_depletion_kernel(
      reaction.data(), mu_realized.data(), fermentation_fraction.data(),
      cells.data(), state.data(), 3, q_consumption, q_maintenance, 1, 1.0,
      kNx, kNy, storage_nx, 0, kNx, 0, nullptr);
  synchronize();
  const auto result = download(reaction, storage_cells);
  const double expected_aerobic = -(q_consumption * mu + q_maintenance);
  const double expected_overflow =
      -(q_consumption * mu * 0.5 + q_maintenance);
  const double expected_fermentative = -q_maintenance;
  assert(std::isfinite(result[1]));
  assert(std::isfinite(result[2]));
  assert(std::isfinite(result[3]));
  assert(result[1] < result[2]);
  assert(result[2] < result[3]);
  assert(close(result[1], expected_aerobic));
  assert(close(result[2], expected_overflow));
  assert(close(result[3], expected_fermentative));
}

void test_metabolism_metabolic_mode() {
  const std::array<double, 2> fractions{0.25, 0.75};
  const MetabolismRun result = run_metabolism(
      0.0, 1.0e-3, false, true, 1.0, 0, 0.0, 0.0, true, 0.0,
      fractions, true, 0.0, 0.0);
  const double aerobic_cost = gutibm::metabolic_mode::interpolate(
      1.0, 4.1, result.fermentation_fraction);
  const double fermentative_cost = gutibm::metabolic_mode::interpolate(
      1.0, 4.1, result.fermentation_fraction_agent1);
  const double aerobic_growth = result.biomass - 1.0;
  const double fermentative_growth = result.biomass_agent1 - 1.0;
  assert(aerobic_growth > 0.0);
  assert(fermentative_growth > 0.0);
  assert(result.fermentation_fraction >= 0.0);
  assert(result.fermentation_fraction <= 1.0);
  assert(result.fermentation_fraction_agent1 >= 0.0);
  assert(result.fermentation_fraction_agent1 <= 1.0);
  const double aerobic_observed_cost =
      -result.carbon_reaction / aerobic_growth / 0.1;
  const double fermentative_observed_cost =
      -result.carbon_reaction_agent1 / fermentative_growth / 0.1;
  assert(std::isfinite(aerobic_observed_cost));
  assert(std::isfinite(fermentative_observed_cost));
  assert(aerobic_observed_cost >= 0.0);
  assert(fermentative_observed_cost >= 0.0);
  assert(close(aerobic_observed_cost, aerobic_cost));
  assert(close(fermentative_observed_cost, fermentative_cost));
  assert(fermentative_observed_cost > aerobic_observed_cost);
  const double aerobic_acid_per_growth =
      result.acetate_reaction / aerobic_growth;
  const double fermentative_acid_per_growth =
      result.acetate_reaction_agent1 / fermentative_growth;
  assert(std::isfinite(aerobic_acid_per_growth));
  assert(std::isfinite(fermentative_acid_per_growth));
  assert(aerobic_acid_per_growth >= 0.0);
  assert(fermentative_acid_per_growth >= 0.0);
  assert(aerobic_acid_per_growth <= 0.1 * aerobic_cost);
  assert(fermentative_acid_per_growth <= 0.1 * fermentative_cost);
  assert(close(aerobic_acid_per_growth,
               result.fermentation_fraction * 0.1 * aerobic_cost));
  assert(close(fermentative_acid_per_growth,
               result.fermentation_fraction_agent1 * 0.1
                   * fermentative_cost));
  assert(fermentative_acid_per_growth > aerobic_acid_per_growth);
}

void test_spatial_hash() {
  constexpr int agents = 6;
  DeviceBuffer<double> x(agents);
  DeviceBuffer<double> y(agents);
  DeviceBuffer<double> z(agents);
  DeviceBuffer<int> state(agents);
  DeviceBuffer<int> keys(agents);
  DeviceBuffer<int> sorted(agents);
  x.upload(std::vector<double>{0.0, 1.0, 3.999, 2.0, -1.0, 2.0});
  y.upload(std::vector<double>{0.0, 1.0, 3.999, 0.0, 5.0, 0.0});
  z.upload(std::vector<double>{0.0, 1.0, 3.999, 3.0, 4.0, 0.0});
  state.upload(std::vector<int>{0, 0, 0, 0, 0, 3});
  gutibm::gpu::launch_spatial_hash_build_kernel(
      x.data(), y.data(), z.data(), state.data(), keys.data(), sorted.data(),
      agents, 0.0, 0.0, 0.0, 1.0, kNx, kNy, kNz, nullptr);
  synchronize();
  const auto result = download(keys, agents);
  assert(result[0] == 0);
  assert(result[1] == 1 + kNx + kNx * kNy);
  assert(result[2] == 3 + 3 * kNx + 3 * kNx * kNy);
  assert(result[3] == 2 + 3 * kNx * kNy);
  assert(result[4] == 3 * kNx + 3 * kNx * kNy);
  // CSR offsets are built outside this wrapper; this test covers only keys.
  assert(result[5] == -1);
  const auto order = download(sorted, agents);
  for (int agent = 0; agent < agents; ++agent) assert(order[agent] == agent);
}

MechanicsLaunchParams mechanics_parameters() {
  MechanicsLaunchParams parameters{};
  parameters.hertzian_enabled = 1;
  parameters.hertz_k = 1.0;
  parameters.viscosity = 1.0;
  parameters.dt = 1.0;
  parameters.max_displacement_fraction = 0.1;
  parameters.cell_size = 1.0;
  parameters.nx_cells = 1;
  parameters.ny_cells = 1;
  parameters.nz_cells = 1;
  parameters.lo0 = 0.0;
  parameters.lo1 = 0.0;
  parameters.lo2 = 0.0;
  parameters.hi0 = 4.0;
  parameters.hi1 = 4.0;
  parameters.hi2 = 4.0;
  return parameters;
}

void test_mechanics_clear() {
  DeviceBuffer<double> dx(2);
  DeviceBuffer<double> dy(2);
  DeviceBuffer<double> dz(2);
  DeviceBuffer<int> clamp(1);
  dx.upload(std::vector<double>{1.0, 2.0});
  dy.upload(std::vector<double>{1.0, 2.0});
  dz.upload(std::vector<double>{1.0, 2.0});
  clamp.upload(std::vector<int>{4});
  gutibm::gpu::launch_mechanics_clear_kernel(
      dx.data(), dy.data(), dz.data(), clamp.data(), 2, nullptr);
  synchronize();
  for (const double value : download(dx, 2)) assert(value == 0.0);
  for (const double value : download(dy, 2)) assert(value == 0.0);
  for (const double value : download(dz, 2)) assert(value == 0.0);
  assert(download(clamp, 1)[0] == 0);
}

void test_mechanics_forces() {
  DeviceBuffer<double> x(2);
  DeviceBuffer<double> y(2);
  DeviceBuffer<double> z(2);
  DeviceBuffer<double> radius(2);
  DeviceBuffer<double> dx(2);
  DeviceBuffer<double> dy(2);
  DeviceBuffer<double> dz(2);
  DeviceBuffer<int> state(2);
  DeviceBuffer<int> offsets(2);
  DeviceBuffer<int> sorted(2);
  x.upload(std::vector<double>{0.5, 1.5});
  y.upload(std::vector<double>{0.5, 0.5});
  z.upload(std::vector<double>{0.5, 0.5});
  radius.upload(std::vector<double>{0.5, 0.5});
  state.upload(std::vector<int>{0, 0});
  offsets.upload(std::vector<int>{0, 2});
  sorted.upload(std::vector<int>{0, 1});
  dx.upload(std::vector<double>(2, 0.0));
  dy.upload(std::vector<double>(2, 0.0));
  dz.upload(std::vector<double>(2, 0.0));
  const MechanicsLaunchParams parameters = mechanics_parameters();
  gutibm::gpu::launch_mechanics_forces_kernel(
      x.data(), y.data(), z.data(), radius.data(), state.data(), offsets.data(),
      sorted.data(), dx.data(), dy.data(), dz.data(), 2, parameters, nullptr);
  synchronize();
  for (const double value : download(dx, 2)) assert(value == 0.0);

  x.upload(std::vector<double>{0.5, 1.500001});
  dx.upload(std::vector<double>(2, 0.0));
  gutibm::gpu::launch_mechanics_forces_kernel(
      x.data(), y.data(), z.data(), radius.data(), state.data(), offsets.data(),
      sorted.data(), dx.data(), dy.data(), dz.data(), 2, parameters, nullptr);
  synchronize();
  for (const double value : download(dx, 2)) assert(value == 0.0);

  x.upload(std::vector<double>{0.5, 0.8});
  gutibm::gpu::launch_mechanics_forces_kernel(
      x.data(), y.data(), z.data(), radius.data(), state.data(), offsets.data(),
      sorted.data(), dx.data(), dy.data(), dz.data(), 2, parameters, nullptr);
  synchronize();
  const auto force = download(dx, 2);
  assert(force[0] < 0.0);
  assert(force[1] > 0.0);
}

void test_mechanics_apply() {
  DeviceBuffer<double> x(2);
  DeviceBuffer<double> y(2);
  DeviceBuffer<double> z(2);
  DeviceBuffer<double> radius(2);
  DeviceBuffer<double> dx(2);
  DeviceBuffer<double> dy(2);
  DeviceBuffer<double> dz(2);
  DeviceBuffer<int> clamp(1);
  x.upload(std::vector<double>{0.5, 0.8});
  y.upload(std::vector<double>{0.5, 0.5});
  z.upload(std::vector<double>{0.5, 0.5});
  radius.upload(std::vector<double>{0.5, 0.5});
  dx.upload(std::vector<double>{-1.0, 1.0});
  dy.upload(std::vector<double>(2, 0.0));
  dz.upload(std::vector<double>(2, 0.0));
  clamp.upload(std::vector<int>{0});
  const MechanicsLaunchParams parameters = mechanics_parameters();
  gutibm::gpu::launch_mechanics_apply_kernel(
      x.data(), y.data(), z.data(), radius.data(), dx.data(), dy.data(),
      dz.data(), clamp.data(), 2, parameters, nullptr);
  synchronize();
  const auto moved = download(x, 2);
  assert(close(moved[0], 0.45));
  assert(close(moved[1], 0.85));
  assert(download(clamp, 1)[0] == 2);
}

void test_receptor_kill_probability() {
  DeviceBuffer<int> cell(1);
  DeviceBuffer<int> state(1);
  DeviceBuffer<double> receptor(8);
  DeviceBuffer<double> ligand(8);
  DeviceBuffer<double> toxin_affinity(8);
  DeviceBuffer<double> immunity(4);
  DeviceBuffer<double> toxin(1);
  DeviceBuffer<double> competitor(1);
  DeviceBuffer<double> zero(1);
  DeviceBuffer<double> kill(1);
  cell.upload(std::vector<int>{0});
  state.upload(std::vector<int>{0});
  receptor.upload(std::vector<double>(8, 1.0));
  ligand.upload(std::vector<double>(8, 1.0));
  toxin_affinity.upload(std::vector<double>(8, 1.0));
  immunity.upload(std::vector<double>(4, 1.0));
  zero.upload(std::vector<double>{0.0});
  competitor.upload(std::vector<double>{0.0});
  toxin.upload(std::vector<double>{1.0});
  gutibm::gpu::launch_receptor_kill_prob_kernel(
      cell.data(), state.data(), receptor.data(), ligand.data(),
      toxin_affinity.data(), immunity.data(), toxin.data(), zero.data(),
      zero.data(), zero.data(), competitor.data(), zero.data(), zero.data(),
      kill.data(), 1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, kNx, kNy, kNx, 0, kNx, 0, nullptr);
  synchronize();
  const double baseline = download(kill, 1)[0];
  assert(baseline > 0.0 && baseline < 1.0);

  competitor.upload(std::vector<double>{10.0});
  gutibm::gpu::launch_receptor_kill_prob_kernel(
      cell.data(), state.data(), receptor.data(), ligand.data(),
      toxin_affinity.data(), immunity.data(), toxin.data(), zero.data(),
      zero.data(), zero.data(), competitor.data(), zero.data(), zero.data(),
      kill.data(), 1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, kNx, kNy, kNx, 0, kNx, 0, nullptr);
  synchronize();
  assert(download(kill, 1)[0] < baseline);

  toxin.upload(std::vector<double>{0.5});
  competitor.upload(std::vector<double>{0.0});
  gutibm::gpu::launch_receptor_kill_prob_kernel(
      cell.data(), state.data(), receptor.data(), ligand.data(),
      toxin_affinity.data(), immunity.data(), toxin.data(), zero.data(),
      zero.data(), zero.data(), competitor.data(), zero.data(), zero.data(),
      kill.data(), 1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, kNx, kNy, kNx, 0, kNx, 0, nullptr);
  synchronize();
  assert(download(kill, 1)[0] < baseline);
}

template <typename Function>
void run_case(const char* name, Function function) {
  function();
  std::cout << "  " << name << ": PASSED\n";
}

#endif

}  // namespace

int main() {
  std::cout << "=== Direct GPU Kernel Unit Tests ===\n";
  const int gpu_status = gutibm::test::require_gpu("gpu_kernel_units");
  if (gpu_status != 0) return gpu_status;
#ifdef GUTIBM_CUDA
  run_case("launch_field_update_kernel", test_field_update);
  run_case("launch_apply_boundaries_kernel", test_apply_boundaries);
  run_case("launch_grid_coupling_kernel", test_grid_coupling);
  run_case("launch_metabolism_kernel", test_metabolism);
  run_case("launch_metabolism_kernel negative mu", test_metabolism_negative_mu);
  run_case("launch_metabolism_kernel uptake limit",
           test_metabolism_uptake_limit);
  run_case("launch_metabolism_kernel carbon maintenance",
           test_metabolism_maintenance);
  run_case("launch_metabolism_kernel carbon maintenance uptake limit",
           test_metabolism_maintenance_uptake_limit);
  run_case("launch_metabolism_kernel Fur", test_metabolism_fur);
  run_case("launch_metabolism_kernel acetate", test_metabolism_acetate);
  run_case("launch_metabolism_kernel metabolic mode",
           test_metabolism_metabolic_mode);
  run_case("launch_diffuse_x_periodic", test_diffuse_x_periodic);
  run_case("launch_diffuse_y_periodic", test_diffuse_y_periodic);
  run_case("launch_diffuse_z_bounded", test_diffuse_z_bounded);
  run_case("launch_set_epithelial_boundary", test_set_epithelial_boundary);
  run_case("launch_set_luminal_neumann", test_set_luminal_neumann);
  run_case("launch_shift_z_gradient", test_shift_z_gradient);
  run_case("launch_clamp_nonneg", test_clamp_nonneg);
  run_case("launch_superpose_kernel", test_superpose);
  run_case("launch_fmm_far_local_kernel", test_fmm_far_local);
  run_case("launch_vbf_coupling_kernel", test_vbf_coupling);
  run_case("launch_vbf_coupling_kernel implicit carbon sink",
           test_vbf_implicit_carbon_sink);
  run_case("launch_o2_depletion_kernel", test_o2_depletion);
  run_case("launch_o2_depletion_kernel metabolic mode",
           test_o2_depletion_metabolic_mode);
  run_case("launch_spatial_hash_build_kernel", test_spatial_hash);
  run_case("launch_mechanics_clear_kernel", test_mechanics_clear);
  run_case("launch_mechanics_forces_kernel", test_mechanics_forces);
  run_case("launch_mechanics_apply_kernel", test_mechanics_apply);
  run_case("launch_receptor_kill_prob_kernel", test_receptor_kill_probability);
  std::cout << "All 22 GPU launch entry points were invoked directly.\n";
  return 0;
#endif
}
