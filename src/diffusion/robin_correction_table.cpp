/* -----------------------------------------------------------------------
   GutIBM – Robin lumen-boundary mode expansion and correction tables
   ----------------------------------------------------------------------- */

#include "robin_correction_table.h"

#include "advection.h"
#include "error.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <string>

namespace gutibm::robin {

namespace {

constexpr double kBisectionTolerance = 1.0e-13;
constexpr double kMinimumRho = 1.0e-12;

double mode_root(int branch, double height, double a, double b) {
  const double pi = std::numbers::pi;
  const double lower = static_cast<double>(branch) * pi / height;
  const double upper = (static_cast<double>(branch) * pi + 0.5 * pi)
      / height;
  // Keep both brackets away from tan() poles by a fixed angular margin.
  // At high mode numbers, forming beta*height otherwise loses enough
  // precision to place the nominal upper endpoint across the pole.
  constexpr double kAngularMargin = 1.0e-6;
  double lo = lower + kAngularMargin / height;
  double hi = upper - kAngularMargin / height;
  auto equation = [height, a, b](double beta) {
    const double denominator = beta * beta + a * b;
    return std::tan(beta * height)
        - (b - a) * beta / denominator;
  };
  double f_lo = equation(lo);
  if (!(f_lo < 0.0 && equation(hi) > 0.0)) {
    throw SimulationError(
        "Robin eigenvalue bracket does not contain a root at branch "
        + std::to_string(branch) + " (f_lo=" + std::to_string(f_lo)
        + ", f_hi=" + std::to_string(equation(hi))
        + ", a=" + std::to_string(a) + ", b=" + std::to_string(b) + ")");
  }
  for (int iteration = 0; iteration < 100; ++iteration) {
    const double mid = 0.5 * (lo + hi);
    const double f_mid = equation(mid);
    if (std::abs(hi - lo) <= kBisectionTolerance
        * std::max(1.0, std::abs(mid))) {
      return mid;
    }
    if (f_mid > 0.0) {
      hi = mid;
    } else {
      lo = mid;
      f_lo = f_mid;
    }
  }
  return 0.5 * (lo + hi);
}

double eigenfunction(double beta, double a, double zeta) {
  if (std::abs(beta * zeta) < 1.0e-8) {
    return 1.0 - a * zeta;
  }
  return std::cos(beta * zeta) - (a / beta) * std::sin(beta * zeta);
}

double normalization(double beta, double a, double height) {
  if (std::abs(beta * height) < 1.0e-8) {
    return height - a * height * height + a * a * height * height * height
        / 3.0;
  }
  const double sin_two = std::sin(2.0 * beta * height);
  const double cos_sq = 0.5 * height + sin_two / (4.0 * beta);
  const double sin_sq = 0.5 * height - sin_two / (4.0 * beta);
  const double cos_sin = std::sin(beta * height)
      * std::sin(beta * height) / (2.0 * beta);
  return cos_sq - 2.0 * a / beta * cos_sin
      + a * a / (beta * beta) * sin_sq;
}

double robin_biot_number(double d_free, double d_eff, double height,
                         double lumen_transfer_length) {
  if (!(lumen_transfer_length > 0.0)
      || !std::isfinite(lumen_transfer_length) || !(d_eff > 0.0)
      || !(height > 0.0)) {
    return 0.0;
  }
  // Keep the D_free/D_eff convention in one helper until the state-variable
  // convention for the lumen flux is finalized.
  return (d_free / lumen_transfer_length) * height / d_eff;
}

double robin_boundary_coefficient(double d_free, double d_eff, double height,
                                  double lumen_transfer_length) {
  return robin_biot_number(
      d_free, d_eff, height, lumen_transfer_length) * d_eff / height;
}

double mode_sum_betas(double z_source, double z_target, double rho,
                      double z_lo, double z_hi, double d_eff,
                      double decay_rate, const std::vector<double>& betas,
                      double flow_x, double flow_y, double flow_z) {
  const double height = z_hi - z_lo;
  const double a = flow_z / (2.0 * d_eff);
  const double mass_rate = decay_rate
      + (flow_x * flow_x + flow_y * flow_y + flow_z * flow_z)
          / (4.0 * d_eff);
  const double lateral_rho = std::max(rho, kMinimumRho);
  double sum = 0.0;
  for (const double beta : betas) {
    const double phi_source =
        eigenfunction(beta, a, z_source - z_lo);
    const double phi_target =
        eigenfunction(beta, a, z_target - z_lo);
    const double kappa = std::sqrt(mass_rate / d_eff + beta * beta);
    const double term = phi_source * phi_target
        * std::cyl_bessel_k(0, kappa * lateral_rho)
        / normalization(beta, a, height);
    sum += term;
  }
  // Return the gauge-transformed radial field. The full drift factor is
  // restored by the caller, which knows the signed lateral displacement.
  return 2.0 * sum;
}

double mode_sum_betas_lookup(
    double z_source, double z_target, double z_lo, double z_hi,
    double d_eff, const std::vector<double>& betas,
    double flow_x, double flow_y, double flow_z, int rho_index,
    const std::vector<double>& bessel_values) {
  const double height = z_hi - z_lo;
  const double a = flow_z / (2.0 * d_eff);
  (void)flow_x;
  (void)flow_y;
  double sum = 0.0;
  for (size_t mode = 0; mode < betas.size(); ++mode) {
    const double beta = betas[mode];
    const double phi_source = eigenfunction(beta, a, z_source - z_lo);
    const double phi_target = eigenfunction(beta, a, z_target - z_lo);
    const double normalization_value = normalization(beta, a, height);
    sum += phi_source * phi_target
        * bessel_values[static_cast<size_t>(rho_index) * betas.size()
                        + mode]
        / normalization_value;
  }
  // The correction table stores the gauge-transformed field only.
  return 2.0 * sum;
}

double mode_sum(double z_source, double z_target, double rho,
                double z_lo, double z_hi, double d_eff, double decay_rate,
                double d_free, double lumen_transfer_length,
                double flow_x, double flow_y, double flow_z,
                bool robin_boundary, int mode_count) {
  const double height = z_hi - z_lo;
  const double a = flow_z / (2.0 * d_eff);
  const double bi = robin_boundary
      ? robin_biot_number(d_free, d_eff, height, lumen_transfer_length) : 0.0;
  const double b = bi / height + flow_z / (2.0 * d_eff);
  std::vector<double> betas;
  betas.reserve(mode_count);
  const bool use_robin_modes = robin_boundary && bi > 0.0;
  for (int n = 0; n < mode_count; ++n) {
    betas.push_back(use_robin_modes
        ? mode_root(n, height, a, b)
        : static_cast<double>(n) * std::numbers::pi / height);
  }
  return mode_sum_betas(z_source, z_target, rho, z_lo, z_hi, d_eff,
                        decay_rate, betas, flow_x, flow_y, flow_z);
}

}  // namespace

bool requires_direct_evaluation(double source_z, double target_z, double rho,
                                double z_lo, double z_hi,
                                double cell_radius) {
  if (!(cell_radius > 0.0) || !(z_hi > z_lo)) return false;
  const double source_wall_distance = std::min(
      source_z - z_lo, z_hi - source_z);
  const double target_wall_distance = std::min(
      target_z - z_lo, z_hi - target_z);
  return rho < cell_radius
      && source_wall_distance < cell_radius
      && target_wall_distance < cell_radius;
}

double normalized_robin_field(
    double z_source, double z_target, double rho, double z_lo, double z_hi,
    double d_eff, double d_free, double decay_rate,
    double lumen_transfer_length,
    double flow_x, double flow_y, double flow_z, int mode_count) {
  const double radial = mode_sum(
      z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate,
      d_free, lumen_transfer_length, flow_x, flow_y, flow_z,
      true, mode_count);
  return radial * std::exp(flow_z * (z_target - z_source)
                            / (2.0 * d_eff));
}

double normalized_correction(
    double z_source, double z_target, double rho, double z_lo, double z_hi,
    double d_eff, double d_free, double decay_rate,
    double lumen_transfer_length,
    double flow_x, double flow_y, double flow_z, int mode_count) {
  if (!(lumen_transfer_length > 0.0)
      || !std::isfinite(lumen_transfer_length)) {
    return 0.0;
  }
  const double robin = mode_sum(
      z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate,
      d_free, lumen_transfer_length, flow_x, flow_y, flow_z,
      true, mode_count);
  const double sealed = mode_sum(
      z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate,
      d_free, lumen_transfer_length, flow_x, flow_y, flow_z,
      false, mode_count);
  return robin - sealed;
}

Table build_table(const AdvectionField& adv, double z_lo, double z_hi,
                  double d_free, double d_eff, double decay_rate,
                  double lumen_transfer_length, double cutoff) {
  if (!(d_eff > 0.0) || !(z_hi > z_lo) || !(cutoff > 0.0)) {
    throw ConfigError("invalid Robin correction-table dimensions");
  }
  Table table;
  table.z_lo = z_lo;
  table.height = z_hi - z_lo;
  table.cutoff = cutoff;
  table.values.resize(kTableValueCount);
  for (int source_index = 0; source_index < kTableNodes; ++source_index) {
    const double z_source = z_lo + table.height * source_index
        / static_cast<double>(kTableNodes - 1);
    const Vec3 source = {0.0, 0.0, z_source};
    const Vec3 flow = adv.velocity(source);
    const double a = flow[2] / (2.0 * d_eff);
  const double kc = robin_boundary_coefficient(
      d_free, d_eff, table.height, lumen_transfer_length);
  const double b = kc / d_eff + flow[2] / (2.0 * d_eff);
    std::vector<double> robin_betas;
    std::vector<double> sealed_betas;
    robin_betas.reserve(kTableModeCount);
    sealed_betas.reserve(kTableModeCount);
    for (int n = 0; n < kTableModeCount; ++n) {
      robin_betas.push_back(mode_root(n, table.height, a, b));
      sealed_betas.push_back(
          static_cast<double>(n) * std::numbers::pi / table.height);
    }
    std::vector<double> robin_bessel(
        static_cast<size_t>(kTableNodes) * kTableModeCount);
    std::vector<double> sealed_bessel(
        static_cast<size_t>(kTableNodes) * kTableModeCount);
    const double mass_rate = decay_rate
        + (flow[0] * flow[0] + flow[1] * flow[1]
           + flow[2] * flow[2]) / (4.0 * d_eff);
    for (int rho_index = 0; rho_index < kTableNodes; ++rho_index) {
      const double rho = cutoff * rho_index
          / static_cast<double>(kTableNodes - 1);
      const double lateral_rho = std::max(rho, kMinimumRho);
      for (int n = 0; n < kTableModeCount; ++n) {
        const double robin_kappa = std::sqrt(
            mass_rate / d_eff + robin_betas[n] * robin_betas[n]);
        const double sealed_kappa = std::sqrt(
            mass_rate / d_eff + sealed_betas[n] * sealed_betas[n]);
        robin_bessel[static_cast<size_t>(rho_index) * kTableModeCount + n] =
            std::cyl_bessel_k(0, robin_kappa * lateral_rho);
        sealed_bessel[static_cast<size_t>(rho_index) * kTableModeCount + n] =
            std::cyl_bessel_k(0, sealed_kappa * lateral_rho);
      }
    }
    for (int target_index = 0; target_index < kTableNodes; ++target_index) {
      const double z_target = z_lo + table.height * target_index
          / static_cast<double>(kTableNodes - 1);
      for (int rho_index = 0; rho_index < kTableNodes; ++rho_index) {
        table.values[table_index(source_index, target_index, rho_index)] =
            mode_sum_betas_lookup(
                z_source, z_target, z_lo, z_hi, d_eff,
                // Bessel values were precomputed for this source-z row.
                robin_betas, flow[0], flow[1], flow[2], rho_index,
                robin_bessel)
            - mode_sum_betas_lookup(
                z_source, z_target, z_lo, z_hi, d_eff,
                sealed_betas, flow[0], flow[1], flow[2], rho_index,
                sealed_bessel);
      }
    }
  }
  return table;
}

}  // namespace gutibm::robin
