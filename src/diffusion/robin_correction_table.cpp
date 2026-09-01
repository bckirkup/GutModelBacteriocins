/* -----------------------------------------------------------------------
   GutIBM – Robin lumen-boundary mode expansion and correction tables
   ----------------------------------------------------------------------- */

#include "robin_correction_table.h"

#include "advection.h"
#include "error.h"
#include "neumann_image_series.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <format>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <string>

namespace gutibm::robin {

namespace {

constexpr uint64_t kFnvOffset = 14695981039346656037ULL;
constexpr uint64_t kFnvPrime = 1099511628211ULL;
constexpr double kBisectionTolerance = 1.0e-13;
constexpr double kMinimumRho = 1.0e-12;
constexpr int kSamplesPerPiInterval = 260;
constexpr int kHistoricalLegacyImageSeriesShells = 3;

void append_hash_bytes(uint64_t& hash, const void* data, size_t size) {
  const auto* bytes = static_cast<const unsigned char*>(data);
  for (size_t index = 0; index < size; ++index) {
    hash ^= bytes[index];
    hash *= kFnvPrime;
  }
}

uint64_t table_identity_hash(const Table& table) {
  uint64_t hash = kFnvOffset;
  for (const int64_t group : table.quantized_key) {
    append_hash_bytes(hash, &group, sizeof(group));
  }
  const int basis = static_cast<int>(table.basis);
  append_hash_bytes(hash, &basis, sizeof(basis));
  append_hash_bytes(hash, &table.z_lo, sizeof(table.z_lo));
  append_hash_bytes(hash, &table.height, sizeof(table.height));
  append_hash_bytes(hash, &table.cutoff, sizeof(table.cutoff));
  for (const double value : table.legacy_values) {
    append_hash_bytes(hash, &value, sizeof(value));
  }
  return hash;
}

}  // namespace

std::string format_identity_hash(uint64_t identity) {
  return std::format("{:016x}", identity);
}

namespace {

Vec3 mean_profile_velocity(const AdvectionField& adv, double z) {
  return adv.mean_velocity({0.0, 0.0, z});
}

int64_t quantize_relative(double value) {
  if (value == 0.0) return std::numeric_limits<int64_t>::min();
  const double magnitude = std::abs(value);
  const double bin = std::round(
      std::log(magnitude) / std::log1p(0.01));
  const auto quantized = static_cast<int64_t>(bin);
  return value < 0.0 ? -quantized : quantized;
}

TableCacheKey make_cache_key(const AdvectionField& adv, double z_lo,
                             double z_hi, double d_free, double d_eff,
                             double decay_rate, double transfer_length,
                             double cutoff, TransferBasis basis,
                             double image_series_relative_tolerance,
                             int image_series_max_shells,
                             bool image_series_max_shells_explicit,
                             bool image_series_legacy_reflections,
                             bool drift_correction) {
  const double height = z_hi - z_lo;
  const Vec3 midpoint = {0.0, 0.0, 0.5 * (z_lo + z_hi)};
  const Vec3 flow = mean_profile_velocity(adv, midpoint[2]);
  const double flow_squared = flow[0] * flow[0] + flow[1] * flow[1]
      + flow[2] * flow[2];
  const double screening_height = std::sqrt(
      (decay_rate + flow_squared / (4.0 * d_eff)) / d_eff) * height;
  const double lower_coefficient_height = flow[2] * height / (2.0 * d_eff);
  const double bi = robin_biot_number(
      d_free, d_eff, height, transfer_length, basis);
  return TableCacheKey{{
      quantize_relative(bi),
      quantize_relative(screening_height),
      quantize_relative(lower_coefficient_height),
      quantize_relative(cutoff / height),
      quantize_relative(image_series_relative_tolerance),
      static_cast<int64_t>(image_series_max_shells),
      image_series_max_shells_explicit ? 1 : 0,
      image_series_legacy_reflections ? 1 : 0,
      drift_correction ? 1 : 0}};
}

void set_table_metadata(Table& table, const AdvectionField& adv, double z_lo,
                        double z_hi, double d_free, double d_eff,
                        double decay_rate, double lumen_transfer_length,
                        double cutoff, TransferBasis basis,
                        double image_series_relative_tolerance,
                        int image_series_max_shells,
                        bool image_series_max_shells_explicit,
                        bool image_series_legacy_reflections, bool drift_correction) {
  table.drift_correction = drift_correction;
  table.biot_number = robin_biot_number(
      d_free, d_eff, table.height, lumen_transfer_length, basis);
  const Vec3 midpoint = {0.0, 0.0, 0.5 * (z_lo + z_hi)};
  const Vec3 midpoint_flow = mean_profile_velocity(adv, midpoint[2]);
  const double flow_squared = midpoint_flow[0] * midpoint_flow[0]
      + midpoint_flow[1] * midpoint_flow[1]
      + midpoint_flow[2] * midpoint_flow[2];
  table.screening_height = std::sqrt(
      (decay_rate + flow_squared / (4.0 * d_eff)) / d_eff)
      * table.height;
  table.lower_coefficient_height =
      midpoint_flow[2] * table.height / (2.0 * d_eff);
  const TableCacheKey key = make_cache_key(
      adv, z_lo, z_hi, d_free, d_eff, decay_rate, lumen_transfer_length,
      cutoff, basis, image_series_relative_tolerance, image_series_max_shells,
      image_series_max_shells_explicit, image_series_legacy_reflections,
      table.drift_correction);
  table.quantized_key = key.groups;
  table.basis = basis;
}

double pole_free_residual(double x, double p, double q) {
  return (x * x + p) * std::sin(x) - q * x * std::cos(x);
}

double imaginary_residual(double y, double p, double q) {
  return (p - y * y) * std::sinh(y) - q * y * std::cosh(y);
}

bool is_sign_change(double first, double second) {
  return (first < 0.0 && second > 0.0)
      || (first > 0.0 && second < 0.0);
}

double bisect_pole_free_root(double left, double right, double f_left,
                             double p, double q) {
  for (int iteration = 0; iteration < 200; ++iteration) {
    const double middle = 0.5 * (left + right);
    const double f_middle = pole_free_residual(middle, p, q);
    if (f_middle == 0.0
        || right - left <= kBisectionTolerance
            * std::max(1.0, std::abs(middle))) {
      return middle;
    }
    if (is_sign_change(f_left, f_middle)) {
      right = middle;
    } else {
      left = middle;
      f_left = f_middle;
    }
  }
  return 0.5 * (left + right);
}

void reject_unrepresented_imaginary_mode(double p, double q, double c_lo,
                                         double c_hi, double height) {
  constexpr int kImaginarySamples = 4096;
  double previous = imaginary_residual(1.0e-12, p, q);
  for (int index = 1; index <= kImaginarySamples; ++index) {
    const double y = 80.0 * index / static_cast<double>(kImaginarySamples);
    const double value = imaginary_residual(y, p, q);
    if (is_sign_change(previous, value)) {
      throw SimulationError(std::format(
          "Robin imaginary mode is not represented by the current kernel "
          "(c_lo={:f}, c_hi={:f}, H={:f})",
          c_lo, c_hi, height));
    }
    previous = value;
  }
}

void validate_root(double x, double p, double q, double c_lo, double c_hi,
                   double height) {
  const double scale = (x * x + std::abs(p)) + std::abs(q) * x;
  const double residual = std::abs(pole_free_residual(x, p, q)) / scale;
  if (!(residual <= 1.0e-9)) {
    throw SimulationError(std::format(
        "Robin eigenvalue residual exceeded tolerance (x={:f}, residual={:f}, "
        "c_lo={:f}, c_hi={:f}, H={:f})",
        x, residual, c_lo, c_hi, height));
  }
}

std::vector<double> robin_mode_roots_impl(double height, double c_lo,
                                          double c_hi, int mode_count) {
  if (!(height > 0.0) || mode_count < 0) {
    throw ConfigError("invalid Robin eigenvalue dimensions");
  }
  const double p = c_lo * c_hi * height * height;
  const double q = (c_hi - c_lo) * height;
  reject_unrepresented_imaginary_mode(p, q, c_lo, c_hi, height);

  std::vector<double> roots;
  roots.reserve(static_cast<size_t>(mode_count));
  const double zero_lhs = c_hi * (1.0 - c_lo * height) - c_lo;
  if (const double zero_scale = std::max(
          {std::abs(c_lo), std::abs(c_hi), 1.0 / height});
      std::abs(zero_lhs) <= 1.0e-12 * zero_scale && mode_count > 0) {
    roots.push_back(0.0);
  }
  if (static_cast<int>(roots.size()) == mode_count) return roots;

  const int intervals = mode_count + 2;
  const int sample_count = intervals * kSamplesPerPiInterval;
  const double step = std::numbers::pi
      / static_cast<double>(kSamplesPerPiInterval);
  double left = step * 1.0e-6;
  double f_left = pole_free_residual(left, p, q);
  for (int index = 1; index <= sample_count; ++index) {
    if (static_cast<int>(roots.size()) >= mode_count) break;
    const double right = step * index;
    const double f_right = pole_free_residual(right, p, q);
    if (f_right == 0.0) {
      roots.push_back(right);
    } else if (is_sign_change(f_left, f_right)) {
      roots.push_back(bisect_pole_free_root(
          left, right, f_left, p, q));
    }
    left = right;
    f_left = f_right;
  }
  if (static_cast<int>(roots.size()) != mode_count) {
    throw SimulationError(std::format(
        "Robin eigenvalue enumeration found {} of {} roots (c_lo={:f}, "
        "c_hi={:f}, H={:f})",
        roots.size(), mode_count, c_lo, c_hi, height));
  }
  for (size_t index = 0; index < roots.size(); ++index) {
    if (index > 0 && !(roots[index] > roots[index - 1])) {
      throw SimulationError("Robin eigenvalue enumeration produced duplicate roots");
    }
    if (roots[index] > 0.0) {
      validate_root(roots[index], p, q, c_lo, c_hi, height);
    }
  }
  return roots;
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

double robin_biot_number_impl(double d_free, double d_eff, double height,
                              double lumen_transfer_length,
                              TransferBasis basis) {
  if (!transfer_enabled(lumen_transfer_length) || !(d_eff > 0.0)
      || !(height > 0.0)) {
    return 0.0;
  }
  const double transfer_coefficient = basis == TransferBasis::Free
      ? d_free / lumen_transfer_length : d_eff / lumen_transfer_length;
  return transfer_coefficient * height / d_eff;
}

double robin_boundary_coefficient(double d_free, double d_eff, double height,
                                  double lumen_transfer_length,
                                  TransferBasis basis) {
  return robin_biot_number_impl(
      d_free, d_eff, height, lumen_transfer_length, basis) * d_eff / height;
}

double mode_sum_betas(double z_source, double z_target, double rho,
                      double z_lo, double z_hi, double d_eff,
                      double decay_rate, const std::vector<double>& betas,
                      double flow_x, double flow_y, double flow_z,
                      double eigen_a) {
  const double height = z_hi - z_lo;
  const double mass_rate = decay_rate
      + (flow_x * flow_x + flow_y * flow_y + flow_z * flow_z)
          / (4.0 * d_eff);
  const double lateral_rho = std::max(rho, kMinimumRho);
  double sum = 0.0;
  for (const double beta : betas) {
    const double phi_source =
        eigenfunction(beta, eigen_a, z_source - z_lo);
    const double phi_target =
        eigenfunction(beta, eigen_a, z_target - z_lo);
    const double kappa = std::sqrt(mass_rate / d_eff + beta * beta);
    const double term = phi_source * phi_target
        * std::cyl_bessel_k(0, kappa * lateral_rho)
        / normalization(beta, eigen_a, height);
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
    const std::vector<double>& bessel_values, double eigen_a) {
  const double height = z_hi - z_lo;
  (void)d_eff;
  (void)flow_x;
  (void)flow_y;
  (void)flow_z;
  double sum = 0.0;
  for (size_t mode = 0; mode < betas.size(); ++mode) {
    const double beta = betas[mode];
    const double phi_source =
        eigenfunction(beta, eigen_a, z_source - z_lo);
    const double phi_target =
        eigenfunction(beta, eigen_a, z_target - z_lo);
    const double normalization_value = normalization(
        beta, eigen_a, height);
    sum += phi_source * phi_target
        * bessel_values[static_cast<size_t>(rho_index) * betas.size()
                        + mode]
        / normalization_value;
  }
  // The correction table stores the gauge-transformed field only.
  return 2.0 * sum;
}

double sealed_zero_mode(double z_source, double z_target, double rho,
                        double z_lo, double height, double d_eff,
                        double decay_rate, double flow_x, double flow_y,
                        double flow_z, SealedReference sealed_reference) {
  const double a = flow_z / (2.0 * d_eff);
  const double eigen_a = sealed_reference == SealedReference::Physical
      ? -a : a;
  const double norm = std::abs(eigen_a) < 1.0e-14
      ? height
      : -std::expm1(-2.0 * eigen_a * height) / (2.0 * eigen_a);
  const double phi_source =
      std::exp(-eigen_a * (z_source - z_lo));
  const double phi_target =
      std::exp(-eigen_a * (z_target - z_lo));
  const double mass_rate = decay_rate
      + (flow_x * flow_x + flow_y * flow_y + flow_z * flow_z)
          / (4.0 * d_eff);
  const double kappa_squared = mass_rate / d_eff - a * a;
  if (!(kappa_squared > 0.0)) return 0.0;
  return phi_source * phi_target
      * std::cyl_bessel_k(0, std::sqrt(kappa_squared)
          * std::max(rho, kMinimumRho)) / norm;
}

double mode_sum(double z_source, double z_target, double rho,
                double z_lo, double z_hi, double d_eff, double decay_rate,
                double d_free, double lumen_transfer_length,
                double flow_x, double flow_y, double flow_z,
                bool robin_boundary, int mode_count, TransferBasis basis,
                SealedReference sealed_reference) {
  const double height = z_hi - z_lo;
  const double a = flow_z / (2.0 * d_eff);
  const double bi = robin_boundary
      ? robin_biot_number_impl(
          d_free, d_eff, height, lumen_transfer_length, basis) : 0.0;
  const double b = bi / height + flow_z / (2.0 * d_eff);
  const double lower = robin_boundary
      ? -a
      : (sealed_reference == SealedReference::Physical ? -a : a);
  std::vector<double> betas;
  betas.reserve(static_cast<size_t>(mode_count));
  const bool use_robin_modes = robin_boundary && bi > 0.0;
  if (use_robin_modes) {
    for (const double root : robin_mode_roots_impl(
             height, lower, b, mode_count)) {
      betas.push_back(root / height);
    }
  } else {
    for (int n = 1; n < mode_count; ++n) {
      betas.push_back(static_cast<double>(n) * std::numbers::pi / height);
    }
  }
  const double positive_modes = mode_sum_betas(
      z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate, betas,
      flow_x, flow_y, flow_z,
      use_robin_modes || sealed_reference == SealedReference::Physical
          ? -a : a);
  if (!use_robin_modes) {
    return positive_modes + 2.0 * sealed_zero_mode(
        z_source, z_target, rho, z_lo, height, d_eff, decay_rate,
        flow_x, flow_y, flow_z, sealed_reference);
  }
  return positive_modes;
}

double normalized_image_value(double image_z, double z_target, double rho,
                              double d_eff, double decay_rate,
                              double flow_x, double flow_z, int reflected) {
  const double dz = z_target - image_z;
  const double reflected_flow_z = reflected != 0 ? -flow_z : flow_z;
  const double radius = std::sqrt(rho * rho + dz * dz);
  if (radius < 1.0e-9) return 0.0;
  const double flow_magnitude = std::sqrt(
      flow_x * flow_x + flow_z * flow_z);
  const double screened_speed = decay_rate <= 0.0
      ? flow_magnitude
      : std::sqrt(flow_magnitude * flow_magnitude
          + 4.0 * d_eff * decay_rate);
  const double exponent = (flow_x * rho + reflected_flow_z * dz
      - screened_speed * radius) / (2.0 * d_eff);
  return std::exp(std::max(exponent, -500.0)) / radius;
}

double normalized_image_series_impl(double z_source, double z_target, double rho,
                               double z_lo, double z_hi, double d_eff,
                               double decay_rate, double flow_x,
                               double flow_z,
                               double image_series_relative_tolerance,
                               int image_series_max_shells,
                               bool image_series_max_shells_explicit,
                               bool image_series_legacy_reflections) {
  const double flow_magnitude = std::sqrt(
      flow_x * flow_x + flow_z * flow_z);
  const auto kernel = [=](double image_z, int reflected) {
    return normalized_image_value(
        image_z, z_target, rho, d_eff, decay_rate, flow_x, flow_z, reflected);
  };
  if (image_series_legacy_reflections) {
    const int shell_limit = image_series_max_shells_explicit
        ? image_series_max_shells : kHistoricalLegacyImageSeriesShells;
    double total = kernel(z_source, 0);
    for (int shell = 1; shell <= shell_limit; ++shell) {
      const double offset = 2.0 * static_cast<double>(shell) * (z_hi - z_lo);
      total += kernel(2.0 * z_lo - z_source - offset, 0)
          + kernel(2.0 * z_hi - z_source + offset, 0)
          + kernel(2.0 * z_lo - z_source + offset, 0)
          + kernel(2.0 * z_hi - z_source - offset, 0);
    }
    return total;
  }
  const auto budget = neumann::image_series_budget(
      d_eff, decay_rate, flow_magnitude, std::abs(flow_z),
      z_hi - z_lo, image_series_relative_tolerance);
  const int shell_limit = image_series_max_shells_explicit
      ? image_series_max_shells : budget.max_shells;
  return neumann::sum_image_series(
      z_source, z_lo, z_hi, kernel, image_series_relative_tolerance,
      shell_limit);
}

}  // namespace

double normalized_image_series(double z_source, double z_target, double rho,
                               double z_lo, double z_hi, double d_eff,
                               double decay_rate, double flow_x, double flow_z,
                               double image_series_relative_tolerance,
                               int image_series_max_shells,
                               bool image_series_max_shells_explicit,
                               bool image_series_legacy_reflections) {
  return normalized_image_series_impl(
      z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate, flow_x, flow_z,
      image_series_relative_tolerance, image_series_max_shells,
      image_series_max_shells_explicit, image_series_legacy_reflections);
}

double robin_biot_number(double d_free, double d_eff, double height,
                         double lumen_transfer_length, TransferBasis basis) {
  return robin_biot_number_impl(
      d_free, d_eff, height, lumen_transfer_length, basis);
}

std::vector<double> robin_mode_roots(double height, double c_lo,
                                     double c_hi, int mode_count) {
  std::vector<double> roots = robin_mode_roots_impl(
      height, c_lo, c_hi, mode_count);
  for (double& root : roots) root /= height;
  return roots;
}

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
    double flow_x, double flow_y, double flow_z, int mode_count,
    TransferBasis basis) {
  const double radial = mode_sum(
      z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate,
      d_free, lumen_transfer_length, flow_x, flow_y, flow_z,
      true, mode_count, basis, SealedReference::Physical);
  return radial * std::exp(
      (flow_x * rho + flow_z * (z_target - z_source)) / (2.0 * d_eff));
}

double normalized_sealed_field(const SealedFieldParams& params) {
  const double radial = mode_sum(
      params.z_source, params.z_target, params.rho, params.z_lo, params.z_hi,
      params.d_eff, params.decay_rate, params.d_eff, kZeroTransferLength,
      params.flow_x, params.flow_y, params.flow_z, false, params.mode_count,
      TransferBasis::Effective, SealedReference::Physical);
  return radial * std::exp(
      (params.flow_x * params.rho
       + params.flow_z * (params.z_target - params.z_source))
      / (2.0 * params.d_eff));
}

double normalized_image_consistent_sealed_field(
    const SealedFieldParams& params) {
  const double radial = mode_sum(
      params.z_source, params.z_target, params.rho, params.z_lo, params.z_hi,
      params.d_eff, params.decay_rate, params.d_eff, kZeroTransferLength,
      params.flow_x, params.flow_y, params.flow_z, false, params.mode_count,
      TransferBasis::Effective, SealedReference::ImageConsistent);
  return radial * std::exp(
      (params.flow_x * params.rho
       + params.flow_z * (params.z_target - params.z_source))
      / (2.0 * params.d_eff));
}

double normalized_correction(
    double z_source, double z_target, double rho, double z_lo, double z_hi,
    double d_eff, double d_free, double decay_rate,
    double lumen_transfer_length,
    double flow_x, double flow_y, double flow_z, int mode_count,
    TransferBasis basis, SealedReference sealed_reference) {
  if (!transfer_enabled(lumen_transfer_length)) {
    return 0.0;
  }
  const double robin = mode_sum(
      z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate,
      d_free, lumen_transfer_length, flow_x, flow_y, flow_z,
      true, mode_count, basis, SealedReference::Physical);
  const double sealed = mode_sum(
      z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate,
      d_free, lumen_transfer_length, flow_x, flow_y, flow_z,
      false, mode_count, basis, sealed_reference);
  return robin - sealed;
}

double interpolate_drift(const Table& table, double source_z,
                         double target_z, double rho) {
  const TableView view{
      table.drift_values.data(), table.z_lo, table.height, table.cutoff};
  return interpolate_log(view, source_z, target_z, rho);
}

double interpolate_legacy(const Table& table, double source_z,
                          double target_z, double rho) {
  const TableView view{
      table.legacy_values.data(), table.z_lo, table.height, table.cutoff, 0.0};
  return interpolate_uniform(view, source_z, target_z, rho);
}

double interpolate_physical_correction(const Table& table, double source_z,
                                       double target_z, double rho) {
  const TableView view{
      table.physical_values.data(), table.z_lo, table.height, table.cutoff};
  return interpolate_log(view, source_z, target_z, rho);
}

Table build_legacy_table(const AdvectionField& adv, double z_lo, double z_hi,
                         double d_free, double d_eff, double decay_rate,
                         double lumen_transfer_length, double cutoff,
                         TransferBasis basis,
                         double image_series_relative_tolerance,
                         int image_series_max_shells,
                         bool image_series_max_shells_explicit,
                         bool image_series_legacy_reflections) {
  Table table;
  table.z_lo = z_lo;
  table.height = z_hi - z_lo;
  table.cutoff = cutoff;
  if (!transfer_enabled(lumen_transfer_length)) {
    set_table_metadata(table, adv, z_lo, z_hi, d_free, d_eff, decay_rate,
                       lumen_transfer_length, cutoff, basis,
                       image_series_relative_tolerance, image_series_max_shells,
                       image_series_max_shells_explicit,
                       image_series_legacy_reflections, false);
    return table;
  }
  std::vector<double> legacy_values(kTableValueCount);
  for (int source_index = 0; source_index < kTableNodes; ++source_index) {
    const double z_source = z_lo + table.height * source_index
        / static_cast<double>(kTableNodes - 1);
    const Vec3 flow = mean_profile_velocity(adv, z_source);
    const double a = flow[2] / (2.0 * d_eff);
    const double kc = robin_boundary_coefficient(
        d_free, d_eff, table.height, lumen_transfer_length, basis);
    const double b = kc / d_eff + flow[2] / (2.0 * d_eff);
    const double lower = -a;
    std::vector<double> robin_betas;
    std::vector<double> sealed_betas;
    robin_betas.reserve(kTableModeCount);
    sealed_betas.reserve(kTableModeCount - 1);
    for (const double root : robin_mode_roots_impl(
             table.height, lower, b, kTableModeCount)) {
      robin_betas.push_back(root / table.height);
    }
    for (int n = 1; n < kTableModeCount; ++n) {
      sealed_betas.push_back(
          static_cast<double>(n) * std::numbers::pi / table.height);
    }
    std::vector<double> robin_bessel(
        static_cast<size_t>(kTableNodes) * kTableModeCount);
    std::vector<double> sealed_bessel(
        static_cast<size_t>(kTableNodes) * (kTableModeCount - 1));
    std::vector<double> sealed_zero_bessel(kTableNodes);
    const double mass_rate = decay_rate
        + (flow[0] * flow[0] + flow[1] * flow[1]
           + flow[2] * flow[2]) / (4.0 * d_eff);
    for (int rho_index = 0; rho_index < kTableNodes; ++rho_index) {
      const double rho = cutoff * rho_index
          / static_cast<double>(kTableNodes - 1);
      const double lateral_rho = std::max(rho, kMinimumRho);
      sealed_zero_bessel[rho_index] = std::cyl_bessel_k(
          0, std::sqrt(mass_rate / d_eff - a * a)
              * lateral_rho);
      for (int n = 0; n < kTableModeCount; ++n) {
        const double robin_kappa = std::sqrt(
            mass_rate / d_eff + robin_betas[n] * robin_betas[n]);
        robin_bessel[static_cast<size_t>(rho_index) * kTableModeCount + n] =
            std::cyl_bessel_k(0, robin_kappa * lateral_rho);
      }
      for (int n = 0; n < kTableModeCount - 1; ++n) {
        const double sealed_kappa = std::sqrt(
            mass_rate / d_eff + sealed_betas[n] * sealed_betas[n]);
        sealed_bessel[static_cast<size_t>(rho_index)
            * (kTableModeCount - 1) + n] =
            std::cyl_bessel_k(0, sealed_kappa * lateral_rho);
      }
    }
    for (int target_index = 0; target_index < kTableNodes; ++target_index) {
      const double z_target = z_lo + table.height * target_index
          / static_cast<double>(kTableNodes - 1);
      for (int rho_index = 0; rho_index < kTableNodes; ++rho_index) {
        legacy_values[table_index(source_index, target_index, rho_index)] =
            mode_sum_betas_lookup(
            z_source, z_target, z_lo, z_hi, d_eff,
            robin_betas, flow[0], flow[1], flow[2], rho_index,
            robin_bessel, -a)
            - mode_sum_betas_lookup(
            z_source, z_target, z_lo, z_hi, d_eff,
            sealed_betas, flow[0], flow[1], flow[2], rho_index,
            sealed_bessel, a);
        const double legacy_zero = 2.0 * std::exp(
            -a * (z_source - z_lo) - a * (z_target - z_lo))
            * sealed_zero_bessel[rho_index]
            / (std::abs(a) < 1.0e-14
                ? table.height
                : -std::expm1(-2.0 * a * table.height)
                    / (2.0 * a));
        legacy_values[table_index(source_index, target_index, rho_index)]
            -= legacy_zero;
      }
    }
  }
  set_table_metadata(table, adv, z_lo, z_hi, d_free, d_eff, decay_rate,
                     lumen_transfer_length, cutoff, basis,
                     image_series_relative_tolerance, image_series_max_shells,
                     image_series_max_shells_explicit,
                     image_series_legacy_reflections, false);
  table.legacy_values = std::move(legacy_values);
  return table;
}

Table build_table(const AdvectionField& adv, double z_lo, double z_hi,
                  double d_free, double d_eff, double decay_rate,
                  double lumen_transfer_length, double cutoff,
                  TransferBasis basis,
                  double image_series_relative_tolerance,
                  int image_series_max_shells,
                  bool image_series_max_shells_explicit,
                  bool image_series_legacy_reflections,
                  bool drift_correction) {
  if (!(d_eff > 0.0) || !(z_hi > z_lo)
      || !(cutoff >= kMinimumTableRho)) {
    throw ConfigError("invalid Robin correction-table dimensions");
  }
  if (!drift_correction) {
    return build_legacy_table(
        adv, z_lo, z_hi, d_free, d_eff, decay_rate, lumen_transfer_length,
        cutoff, basis, image_series_relative_tolerance, image_series_max_shells,
        image_series_max_shells_explicit, image_series_legacy_reflections);
  }
  Table table;
  table.z_lo = z_lo;
  table.height = z_hi - z_lo;
  table.cutoff = cutoff;
  table.drift_correction = drift_correction;
  const bool transfer = transfer_enabled(lumen_transfer_length);
  table.drift_values.resize(kTableValueCount);
  if (transfer) {
    table.physical_values.resize(kTableValueCount);
  }
  for (int source_index = 0; source_index < kTableNodes; ++source_index) {
    const double z_source = z_lo + table.height * source_index
        / static_cast<double>(kTableNodes - 1);
    const Vec3 flow = mean_profile_velocity(adv, z_source);
    const double a = flow[2] / (2.0 * d_eff);
    const double kc = robin_boundary_coefficient(
        d_free, d_eff, table.height, lumen_transfer_length, basis);
    const double b = kc / d_eff + flow[2] / (2.0 * d_eff);
    std::vector<double> robin_betas;
    std::vector<double> sealed_betas;
    robin_betas.reserve(kTableModeCount);
    sealed_betas.reserve(kTableModeCount - 1);
    if (transfer) {
      for (const double root : robin_mode_roots_impl(
               table.height, -a, b, kTableModeCount)) {
        robin_betas.push_back(root / table.height);
      }
    }
    for (int n = 1; n < kTableModeCount; ++n) {
      sealed_betas.push_back(
          static_cast<double>(n) * std::numbers::pi / table.height);
    }
    std::vector<double> robin_bessel;
    std::vector<double> sealed_bessel;
    std::vector<double> sealed_zero_bessel;
    robin_bessel.resize(static_cast<size_t>(kTableNodes)
                        * robin_betas.size());
    sealed_bessel.resize(static_cast<size_t>(kTableNodes)
                         * sealed_betas.size());
    sealed_zero_bessel.resize(kTableNodes);
    const double mass_rate = decay_rate
        + (flow[0] * flow[0] + flow[1] * flow[1]
           + flow[2] * flow[2]) / (4.0 * d_eff);
    for (int rho_index = 0; rho_index < kTableNodes; ++rho_index) {
      const double log_rho = kMinimumTableRho * std::pow(
          cutoff / kMinimumTableRho,
          static_cast<double>(rho_index)
              / static_cast<double>(kTableNodes - 1));
      const double lateral_rho = std::max(log_rho, kMinimumRho);
      sealed_zero_bessel[rho_index] = std::cyl_bessel_k(
          0, std::sqrt(mass_rate / d_eff - a * a) * lateral_rho);
      for (size_t n = 0; n < robin_betas.size(); ++n) {
        const double robin_kappa = std::sqrt(
            mass_rate / d_eff + robin_betas[n] * robin_betas[n]);
        robin_bessel[rho_index * robin_betas.size() + n] =
            std::cyl_bessel_k(0, robin_kappa * lateral_rho);
      }
      for (size_t n = 0; n < sealed_betas.size(); ++n) {
        const double sealed_kappa = std::sqrt(
            mass_rate / d_eff + sealed_betas[n] * sealed_betas[n]);
        sealed_bessel[rho_index * sealed_betas.size() + n] =
            std::cyl_bessel_k(0, sealed_kappa * lateral_rho);
      }
    }
    for (int target_index = 0; target_index < kTableNodes; ++target_index) {
      const double z_target = z_lo + table.height * target_index
          / static_cast<double>(kTableNodes - 1);
      for (int rho_index = 0; rho_index < kTableNodes; ++rho_index) {
        const size_t index = static_cast<size_t>(
            table_index(source_index, target_index, rho_index));
        const double rho = kMinimumTableRho * std::pow(
            cutoff / kMinimumTableRho,
            static_cast<double>(rho_index)
                / static_cast<double>(kTableNodes - 1));
        const double robin = robin_betas.empty() ? 0.0
            : mode_sum_betas_lookup(
                z_source, z_target, z_lo, z_hi, d_eff, robin_betas,
                flow[0], flow[1], flow[2], rho_index, robin_bessel, -a);
        const double physical_sealed = mode_sum_betas_lookup(
            z_source, z_target, z_lo, z_hi, d_eff, sealed_betas,
            flow[0], flow[1], flow[2], rho_index, sealed_bessel, -a)
            + 2.0 * std::exp(
                a * (z_source - z_lo) + a * (z_target - z_lo))
                * sealed_zero_bessel[rho_index]
                / (std::abs(a) < 1.0e-14
                    ? table.height
                    : std::expm1(2.0 * a * table.height)
                        / (2.0 * a));
        const double gauge = std::exp(
            (flow[0] * rho + flow[2] * (z_target - z_source))
            / (2.0 * d_eff));
        const double image = normalized_image_series(
            z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate,
            flow[0], flow[2], image_series_relative_tolerance,
            image_series_max_shells, image_series_max_shells_explicit,
            image_series_legacy_reflections) / gauge;
        table.drift_values[index] = std::abs(a) < 1.0e-14
            ? 0.0 : physical_sealed - image;
        if (transfer) {
          table.physical_values[index] = robin - physical_sealed;
        }
      }
    }
  }
  set_table_metadata(table, adv, z_lo, z_hi, d_free, d_eff, decay_rate,
                     lumen_transfer_length, cutoff, basis,
                     image_series_relative_tolerance, image_series_max_shells,
                     image_series_max_shells_explicit,
                     image_series_legacy_reflections, drift_correction);
  return table;
}

std::shared_ptr<const Table> TableCache::get(
    const AdvectionField& adv, double z_lo, double z_hi,
    double d_free, double d_eff, double decay_rate,
    double lumen_transfer_length, double cutoff, TransferBasis basis,
    double image_series_relative_tolerance, int image_series_max_shells,
    bool image_series_max_shells_explicit,
    bool image_series_legacy_reflections, bool drift_correction) {
  const TableCacheKey key = make_cache_key(
      adv, z_lo, z_hi, d_free, d_eff, decay_rate, lumen_transfer_length,
      cutoff, basis, image_series_relative_tolerance, image_series_max_shells,
      image_series_max_shells_explicit, image_series_legacy_reflections,
      drift_correction);
  std::scoped_lock lock(mutex_);
  if (const auto found = entries_.find(key); found != entries_.end()) {
    lru_.splice(lru_.begin(), lru_, found->second);
    return found->second->table;
  }
  auto table = std::make_shared<Table>(build_table(
      adv, z_lo, z_hi, d_free, d_eff, decay_rate, lumen_transfer_length,
      cutoff, basis, image_series_relative_tolerance, image_series_max_shells,
      image_series_max_shells_explicit, image_series_legacy_reflections,
      drift_correction));
  if (lru_.size() >= kMaximumTables) {
    entries_.erase(lru_.back().key);
    lru_.pop_back();
    ++table_evictions_;
  }
  lru_.emplace_front(key, table);
  entries_[key] = lru_.begin();
  built_identity_ ^= table_identity_hash(*table);
  ++tables_built_;
  ++generation_;
  return table;
}

TableCacheSnapshot TableCache::snapshot() const {
  std::scoped_lock lock(mutex_);
  TableCacheSnapshot result;
  result.generation = generation_;
  result.tables.reserve(lru_.size());
  for (const auto& [key, iterator] : entries_) {
    (void)key;
    result.tables.push_back(iterator->table);
  }
  return result;
}

size_t TableCache::size() const {
  std::scoped_lock lock(mutex_);
  return lru_.size();
}

uint64_t TableCache::tables_built() const {
  std::scoped_lock lock(mutex_);
  return tables_built_;
}

uint64_t TableCache::table_evictions() const {
  std::scoped_lock lock(mutex_);
  return table_evictions_;
}

uint64_t TableCache::built_identity() const {
  std::scoped_lock lock(mutex_);
  return built_identity_;
}

std::string TableCache::metadata() const {
  std::scoped_lock lock(mutex_);
  std::ostringstream output;
  output << "capacity=" << kMaximumTables
         << ";tables_built=" << tables_built_
         << ";table_evictions=" << table_evictions_;
  for (const auto& [key, iterator] : entries_) {
    const Table& table = *iterator->table;
    output << ";table{Bi=" << std::format("{:.17g}", table.biot_number)
           << ",kH=" << std::format("{:.17g}", table.screening_height)
           << ",aH=" << std::format("{:.17g}", table.lower_coefficient_height)
           << ",basis="
           << (table.basis == TransferBasis::Free ? "free" : "effective")
           << ",key=" << key.groups[0] << ',' << key.groups[1] << ','
           << key.groups[2] << ',' << key.groups[3] << '}';
  }
  return output.str();
}

std::string TableCache::built_identity_hash() const {
  std::scoped_lock lock(mutex_);
  return format_identity_hash(built_identity_);
}

TableCache& global_table_cache() {
  static TableCache cache;
  return cache;
}

}  // namespace gutibm::robin
