/* -----------------------------------------------------------------------
   GutIBM – Robin lumen-boundary mode expansion and correction tables
   ----------------------------------------------------------------------- */

#include "robin_correction_table.h"

#include "advection.h"
#include "error.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <string>

namespace gutibm::robin {

namespace {

constexpr double kBisectionTolerance = 1.0e-13;
constexpr double kMinimumRho = 1.0e-12;

Vec3 mean_profile_velocity(const AdvectionField& adv, double z) {
  return adv.mean_velocity({0.0, 0.0, z});
}

int64_t quantize_relative(double value) {
  if (value == 0.0) return std::numeric_limits<int64_t>::min();
  const double magnitude = std::abs(value);
  const double bin = std::round(
      std::log(magnitude) / std::log1p(0.01));
  const int64_t quantized = static_cast<int64_t>(bin);
  return value < 0.0 ? -quantized : quantized;
}

TableCacheKey make_cache_key(const AdvectionField& adv, double z_lo,
                             double z_hi, double d_free, double d_eff,
                             double decay_rate, double transfer_length,
                             double cutoff, TransferBasis basis) {
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
      quantize_relative(cutoff / height)}};
}

double mode_root(int branch, double height, double c_lo, double c_hi,
                 double bi, double a_height) {
  const double pi = std::numbers::pi;
  const double branch_lo = static_cast<double>(branch) * pi / height;
  const double branch_hi = (static_cast<double>(branch) * pi + 0.5 * pi)
      / height;
  // Keep both brackets away from tan() poles by a fixed angular margin.
  // At high mode numbers, forming beta*height otherwise loses enough
  // precision to place the nominal upper endpoint across the pole.
  constexpr double kAngularMargin = 1.0e-6;
  double lo = branch_lo + kAngularMargin / height;
  double hi = branch_hi - kAngularMargin / height;
  const double coefficient_product = c_lo * c_hi;
  if (coefficient_product < 0.0) {
    const double rational_pole = std::sqrt(-coefficient_product);
    lo = std::max(lo, rational_pole + kAngularMargin / height);
  }
  auto equation = [height, c_lo, c_hi](double beta) {
    const double denominator = beta * beta + c_lo * c_hi;
    return std::tan(beta * height)
        - (c_hi - c_lo) * beta / denominator;
  };
  double f_lo = equation(lo);
  if (!(f_lo < 0.0 && equation(hi) > 0.0)) {
    throw SimulationError(
        "Robin eigenvalue bracket does not contain a root at branch "
        + std::to_string(branch) + " (f_lo=" + std::to_string(f_lo)
        + ", f_hi=" + std::to_string(equation(hi))
        + ", c_lo=" + std::to_string(c_lo)
        + ", c_hi=" + std::to_string(c_hi)
        + ", Bi=" + std::to_string(bi)
        + ", aH=" + std::to_string(a_height) + ")");
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

double robin_biot_number_impl(double d_free, double d_eff, double height,
                              double lumen_transfer_length,
                              TransferBasis basis) {
  if (!(lumen_transfer_length > 0.0)
      || !std::isfinite(lumen_transfer_length) || !(d_eff > 0.0)
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
                        double flow_z) {
  const double a = flow_z / (2.0 * d_eff);
  const double norm = std::abs(a) < 1.0e-14
      ? height : -std::expm1(-2.0 * a * height) / (2.0 * a);
  const double phi_source = std::exp(-a * (z_source - z_lo));
  const double phi_target = std::exp(-a * (z_target - z_lo));
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
                bool robin_boundary, int mode_count, TransferBasis basis) {
  const double height = z_hi - z_lo;
  const double a = flow_z / (2.0 * d_eff);
  const double bi = robin_boundary
      ? robin_biot_number_impl(
          d_free, d_eff, height, lumen_transfer_length, basis) : 0.0;
  const double b = bi / height + flow_z / (2.0 * d_eff);
  const double lower = robin_boundary ? -a : a;
  std::vector<double> betas;
  betas.reserve(static_cast<size_t>(mode_count));
  const bool use_robin_modes = robin_boundary && bi > 0.0;
  const int first_mode = use_robin_modes ? 0 : 1;
  for (int n = first_mode; n < mode_count; ++n) {
    betas.push_back(use_robin_modes
        ? mode_root(n, height, lower, b, bi, a * height)
        : static_cast<double>(n) * std::numbers::pi / height);
  }
  const double positive_modes = mode_sum_betas(
      z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate, betas,
      flow_x, flow_y, flow_z, use_robin_modes ? -a : a);
  if (!use_robin_modes) {
    return positive_modes + 2.0 * sealed_zero_mode(
        z_source, z_target, rho, z_lo, height, d_eff, decay_rate,
        flow_x, flow_y, flow_z);
  }
  return positive_modes;
}

}  // namespace

double robin_biot_number(double d_free, double d_eff, double height,
                         double lumen_transfer_length, TransferBasis basis) {
  return robin_biot_number_impl(
      d_free, d_eff, height, lumen_transfer_length, basis);
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
      true, mode_count, basis);
  return radial * std::exp(
      (flow_x * rho + flow_z * (z_target - z_source)) / (2.0 * d_eff));
}

double normalized_correction(
    double z_source, double z_target, double rho, double z_lo, double z_hi,
    double d_eff, double d_free, double decay_rate,
    double lumen_transfer_length,
    double flow_x, double flow_y, double flow_z, int mode_count,
    TransferBasis basis) {
  if (!(lumen_transfer_length > 0.0)
      || !std::isfinite(lumen_transfer_length)) {
    return 0.0;
  }
  const double robin = mode_sum(
      z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate,
      d_free, lumen_transfer_length, flow_x, flow_y, flow_z,
      true, mode_count, basis);
  const double sealed = mode_sum(
      z_source, z_target, rho, z_lo, z_hi, d_eff, decay_rate,
      d_free, lumen_transfer_length, flow_x, flow_y, flow_z,
      false, mode_count, basis);
  return robin - sealed;
}

Table build_table(const AdvectionField& adv, double z_lo, double z_hi,
                  double d_free, double d_eff, double decay_rate,
                  double lumen_transfer_length, double cutoff,
                  TransferBasis basis) {
  if (!(d_eff > 0.0) || !(z_hi > z_lo) || !(cutoff > 0.0)) {
    throw ConfigError("invalid Robin correction-table dimensions");
  }
  Table table;
  table.z_lo = z_lo;
  table.height = z_hi - z_lo;
  table.cutoff = cutoff;
  table.values.resize(kTableValueCount);
  const auto set_metadata = [&table, &adv, z_lo, z_hi, d_free, d_eff,
                             decay_rate, lumen_transfer_length, cutoff,
                             basis]() {
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
        cutoff, basis);
    table.quantized_key = key.groups;
    table.basis = basis;
  };
  if (!(lumen_transfer_length > 0.0)
      || !std::isfinite(lumen_transfer_length)) {
    set_metadata();
    return table;
  }
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
    for (int n = 0; n < kTableModeCount; ++n) {
      robin_betas.push_back(mode_root(
          n, table.height, lower, b,
          robin_biot_number(d_free, d_eff, table.height,
                            lumen_transfer_length, basis),
          a * table.height));
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
        table.values[table_index(source_index, target_index, rho_index)] =
            mode_sum_betas_lookup(
                z_source, z_target, z_lo, z_hi, d_eff,
                // Bessel values were precomputed for this source-z row.
                robin_betas, flow[0], flow[1], flow[2], rho_index,
                robin_bessel, -a)
            - mode_sum_betas_lookup(
                z_source, z_target, z_lo, z_hi, d_eff,
                sealed_betas, flow[0], flow[1], flow[2], rho_index,
                sealed_bessel, a);
            const double sealed_zero = 2.0 * std::exp(
                -a * (z_source - z_lo) - a * (z_target - z_lo))
                * sealed_zero_bessel[rho_index]
                / (std::abs(a) < 1.0e-14
                    ? table.height
                    : -std::expm1(-2.0 * a * table.height) / (2.0 * a));
            table.values[table_index(source_index, target_index, rho_index)]
                -= sealed_zero;
      }
    }
  }
  set_metadata();
  return table;
}

std::shared_ptr<const Table> TableCache::get(
    const AdvectionField& adv, double z_lo, double z_hi,
    double d_free, double d_eff, double decay_rate,
    double lumen_transfer_length, double cutoff, TransferBasis basis) {
  const TableCacheKey key = make_cache_key(
      adv, z_lo, z_hi, d_free, d_eff, decay_rate, lumen_transfer_length,
      cutoff, basis);
  std::lock_guard<std::mutex> lock(mutex_);
  if (const auto found = entries_.find(key); found != entries_.end()) {
    lru_.splice(lru_.begin(), lru_, found->second);
    return found->second->table;
  }
  auto table = std::make_shared<Table>(build_table(
      adv, z_lo, z_hi, d_free, d_eff, decay_rate, lumen_transfer_length,
      cutoff, basis));
  if (lru_.size() >= kMaximumTables) {
    entries_.erase(lru_.back().key);
    lru_.pop_back();
    ++table_evictions_;
  }
  lru_.push_front(Entry{key, table});
  entries_[key] = lru_.begin();
  ++tables_built_;
  ++generation_;
  return table;
}

TableCacheSnapshot TableCache::snapshot() const {
  std::lock_guard<std::mutex> lock(mutex_);
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
  std::lock_guard<std::mutex> lock(mutex_);
  return lru_.size();
}

uint64_t TableCache::tables_built() const {
  std::lock_guard<std::mutex> lock(mutex_);
  return tables_built_;
}

uint64_t TableCache::table_evictions() const {
  std::lock_guard<std::mutex> lock(mutex_);
  return table_evictions_;
}

std::string TableCache::metadata() const {
  std::lock_guard<std::mutex> lock(mutex_);
  std::ostringstream output;
  output << "capacity=" << kMaximumTables
         << ";tables_built=" << tables_built_
         << ";table_evictions=" << table_evictions_;
  for (const auto& [key, iterator] : entries_) {
    const Table& table = *iterator->table;
    output << ";table{Bi=" << std::setprecision(17) << table.biot_number
           << ",kH=" << table.screening_height
           << ",aH=" << table.lower_coefficient_height
           << ",basis="
           << (table.basis == TransferBasis::Free ? "free" : "effective")
           << ",key=" << key.groups[0] << ',' << key.groups[1] << ','
           << key.groups[2] << ',' << key.groups[3] << '}';
  }
  return output.str();
}

std::string TableCache::values_hash() const {
  std::lock_guard<std::mutex> lock(mutex_);
  uint64_t hash = 14695981039346656037ULL;
  for (const auto& [key, iterator] : entries_) {
    (void)key;
    for (const double value : iterator->table->values) {
      const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
      for (size_t byte = 0; byte < sizeof(value); ++byte) {
        hash ^= bytes[byte];
        hash *= 1099511628211ULL;
      }
    }
  }
  std::ostringstream output;
  output << std::hex << std::setfill('0') << std::setw(16) << hash;
  return output.str();
}

TableCache& global_table_cache() {
  static TableCache cache;
  return cache;
}

}  // namespace gutibm::robin
