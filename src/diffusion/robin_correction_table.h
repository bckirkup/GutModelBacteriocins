/* -----------------------------------------------------------------------
   GutIBM – Shared Robin lumen-boundary correction table interpolation
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_ROBIN_CORRECTION_TABLE_H
#define GUTIBM_ROBIN_CORRECTION_TABLE_H

#include <cstddef>
#include <cstdint>
#include <array>
#include <cmath>
#include <list>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <tuple>
#include <vector>

#ifdef __CUDACC__
#define GUTIBM_ROBIN_HOST_DEVICE __host__ __device__
#else
#define GUTIBM_ROBIN_HOST_DEVICE
#endif

namespace gutibm {

class AdvectionField;

namespace robin {

enum class TransferBasis {
  Effective,
  Free
};

constexpr int kTableNodes = 33;
constexpr int kTableValueCount =
    kTableNodes * kTableNodes * kTableNodes;
constexpr int kTableModeCount = 2048;
constexpr double kTableRelativeTolerance = 1.0e-6;
constexpr double kDefaultCutoff = 200.0e-6;
constexpr double kMinimumTableRho = 0.25e-6;
constexpr double kZeroTransferLength =
    std::numeric_limits<double>::infinity();

GUTIBM_ROBIN_HOST_DEVICE inline bool transfer_enabled(double transfer_length) {
  return transfer_length > 0.0 && std::isfinite(transfer_length);
}

struct TableView {
  const double* values;
  double z_lo;
  double height;
  double cutoff;
  double rho_min = kMinimumTableRho;
};

GUTIBM_ROBIN_HOST_DEVICE inline double clamp_table_coordinate(
    double coordinate) {
  if (coordinate < 0.0) return 0.0;
  if (coordinate > static_cast<double>(kTableNodes - 1)) {
    return static_cast<double>(kTableNodes - 1);
  }
  return coordinate;
}

GUTIBM_ROBIN_HOST_DEVICE inline int table_index(
    int source_index, int target_index, int rho_index) {
  return (source_index * kTableNodes + target_index) * kTableNodes
      + rho_index;
}

GUTIBM_ROBIN_HOST_DEVICE inline double table_value(
    const TableView& table, int source_index, int target_index, int rho_index) {
  return table.values[table_index(source_index, target_index, rho_index)];
}

GUTIBM_ROBIN_HOST_DEVICE inline double interpolate_impl(
    const TableView& table, double source_z, double target_z, double rho,
    bool logarithmic_rho) {
  const double source_coordinate =
      (source_z - table.z_lo) / table.height
      * static_cast<double>(kTableNodes - 1);
  const double target_coordinate =
      (target_z - table.z_lo) / table.height
      * static_cast<double>(kTableNodes - 1);
  const double rho_coordinate = logarithmic_rho
      ? (rho <= table.rho_min
          ? 0.0
          : std::log(rho / table.rho_min)
              / std::log(table.cutoff / table.rho_min)
              * static_cast<double>(kTableNodes - 1))
      : rho / table.cutoff * static_cast<double>(kTableNodes - 1);
  const double source_clamped = clamp_table_coordinate(source_coordinate);
  const double target_clamped = clamp_table_coordinate(target_coordinate);
  const double rho_clamped = clamp_table_coordinate(rho_coordinate);

  const auto source_low = static_cast<int>(source_clamped);
  const auto target_low = static_cast<int>(target_clamped);
  const auto rho_low = static_cast<int>(rho_clamped);
  const int source_high = source_low == kTableNodes - 1
      ? source_low : source_low + 1;
  const int target_high = target_low == kTableNodes - 1
      ? target_low : target_low + 1;
  const int rho_high = rho_low == kTableNodes - 1
      ? rho_low : rho_low + 1;
  const double source_weight = source_clamped - source_low;
  const double target_weight = target_clamped - target_low;
  const double rho_weight = rho_clamped - rho_low;

  const double c000 = table_value(
      table, source_low, target_low, rho_low);
  const double c001 = table_value(
      table, source_low, target_low, rho_high);
  const double c010 = table_value(
      table, source_low, target_high, rho_low);
  const double c011 = table_value(
      table, source_low, target_high, rho_high);
  const double c100 = table_value(
      table, source_high, target_low, rho_low);
  const double c101 = table_value(
      table, source_high, target_low, rho_high);
  const double c110 = table_value(
      table, source_high, target_high, rho_low);
  const double c111 = table_value(
      table, source_high, target_high, rho_high);

  const double c00 = c000 + rho_weight * (c001 - c000);
  const double c01 = c010 + rho_weight * (c011 - c010);
  const double c10 = c100 + rho_weight * (c101 - c100);
  const double c11 = c110 + rho_weight * (c111 - c110);
  const double c0 = c00 + target_weight * (c01 - c00);
  const double c1 = c10 + target_weight * (c11 - c10);
  return c0 + source_weight * (c1 - c0);
}

GUTIBM_ROBIN_HOST_DEVICE inline double interpolate_log(
    const TableView& table, double source_z, double target_z, double rho) {
  return interpolate_impl(table, source_z, target_z, rho, true);
}

GUTIBM_ROBIN_HOST_DEVICE inline double interpolate_uniform(
    const TableView& table, double source_z, double target_z, double rho) {
  return interpolate_impl(table, source_z, target_z, rho, false);
}

struct Table {
  std::vector<double> legacy_values;
  std::vector<double> physical_values;
  std::vector<double> drift_values;
  double z_lo = 0.0;
  double height = 0.0;
  double cutoff = kDefaultCutoff;
  double biot_number = 0.0;
  double screening_height = 0.0;
  double lower_coefficient_height = 0.0;
  bool drift_correction = false;
  std::array<int64_t, 9> quantized_key{};
  TransferBasis basis = TransferBasis::Effective;
};

struct TableCacheKey {
  std::array<int64_t, 9> groups{};
  bool operator<(const TableCacheKey& other) const {
    return groups < other.groups;
  }
};

struct TableCacheSnapshot {
  uint64_t generation = 0;
  std::vector<std::shared_ptr<const Table>> tables;
};

struct SealedFieldParams {
  double z_source;
  double z_target;
  double rho;
  double z_lo;
  double z_hi;
  double d_eff;
  double decay_rate;
  double flow_x;
  double flow_y;
  double flow_z;
  int mode_count;
};

enum class SealedReference {
  Physical,
  ImageConsistent
};

class TableCache {
 public:
  static constexpr size_t kMaximumTables = 64;

  std::shared_ptr<const Table> get(
      const AdvectionField& adv, double z_lo, double z_hi,
      double d_free, double d_eff, double decay_rate,
      double lumen_transfer_length, double cutoff, TransferBasis basis,
      double image_series_relative_tolerance = 1.0e-10,
      int image_series_max_shells = 512,
      bool image_series_max_shells_explicit = false,
      bool image_series_legacy_reflections = false,
      bool drift_correction = false);
  TableCacheSnapshot snapshot() const;
  size_t size() const;
  uint64_t tables_built() const;
  uint64_t table_evictions() const;
  uint64_t built_identity() const;
  std::string metadata() const;
  std::string built_identity_hash() const;

 private:
  struct Entry {
    TableCacheKey key;
    std::shared_ptr<const Table> table;
  };

  mutable std::mutex mutex_;
  std::list<Entry> lru_;
  std::map<TableCacheKey, std::list<Entry>::iterator> entries_;
  uint64_t generation_ = 0;
  uint64_t tables_built_ = 0;
  uint64_t table_evictions_ = 0;
  uint64_t built_identity_ = 0;
};

TableCache& global_table_cache();

std::string format_identity_hash(uint64_t identity);

double robin_biot_number(double d_free, double d_eff, double height,
                         double lumen_transfer_length, TransferBasis basis);

std::vector<double> robin_mode_roots(double height, double c_lo,
                                     double c_hi, int mode_count);

bool requires_direct_evaluation(double source_z, double target_z, double rho,
                                double z_lo, double z_hi,
                                double cell_radius);

Table build_table(const AdvectionField& adv, double z_lo, double z_hi,
                  double d_free, double d_eff, double decay_rate,
                  double lumen_transfer_length, double cutoff,
                  TransferBasis basis = TransferBasis::Effective,
                  double image_series_relative_tolerance = 1.0e-10,
                  int image_series_max_shells = 512,
                  bool image_series_max_shells_explicit = false,
                  bool image_series_legacy_reflections = false,
                  bool drift_correction = false);

double normalized_robin_field(double z_source, double z_target, double rho,
                              double z_lo, double z_hi, double d_eff,
                              double d_free, double decay_rate,
                              double lumen_transfer_length,
                              double flow_x, double flow_y, double flow_z,
                              int mode_count,
                              TransferBasis basis = TransferBasis::Effective);

double normalized_sealed_field(const SealedFieldParams& params);
double normalized_image_consistent_sealed_field(
    const SealedFieldParams& params);

double normalized_image_series(double z_source, double z_target, double rho,
                               double z_lo, double z_hi, double d_eff,
                               double decay_rate, double flow_x, double flow_z,
                               double image_series_relative_tolerance,
                               int image_series_max_shells,
                               bool image_series_max_shells_explicit,
                               bool image_series_legacy_reflections);

double normalized_correction(double z_source, double z_target, double rho,
                             double z_lo, double z_hi, double d_eff,
                             double d_free, double decay_rate,
                             double lumen_transfer_length,
                             double flow_x, double flow_y, double flow_z,
                             int mode_count,
                             TransferBasis basis,
                             SealedReference sealed_reference);

double interpolate_drift(const Table& table, double source_z,
                         double target_z, double rho);
double interpolate_legacy(const Table& table, double source_z, double target_z,
                          double rho);
double interpolate_physical_correction(const Table& table, double source_z,
                                       double target_z, double rho);

}  // namespace robin
}  // namespace gutibm

#undef GUTIBM_ROBIN_HOST_DEVICE

#endif  // GUTIBM_ROBIN_CORRECTION_TABLE_H
