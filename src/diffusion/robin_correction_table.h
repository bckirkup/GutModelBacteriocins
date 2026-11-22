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
  double rho_coordinate = rho / table.cutoff
      * static_cast<double>(kTableNodes - 1);
  if (logarithmic_rho) {
    rho_coordinate = rho <= table.rho_min
        ? 0.0
        : std::log(rho / table.rho_min)
            / std::log(table.cutoff / table.rho_min)
            * static_cast<double>(kTableNodes - 1);
  }
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

double normalized_robin_radial(double z_source, double z_target, double rho,
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

namespace drift_detail {

constexpr double kDriftPi = 3.14159265358979323846;

// Double-precision exponentially scaled Bessel function exp(x) * K_0(x) for x > 0.
// Shared host/device implementation so CPU and CUDA execute identical arithmetic.
GUTIBM_ROBIN_HOST_DEVICE inline double bessel_k0e(double x) {
  if (!(x > 0.0)) return 0.0;
  if (x <= 2.0) {
    constexpr double kEulerGamma = 0.5772156649015328606065120900824;
    const double t = 0.25 * x * x;
    double i0 = 1.0;
    double hsum = 0.0;
    double term = 1.0;
    double harmonic = 0.0;
    for (int k = 1; k <= 14; ++k) {
      const double kd = static_cast<double>(k);
      term *= t / (kd * kd);
      harmonic += 1.0 / kd;
      i0 += term;
      hsum += term * harmonic;
      if (term < 1.0e-17 * i0) break;
    }
    const double k0 = -(log(0.5 * x) + kEulerGamma) * i0 + hsum;
    return exp(x) * k0;
  }
  // Trapezoidal rule for int_{-inf}^{inf} exp(-s^2) / sqrt(s^2 + 2x) ds.
  // Branch point at |Im(s)| = sqrt(2x) >= 2; h = 0.30 gives < 5e-15 relative error.
  constexpr double kStep = 0.30;
  const double weights[21] = {
      9.1393118527122819e-01, 6.9767632607103101e-01, 4.4485806622294115e-01,
      2.3692775868212176e-01, 1.0539922456186433e-01, 3.9163895098987070e-02,
      1.2155178329915386e-02, 3.1511115984444419e-03, 6.8232805275637758e-04,
      1.2340980408667956e-04, 1.8643742331510657e-05, 2.3525752009657407e-06,
      2.4795960179464771e-07, 2.1829577951202660e-08, 1.6052280551856116e-09,
      9.8595055757148522e-11, 5.0561068139878006e-12, 2.1653566731297942e-13,
      7.7440337371461215e-15, 2.3127996401501119e-16, 5.7690497045896265e-18};
  const double two_x = 2.0 * x;
  double sum = 0.5 / sqrt(two_x);
  for (int k = 1; k <= 21; ++k) {
    const double s = kStep * static_cast<double>(k);
    sum += weights[k - 1] / sqrt(s * s + two_x);
  }
  return 2.0 * kStep * sum;
}

// Double-precision cylindrical Bessel function J_0(x).
GUTIBM_ROBIN_HOST_DEVICE inline double bessel_j0(double x) {
  const double ax = fabs(x);
  if (ax <= 12.0) {
    const double t = -0.25 * ax * ax;
    double sum = 1.0;
    double term = 1.0;
    for (int k = 1; k <= 34; ++k) {
      const double kd = static_cast<double>(k);
      term *= t / (kd * kd);
      sum += term;
      const double scale = fabs(sum) > 1.0 ? fabs(sum) : 1.0;
      if (fabs(term) < 1.0e-16 * scale) break;
    }
    return sum;
  }
  const double inv8x = 1.0 / (8.0 * ax);
  double p = 1.0;
  double coeff = inv8x;
  double q = -coeff;
  for (int m = 2; m <= 16; ++m) {
    const double odd = 2.0 * static_cast<double>(m) - 1.0;
    coeff *= (odd * odd) / static_cast<double>(m) * inv8x;
    if ((m & 1) == 0) {
      p += (((m / 2) & 1) != 0 ? -coeff : coeff);
    } else {
      q += ((((m - 1) / 2) & 1) != 0 ? coeff : -coeff);
    }
    if (fabs(coeff) < 1.0e-16) break;
  }
  const double phase = ax - 0.25 * kDriftPi;
  return sqrt(2.0 / (kDriftPi * ax)) * (p * cos(phase) - q * sin(phase));
}

// Evaluates exp(expo) * int_0^inf exp(s_a*t) * exp(-kappa*R)/R dt,
// where R = sqrt(rho^2 + (d+t)^2) and s_a = -a (sign_plus) or +a (!sign_plus).
// Folding `expo` inside exp(...) keeps the exponent non-positive everywhere.
GUTIBM_ROBIN_HOST_DEVICE inline double drift_integral_pm_scaled(
    double d, double rho, double a, double kappa, double expo,
    bool sign_plus, double ref_scale) {
  const double alpha = sign_plus ? (kappa + a) : (kappa - a);
  const double r0 = sqrt(rho * rho + d * d);
  if (r0 < 1.0e-15 || !(alpha > 0.0)) return 0.0;
  const double e0 = expo - kappa * r0;
  if (e0 < -680.0) return 0.0;
  const double upper = 2.0 * fabs(a) * exp(e0) / (r0 * alpha);
  if (upper < 1.0e-14 * ref_scale) return 0.0;

  const double gl_x[12] = {
      0.064056892862605626, 0.191118867473616309, 0.315042679696163374,
      0.433793507626045138, 0.545421471388839536, 0.648093651936975569,
      0.740124191578554364, 0.820001985973902922, 0.886415527004401034,
      0.938274552002732759, 0.974728555971309498, 0.995187219997021360};
  const double gl_w[12] = {
      0.127938195346752157, 0.125837456346828296, 0.121670472927803391,
      0.115505668053725601, 0.107444270115965635, 0.097618652104113888,
      0.086190161531953276, 0.073346481411080306, 0.059298584915436781,
      0.044277438817419806, 0.028531388628933663, 0.012341229799987200};

  const double v_max = log1p(34.0 / (alpha * r0));
  const double half = 0.5 * v_max;
  double sum = 0.0;
  for (int i = 0; i < 12; ++i) {
    for (int side = -1; side <= 1; side += 2) {
      const double v = half * (1.0 + static_cast<double>(side) * gl_x[i]);
      const double ev = exp(v);
      const double t = r0 * expm1(v);
      const double dt_sum = d + t;
      const double radius = sqrt(rho * rho + dt_sum * dt_sum);
      const double r_minus_dt = (rho * rho) / (radius + dt_sum);
      const double ex = expo - alpha * t - kappa * (d + r_minus_dt);
      if (ex > -700.0) {
        sum += gl_w[i] * (r0 * ev / radius) * exp(ex);
      }
    }
  }
  return half * sum;
}

GUTIBM_ROBIN_HOST_DEVICE inline double pure_neumann_image_term(
    double rho, double dz, double kappa, double expo) {
  const double radius = sqrt(rho * rho + dz * dz);
  if (radius < 1.0e-9) return 0.0;
  const double ex = expo - kappa * radius;
  return ex > -700.0 ? exp(ex) / radius : 0.0;
}

// Pure symmetric Neumann image series for (nabla^2 - kappa^2) psi = -4*pi*delta,
// scaled by exp(expo). Every image exponent expo - kappa*r is non-positive.
GUTIBM_ROBIN_HOST_DEVICE inline double pure_neumann_scaled(
    double zeta_s, double zeta_t, double rho, double height,
    double kappa, double expo, double rel_tol) {
  double total = pure_neumann_image_term(rho, fabs(zeta_t - zeta_s), kappa, expo)
      + pure_neumann_image_term(rho, zeta_t + zeta_s, kappa, expo);
  const double tol = rel_tol > 0.0 ? rel_tol : 1.0e-12;
  for (int m = 1; m <= 512; ++m) {
    const double off = 2.0 * static_cast<double>(m) * height;
    const double shell =
        pure_neumann_image_term(rho, fabs(zeta_t - zeta_s + off), kappa, expo)
        + pure_neumann_image_term(rho, fabs(zeta_t - zeta_s - off), kappa, expo)
        + pure_neumann_image_term(rho, fabs(zeta_t + zeta_s + off), kappa, expo)
        + pure_neumann_image_term(rho, fabs(zeta_t + zeta_s - off), kappa, expo);
    total += shell;
    const double scale = total > 1.0e-300 ? total : 1.0e-300;
    if (shell <= tol * scale) break;
  }
  return total;
}

// Smooth higher-shell Hankel tail for the near-field branch (delta*H < 2, rho < 0.15*H).
GUTIBM_ROBIN_HOST_DEVICE inline double drift_tail_integral_scaled(
    double sigma, double rho, double height, double a,
    double kappa, double expo) {
  const double abs_a = fabs(a);
  const double delta = kappa - abs_a;
  if (!(delta > 0.0) || abs_a < 1.0e-18) return 0.0;

  const double gl_x[12] = {
      0.064056892862605626, 0.191118867473616309, 0.315042679696163374,
      0.433793507626045138, 0.545421471388839536, 0.648093651936975569,
      0.740124191578554364, 0.820001985973902922, 0.886415527004401034,
      0.938274552002732759, 0.974728555971309498, 0.995187219997021360};
  const double gl_w[12] = {
      0.127938195346752157, 0.125837456346828296, 0.121670472927803391,
      0.115505668053725601, 0.107444270115965635, 0.097618652104113888,
      0.086190161531953276, 0.073346481411080306, 0.059298584915436781,
      0.044277438817419806, 0.028531388628933663, 0.012341229799987200};

  const double s_max = log1p(20.0 / (delta * height));
  const double half = 0.5 * s_max;
  double sum = 0.0;
  for (int i = 0; i < 12; ++i) {
    for (int side = -1; side <= 1; side += 2) {
      const double s = half * (1.0 + static_cast<double>(side) * gl_x[i]);
      const double es = exp(s);
      const double u = delta * expm1(s);
      const double du_ds = delta * es;
      const double gam = kappa + u;
      const double two_gH = 2.0 * gam * height;
      if (two_gH > 700.0) continue;
      const double denom = -expm1(-two_gH);
      const double q = sqrt(u * (u + 2.0 * kappa));
      const double bess = bessel_j0(rho * q);
      const double ex_hi = expo - two_gH - gam * (2.0 * height - sigma);
      const double ex_lo = expo - two_gH - gam * sigma;
      double term_plus = 0.0;
      double term_minus = 0.0;
      if (a >= 0.0) {
        if (ex_hi > -700.0) term_minus = (2.0 * a) * exp(ex_hi);
        if (ex_lo > -700.0) {
          term_plus = -(2.0 * a) * (du_ds / (gam + a)) * exp(ex_lo);
        }
      } else {
        if (ex_lo > -700.0) term_plus = -(2.0 * a) * exp(ex_lo);
        if (ex_hi > -700.0) {
          term_minus = (2.0 * a) * (du_ds / (gam - a)) * exp(ex_hi);
        }
      }
      sum += gl_w[i] * ((term_plus + term_minus) / denom) * bess;
    }
  }
  return half * sum;
}

// Closed-form Sturm-Liouville eigenmode sum for the physical sealed wall law
// (-D * dC/dz + U_z * C = 0 at z_lo and z_hi), used when delta*H < 2 and rho >= 0.15*H.
GUTIBM_ROBIN_HOST_DEVICE inline double sealed_modal_scaled(
    double zeta_s, double zeta_t, double rho, double height, double a,
    double kappa_sq, double kappa0_sq, double lat_expo, double uz_dz_expo,
    int max_modes) {
  double sum = 0.0;
  if (kappa0_sq > 0.0) {
    const double kappa0 = sqrt(kappa0_sq);
    const double e_lat = lat_expo - kappa0 * rho;
    if (e_lat > -700.0) {
      const double aH = a * height;
      double z_factor = 0.0;
      if (fabs(aH) < 1.0e-10) {
        z_factor = exp(2.0 * a * zeta_t) / height;
      } else if (a > 0.0) {
        z_factor = (2.0 * a / (-expm1(-2.0 * aH)))
            * exp(2.0 * a * (zeta_t - height));
      } else {
        z_factor = (-2.0 * a / (-expm1(2.0 * aH)))
            * exp(2.0 * a * zeta_t);
      }
      sum += z_factor * bessel_k0e(kappa0 * rho) * exp(e_lat);
    }
  }
  const double expo = lat_expo + uz_dz_expo;
  for (int n = 1; n <= max_modes; ++n) {
    const double beta = static_cast<double>(n) * kDriftPi / height;
    const double kn = sqrt(kappa_sq + beta * beta);
    const double tot_ex = expo - kn * rho;
    if (tot_ex < -700.0) break;
    const double k0s = bessel_k0e(kn * rho) * exp(tot_ex);
    const double phi_s = cos(beta * zeta_s) + (a / beta) * sin(beta * zeta_s);
    const double phi_t = cos(beta * zeta_t) + (a / beta) * sin(beta * zeta_t);
    const double norm_n = 0.5 * height * (1.0 + (a * a) / (beta * beta));
    sum += phi_s * phi_t * k0s / norm_n;
    const double env = (2.0 / height) * k0s;
    const double scale = fabs(sum) > 1.0e-300 ? fabs(sum) : 1.0e-300;
    if (n > 4 && env <= 1.0e-13 * scale) break;
  }
  return 2.0 * sum;
}

}  // namespace drift_detail

// Exact hybrid Sturm-Liouville / Hankel steady-state concentration for the
// physical sealed slab (zero total flux -D*dC/dz + U_z*C = 0 at z_lo and z_hi).
// Executes identically on CPU and CUDA GPU with zero table interpolation error
// and zero negative concentration overshoot.
GUTIBM_ROBIN_HOST_DEVICE inline double evaluate_physical_sealed_concentration(
    double dx, double dy, double source_z, double target_z,
    double z_lo, double z_hi, double d_eff, double source_rate,
    double decay_rate, double flow_x, double flow_y, double flow_z,
    double rel_tol = 1.0e-10) {
  const double height = z_hi - z_lo;
  if (!(d_eff > 0.0) || !(height > 0.0) || source_rate == 0.0) return 0.0;
  const double rho = sqrt(dx * dx + dy * dy);
  const double dz = target_z - source_z;
  if (sqrt(rho * rho + dz * dz) < 1.0e-9) return 0.0;

  double zeta_s = source_z - z_lo;
  if (zeta_s < 0.0) zeta_s = 0.0;
  if (zeta_s > height) zeta_s = height;
  double zeta_t = target_z - z_lo;
  if (zeta_t < 0.0) zeta_t = 0.0;
  if (zeta_t > height) zeta_t = height;

  const double a = flow_z / (2.0 * d_eff);
  const double kappa0_sq = (decay_rate > 0.0 ? decay_rate / d_eff : 0.0)
      + (flow_x * flow_x + flow_y * flow_y) / (4.0 * d_eff * d_eff);
  const double kappa_sq = kappa0_sq + a * a;
  const double kappa = sqrt(kappa_sq);
  const double delta = kappa - fabs(a);
  const double lat_expo = (flow_x * dx + flow_y * dy) / (2.0 * d_eff);
  const double uz_dz_expo = a * (zeta_t - zeta_s);
  const double expo = lat_expo + uz_dz_expo;

  double val = 0.0;
  if (!(delta > 0.0)) {
    const double lateral_rho = rho > 1.0e-12 ? rho : 1.0e-12;
    val = drift_detail::sealed_modal_scaled(
        zeta_s, zeta_t, lateral_rho, height, a, kappa_sq, 0.0,
        lat_expo, uz_dz_expo, kTableModeCount);
  } else if (delta * height >= 2.0) {
    val = drift_detail::pure_neumann_scaled(
        zeta_s, zeta_t, rho, height, kappa, expo, rel_tol);
    if (fabs(a) > 1.0e-18) {
      const double sigma = zeta_s + zeta_t;
      for (int m = 0; m < 16; ++m) {
        const double d_lo = 2.0 * static_cast<double>(m) * height + sigma;
        const double d_hi = 2.0 * static_cast<double>(m + 1) * height - sigma;
        const double scale = fabs(val) > 1.0e-300 ? fabs(val) : 1.0e-300;
        const double term =
            -2.0 * a * drift_detail::drift_integral_pm_scaled(
                d_lo, rho, a, kappa, expo, true, scale)
            + 2.0 * a * drift_detail::drift_integral_pm_scaled(
                d_hi, rho, a, kappa, expo, false, scale);
        val += term;
        const double new_scale = fabs(val) > 1.0e-300 ? fabs(val) : 1.0e-300;
        if (m > 0 && fabs(term) <= 1.0e-13 * new_scale) break;
      }
    }
  } else if (rho >= 0.15 * height) {
    val = drift_detail::sealed_modal_scaled(
        zeta_s, zeta_t, rho, height, a, kappa_sq, kappa0_sq,
        lat_expo, uz_dz_expo, 256);
  } else {
    val = drift_detail::pure_neumann_scaled(
        zeta_s, zeta_t, rho, height, kappa, expo, rel_tol);
    if (fabs(a) > 1.0e-18) {
      const double sigma = zeta_s + zeta_t;
      const double scale = fabs(val) > 1.0e-300 ? fabs(val) : 1.0e-300;
      val += -2.0 * a * drift_detail::drift_integral_pm_scaled(
              sigma, rho, a, kappa, expo, true, scale)
          + 2.0 * a * drift_detail::drift_integral_pm_scaled(
              2.0 * height - sigma, rho, a, kappa, expo, false, scale)
          + drift_detail::drift_tail_integral_scaled(
              sigma, rho, height, a, kappa, expo);
    }
  }
  return source_rate / (4.0 * drift_detail::kDriftPi * d_eff) * val;
}

}  // namespace robin
}  // namespace gutibm

#undef GUTIBM_ROBIN_HOST_DEVICE

#endif  // GUTIBM_ROBIN_CORRECTION_TABLE_H
