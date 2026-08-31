/* -----------------------------------------------------------------------
   GutIBM – Chemical field implementation
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "domain.h"
#include "delivery_support.h"
#include "error.h"
#include "species_names.h"
#include "tridiagonal_factorization.h"
#include "chemical_field_gpu.h"
#include "diffusion_gpu.h"
#include "dispatch.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <ranges>
#include <string>
#include <vector>

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {

namespace {

Real z_gradient_reference(const ChemicalSpec& spec,
                          const Domain& domain,
                          Int iz);
Real owned_content(const std::vector<Real>& concentration,
                   const Domain& domain, Real cell_volume);
Real owned_content_slab(const std::vector<Real>& concentration,
                        const Domain& domain, Int storage_nx,
                        Int halo_width, Real cell_volume);

bool nutrient_debug_enabled() {
  static const bool enabled = [] {
    const char* value = std::getenv("GUTIBM_DEBUG_NUTRIENT_LEDGER");
    return value != nullptr && value[0] != '\0' && value[0] != '0';
  }();
  return enabled;
}

bool delivery_axis_cap_enabled() {
  static const bool enabled = [] {
    const char* value = std::getenv("GUTIBM_DELIVERY_AXIS_CAP");
    return value != nullptr && value[0] != '\0' && value[0] != '0';
  }();
  return enabled;
}

Int nutrient_debug_transport_step() {
  static const Int step = [] {
    const char* value =
        std::getenv("GUTIBM_DEBUG_NUTRIENT_TRANSPORT_STEP");
    return value != nullptr && value[0] != '\0' ? std::atoi(value) : 1;
  }();
  return step;
}

Int& nutrient_debug_step_counter() {
  static Int step = 0;
  return step;
}

std::string debug_real(Real value) {
  std::array<char, 64> buffer{};
  const auto result = std::to_chars(
      buffer.data(), buffer.data() + buffer.size(), value,
      std::chars_format::general, 17);
  return result.ec == std::errc{}
      ? std::string(buffer.data(), result.ptr) : std::string{};
}

void solve_tridiagonal_with_diagonal(
    std::vector<Real>& values, const std::vector<Real>& diagonal,
    Real alpha, Real first_source = 0.0) {
  if (values.empty()) return;
  std::vector<Real> diag = diagonal;
  std::vector upper(diag.size() > 1 ? diag.size() - 1 : 0, -alpha);
  values.front() += first_source;
  for (size_t i = 1; i < diag.size(); ++i) {
    const Real multiplier = -alpha / diag[i - 1];
    diag[i] -= multiplier * upper[i - 1];
    values[i] -= multiplier * values[i - 1];
  }
  values.back() /= diag.back();
  for (size_t i = values.size() - 1; i > 0; --i) {
    values[i - 1] =
        (values[i - 1] - upper[i - 1] * values[i]) / diag[i - 1];
  }
}

void solve_periodic_with_sink(
    std::vector<Real>& values, const std::vector<Real>& sink, Real alpha,
    const std::vector<Real>* gradient_profile = nullptr,
    const std::vector<Real>* prescribed = nullptr) {
  const size_t n = values.size();
  if (n == 0) return;
  if (gradient_profile != nullptr) {
    for (size_t i = 0; i < n; ++i) {
      values[i] -= sink[i] * (*gradient_profile)[i];
    }
  }
  if (prescribed != nullptr) {
    for (size_t i = 0; i < n; ++i) values[i] -= (*prescribed)[i];
  }
  if (n == 1) {
    values[0] /= 1.0 + sink[0];
    return;
  }
  if (n == 2) {
    const Real d0 = 1.0 + 2.0 * alpha + sink[0];
    const Real d1 = 1.0 + 2.0 * alpha + sink[1];
    const Real off = 2.0 * alpha;
    const Real determinant = d0 * d1 - off * off;
    const Real first = (d1 * values[0] + off * values[1]) / determinant;
    const Real second = (off * values[0] + d0 * values[1]) / determinant;
    values[0] = first;
    values[1] = second;
    return;
  }
  const Real corner = -alpha;
  const Real gamma = -(1.0 + 2.0 * alpha + sink[0]);
  std::vector diagonal(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    diagonal[i] = 1.0 + 2.0 * alpha + sink[i];
  }
  diagonal[0] -= gamma;
  diagonal[n - 1] -= corner * corner / gamma;
  TridiagonalFactorization factorization;
  std::vector lower(n - 1, -alpha);
  std::vector upper(n - 1, -alpha);
  factorization.factorize(lower, diagonal, upper);
  std::vector correction(n, 0.0);
  correction[0] = gamma;
  correction[n - 1] = corner;
  factorization.solve_in_place(correction);
  const Real denominator = 1.0 + correction[0]
      + corner * correction[n - 1] / gamma;
  factorization.solve_in_place(values);
  const Real adjustment =
      (values[0] + corner * values[n - 1] / gamma) / denominator;
  for (size_t i = 0; i < n; ++i) values[i] -= adjustment * correction[i];
}

class PeriodicLineSolver {
 public:
  PeriodicLineSolver(Int size, Real alpha) : size_(size), alpha_(alpha) {
    if (size_ < 3) return;

    const Real diagonal_value = 1.0 + 2.0 * alpha_;
    const Real gamma = -diagonal_value;
    const Real corner = -alpha_;
    std::vector lower(static_cast<size_t>(size_ - 1), -alpha_);
    std::vector upper(static_cast<size_t>(size_ - 1), -alpha_);
    std::vector diagonal(static_cast<size_t>(size_), diagonal_value);
    diagonal.front() -= gamma;
    diagonal.back() -= corner * corner / gamma;
    factorization_.factorize(lower, diagonal, upper);

    correction_.assign(static_cast<size_t>(size_), 0.0);
    correction_.front() = gamma;
    correction_.back() = corner;
    factorization_.solve_in_place(correction_);
    denominator_ = 1.0 + correction_.front()
        + corner * correction_.back() / gamma;
  }

  void solve(std::vector<Real>& values) const {
    if (size_ <= 1) return;
    if (size_ == 2) {
      const Real diagonal = 1.0 + 2.0 * alpha_;
      const Real off_diagonal = 2.0 * alpha_;
      const Real determinant = 1.0 + 4.0 * alpha_;
      const Real first =
          (diagonal * values[0] + off_diagonal * values[1]) / determinant;
      const Real second =
          (off_diagonal * values[0] + diagonal * values[1]) / determinant;
      values[0] = first;
      values[1] = second;
      return;
    }

    const Real diagonal_value = 1.0 + 2.0 * alpha_;
    const Real gamma = -diagonal_value;
    const Real corner = -alpha_;
    factorization_.solve_in_place(values);
    const Real numerator = values.front() + corner * values.back() / gamma;
    const Real adjustment = numerator / denominator_;
    for (size_t i = 0; i < values.size(); ++i) {
      values[i] -= adjustment * correction_[i];
    }
  }

 private:
  Int size_ = 0;
  Real alpha_ = 0.0;
  Real denominator_ = 1.0;
  TridiagonalFactorization factorization_;
  std::vector<Real> correction_;
};

class NeumannTopLineSolver {
 public:
  NeumannTopLineSolver(Int size, Real alpha) : size_(size), alpha_(alpha) {
    if (size_ <= 0) return;
    std::vector lower(static_cast<size_t>(std::max(size_ - 1, 0)), -alpha_);
    std::vector upper(static_cast<size_t>(std::max(size_ - 1, 0)), -alpha_);
    std::vector diagonal(static_cast<size_t>(size_), 1.0 + 2.0 * alpha_);
    diagonal.back() = 1.0 + alpha_;
    factorization_.factorize(lower, diagonal, upper);
  }

  void solve(std::vector<Real>& values, Real boundary_conc) const {
    if (size_ <= 0) return;
    values.front() += alpha_ * boundary_conc;
    factorization_.solve_in_place(values);
  }

 private:
  Int size_ = 0;
  Real alpha_ = 0.0;
  TridiagonalFactorization factorization_;
};

void diffuse_periodic_x(std::vector<Real>& concentration,
                        const Domain& domain,
                        Real alpha) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const PeriodicLineSolver solver(nx, alpha);

  #ifdef GUTIBM_OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<Real> line(static_cast<size_t>(nx));
    #ifdef GUTIBM_OPENMP
    #pragma omp for collapse(2) schedule(static)
    #endif
    for (Int iz = 0; iz < nz; ++iz) {
      for (Int iy = 0; iy < ny; ++iy) {
        for (Int ix = 0; ix < nx; ++ix) {
          line[static_cast<size_t>(ix)] =
              concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))];
        }
        solver.solve(line);
        for (Int ix = 0; ix < nx; ++ix) {
          concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
              line[static_cast<size_t>(ix)];
        }
      }
    }
  }
}

void diffuse_periodic_y(std::vector<Real>& concentration,
                        const Domain& domain,
                        Real alpha) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const PeriodicLineSolver solver(ny, alpha);

  #ifdef GUTIBM_OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<Real> line(static_cast<size_t>(ny));
    #ifdef GUTIBM_OPENMP
    #pragma omp for collapse(2) schedule(static)
    #endif
    for (Int iz = 0; iz < nz; ++iz) {
      for (Int ix = 0; ix < nx; ++ix) {
        for (Int iy = 0; iy < ny; ++iy) {
          line[static_cast<size_t>(iy)] =
              concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))];
        }
        solver.solve(line);
        for (Int iy = 0; iy < ny; ++iy) {
          concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
              line[static_cast<size_t>(iy)];
        }
      }
    }
  }
}

struct DeliverySinkParameters {
  const std::vector<Real>& sink_rate;
  std::vector<Real>& realized;
  Real sink_dt = 0.0;
  Real cell_volume = 0.0;
  const ChemicalSpec* gradient_spec = nullptr;
  const std::vector<Real>* prescribed_mass = nullptr;
  NutrientFluxAccounting* flux = nullptr;
  Int spec = -1;
};

void record_negative_delivery_excursion(
    NutrientFluxAccounting* flux, Int spec, Real concentration,
    Real sink_amount, Real cell_volume, Int axis) {
  if (flux == nullptr || concentration >= 0.0 || sink_amount <= 0.0) {
    return;
  }
  flux->add_negative_delivery_excursion(
      spec, 1.0, -sink_amount * concentration * cell_volume);
  flux->add_negative_delivery_axis(spec, axis, 1.0);
  flux->add_negative_delivery_min(spec, concentration);
}

void record_negative_delivery_excursion(
    NutrientFluxAccounting* flux, Int spec, Real created_mass, Int axis) {
  if (flux == nullptr || created_mass <= 0.0) {
    return;
  }
  flux->add_negative_delivery_excursion(spec, 1.0, created_mass);
  flux->add_negative_delivery_axis(spec, axis, 1.0);
}

template <typename OwnsCell, typename StorageCell>
bool has_negative_owned_cell(
    const std::vector<Real>& concentration, Int global_ncells,
    OwnsCell owns_cell, StorageCell storage_cell) {
  for (Int cell = 0; cell < global_ncells; ++cell) {
    if (!owns_cell(cell)) continue;
    const Int storage = storage_cell(cell);
    if (storage >= 0
        && concentration[static_cast<size_t>(storage)] < 0.0) {
      return true;
    }
  }
  return false;
}

template <typename OwnsCell, typename StorageCell>
Real owned_negative_fraction(
    const std::vector<Real>& concentration, Int global_ncells,
    OwnsCell owns_cell, StorageCell storage_cell) {
  Real negative = 0.0;
  Real owned = 0.0;
  for (Int cell = 0; cell < global_ncells; ++cell) {
    if (!owns_cell(cell)) continue;
    const Int storage = storage_cell(cell);
    if (storage < 0) continue;
    owned += 1.0;
    negative += concentration[static_cast<size_t>(storage)] < 0.0 ? 1.0 : 0.0;
  }
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized && !finalized) {
    std::array<Real, 2> totals = {negative, owned};
    MPI_Allreduce(MPI_IN_PLACE, totals.data(), 2, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    negative = totals[0];
    owned = totals[1];
  }
#endif
  return owned > 0.0 ? negative / owned : 0.0;
}

template <typename OwnsCell, typename StorageCell>
void mark_affected_storage(
    Int affected_cell, std::vector<Real>& prescribed,
    std::vector<char>& affected, std::vector<Int>& affected_cells,
    OwnsCell owns_cell, StorageCell storage_cell) {
  if (!owns_cell(affected_cell)) return;
  const Int affected_storage = storage_cell(affected_cell);
  if (affected_storage < 0
      || static_cast<size_t>(affected_storage) >= prescribed.size()) {
    return;
  }
  const auto index = static_cast<size_t>(affected_storage);
  if (affected[index] != 0) return;
  affected[index] = 1;
  affected_cells.push_back(affected_storage);
}

template <typename OwnsCell, typename StorageCell>
void mark_delivery_support_cells(
    Int cell, const Domain& domain, const DeliverySupportStencil& stencil,
    std::vector<Real>& prescribed, std::vector<char>& affected,
    std::vector<Int>& affected_cells, OwnsCell owns_cell,
    StorageCell storage_cell) {
  const Int ix = cell % domain.nx();
  const Int yz = cell / domain.nx();
  const Int iz = yz / domain.ny();
  const Int iy = yz % domain.ny();
  if (stencil.radius <= 0.0) {
    mark_affected_storage(
        cell, prescribed, affected, affected_cells, owns_cell, storage_cell);
    return;
  }
  for (const auto& offset : stencil.offsets) {
    if (!offset.within_radius) continue;
    Int affected_cell = -1;
    if (delivery_support_target(
            domain, ix, iy, iz, offset, affected_cell)) {
      mark_affected_storage(
          affected_cell, prescribed, affected, affected_cells,
          owns_cell, storage_cell);
    }
  }
}

template <typename OwnsCell, typename StorageCell>
Real reduce_prescribed_near_negative_cells(
    const std::vector<Real>& concentration,
    std::vector<Real>& prescribed, const Domain& domain,
    const DeliverySupportStencil& stencil, std::vector<char>& affected,
    std::vector<Int>& affected_cells, Int global_ncells,
    OwnsCell owns_cell, StorageCell storage_cell) {
  Real reduced = 0.0;
  for (const Int affected_storage : affected_cells) {
    affected[static_cast<size_t>(affected_storage)] = 0;
  }
  affected_cells.clear();
  for (Int cell = 0; cell < global_ncells; ++cell) {
    if (!owns_cell(cell)) continue;
    if (const Int storage = storage_cell(cell);
        storage < 0
        || concentration[static_cast<size_t>(storage)] >= 0.0) {
      continue;
    }
    mark_delivery_support_cells(
        cell, domain, stencil, prescribed, affected, affected_cells,
        owns_cell, storage_cell);
  }
  for (const Int affected_storage : affected_cells) {
    const auto index = static_cast<size_t>(affected_storage);
    const auto old = prescribed[index];
    prescribed[index] = old * 0.5;
    reduced += old - prescribed[index];
  }
  return reduced;
}

constexpr Int kMaxDeliveryLocalRetries = 4;  // Bounded local feasibility pass.
// Above 25%, dilation covers the domain and only adds global-equivalent work.
constexpr Real kMaxNegativeFractionForLocalRationing = 0.25;
constexpr Int kDeliveryBisectionIterations = 12;  // Fixed feasibility search.

struct DeliveryRationingCallbacks {
  std::function<void()> restore;
  std::function<void()> restore_original;
  std::function<void()> solve;
  std::function<bool()> has_negative;
  std::function<Real()> negative_fraction;
  std::function<Real()> reduce;
  std::function<std::vector<Real>&()> prescribed;
  std::function<Real(const std::vector<Real>&)> owned_sum;
  std::function<Real(const std::vector<Real>&, const std::vector<Real>&)>
      ratio;
};

bool collective_positive(Real local_value) {
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized && !finalized) {
    Real global_value = 0.0;
    MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    return global_value > 0.0;
  }
#else
  (void)local_value;
#endif
  return local_value > 0.0;
}

DeliveryRetryResult run_delivery_rationing(
    const std::vector<Real>& original,
    const DeliveryRationingCallbacks& callbacks) {
  DeliveryRetryResult result;
  callbacks.solve();
  result.negative_after_solve = callbacks.has_negative();
  if (const auto negative_fraction_value = result.negative_after_solve
          ? callbacks.negative_fraction() : 0.0;
      negative_fraction_value <= kMaxNegativeFractionForLocalRationing) {
    for (Int attempt = 0;
         result.negative_after_solve && attempt < kMaxDeliveryLocalRetries;
         ++attempt) {
      const auto reduced = callbacks.reduce();
      if (!collective_positive(reduced)) break;
      callbacks.restore();
      callbacks.solve();
      result.retry_events += 1.0;
      result.negative_after_solve = callbacks.has_negative();
    }
  }
  if (result.negative_after_solve) {
    callbacks.restore_original();
    Real lo = 0.0;
    Real hi = 1.0;
    Real best = 0.0;
    for (Int iteration = 0;
         iteration < kDeliveryBisectionIterations; ++iteration) {
      const Real mid = (lo + hi) * 0.5;
      callbacks.restore_original();
      auto& current = callbacks.prescribed();
      std::ranges::transform(
          original, current.begin(),
          [mid](const Real value) { return value * mid; });
      callbacks.solve();
      result.retry_events += 1.0;
      if (callbacks.has_negative()) {
        hi = mid;
      } else {
        best = mid;
        lo = mid;
      }
    }
    callbacks.restore_original();
    auto& current = callbacks.prescribed();
    std::ranges::transform(
        original, current.begin(),
        [best](const Real value) { return value * best; });
    callbacks.solve();
    result.retry_events += 1.0;
    result.negative_after_solve = callbacks.has_negative();
  }
  result.delivery_reduction =
      callbacks.owned_sum(original) - callbacks.owned_sum(callbacks.prescribed());
  result.rationing_factor =
      callbacks.ratio(original, callbacks.prescribed());
  return result;
}

bool collective_negative(bool local_negative) {
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized && !finalized) {
    int local = local_negative ? 1 : 0;
    int global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    return global != 0;
  }
#else
  (void)local_negative;
#endif
  return local_negative;
}

template <typename OwnsCell, typename StorageCell>
Real owned_prescribed_sum(
    const std::vector<Real>& prescribed, Int global_ncells,
    OwnsCell owns_cell, StorageCell storage_cell) {
  Real total = 0.0;
  for (Int cell = 0; cell < global_ncells; ++cell) {
    if (!owns_cell(cell)) continue;
    const Int storage = storage_cell(cell);
    if (storage >= 0) total += prescribed[static_cast<size_t>(storage)];
  }
  return total;
}

template <typename OwnsCell, typename StorageCell>
Real minimum_prescribed_ratio(
    const std::vector<Real>& original, const std::vector<Real>& final_values,
    Int global_ncells, OwnsCell owns_cell, StorageCell storage_cell) {
  Real minimum = 1.0;
  for (Int cell = 0; cell < global_ncells; ++cell) {
    if (!owns_cell(cell)) continue;
    const Int storage = storage_cell(cell);
    if (storage < 0) continue;
    const auto index = static_cast<size_t>(storage);
    if (original[index] > 0.0) {
      minimum = std::min(minimum, final_values[index] / original[index]);
    }
  }
  return minimum;
}

template <typename OwnsCell, typename StorageCell>
std::pair<Real, Real> negative_diagnostics(
    const std::vector<Real>& concentration, Int global_ncells,
    OwnsCell owns_cell, StorageCell storage_cell, bool replicated) {
  Real minimum = 0.0;
  Real count = 0.0;
  bool found = false;
  for (Int cell = 0; cell < global_ncells; ++cell) {
    if (!owns_cell(cell)) continue;
    const Int storage = storage_cell(cell);
    if (storage < 0) continue;
    const Real value = concentration[static_cast<size_t>(storage)];
    if (value < 0.0) {
      minimum = found ? std::min(minimum, value) : value;
      found = true;
      count += 1.0;
    }
  }
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized && !finalized) {
    Real global_minimum = 0.0;
    MPI_Allreduce(&minimum, &global_minimum, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_DOUBLE,
                  replicated ? MPI_MAX : MPI_SUM,
                  MPI_COMM_WORLD);
    minimum = global_minimum;
  }
#endif
#ifndef GUTIBM_MPI
  (void)replicated;
#endif
  return {minimum, count};
}

void apply_gradient_sink(std::vector<Real>& line,
                         const std::vector<Real>& sink,
                         std::vector<Real>& gradient,
                         const ChemicalSpec* gradient_spec,
                         const Domain& domain,
                         Int first_iz) {
  if (gradient_spec == nullptr) return;
  for (size_t index = 0; index < line.size(); ++index) {
    const Int iz = first_iz + static_cast<Int>(index);
    gradient[index] = z_gradient_reference(*gradient_spec, domain, iz);
    line[index] -= sink[index] * gradient[index];
  }
}

void cap_delivery_axis_prescribed(
    const std::vector<Real>& line, const std::vector<Real>& sink,
    const std::vector<Real>& gradient, Real sink_dt,
    std::vector<Real>& prescribed, NutrientFluxAccounting* flux, Int spec,
    Real cell_volume) {
  if (!delivery_axis_cap_enabled()) return;
  for (size_t index = 0; index < prescribed.size(); ++index) {
    const Real total = line[index] + gradient[index];
    const Real available =
        std::max(0.0, total / (1.0 + sink[index] * sink_dt));
    if (prescribed[index] <= available) continue;
    if (flux != nullptr) {
      flux->add_delivery_axis_deferred_mass(
          spec, (prescribed[index] - available) * cell_volume);
    }
    prescribed[index] = available;
  }
  if (flux != nullptr) {
    for (const Real value : prescribed) {
      flux->add_delivery_axis_prescribed_mass(spec, value * cell_volume);
    }
  }
}

void fill_gradient_profile(std::vector<Real>& gradient,
                           const ChemicalSpec* gradient_spec,
                           const Domain& domain,
                           Int iz) {
  if (gradient_spec == nullptr) return;
  const Real value = z_gradient_reference(*gradient_spec, domain, iz);
  std::ranges::fill(gradient, value);
}

void fill_gradient_line(std::vector<Real>& gradient,
                        const ChemicalSpec* gradient_spec,
                        const Domain& domain, Int first_iz) {
  if (gradient_spec == nullptr) return;
  for (size_t index = 0; index < gradient.size(); ++index) {
    gradient[index] = z_gradient_reference(
        *gradient_spec, domain, first_iz + static_cast<Int>(index));
  }
}

void load_periodic_x_delivery_line(
    std::vector<Real>& line, std::vector<Real>& sink,
    std::vector<Real>& prescribed, const std::vector<Real>& concentration,
    const Domain& domain, Int iy, Int iz,
    const DeliverySinkParameters& sink_params) {
  for (Int ix = 0; ix < domain.nx(); ++ix) {
    const Int cell = domain.cell_index(ix, iy, iz);
    const auto index = static_cast<size_t>(ix);
    line[index] = concentration[static_cast<size_t>(cell)];
    sink[index] = sink_params.sink_rate[static_cast<size_t>(cell)]
        * sink_params.sink_dt;
    if (sink_params.prescribed_mass != nullptr) {
      prescribed[index] =
          (*sink_params.prescribed_mass)[static_cast<size_t>(cell)]
          / (3.0 * sink_params.cell_volume);
    }
  }
}

void load_periodic_y_delivery_line(
    std::vector<Real>& line, std::vector<Real>& sink,
    std::vector<Real>& prescribed, const std::vector<Real>& concentration,
    const Domain& domain, Int ix, Int iz,
    const DeliverySinkParameters& sink_params) {
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    const Int cell = domain.cell_index(ix, iy, iz);
    const auto index = static_cast<size_t>(iy);
    line[index] = concentration[static_cast<size_t>(cell)];
    sink[index] = sink_params.sink_rate[static_cast<size_t>(cell)]
        * sink_params.sink_dt;
    if (sink_params.prescribed_mass != nullptr) {
      prescribed[index] =
          (*sink_params.prescribed_mass)[static_cast<size_t>(cell)]
          / (3.0 * sink_params.cell_volume);
    }
  }
}

void diffuse_periodic_x_delivery(
    std::vector<Real>& concentration, const Domain& domain, Real alpha,
    const DeliverySinkParameters& sink_params) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int iy = 0; iy < ny; ++iy) {
      std::vector line(static_cast<size_t>(nx), 0.0);
      std::vector sink(static_cast<size_t>(nx), 0.0);
      std::vector prescribed(static_cast<size_t>(nx), 0.0);
      load_periodic_x_delivery_line(
          line, sink, prescribed, concentration, domain, iy, iz,
          sink_params);
      std::vector gradient(static_cast<size_t>(nx), 0.0);
      fill_gradient_profile(gradient, sink_params.gradient_spec, domain, iz);
      cap_delivery_axis_prescribed(
          line, sink, gradient, 1.0, prescribed, sink_params.flux,
          sink_params.spec, sink_params.cell_volume);
      solve_periodic_with_sink(
          line, sink, alpha,
          sink_params.gradient_spec != nullptr ? &gradient : nullptr,
          sink_params.prescribed_mass != nullptr ? &prescribed : nullptr);
      for (Int ix = 0; ix < nx; ++ix) {
        const Int cell = domain.cell_index(ix, iy, iz);
        concentration[static_cast<size_t>(cell)] = line[static_cast<size_t>(ix)];
        const Real total = sink_params.gradient_spec != nullptr
            ? line[static_cast<size_t>(ix)]
                + gradient[static_cast<size_t>(ix)]
            : line[static_cast<size_t>(ix)];
        sink_params.realized[static_cast<size_t>(cell)] += sink[
            static_cast<size_t>(ix)] * total
            * sink_params.cell_volume;
        record_negative_delivery_excursion(
            sink_params.flux, sink_params.spec, total,
            sink[static_cast<size_t>(ix)], sink_params.cell_volume, 0);
      }
    }
  }
}

void diffuse_periodic_y_delivery(
    std::vector<Real>& concentration, const Domain& domain, Real alpha,
    const DeliverySinkParameters& sink_params) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int ix = 0; ix < nx; ++ix) {
      std::vector line(static_cast<size_t>(ny), 0.0);
      std::vector sink(static_cast<size_t>(ny), 0.0);
      std::vector prescribed(static_cast<size_t>(ny), 0.0);
      load_periodic_y_delivery_line(
          line, sink, prescribed, concentration, domain, ix, iz,
          sink_params);
      std::vector gradient(static_cast<size_t>(ny), 0.0);
      fill_gradient_profile(gradient, sink_params.gradient_spec, domain, iz);
      cap_delivery_axis_prescribed(
          line, sink, gradient, 1.0, prescribed, sink_params.flux,
          sink_params.spec, sink_params.cell_volume);
      solve_periodic_with_sink(
          line, sink, alpha,
          sink_params.gradient_spec != nullptr ? &gradient : nullptr,
          sink_params.prescribed_mass != nullptr ? &prescribed : nullptr);
      for (Int iy = 0; iy < ny; ++iy) {
        const Int cell = domain.cell_index(ix, iy, iz);
        concentration[static_cast<size_t>(cell)] = line[static_cast<size_t>(iy)];
        const Real total = sink_params.gradient_spec != nullptr
            ? line[static_cast<size_t>(iy)]
                + gradient[static_cast<size_t>(iy)]
            : line[static_cast<size_t>(iy)];
        sink_params.realized[static_cast<size_t>(cell)] += sink[
            static_cast<size_t>(iy)] * total
            * sink_params.cell_volume;
        record_negative_delivery_excursion(
            sink_params.flux, sink_params.spec, total,
            sink[static_cast<size_t>(iy)], sink_params.cell_volume, 1);
      }
    }
  }
}

Real diffuse_bounded_z(std::vector<Real>& concentration,
                       const Domain& domain,
                       Real alpha,
                       Real boundary_conc,
                       Real cell_volume) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  if (nz <= 1) return 0.0;

  const NeumannTopLineSolver solver(nz - 1, alpha);
  Real face_exchange = 0.0;
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<Real> line(static_cast<size_t>(nz - 1));
    #ifdef GUTIBM_OPENMP
    #pragma omp for collapse(2) schedule(static)
    #endif
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        for (Int iz = 1; iz < nz; ++iz) {
          line[static_cast<size_t>(iz - 1)] =
              concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))];
        }
        solver.solve(line, boundary_conc);
        #ifdef GUTIBM_OPENMP
        #pragma omp atomic
        #endif
        face_exchange += alpha * (boundary_conc - line.front()) * cell_volume;
        for (Int iz = 1; iz < nz; ++iz) {
          concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
              line[static_cast<size_t>(iz - 1)];
        }
      }
    }
  }
  return face_exchange;
}

struct DeliveryBoundaryParameters {
  Real boundary_conc = 0.0;
  Real beta = 0.0;
  Real flux_source = 0.0;
  EpithelialBoundaryMode mode = EpithelialBoundaryMode::Flux;
  Real cell_volume = 0.0;
};

struct ReplicatedDeliveryLineContext {
  std::vector<Real>& concentration;
  const Domain& domain;
  Real alpha = 0.0;
  const DeliveryBoundaryParameters& boundary;
  const DeliverySinkParameters& sink;
};

Real delivery_face_exchange(
    Real alpha, const DeliveryBoundaryParameters& params,
    const std::vector<Real>& line) {
  if (params.mode == EpithelialBoundaryMode::Robin) {
    return params.beta * (params.boundary_conc - line.front())
        * params.cell_volume;
  }
  if (params.mode == EpithelialBoundaryMode::Flux) {
    return params.flux_source * params.cell_volume;
  }
  return alpha * (params.boundary_conc - line.front())
      * params.cell_volume;
}

Real delivery_boundary_source(
    Real alpha, const DeliveryBoundaryParameters& params) {
  if (params.mode == EpithelialBoundaryMode::Robin) {
    return params.beta * params.boundary_conc;
  }
  if (params.mode == EpithelialBoundaryMode::Flux) {
    return params.flux_source;
  }
  return alpha * params.boundary_conc;
}

void solve_replicated_delivery_z_line(
    Int ix, Int iy, ReplicatedDeliveryLineContext& context,
    Real& face_exchange) {
  const Int nz = context.domain.nz();
  std::vector<Real> line(static_cast<size_t>(nz - 1));
  std::vector<Real> sink(static_cast<size_t>(nz - 1));
  std::vector<Real> prescribed(static_cast<size_t>(nz - 1), 0.0);
  for (Int iz = 1; iz < nz; ++iz) {
    const Int cell = context.domain.cell_index(ix, iy, iz);
    const auto index = static_cast<size_t>(iz - 1);
    line[index] = context.concentration[static_cast<size_t>(cell)];
    sink[index] = context.sink.sink_rate[static_cast<size_t>(cell)]
        * context.sink.sink_dt;
    if (context.sink.prescribed_mass != nullptr) {
      prescribed[index] =
          (*context.sink.prescribed_mass)[static_cast<size_t>(cell)]
          / (3.0 * context.sink.cell_volume);
    }
  }
  std::vector gradient(static_cast<size_t>(nz - 1), 0.0);
  fill_gradient_line(
      gradient, context.sink.gradient_spec, context.domain, 1);
  cap_delivery_axis_prescribed(
      line, sink, gradient, 1.0, prescribed, context.sink.flux,
      context.sink.spec, context.sink.cell_volume);
  apply_gradient_sink(
      line, sink, gradient, context.sink.gradient_spec,
      context.domain, 1);
  for (size_t index = 0; index < line.size(); ++index) {
    line[index] -= prescribed[index];
  }
  std::vector diagonal(static_cast<size_t>(nz - 1), 0.0);
  for (Int iz = 1; iz < nz; ++iz) {
    diagonal[static_cast<size_t>(iz - 1)] =
        1.0 + 2.0 * context.alpha + sink[static_cast<size_t>(iz - 1)];
  }
  diagonal.back() = 1.0 + context.alpha + sink.back();
  solve_tridiagonal_with_diagonal(
      line, diagonal, context.alpha,
      context.alpha * context.boundary.boundary_conc);
  face_exchange += delivery_face_exchange(
      context.alpha, context.boundary, line);
  for (Int iz = 1; iz < nz; ++iz) {
    const Int cell = context.domain.cell_index(ix, iy, iz);
    const auto index = static_cast<size_t>(iz - 1);
    context.concentration[static_cast<size_t>(cell)] = line[index];
    const Real total = context.sink.gradient_spec != nullptr
        ? line[index] + gradient[index] : line[index];
    context.sink.realized[static_cast<size_t>(cell)] +=
        sink[index] * total * context.sink.cell_volume;
    record_negative_delivery_excursion(
        context.sink.flux, context.sink.spec, total, sink[index],
        context.sink.cell_volume, 2);
  }
}

Real diffuse_bounded_z_delivery(
    std::vector<Real>& concentration, const Domain& domain, Real alpha,
    const DeliveryBoundaryParameters& params,
    const DeliverySinkParameters& sink_params) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  if (const Int nz = domain.nz(); nz <= 1) return 0.0;
  ReplicatedDeliveryLineContext line_context{
      concentration, domain, alpha, params, sink_params};
  Real face_exchange = 0.0;
  for (Int iy = 0; iy < ny; ++iy) {
    for (Int ix = 0; ix < nx; ++ix) {
      solve_replicated_delivery_z_line(
          ix, iy, line_context, face_exchange);
    }
  }
  return face_exchange;
}

class DeliveryBottomLineSolver {
 public:
  DeliveryBottomLineSolver(Int size, Real alpha,
                           EpithelialBoundaryMode mode, Real beta)
      : size_(size), alpha_(alpha), beta_(beta), mode_(mode) {
    if (size_ <= 0) return;
    std::vector lower(static_cast<size_t>(std::max(size_ - 1, 0)), -alpha_);
    std::vector upper(static_cast<size_t>(std::max(size_ - 1, 0)), -alpha_);
    std::vector diagonal(static_cast<size_t>(size_), 1.0 + 2.0 * alpha_);
    diagonal.front() = mode_ == EpithelialBoundaryMode::Robin
        ? 1.0 + alpha_ + beta_ : 1.0 + alpha_;
    diagonal.back() = 1.0 + alpha_;
    factorization_.factorize(lower, diagonal, upper);
  }

  void solve(std::vector<Real>& values, Real source) const {
    if (size_ <= 0) return;
    values.front() += source;
    factorization_.solve_in_place(values);
  }

 private:
  Int size_ = 0;
  Real alpha_ = 0.0;
  Real beta_ = 0.0;
  EpithelialBoundaryMode mode_ = EpithelialBoundaryMode::Flux;
  TridiagonalFactorization factorization_;
};

struct DeliveryGridParameters {
  Int nx = 0;
  Int ny = 0;
  Int nz = 0;
  Real alpha = 0.0;
};

template <typename LoadLine, typename StoreLine, typename LoadSink,
          typename LoadPrescribed, typename LoadProfile, typename AddRealized>
struct DeliveryLineOperations {
  LoadLine load_line;
  StoreLine store_line;
  LoadSink load_sink;
  LoadPrescribed load_prescribed;
  LoadProfile load_profile;
  AddRealized add_realized;
  Real sink_dt = 0.0;
  NutrientFluxAccounting* flux = nullptr;
  Int spec = -1;
};

template <typename LoadLine, typename StoreLine, typename LoadSink,
          typename LoadPrescribed, typename LoadProfile, typename AddRealized>
void solve_delivery_z_line(
    Int ix, Int iy, const DeliveryGridParameters& grid,
    const DeliveryBoundaryParameters& params,
    const DeliveryLineOperations<
        LoadLine, StoreLine, LoadSink, LoadPrescribed, LoadProfile,
        AddRealized>& operations,
    Real& face_exchange) {
  std::vector<Real> line(static_cast<size_t>(grid.nz));
  std::vector<Real> sink(static_cast<size_t>(grid.nz));
  std::vector gradient(static_cast<size_t>(grid.nz), 0.0);
  std::vector diagonal(static_cast<size_t>(grid.nz), 0.0);
  operations.load_line(ix, iy, line);
  operations.load_sink(ix, iy, sink);
  std::vector prescribed(static_cast<size_t>(grid.nz), 0.0);
  operations.load_prescribed(ix, iy, prescribed);
  const bool has_gradient = operations.load_profile(ix, iy, gradient);
  cap_delivery_axis_prescribed(
      line, sink, gradient, 1.0, prescribed,
      operations.flux, operations.spec, params.cell_volume);
  if (has_gradient) {
    for (Int iz = 0; iz < grid.nz; ++iz) {
    const auto index = static_cast<size_t>(iz);
      line[index] -= sink[index] * operations.sink_dt * gradient[index];
    }
  }
  for (Int iz = 0; iz < grid.nz; ++iz) {
    line[static_cast<size_t>(iz)] -= prescribed[static_cast<size_t>(iz)];
  }
  for (Int iz = 0; iz < grid.nz; ++iz) {
    diagonal[static_cast<size_t>(iz)] =
        1.0 + 2.0 * grid.alpha
        + sink[static_cast<size_t>(iz)] * operations.sink_dt;
  }
  diagonal.front() += params.mode == EpithelialBoundaryMode::Robin
      ? -grid.alpha + params.beta : -grid.alpha;
  diagonal.back() -= grid.alpha;
  const Real source = delivery_boundary_source(
      grid.alpha, params);
  solve_tridiagonal_with_diagonal(line, diagonal, grid.alpha, source);
  const Real boundary_realized = delivery_face_exchange(
      grid.alpha, params, line);
  #ifdef GUTIBM_OPENMP
  #pragma omp atomic
  #endif
  face_exchange += boundary_realized;
  for (Int iz = 0; iz < grid.nz; ++iz) {
    const auto index = static_cast<size_t>(iz);
    operations.store_line(ix, iy, iz, line[index]);
    const Real total = has_gradient
        ? line[index] + gradient[index] : line[index];
    const Real amount = sink[index] * operations.sink_dt * total
        * params.cell_volume;
    operations.add_realized(ix, iy, iz, amount, total);
  }
}

template <typename LoadLine, typename StoreLine>
Real diffuse_bounded_z_delivery_impl(
    Int nx, Int ny, Int nz, Real alpha,
    const DeliveryBoundaryParameters& params,
    LoadLine load_line, StoreLine store_line) {
  if (nz <= 0) return 0.0;

  const DeliveryBottomLineSolver solver(nz, alpha, params.mode, params.beta);
  Real face_exchange = 0.0;
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<Real> line(static_cast<size_t>(nz));
    #ifdef GUTIBM_OPENMP
    #pragma omp for collapse(2) schedule(static)
    #endif
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        load_line(ix, iy, line);
        const Real source = delivery_boundary_source(alpha, params);
        solver.solve(line, source);
        const Real realized = delivery_face_exchange(
            alpha, params, line);
        #ifdef GUTIBM_OPENMP
        #pragma omp atomic
        #endif
        face_exchange += realized;
        store_line(ix, iy, line);
      }
    }
  }
  return face_exchange;
}

Real diffuse_bounded_z_delivery(
    std::vector<Real>& concentration, const Domain& domain, Real alpha,
    const DeliveryBoundaryParameters& params) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const auto load_line = [&concentration, &domain, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] =
          concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))];
    }
  };
  const auto store_line = [&concentration, &domain, nz](
                              Int ix, Int iy,
                              const std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
          line[static_cast<size_t>(iz)];
    }
  };
  return diffuse_bounded_z_delivery_impl(
      nx, ny, nz, alpha, params, load_line, store_line);
}

template <typename LoadLine, typename StoreLine, typename LoadSink,
          typename LoadPrescribed, typename LoadProfile, typename AddRealized>
Real diffuse_bounded_z_delivery_with_sink_impl(
    const DeliveryGridParameters& grid,
    const DeliveryBoundaryParameters& params,
    const DeliveryLineOperations<
        LoadLine, StoreLine, LoadSink, LoadPrescribed, LoadProfile,
        AddRealized>& operations) {
  const Int nx = grid.nx;
  const Int ny = grid.ny;
  if (const Int nz = grid.nz; nz <= 0) return 0.0;

  Real face_exchange = 0.0;
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef GUTIBM_OPENMP
    #pragma omp for collapse(2) schedule(static)
    #endif
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        solve_delivery_z_line(
            ix, iy, grid, params, operations, face_exchange);
      }
    }
  }
  return face_exchange;
}

Real diffuse_bounded_z_delivery_with_sink(
    std::vector<Real>& concentration, const Domain& domain, Real alpha,
    const DeliveryBoundaryParameters& params,
    const DeliverySinkParameters& sink_params) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const auto load_line = [&concentration, &domain, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] =
          concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))];
    }
  };
  const auto store_line = [&concentration, &domain](
                              Int ix, Int iy, Int iz, Real value) {
    concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] =
        value;
  };
  const auto load_sink = [&sink_params, &domain, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] =
          sink_params.sink_rate[
              static_cast<size_t>(domain.cell_index(ix, iy, iz))];
    }
  };
  const auto load_prescribed = [&sink_params, &domain, nz](
                                   Int ix, Int iy,
                                   std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      const Int cell = domain.cell_index(ix, iy, iz);
      line[static_cast<size_t>(iz)] =
          sink_params.prescribed_mass == nullptr
              ? 0.0
              : (*sink_params.prescribed_mass)[static_cast<size_t>(cell)]
                    / (3.0 * sink_params.cell_volume);
    }
  };
  const auto load_profile = [gradient_spec = sink_params.gradient_spec,
                             &domain, nz](
                                 Int, Int, std::vector<Real>& line) {
    if (gradient_spec == nullptr) return false;
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] =
          z_gradient_reference(*gradient_spec, domain, iz);
    }
    return true;
  };
  const auto add_realized = [
      &realized = sink_params.realized, &domain, &sink_params](
      Int ix, Int iy, Int iz, Real amount, Real total) {
    realized[static_cast<size_t>(domain.cell_index(ix, iy, iz))] += amount;
    record_negative_delivery_excursion(
        sink_params.flux, sink_params.spec, -amount, 2);
    if (sink_params.flux != nullptr && total < 0.0) {
      sink_params.flux->add_negative_delivery_min(
          sink_params.spec, total);
    }
  };
  const DeliveryLineOperations<
      decltype(load_line), decltype(store_line), decltype(load_sink),
      decltype(load_prescribed), decltype(load_profile),
      decltype(add_realized)>
      operations{load_line, store_line, load_sink, load_prescribed,
                 load_profile, add_realized, sink_params.sink_dt,
                 sink_params.flux, sink_params.spec};
  return diffuse_bounded_z_delivery_with_sink_impl(
      {nx, ny, nz, alpha}, params, operations);
}

Real set_epithelial_boundary(std::vector<Real>& concentration,
                             const Domain& domain,
                             Real boundary_conc,
                             Real cell_volume) {
  Real amount = 0.0;
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.nx(); ++ix) {
      const auto index = static_cast<size_t>(
          domain.cell_index(ix, iy, 0));
      amount += (boundary_conc - concentration[index]) * cell_volume;
      concentration[index] = boundary_conc;
    }
  }
  return amount;
}

void set_luminal_neumann_boundary(std::vector<Real>& concentration,
                                  const Domain& domain) {
  if (domain.nz() < 2) return;
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.nx(); ++ix) {
      const auto top = static_cast<size_t>(
          domain.cell_index(ix, iy, domain.nz() - 1));
      const auto below = static_cast<size_t>(
          domain.cell_index(ix, iy, domain.nz() - 2));
      concentration[top] = concentration[below];
    }
  }
}

Real z_gradient_reference(const ChemicalSpec& spec,
                          const Domain& domain,
                          Int iz) {
  if (iz == 0) return spec.boundary_conc;
  const Int profile_iz = (domain.nz() >= 2 && iz == domain.nz() - 1)
      ? domain.nz() - 2 : iz;
    const Real z_rel = (profile_iz + 0.5) * domain.dx_z();
  return spec.initial_conc * std::exp(-z_rel / spec.z_gradient_lambda);
}

void shift_z_gradient(std::vector<Real>& concentration,
                      const ChemicalSpec& spec,
                      const Domain& domain,
                      Real scale) {
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    const Real shift = scale * z_gradient_reference(spec, domain, iz);
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.nx(); ++ix) {
        concentration[static_cast<size_t>(domain.cell_index(ix, iy, iz))] += shift;
      }
    }
  }
}

void clamp_nonnegative(std::vector<Real>& concentration) {
  const auto size = static_cast<Int>(concentration.size());
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (Int cell = 0; cell < size; ++cell) {
    concentration[static_cast<size_t>(cell)] =
        std::max(concentration[static_cast<size_t>(cell)], 0.0);
  }
}

void validate_chemical_decomposition(const Domain& domain,
                                    std::string_view decomposition) {
  if (decomposition != "replicated" && decomposition != "slab") {
    throw ConfigError("invalid chemistry decomposition mode");
  }
  if (decomposition != "slab") return;

  if (const auto required_halo = static_cast<Int>(
          std::ceil(domain.ghost_width() / domain.dx_x()));
      domain.grid_halo_width() < required_halo) {
    throw ConfigError(
        "slab chemistry requires grid_halo_width >= ceil(ghost_width / dx)");
  }
  if (domain.local_grid_nx() < domain.grid_halo_width()) {
    throw ConfigError(
        "slab chemistry requires each owned x-slab to be at least "
        "grid_halo_width cells wide");
  }
}

void validate_epithelial_delivery(const std::vector<ChemicalSpec>& specs) {
  for (const auto& spec : specs) {
    if (spec.epithelial_transfer_coeff < 0.0
        || spec.epithelial_flux < 0.0) {
      throw ConfigError(
          "epithelial delivery coefficients and fluxes must be nonnegative");
    }
    if (spec.epithelial_boundary_mode != EpithelialBoundaryMode::Dirichlet
        && spec.z_gradient_enabled) {
      throw ConfigError(
          "z_gradient_enabled cannot be combined with Robin or flux "
          "epithelial boundary modes");
    }
  }
}

}  // namespace

void ChemicalField::init(const Domain& domain,
                          const std::vector<ChemicalSpec>& specs,
                          std::string_view decomposition,
                          Real delivery_far_field_radius) {
  validate_chemical_decomposition(domain, decomposition);
  validate_epithelial_delivery(specs);
  delivery_far_field_radius_ = delivery_far_field_radius;

  domain_ = &domain;
  mode_ = decomposition == "slab"
      ? DecompositionMode::Slab : DecompositionMode::Replicated;
  global_nx_ = domain.nx();
  global_ny_ = domain.ny();
  global_nz_ = domain.nz();
  global_ncells_ = domain.ncells();
  owned_x_begin_ = decomposition == "slab"
      ? domain.local_grid_x_begin() : 0;
  owned_x_end_ = decomposition == "slab"
      ? domain.local_grid_x_end() : domain.nx();
  halo_width_ = decomposition == "slab" ? domain.grid_halo_width() : 0;
  storage_nx_ = decomposition == "slab"
      ? domain.local_grid_storage_nx() : domain.nx();
  specs_  = specs;
  nspec_  = static_cast<Int>(specs.size());
  ncells_ = storage_nx_ * global_ny_ * global_nz_;
  flux_accounting_.init(specs.size());
  delivery_support_stencil_ =
      make_delivery_support_stencil(domain, delivery_far_field_radius_);
  delivery_affected_mask_.assign(static_cast<size_t>(ncells_), 0);
  delivery_affected_cells_.reserve(static_cast<size_t>(ncells_));

  conc_.resize(nspec_);
  reac_.resize(nspec_);
  sink_rate_.assign(static_cast<size_t>(nspec_),
                    std::vector<Real>(static_cast<size_t>(ncells_), 0.0));
  vbf_sink_rate_.assign(
      static_cast<size_t>(nspec_),
      std::vector<Real>(static_cast<size_t>(ncells_), 0.0));
  sink_realized_.assign(static_cast<size_t>(nspec_),
                        std::vector<Real>(
                            static_cast<size_t>(ncells_), 0.0));
  total_sink_realized_.assign(
      static_cast<size_t>(nspec_),
      std::vector<Real>(static_cast<size_t>(ncells_), 0.0));
  vbf_sink_realized_.assign(
      static_cast<size_t>(nspec_),
      std::vector<Real>(static_cast<size_t>(ncells_), 0.0));
  prescribed_sink_.assign(
      static_cast<size_t>(nspec_),
      std::vector<Real>(static_cast<size_t>(ncells_), 0.0));
  prescribed_active_.assign(static_cast<size_t>(nspec_), false);
  for (Int s = 0; s < nspec_; ++s) {
    conc_[s].assign(ncells_, specs_[s].initial_conc);
    reac_[s].assign(ncells_, 0.0);

    if (specs_[s].z_gradient_enabled) {
      for (Int storage_cell = 0; storage_cell < ncells_; ++storage_cell) {
        const Int global_cell = storage_to_global_cell(storage_cell);
        if (global_cell < 0) continue;
        const Int iz = global_cell / (global_nx_ * global_ny_);
        const Real z_rel = (iz + 0.5) * domain.dx_z();
        conc_[s][static_cast<size_t>(storage_cell)] =
            specs_[s].initial_conc
            * std::exp(-z_rel / specs_[s].z_gradient_lambda);
      }
    }
  }
}

Int ChemicalField::global_to_storage_cell(Int global_cell) const {
  assert(global_cell >= 0 && global_cell < global_ncells_);
  if (mode_ == DecompositionMode::Replicated) return global_cell;

  const Int ix = global_cell % global_nx_;
  const Int yz = global_cell / global_nx_;
  const Int local_ix = domain_->global_to_local_grid_x(ix);
  if (local_ix < 0) return -1;
  return yz * storage_nx_ + local_ix;
}

Int ChemicalField::storage_to_global_cell(Int storage_cell) const {
  assert(storage_cell >= 0 && storage_cell < ncells_);
  if (mode_ == DecompositionMode::Replicated) return storage_cell;

  const Int local_ix = storage_cell % storage_nx_;
  const Int yz = storage_cell / storage_nx_;
  const Int global_ix = domain_->local_to_global_grid_x(local_ix);
  if (global_ix < 0) return -1;
  return yz * global_nx_ + global_ix;
}

bool ChemicalField::owns_global_cell(Int global_cell) const {
  if (global_cell < 0 || global_cell >= global_ncells_) return false;
  const Int ix = global_cell % global_nx_;
  return mode_ != DecompositionMode::Slab
      || (ix >= owned_x_begin_ && ix < owned_x_end_);
}

bool ChemicalField::global_cell_in_halo(Int global_cell) const {
  if (global_cell < 0 || global_cell >= global_ncells_) return false;
  return mode_ == DecompositionMode::Slab
      && !owns_global_cell(global_cell)
      && global_to_storage_cell(global_cell) >= 0;
}

void ChemicalField::exchange_concentration_halos() {
  if (mode_ != DecompositionMode::Slab) return;

  const Int plane_cells = global_ny_ * global_nz_;
  const Int halo_cells = halo_width_ * plane_cells;
  const Int local_nx = owned_x_end_ - owned_x_begin_;
  if (halo_cells <= 0 || local_nx <= 0) return;

  if (domain_->nprocs() <= 1) {
    for (Int s = 0; s < nspec_; ++s) {
      for (Int iz = 0; iz < global_nz_; ++iz) {
        for (Int iy = 0; iy < global_ny_; ++iy) {
          for (Int h = 0; h < halo_width_; ++h) {
            const Int lo_storage = iz * storage_nx_ * global_ny_
                + iy * storage_nx_ + h;
            const Int hi_storage = iz * storage_nx_ * global_ny_
                + iy * storage_nx_ + halo_width_ + local_nx + h;
            const Int lo_source = iz * storage_nx_ * global_ny_
                + iy * storage_nx_ + halo_width_ + local_nx
                - halo_width_ + h;
            const Int hi_source = iz * storage_nx_ * global_ny_
                + iy * storage_nx_ + halo_width_ + h;
            if (!domain_->config().periodic[0]) continue;
            conc_[static_cast<size_t>(s)][static_cast<size_t>(lo_storage)] =
                conc_[static_cast<size_t>(s)][static_cast<size_t>(lo_source)];
            conc_[static_cast<size_t>(s)][static_cast<size_t>(hi_storage)] =
                conc_[static_cast<size_t>(s)][static_cast<size_t>(hi_source)];
          }
        }
      }
    }
    return;
  }

#ifdef GUTIBM_MPI
  const int rank_lo = domain_->rank_lo() >= 0
      ? domain_->rank_lo() : MPI_PROC_NULL;
  const int rank_hi = domain_->rank_hi() >= 0
      ? domain_->rank_hi() : MPI_PROC_NULL;
  for (Int s = 0; s < nspec_; ++s) {
    std::vector<Real> send_lo(static_cast<size_t>(halo_cells));
    std::vector<Real> send_hi(static_cast<size_t>(halo_cells));
    std::vector<Real> recv_lo(static_cast<size_t>(halo_cells));
    std::vector<Real> recv_hi(static_cast<size_t>(halo_cells));
    for (Int iz = 0; iz < global_nz_; ++iz) {
      for (Int iy = 0; iy < global_ny_; ++iy) {
        for (Int h = 0; h < halo_width_; ++h) {
          const Int lo_storage = iz * storage_nx_ * global_ny_
              + iy * storage_nx_ + halo_width_ + h;
          const Int hi_storage = iz * storage_nx_ * global_ny_
              + iy * storage_nx_ + halo_width_ + local_nx
              - halo_width_ + h;
          const Int offset = (iz * global_ny_ + iy) * halo_width_ + h;
          send_lo[static_cast<size_t>(offset)] =
              conc_[static_cast<size_t>(s)][static_cast<size_t>(lo_storage)];
          send_hi[static_cast<size_t>(offset)] =
              conc_[static_cast<size_t>(s)][static_cast<size_t>(hi_storage)];
        }
      }
    }

    const int tag = 700 + 2 * static_cast<int>(s);
    MPI_Sendrecv(send_lo.data(), halo_cells, MPI_DOUBLE, rank_lo, tag,
                 recv_hi.data(), halo_cells, MPI_DOUBLE, rank_hi, tag,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(send_hi.data(), halo_cells, MPI_DOUBLE, rank_hi, tag + 1,
                 recv_lo.data(), halo_cells, MPI_DOUBLE, rank_lo, tag + 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (Int iz = 0; iz < global_nz_; ++iz) {
      for (Int iy = 0; iy < global_ny_; ++iy) {
        for (Int h = 0; h < halo_width_; ++h) {
          const Int lo_storage = iz * storage_nx_ * global_ny_
              + iy * storage_nx_ + h;
          const Int hi_storage = iz * storage_nx_ * global_ny_
              + iy * storage_nx_ + halo_width_ + local_nx + h;
          const Int offset = (iz * global_ny_ + iy) * halo_width_ + h;
          conc_[static_cast<size_t>(s)][static_cast<size_t>(lo_storage)] =
              recv_lo[static_cast<size_t>(offset)];
          conc_[static_cast<size_t>(s)][static_cast<size_t>(hi_storage)] =
              recv_hi[static_cast<size_t>(offset)];
        }
      }
    }
  }
#endif
}

void ChemicalField::zero_reactions() {
  for (Int s = 0; s < nspec_; ++s) {
    std::ranges::fill(reac_[s], 0.0);
    std::ranges::fill(sink_rate_[s], 0.0);
    std::ranges::fill(vbf_sink_rate_[s], 0.0);
    std::ranges::fill(sink_realized_[s], 0.0);
    std::ranges::fill(total_sink_realized_[s], 0.0);
    std::ranges::fill(vbf_sink_realized_[s], 0.0);
    std::ranges::fill(prescribed_sink_[s], 0.0);
    flux_accounting_.negative_delivery_events_step[
        static_cast<size_t>(s)] = 0.0;
    flux_accounting_.negative_delivery_mass_step[
        static_cast<size_t>(s)] = 0.0;
    prescribed_active_[static_cast<size_t>(s)] = false;
  }
  if (!nutrient_debug_enabled() || domain_ == nullptr) return;
  ++nutrient_debug_step_counter();
  debug_initial_content_.resize(static_cast<size_t>(nspec_));
  for (Int s = 0; s < nspec_; ++s) {
    debug_initial_content_[static_cast<size_t>(s)] =
        mode_ == DecompositionMode::Slab
        ? owned_content_slab(
              conc_[static_cast<size_t>(s)], *domain_, storage_nx_,
              halo_width_, domain_->cell_volume())
        : owned_content(
              conc_[static_cast<size_t>(s)], *domain_, domain_->cell_volume());
  }
}

void ChemicalField::add_sink_rate_global(Int spec, Int cell, Real rate) {
  const Int storage_cell = global_to_storage_cell(cell);
  if (spec < 0 || spec >= nspec_ || storage_cell < 0 || rate <= 0.0) {
    return;
  }
  #ifdef GUTIBM_OPENMP
  #pragma omp atomic
  #endif
  sink_rate_[static_cast<size_t>(spec)][static_cast<size_t>(storage_cell)]
      += rate;
}

void ChemicalField::add_prescribed_sink_global(
    Int spec, Int cell, Real amount) {
  const Int storage_cell = global_to_storage_cell(cell);
  if (spec < 0 || spec >= nspec_ || storage_cell < 0
      || !owns_global_cell(cell) || amount <= 0.0) {
    return;
  }
  #ifdef GUTIBM_OPENMP
  #pragma omp atomic
  #endif
  prescribed_sink_[static_cast<size_t>(spec)]
      [static_cast<size_t>(storage_cell)] += amount;
  prescribed_active_[static_cast<size_t>(spec)] = true;
}

Real ChemicalField::prescribed_sink_global(Int spec, Int cell) const {
  const Int storage_cell = global_to_storage_cell(cell);
  if (spec < 0 || spec >= nspec_ || storage_cell < 0) return 0.0;
  return prescribed_sink_[static_cast<size_t>(spec)]
      [static_cast<size_t>(storage_cell)];
}

void ChemicalField::add_delivery_reduction(Int spec, Real amount) {
  if (spec < 0 || spec >= nspec_ || amount <= 0.0) return;
  #ifdef GUTIBM_OPENMP
  #pragma omp atomic
  #endif
  flux_accounting_.delivery_reduction_step[
      static_cast<size_t>(spec)] += amount;
}

void ChemicalField::add_vbf_sink_rate_global(Int spec, Int cell, Real rate) {
  const Int storage_cell = global_to_storage_cell(cell);
  if (spec < 0 || spec >= nspec_ || storage_cell < 0 || rate <= 0.0) {
    return;
  }
  #ifdef GUTIBM_OPENMP
  #pragma omp atomic
  #endif
  vbf_sink_rate_[static_cast<size_t>(spec)]
                [static_cast<size_t>(storage_cell)] += rate;
  #ifdef GUTIBM_OPENMP
  #pragma omp atomic
  #endif
  sink_rate_[static_cast<size_t>(spec)]
            [static_cast<size_t>(storage_cell)] += rate;
}

void ChemicalField::add_vbf_sink_rates(
    Int spec, const std::vector<Real>& rates) {
  if (spec < 0 || spec >= nspec_) return;
  auto& vbf_rates = vbf_sink_rate_[static_cast<size_t>(spec)];
  auto& total_rates = sink_rate_[static_cast<size_t>(spec)];
  if (rates.size() != vbf_rates.size()) {
    throw Error("VBF sink-rate buffer size mismatch");
  }
  for (size_t i = 0; i < rates.size(); ++i) {
    const Real rate = rates[i];
    if (rate <= 0.0) continue;
    vbf_rates[i] += rate;
    total_rates[i] += rate;
  }
}

void ChemicalField::split_delivery_sink_realized(Int spec) {
  if (spec < 0 || spec >= nspec_) return;
  for (Int cell = 0; cell < global_ncells_; ++cell) {
    if (!owns_global_cell(cell)) continue;
    const Int storage_cell = global_to_storage_cell(cell);
    if (storage_cell < 0) continue;
    const auto index = static_cast<size_t>(storage_cell);
    const Real total_rate = sink_rate_[static_cast<size_t>(spec)][index];
    const Real vbf_rate = vbf_sink_rate_[static_cast<size_t>(spec)][index];
    const Real total = total_sink_realized_[static_cast<size_t>(spec)][index];
    const Real vbf_share = total_rate > 0.0
        ? total * std::clamp(vbf_rate / total_rate, 0.0, 1.0) : 0.0;
    vbf_sink_realized_[static_cast<size_t>(spec)][index] = vbf_share;
    sink_realized_[static_cast<size_t>(spec)][index] = total - vbf_share;
  }
}

void ChemicalField::add_sink_rate_global(Int cell, Real rate) {
  const Int carbon = find(species::CARBON);
  if (carbon >= 0) add_sink_rate_global(carbon, cell, rate);
}

Real ChemicalField::sink_realized_global(Int spec, Int cell) const {
  const Int storage_cell = global_to_storage_cell(cell);
  if (spec < 0 || spec >= nspec_ || storage_cell < 0) return 0.0;
  return sink_realized_[static_cast<size_t>(spec)]
      [static_cast<size_t>(storage_cell)];
}

Real ChemicalField::sink_realized_global(Int cell) const {
  const Int carbon = find(species::CARBON);
  return carbon >= 0 ? sink_realized_global(carbon, cell) : 0.0;
}

Real ChemicalField::total_sink_realized_global(Int spec, Int cell) const {
  const Int storage_cell = global_to_storage_cell(cell);
  if (spec < 0 || spec >= nspec_ || storage_cell < 0) return 0.0;
  return total_sink_realized_[static_cast<size_t>(spec)]
      [static_cast<size_t>(storage_cell)];
}

Real ChemicalField::vbf_sink_realized_global(Int spec, Int cell) const {
  const Int storage_cell = global_to_storage_cell(cell);
  if (spec < 0 || spec >= nspec_ || storage_cell < 0) return 0.0;
  return vbf_sink_realized_[static_cast<size_t>(spec)]
      [static_cast<size_t>(storage_cell)];
}

Real ChemicalField::vbf_sink_realized(Int spec) const {
  if (spec < 0 || spec >= nspec_) return 0.0;
  Real total = 0.0;
  for (Int cell = 0; cell < global_ncells_; ++cell) {
    if (owns_global_cell(cell)) {
      total += vbf_sink_realized_global(spec, cell);
    }
  }
  return total;
}

bool ChemicalField::has_sink_rate(Int spec) const {
  if (spec < 0 || spec >= nspec_) return false;
  const bool local = std::ranges::any_of(
      sink_rate_[static_cast<size_t>(spec)],
      [](Real value) { return value > 0.0; });
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized && !finalized) {
    int local_value = local ? 1 : 0;
    int global_value = 0;
    MPI_Allreduce(&local_value, &global_value, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    return global_value != 0;
  }
#endif
  return local;
}

bool ChemicalField::has_sink_rate() const {
  const Int carbon = find(species::CARBON);
  return carbon >= 0 && has_sink_rate(carbon);
}

void ChemicalField::sum_reactions_across_ranks() {
#ifdef GUTIBM_MPI
  if (mode_ == DecompositionMode::Slab) return;
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return;

  int ranks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  if (ranks <= 1) return;
  for (auto& reaction : reac_) {
    MPI_Allreduce(MPI_IN_PLACE, reaction.data(), ncells_, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
  }
  for (auto& sink : sink_rate_) {
    MPI_Allreduce(MPI_IN_PLACE, sink.data(), ncells_, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
  }
#endif
}

void ChemicalField::sum_agent_uptake_across_ranks() {
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return;
  int ranks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  if (ranks <= 1) return;
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.agent_uptake_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.maintenance_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.maintenance_shortfall_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.maintenance_limited_agents_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.uptake_demand_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.uptake_shortfall_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,
                flux_accounting_.uptake_limited_step.data(), nspec_,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}

void ChemicalField::sum_prescribed_sinks_across_ranks() {
#ifdef GUTIBM_MPI
  if (mode_ == DecompositionMode::Slab) return;
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return;
  int ranks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  if (ranks <= 1) return;
  for (auto& prescribed : prescribed_sink_) {
    MPI_Allreduce(MPI_IN_PLACE, prescribed.data(), ncells_, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
  }
  std::vector active(static_cast<size_t>(nspec_), 0);
  for (Int s = 0; s < nspec_; ++s) {
    active[static_cast<size_t>(s)] =
        prescribed_active_[static_cast<size_t>(s)] ? 1 : 0;
  }
  MPI_Allreduce(
      MPI_IN_PLACE, active.data(), nspec_, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  for (Int s = 0; s < nspec_; ++s) {
    prescribed_active_[static_cast<size_t>(s)] =
        active[static_cast<size_t>(s)] != 0;
  }
#endif
}

void ChemicalField::sum_values_across_ranks(
    std::vector<Real>& values) const {
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return;
  int ranks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  if (ranks <= 1 || values.empty()) return;
  MPI_Allreduce(
      MPI_IN_PLACE, values.data(), static_cast<int>(values.size()),
      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  (void)values;
#endif
}

void ChemicalField::sum_accounting_across_ranks() {
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return;
  const auto count = static_cast<int>(nspec_);
  if (mode_ != DecompositionMode::Slab) {
    MPI_Allreduce(
        MPI_IN_PLACE, flux_accounting_.delivery_reduction_step.data(),
        count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(
        MPI_IN_PLACE, flux_accounting_.delivery_retry_events_step.data(),
        count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(
        MPI_IN_PLACE, flux_accounting_.delivery_rationing_factor_step.data(),
        count, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(
        MPI_IN_PLACE, flux_accounting_.delivery_infeasible_step.data(),
        count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(
        MPI_IN_PLACE,
        flux_accounting_.delivery_axis_deferred_mass_step.data(),
        count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(
        MPI_IN_PLACE,
        flux_accounting_.delivery_axis_prescribed_mass_step.data(),
        count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    for (auto* values : {
             &flux_accounting_.negative_delivery_events_x_step,
             &flux_accounting_.negative_delivery_events_y_step,
             &flux_accounting_.negative_delivery_events_z_step}) {
      MPI_Allreduce(
          MPI_IN_PLACE, values->data(), count, MPI_DOUBLE, MPI_MAX,
          MPI_COMM_WORLD);
    }
    MPI_Allreduce(
        MPI_IN_PLACE, flux_accounting_.negative_delivery_min_step.data(),
        count, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return;
  }
  auto reduce = [count](std::vector<Real>& values) {
    MPI_Allreduce(MPI_IN_PLACE, values.data(), count, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
  };
  reduce(flux_accounting_.boundary_step);
  reduce(flux_accounting_.gradient_source_step);
  reduce(flux_accounting_.reaction_clip_step);
  reduce(flux_accounting_.negative_delivery_events_step);
  reduce(flux_accounting_.negative_delivery_mass_step);
  reduce(flux_accounting_.delivery_axis_deferred_mass_step);
  reduce(flux_accounting_.delivery_axis_prescribed_mass_step);
  reduce(flux_accounting_.negative_delivery_events_x_step);
  reduce(flux_accounting_.negative_delivery_events_y_step);
  reduce(flux_accounting_.negative_delivery_events_z_step);
  MPI_Allreduce(
      MPI_IN_PLACE, flux_accounting_.negative_delivery_min_step.data(),
      count, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  reduce(flux_accounting_.delivery_reduction_step);
  MPI_Allreduce(
      MPI_IN_PLACE, flux_accounting_.delivery_retry_events_step.data(),
      count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(
      MPI_IN_PLACE, flux_accounting_.delivery_rationing_factor_step.data(),
      count, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(
      MPI_IN_PLACE, flux_accounting_.delivery_infeasible_step.data(),
      count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}

namespace {

Int slab_storage_index(
    Int local_ix, Int iy, Int iz, Int storage_nx, Int ny) {
  return iz * storage_nx * ny + iy * storage_nx + local_ix;
}

Real owned_content(const std::vector<Real>& concentration,
                   const Domain& domain, Real cell_volume) {
  Real content = 0.0;
  for (Int cell = 0; cell < domain.ncells(); ++cell) {
    content += concentration[static_cast<size_t>(cell)] * cell_volume;
  }
  return content;
}

Real owned_content_slab(const std::vector<Real>& concentration,
                        const Domain& domain, Int storage_nx,
                        Int halo_width, Real cell_volume) {
  Real content = 0.0;
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.local_grid_nx(); ++ix) {
        const Int cell = slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, domain.ny());
        content += concentration[static_cast<size_t>(cell)] * cell_volume;
      }
    }
  }
  return content;
}

struct SlabTransportContext {
  std::vector<Real>& concentration;
  const Domain& domain;
  Int storage_nx;
  Int halo_width;
  Real alpha;
  const std::vector<Real>* sink_rate = nullptr;
  std::vector<Real>* realized = nullptr;
  const ChemicalSpec* gradient_spec = nullptr;
  Real sink_dt = 0.0;
  Real cell_volume = 0.0;
  const std::vector<Real>* prescribed_mass = nullptr;
  NutrientFluxAccounting* flux = nullptr;
  Int spec = -1;
};

struct SlabDeliveryLineContext {
  std::vector<Real>& concentration;
  std::vector<Real>& realized;
  const std::vector<Real>& sink_rate;
  const std::vector<Real>* prescribed_mass = nullptr;
  const Domain& domain;
  Int storage_nx = 0;
  Int halo_width = 0;
  Real alpha = 0.0;
  Real boundary_conc = 0.0;
  Real cell_volume = 0.0;
  Real sink_dt = 0.0;
  const ChemicalSpec* gradient_spec = nullptr;
  NutrientFluxAccounting* flux = nullptr;
  Int spec = -1;
};

void load_slab_periodic_x_delivery_line(
    std::vector<Real>& line, std::vector<Real>& sink,
    std::vector<Real>& prescribed, const SlabTransportContext& context,
    Int iy, Int iz) {
  for (Int ix = 0; ix < context.domain.nx(); ++ix) {
    const Int cell = slab_storage_index(
        context.halo_width + ix, iy, iz, context.storage_nx,
        context.domain.ny());
    const auto index = static_cast<size_t>(ix);
    line[index] = context.concentration[static_cast<size_t>(cell)];
    if (context.sink_rate != nullptr) {
      sink[index] = (*context.sink_rate)[static_cast<size_t>(cell)]
          * context.sink_dt;
    }
    if (context.prescribed_mass != nullptr) {
      prescribed[index] =
          (*context.prescribed_mass)[static_cast<size_t>(cell)]
          / (3.0 * context.cell_volume);
    }
  }
}

void solve_slab_delivery_line(
    const SlabDeliveryLineContext& context, Int ix, Int iy,
    Real& face_exchange) {
  const Int nz = context.domain.nz();
  const Int ny = context.domain.ny();
  std::vector<Real> line(static_cast<size_t>(nz - 1));
  std::vector<Real> sink(static_cast<size_t>(nz - 1));
  std::vector<Real> prescribed(static_cast<size_t>(nz - 1), 0.0);
  for (Int iz = 1; iz < nz; ++iz) {
    const Int cell = slab_storage_index(
        context.halo_width + ix, iy, iz, context.storage_nx, ny);
    const auto index = static_cast<size_t>(iz - 1);
    line[index] = context.concentration[static_cast<size_t>(cell)];
    sink[index] = context.sink_rate[static_cast<size_t>(cell)]
        * context.sink_dt;
    if (context.prescribed_mass != nullptr) {
      prescribed[index] =
          (*context.prescribed_mass)[static_cast<size_t>(cell)]
          / (3.0 * context.cell_volume);
    }
  }
  std::vector gradient(static_cast<size_t>(nz - 1), 0.0);
  fill_gradient_line(
      gradient, context.gradient_spec, context.domain, 1);
  cap_delivery_axis_prescribed(
      line, sink, gradient, 1.0, prescribed, context.flux,
      context.spec, context.cell_volume);
  apply_gradient_sink(
      line, sink, gradient, context.gradient_spec, context.domain, 1);
  for (size_t index = 0; index < line.size(); ++index) {
    line[index] -= prescribed[index];
  }
  std::vector diagonal(static_cast<size_t>(nz - 1), 0.0);
  for (Int iz = 1; iz < nz; ++iz) {
    diagonal[static_cast<size_t>(iz - 1)] =
        1.0 + 2.0 * context.alpha + sink[static_cast<size_t>(iz - 1)];
  }
  diagonal.back() = 1.0 + context.alpha + sink.back();
  solve_tridiagonal_with_diagonal(
      line, diagonal, context.alpha,
      context.alpha * context.boundary_conc);
  face_exchange += context.alpha
      * (context.boundary_conc - line.front()) * context.cell_volume;
  for (Int iz = 1; iz < nz; ++iz) {
    const Int cell = slab_storage_index(
        context.halo_width + ix, iy, iz, context.storage_nx, ny);
    const auto index = static_cast<size_t>(iz - 1);
    context.concentration[static_cast<size_t>(cell)] = line[index];
    const Real total = context.gradient_spec != nullptr
        ? line[index] + gradient[index] : line[index];
    context.realized[static_cast<size_t>(cell)] +=
        sink[index] * total * context.cell_volume;
    record_negative_delivery_excursion(
        context.flux, context.spec, total, sink[index],
        context.cell_volume, 2);
  }
}

void diffuse_periodic_x_slab_single(
    const SlabTransportContext& context, const PeriodicLineSolver& solver) {
  auto& concentration = context.concentration;
  const auto& domain = context.domain;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  const Real alpha = context.alpha;
  const auto* sink_rate = context.sink_rate;
  const auto* gradient_spec = context.gradient_spec;
  auto* realized = context.realized;
  const Real cell_volume = context.cell_volume;
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int line_count = ny * domain.nz();
  for (Int line_id = 0; line_id < line_count; ++line_id) {
    const Int iy = line_id % ny;
    const Int iz = line_id / ny;
    std::vector line(static_cast<size_t>(nx), 0.0);
    std::vector sink(static_cast<size_t>(nx), 0.0);
    std::vector prescribed(static_cast<size_t>(nx), 0.0);
    load_slab_periodic_x_delivery_line(
        line, sink, prescribed, context, iy, iz);
    if (sink_rate != nullptr || context.prescribed_mass != nullptr) {
      std::vector<Real> gradient;
      if (gradient_spec != nullptr) {
        gradient.assign(
            static_cast<size_t>(nx),
            z_gradient_reference(*gradient_spec, domain, iz));
      }
      cap_delivery_axis_prescribed(
          line, sink, gradient, 1.0, prescribed, context.flux,
          context.spec, cell_volume);
      solve_periodic_with_sink(
          line, sink, alpha,
          gradient_spec != nullptr ? &gradient : nullptr,
          context.prescribed_mass != nullptr ? &prescribed : nullptr);
    } else {
      solver.solve(line);
    }
    for (Int ix = 0; ix < nx; ++ix) {
      const Int cell = slab_storage_index(
          halo_width + ix, iy, iz, storage_nx, ny);
      concentration[static_cast<size_t>(cell)] = line[static_cast<size_t>(ix)];
      if (realized != nullptr) {
        const Real amount = gradient_spec != nullptr
            ? sink[static_cast<size_t>(ix)]
                * std::max(
                    line[static_cast<size_t>(ix)]
                        + z_gradient_reference(*gradient_spec, domain, iz),
                    0.0) * cell_volume
            : sink[static_cast<size_t>(ix)]
                * line[static_cast<size_t>(ix)] * cell_volume;
        (*realized)[static_cast<size_t>(cell)] += amount;
      }
    }
  }
}

#ifdef GUTIBM_MPI
struct SlabPeriodicXLayout {
  Int nx = 0;
  Int ny = 0;
  Int local_nx = 0;
  Int line_count = 0;
  Int process_count = 0;
  Int local_rank = 0;
  std::vector<Int> line_counts;
  std::vector<Int> line_displacements;
  std::vector<Int> x_counts;
  std::vector<Int> x_displacements;
  std::vector<Int> send_counts;
  std::vector<Int> send_displacements;
  std::vector<Int> recv_counts;
  std::vector<Int> recv_displacements;
};

struct SlabPeriodicXBuffers {
  std::vector<Real> send;
  std::vector<Real> recv;
  std::vector<Real> solved;
  std::vector<Real> sink_send;
  std::vector<Real> sink_recv;
  std::vector<Real> prescribed_send;
  std::vector<Real> prescribed_recv;
};

SlabPeriodicXLayout make_slab_periodic_x_layout(
    const Domain& domain, Int nx) {
  SlabPeriodicXLayout layout;
  layout.nx = nx;
  layout.ny = domain.ny();
  layout.local_nx = domain.local_grid_nx();
  layout.line_count = layout.ny * domain.nz();
  layout.process_count = domain.nprocs();
  layout.local_rank = domain.rank();
  const auto count = static_cast<size_t>(layout.process_count);
  layout.line_counts.assign(count, 0);
  layout.line_displacements.assign(count, 0);
  layout.x_counts.assign(count, 0);
  layout.x_displacements.assign(count, 0);
  for (Int rank = 0; rank < layout.process_count; ++rank) {
    const auto [begin, end] =
        Domain::grid_x_range_for_rank(nx, layout.process_count, rank);
    layout.x_counts[static_cast<size_t>(rank)] = end - begin;
    layout.x_displacements[static_cast<size_t>(rank)] = begin;
    layout.line_counts[static_cast<size_t>(rank)] =
        (layout.line_count + layout.process_count - rank - 1)
        / layout.process_count;
    if (rank > 0) {
      layout.line_displacements[static_cast<size_t>(rank)] =
          layout.line_displacements[static_cast<size_t>(rank - 1)]
          + layout.line_counts[static_cast<size_t>(rank - 1)];
    }
  }
  layout.send_counts.assign(count, 0);
  layout.send_displacements.assign(count, 0);
  layout.recv_counts.assign(count, 0);
  layout.recv_displacements.assign(count, 0);
  for (Int rank = 0; rank < layout.process_count; ++rank) {
    layout.send_counts[static_cast<size_t>(rank)] =
        layout.line_counts[static_cast<size_t>(rank)] * layout.local_nx;
    layout.recv_counts[static_cast<size_t>(rank)] =
        layout.line_counts[static_cast<size_t>(layout.local_rank)]
        * layout.x_counts[static_cast<size_t>(rank)];
    if (rank > 0) {
      layout.send_displacements[static_cast<size_t>(rank)] =
          layout.send_displacements[static_cast<size_t>(rank - 1)]
          + layout.send_counts[static_cast<size_t>(rank - 1)];
      layout.recv_displacements[static_cast<size_t>(rank)] =
          layout.recv_displacements[static_cast<size_t>(rank - 1)]
          + layout.recv_counts[static_cast<size_t>(rank - 1)];
    }
  }
  return layout;
}

void pack_slab_periodic_x_line(
    const SlabTransportContext& context,
    const SlabPeriodicXLayout& layout, Int iy, Int iz, Int offset,
    SlabPeriodicXBuffers& buffers) {
  const auto& concentration = context.concentration;
  const auto* sink_rate = context.sink_rate;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  for (Int ix = 0; ix < layout.local_nx; ++ix) {
    const Int cell = slab_storage_index(
        halo_width + ix, iy, iz, storage_nx, layout.ny);
    buffers.send[static_cast<size_t>(offset + ix)] =
        concentration[static_cast<size_t>(cell)];
    if (sink_rate != nullptr) {
      buffers.sink_send[static_cast<size_t>(offset + ix)] =
          (*sink_rate)[static_cast<size_t>(cell)] * context.sink_dt;
    }
    if (context.prescribed_mass != nullptr
        && delivery_axis_cap_enabled()) {
      buffers.prescribed_send[static_cast<size_t>(offset + ix)] =
          (*context.prescribed_mass)[static_cast<size_t>(cell)]
          / (3.0 * context.cell_volume);
    }
  }
}

void pack_slab_periodic_x(
    const SlabTransportContext& context,
    const SlabPeriodicXLayout& layout, SlabPeriodicXBuffers& buffers) {
  for (Int destination = 0; destination < layout.process_count;
       ++destination) {
    Int offset = layout.send_displacements[static_cast<size_t>(destination)];
    for (Int line_id = destination; line_id < layout.line_count;
         line_id += layout.process_count) {
      pack_slab_periodic_x_line(
          context, layout, line_id % layout.ny, line_id / layout.ny, offset,
          buffers);
      offset += layout.local_nx;
    }
  }
}

void solve_slab_periodic_x_lines(
    const SlabTransportContext& context,
    const SlabPeriodicXLayout& layout, const PeriodicLineSolver& solver,
    SlabPeriodicXBuffers& buffers) {
  const auto* sink_rate = context.sink_rate;
  const auto* gradient_spec = context.gradient_spec;
  const bool axis_cap_enabled = delivery_axis_cap_enabled();
  std::vector sink(static_cast<size_t>(layout.nx), 0.0);
  const auto& gathered_displacements = layout.recv_displacements;
  for (Int line_index = 0;
       line_index < layout.line_counts[static_cast<size_t>(layout.local_rank)];
       ++line_index) {
    std::vector line(static_cast<size_t>(layout.nx), 0.0);
    for (Int source = 0; source < layout.process_count; ++source) {
      const Int segment_offset =
          gathered_displacements[static_cast<size_t>(source)]
          + line_index * layout.x_counts[static_cast<size_t>(source)];
      std::copy_n(
          buffers.recv.begin() + segment_offset,
          layout.x_counts[static_cast<size_t>(source)],
          line.begin() + layout.x_displacements[static_cast<size_t>(source)]);
      if (sink_rate != nullptr) {
        std::copy_n(
            buffers.sink_recv.begin() + segment_offset,
            layout.x_counts[static_cast<size_t>(source)],
            sink.begin() + layout.x_displacements[static_cast<size_t>(source)]);
      }
    }
    if (sink_rate != nullptr) {
      std::vector<Real> gradient(static_cast<size_t>(layout.nx), 0.0);
      std::vector<Real> prescribed;
      if (gradient_spec != nullptr) {
        const Int line_id = layout.local_rank
            + line_index * layout.process_count;
        const Int iz = line_id / layout.ny;
        std::ranges::fill(
            gradient,
            z_gradient_reference(*gradient_spec, context.domain, iz));
      }
      if (context.prescribed_mass != nullptr && axis_cap_enabled) {
        prescribed.assign(static_cast<size_t>(layout.nx), 0.0);
        for (Int source = 0; source < layout.process_count; ++source) {
          const Int segment_offset =
              gathered_displacements[static_cast<size_t>(source)]
              + line_index * layout.x_counts[static_cast<size_t>(source)];
          std::copy_n(
              buffers.prescribed_recv.begin() + segment_offset,
              layout.x_counts[static_cast<size_t>(source)],
              prescribed.begin()
                  + layout.x_displacements[static_cast<size_t>(source)]);
        }
        cap_delivery_axis_prescribed(
            line, sink, gradient, 1.0, prescribed, context.flux,
            context.spec, context.cell_volume);
      }
      solve_periodic_with_sink(
          line, sink, context.alpha,
          gradient_spec != nullptr ? &gradient : nullptr,
          context.prescribed_mass != nullptr && axis_cap_enabled
              ? &prescribed : nullptr);
    } else {
      solver.solve(line);
    }
    for (Int destination = 0; destination < layout.process_count;
         ++destination) {
      const Int segment_offset =
          gathered_displacements[static_cast<size_t>(destination)]
          + line_index * layout.x_counts[static_cast<size_t>(destination)];
      std::copy_n(
          line.begin() + layout.x_displacements[static_cast<size_t>(destination)],
          layout.x_counts[static_cast<size_t>(destination)],
          buffers.solved.begin() + segment_offset);
    }
  }
}

void store_slab_periodic_x_line(
    const SlabTransportContext& context,
    const SlabPeriodicXLayout& layout, Int iy, Int iz, Int offset,
    const SlabPeriodicXBuffers& buffers) {
  auto& concentration = context.concentration;
  auto* realized = context.realized;
  const auto* sink_rate = context.sink_rate;
  const auto* gradient_spec = context.gradient_spec;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  for (Int ix = 0; ix < layout.local_nx; ++ix) {
    const Int cell = slab_storage_index(
        halo_width + ix, iy, iz, storage_nx, layout.ny);
    concentration[static_cast<size_t>(cell)] =
        buffers.recv[static_cast<size_t>(offset + ix)];
    if (realized != nullptr && sink_rate != nullptr) {
      const Real total = gradient_spec != nullptr
          ? concentration[static_cast<size_t>(cell)]
              + z_gradient_reference(*gradient_spec, context.domain, iz)
          : concentration[static_cast<size_t>(cell)];
      (*realized)[static_cast<size_t>(cell)] +=
          (*sink_rate)[static_cast<size_t>(cell)] * context.sink_dt
          * total
          * context.cell_volume;
      record_negative_delivery_excursion(
          context.flux, context.spec, total,
          (*sink_rate)[static_cast<size_t>(cell)] * context.sink_dt,
          context.cell_volume, 0);
    }
  }
}

void diffuse_periodic_x_slab_mpi(
    const SlabTransportContext& context, const PeriodicLineSolver& solver) {
  const auto& domain = context.domain;
  const SlabPeriodicXLayout layout =
      make_slab_periodic_x_layout(domain, domain.nx());
  SlabPeriodicXBuffers buffers;
  const auto solved_total = std::accumulate(
      layout.recv_counts.begin(), layout.recv_counts.end(), 0);
  const auto output_total = std::accumulate(
      layout.send_counts.begin(), layout.send_counts.end(), 0);
  buffers.send.resize(static_cast<size_t>(output_total));
  buffers.recv.resize(static_cast<size_t>(solved_total));
  buffers.solved.resize(static_cast<size_t>(solved_total));
  if (context.sink_rate != nullptr) {
    buffers.sink_send.resize(static_cast<size_t>(output_total));
    buffers.sink_recv.resize(static_cast<size_t>(solved_total));
  }
  if (context.prescribed_mass != nullptr
      && delivery_axis_cap_enabled()) {
    buffers.prescribed_send.resize(static_cast<size_t>(output_total));
    buffers.prescribed_recv.resize(static_cast<size_t>(solved_total));
  }
  pack_slab_periodic_x(context, layout, buffers);
  MPI_Alltoallv(
      buffers.send.data(), layout.send_counts.data(),
      layout.send_displacements.data(), MPI_DOUBLE, buffers.recv.data(),
      layout.recv_counts.data(), layout.recv_displacements.data(), MPI_DOUBLE,
      MPI_COMM_WORLD);
  if (context.sink_rate != nullptr) {
    MPI_Alltoallv(
        buffers.sink_send.data(), layout.send_counts.data(),
        layout.send_displacements.data(), MPI_DOUBLE, buffers.sink_recv.data(),
        layout.recv_counts.data(), layout.recv_displacements.data(), MPI_DOUBLE,
        MPI_COMM_WORLD);
  }
  if (context.prescribed_mass != nullptr
      && delivery_axis_cap_enabled()) {
    MPI_Alltoallv(
        buffers.prescribed_send.data(), layout.send_counts.data(),
        layout.send_displacements.data(), MPI_DOUBLE,
        buffers.prescribed_recv.data(), layout.recv_counts.data(),
        layout.recv_displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  }
  solve_slab_periodic_x_lines(context, layout, solver, buffers);
  const auto& output_counts = layout.recv_counts;
  const auto& output_displacements = layout.recv_displacements;
  const auto& solved_counts = layout.send_counts;
  const auto& solved_displacements = layout.send_displacements;
  const auto& gathered_displacements = layout.recv_displacements;
  buffers.send.assign(
      static_cast<size_t>(std::accumulate(
          output_counts.begin(), output_counts.end(), 0)), 0.0);
  buffers.recv.assign(
      static_cast<size_t>(std::accumulate(
          solved_counts.begin(), solved_counts.end(), 0)), 0.0);
  for (Int destination = 0; destination < layout.process_count;
       ++destination) {
    const Int segment_length =
        layout.x_counts[static_cast<size_t>(destination)];
    Int offset = output_displacements[static_cast<size_t>(destination)];
    for (Int line_index = 0;
         line_index < layout.line_counts[static_cast<size_t>(layout.local_rank)];
         ++line_index) {
      const Int source_offset =
          gathered_displacements[static_cast<size_t>(destination)]
          + line_index * segment_length;
      std::copy_n(
          buffers.solved.begin() + source_offset, segment_length,
          buffers.send.begin() + offset);
      offset += segment_length;
    }
  }
  MPI_Alltoallv(
      buffers.send.data(), output_counts.data(), output_displacements.data(),
      MPI_DOUBLE, buffers.recv.data(), solved_counts.data(),
      solved_displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  for (Int source = 0; source < layout.process_count; ++source) {
    Int offset = layout.send_displacements[static_cast<size_t>(source)];
    for (Int line_id = source; line_id < layout.line_count;
         line_id += layout.process_count) {
      store_slab_periodic_x_line(
          context, layout, line_id % layout.ny, line_id / layout.ny, offset,
          buffers);
      offset += layout.local_nx;
    }
  }
}
#endif

void diffuse_periodic_x_slab(
    const SlabTransportContext& context) {
  const PeriodicLineSolver solver(
      context.domain.nx(), context.alpha);
  if (context.domain.nprocs() == 1) {
    diffuse_periodic_x_slab_single(context, solver);
    return;
  }
#ifdef GUTIBM_MPI
  diffuse_periodic_x_slab_mpi(context, solver);
#else
  (void)solver;
#endif
}

void diffuse_periodic_y_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width, Real alpha) {
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const Int local_nx = domain.local_grid_nx();
  const PeriodicLineSolver solver(ny, alpha);
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int ix = 0; ix < local_nx; ++ix) {
      std::vector<Real> line(static_cast<size_t>(ny));
      for (Int iy = 0; iy < ny; ++iy) {
        line[static_cast<size_t>(iy)] = concentration[
            static_cast<size_t>(slab_storage_index(
                halo_width + ix, iy, iz, storage_nx, ny))];
      }
      solver.solve(line);
      for (Int iy = 0; iy < ny; ++iy) {
        concentration[static_cast<size_t>(slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny))] =
            line[static_cast<size_t>(iy)];
      }
    }
  }
}

void diffuse_periodic_y_slab_delivery(
    const SlabTransportContext& context) {
  auto& concentration = context.concentration;
  const auto& domain = context.domain;
  const auto& sink_rate = *context.sink_rate;
  auto& realized = *context.realized;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  const Real alpha = context.alpha;
  const Real sink_dt = context.sink_dt;
  const Real cell_volume = context.cell_volume;
  const auto* gradient_spec = context.gradient_spec;
  const Int nx = domain.local_grid_nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int ix = 0; ix < nx; ++ix) {
      std::vector<Real> line(static_cast<size_t>(ny));
      std::vector<Real> sink(static_cast<size_t>(ny));
      std::vector<Real> prescribed(static_cast<size_t>(ny), 0.0);
      for (Int iy = 0; iy < ny; ++iy) {
        const Int cell = slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny);
        line[static_cast<size_t>(iy)] = concentration[static_cast<size_t>(cell)];
        sink[static_cast<size_t>(iy)] =
            sink_rate[static_cast<size_t>(cell)] * sink_dt;
        if (context.prescribed_mass != nullptr) {
          prescribed[static_cast<size_t>(iy)] =
              (*context.prescribed_mass)[static_cast<size_t>(cell)]
              / (3.0 * cell_volume);
        }
      }
      std::vector<Real> gradient(static_cast<size_t>(ny), 0.0);
      if (gradient_spec != nullptr) {
        std::ranges::fill(
            gradient, z_gradient_reference(*gradient_spec, domain, iz));
      }
      cap_delivery_axis_prescribed(
          line, sink, gradient, 1.0, prescribed, context.flux,
          context.spec, cell_volume);
      solve_periodic_with_sink(
          line, sink, alpha,
          gradient_spec != nullptr ? &gradient : nullptr,
          context.prescribed_mass != nullptr ? &prescribed : nullptr);
      for (Int iy = 0; iy < ny; ++iy) {
        const Int cell = slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny);
        concentration[static_cast<size_t>(cell)] = line[static_cast<size_t>(iy)];
        const Real total = gradient_spec != nullptr
            ? line[static_cast<size_t>(iy)]
                + gradient[static_cast<size_t>(iy)]
            : line[static_cast<size_t>(iy)];
        realized[static_cast<size_t>(cell)] +=
            sink[static_cast<size_t>(iy)] * total * cell_volume;
        record_negative_delivery_excursion(
            context.flux, context.spec, total,
            sink[static_cast<size_t>(iy)], cell_volume, 1);
      }
    }
  }
}

Real diffuse_bounded_z_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width, Real alpha, Real boundary_conc,
    Real cell_volume) {
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const Int local_nx = domain.local_grid_nx();
  if (nz <= 1) return 0.0;
  const NeumannTopLineSolver solver(nz - 1, alpha);
  Real face_exchange = 0.0;
  for (Int iy = 0; iy < ny; ++iy) {
    for (Int ix = 0; ix < local_nx; ++ix) {
      std::vector<Real> line(static_cast<size_t>(nz - 1));
      for (Int iz = 1; iz < nz; ++iz) {
        line[static_cast<size_t>(iz - 1)] = concentration[
            static_cast<size_t>(slab_storage_index(
                halo_width + ix, iy, iz, storage_nx, ny))];
      }
      solver.solve(line, boundary_conc);
      face_exchange += alpha * (boundary_conc - line.front()) * cell_volume;
      for (Int iz = 1; iz < nz; ++iz) {
        concentration[static_cast<size_t>(slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, ny))] =
            line[static_cast<size_t>(iz - 1)];
      }
    }
  }
  return face_exchange;
}

Real diffuse_bounded_z_slab_delivery(
    const SlabTransportContext& context, Real boundary_conc) {
  auto& concentration = context.concentration;
  const auto& domain = context.domain;
  const auto& sink_rate = *context.sink_rate;
  auto& realized = *context.realized;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  const Real alpha = context.alpha;
  const Real cell_volume = context.cell_volume;
  const Real sink_dt = context.sink_dt;
  const auto* gradient_spec = context.gradient_spec;
  const Int nx = domain.local_grid_nx();
  const Int ny = domain.ny();
  if (const Int nz = domain.nz(); nz <= 1) return 0.0;
  Real face_exchange = 0.0;
  const SlabDeliveryLineContext line_context{
      concentration, realized, sink_rate, context.prescribed_mass, domain,
      storage_nx, halo_width, alpha, boundary_conc, cell_volume, sink_dt,
      gradient_spec, context.flux, context.spec};
  for (Int iy = 0; iy < ny; ++iy) {
    for (Int ix = 0; ix < nx; ++ix) {
      solve_slab_delivery_line(line_context, ix, iy, face_exchange);
    }
  }
  return face_exchange;
}

Real diffuse_bounded_z_delivery_slab(
    std::vector<Real>& concentration, const Domain& domain, Int storage_nx,
    Int halo_width, Real alpha,
    const DeliveryBoundaryParameters& params) {
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const Int local_nx = domain.local_grid_nx();
  const auto load_line = [&concentration, storage_nx, halo_width, ny, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] = concentration[
          static_cast<size_t>(slab_storage_index(
              halo_width + ix, iy, iz, storage_nx, ny))];
    }
  };
  const auto store_line = [&concentration, storage_nx, halo_width, ny, nz](
                              Int ix, Int iy,
                              const std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      concentration[static_cast<size_t>(slab_storage_index(
          halo_width + ix, iy, iz, storage_nx, ny))] =
          line[static_cast<size_t>(iz)];
    }
  };
  return diffuse_bounded_z_delivery_impl(
      local_nx, ny, nz, alpha, params, load_line, store_line);
}

Real diffuse_bounded_z_delivery_with_sink_slab(
    const SlabTransportContext& context,
    const DeliveryBoundaryParameters& params) {
  auto& concentration = context.concentration;
  const auto& sink_rate = *context.sink_rate;
  auto& realized = *context.realized;
  const auto& domain = context.domain;
  const Int storage_nx = context.storage_nx;
  const Int halo_width = context.halo_width;
  const Real alpha = context.alpha;
  const Real sink_dt = context.sink_dt;
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const Int local_nx = domain.local_grid_nx();
  const auto load_line = [&concentration, storage_nx, halo_width, ny, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] = concentration[
          static_cast<size_t>(slab_storage_index(
              halo_width + ix, iy, iz, storage_nx, ny))];
    }
  };
  const auto store_line = [&concentration, storage_nx, halo_width, ny](
                              Int ix, Int iy, Int iz, Real value) {
    concentration[static_cast<size_t>(slab_storage_index(
        halo_width + ix, iy, iz, storage_nx, ny))] = value;
  };
  const auto load_sink = [&sink_rate, storage_nx, halo_width, ny, nz](
                             Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] = sink_rate[
          static_cast<size_t>(slab_storage_index(
              halo_width + ix, iy, iz, storage_nx, ny))];
    }
  };
  const auto load_prescribed = [
      &context, storage_nx, halo_width, ny, nz](
      Int ix, Int iy, std::vector<Real>& line) {
    for (Int iz = 0; iz < nz; ++iz) {
      const Int cell = slab_storage_index(
          halo_width + ix, iy, iz, storage_nx, ny);
      line[static_cast<size_t>(iz)] =
          context.prescribed_mass == nullptr
              ? 0.0
              : (*context.prescribed_mass)[static_cast<size_t>(cell)]
                    / (3.0 * context.cell_volume);
    }
  };
  const auto* gradient_spec = context.gradient_spec;
  const auto load_profile = [gradient_spec, &domain, nz](
                                 Int, Int, std::vector<Real>& line) {
    if (gradient_spec == nullptr) return false;
    for (Int iz = 0; iz < nz; ++iz) {
      line[static_cast<size_t>(iz)] =
          z_gradient_reference(*gradient_spec, domain, iz);
    }
    return true;
  };
  const auto add_realized = [
      &realized, storage_nx, halo_width, ny, &context](
      Int ix, Int iy, Int iz, Real amount, Real total) {
    realized[static_cast<size_t>(slab_storage_index(
        halo_width + ix, iy, iz, storage_nx, ny))] += amount;
    record_negative_delivery_excursion(
        context.flux, context.spec, -amount, 2);
    if (context.flux != nullptr && total < 0.0) {
      context.flux->add_negative_delivery_min(context.spec, total);
    }
  };
  const DeliveryLineOperations<
      decltype(load_line), decltype(store_line), decltype(load_sink),
      decltype(load_prescribed), decltype(load_profile),
      decltype(add_realized)>
      operations{load_line, store_line, load_sink, load_prescribed,
                 load_profile, add_realized, sink_dt,
                 context.flux, context.spec};
  return diffuse_bounded_z_delivery_with_sink_impl(
      {local_nx, ny, nz, alpha}, params, operations);
}

Real set_epithelial_boundary_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width, Real boundary_conc, Real cell_volume) {
  Real amount = 0.0;
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.local_grid_nx(); ++ix) {
      const Int index = slab_storage_index(
          halo_width + ix, iy, 0, storage_nx, domain.ny());
      amount += (boundary_conc - concentration[static_cast<size_t>(index)])
          * cell_volume;
      concentration[static_cast<size_t>(index)] = boundary_conc;
    }
  }
  return amount;
}

void set_luminal_neumann_boundary_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width) {
  if (domain.nz() < 2) return;
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.local_grid_nx(); ++ix) {
      const Int top = slab_storage_index(
          halo_width + ix, iy, domain.nz() - 1, storage_nx, domain.ny());
      const Int below = slab_storage_index(
          halo_width + ix, iy, domain.nz() - 2, storage_nx, domain.ny());
      concentration[static_cast<size_t>(top)] =
          concentration[static_cast<size_t>(below)];
    }
  }
}

void shift_z_gradient_slab(
    std::vector<Real>& concentration, const ChemicalSpec& spec,
    const Domain& domain, Int storage_nx, Int halo_width, Real scale) {
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    const Int profile_iz = domain.nz() >= 2 && iz == domain.nz() - 1
        ? domain.nz() - 2 : iz;
  const Real z_rel = (profile_iz + 0.5) * domain.dx_z();
    const Real shift = scale * spec.initial_conc
        * std::exp(-z_rel / spec.z_gradient_lambda);
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.local_grid_nx(); ++ix) {
        concentration[static_cast<size_t>(slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, domain.ny()))] += shift;
      }
    }
  }
}

void clamp_nonnegative_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width) {
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.local_grid_nx(); ++ix) {
        const Int index = slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, domain.ny());
        concentration[static_cast<size_t>(index)] =
            std::max(concentration[static_cast<size_t>(index)], 0.0);
      }
    }
  }
}

void clamp_delivery_total(
    std::vector<Real>& concentration, Real cell_volume,
    NutrientFluxAccounting& flux, Int spec) {
  Real clipped = 0.0;
  for (Real& value : concentration) {
    if (value < 0.0) {
      clipped -= value * cell_volume;
      value = 0.0;
    }
  }
  if (clipped > 0.0) flux.add_reaction_clip(spec, clipped);
}

Real clamp_delivery_value(Real& value, Real cell_volume) {
  if (value >= 0.0) return 0.0;
  const Real clipped = -value * cell_volume;
  value = 0.0;
  return clipped;
}

void clamp_delivery_total_slab(
    std::vector<Real>& concentration, const Domain& domain,
    Int storage_nx, Int halo_width, Real cell_volume,
    NutrientFluxAccounting& flux, Int spec) {
  Real clipped = 0.0;
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.local_grid_nx(); ++ix) {
        const Int index = slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, domain.ny());
        Real& value = concentration[static_cast<size_t>(index)];
        clipped += clamp_delivery_value(value, cell_volume);
      }
    }
  }
  if (clipped > 0.0) flux.add_reaction_clip(spec, clipped);
}

struct ReplicatedDiffusionContext {
  std::vector<Real>& concentration;
  const Domain& domain;
  const ChemicalSpec& chemical;
  NutrientFluxAccounting& flux;
  std::vector<Real>& sink_rate;
  std::vector<Real>& sink_realized;
  Int spec;
  Real dt;
  Real alpha_x;
  Real alpha_y;
  Real alpha_z;
  Real cell_volume;
  Real diffusion_boundary;
  bool preserve_gradient;
  bool delivery;
  const ChemicalSpec* gradient_spec;
  const std::vector<Real>* prescribed_mass;
  bool prescribed_active;
};

struct DebugStageSnapshot {
  Real content = 0.0;
  Real sink_realized = 0.0;
  Real boundary = 0.0;
};

bool debug_stage_enabled(const ChemicalSpec& chemical) {
  return nutrient_debug_enabled()
      && nutrient_debug_step_counter() == nutrient_debug_transport_step()
      && (chemical.name == species::OXYGEN
          || chemical.name == species::CARBON);
}

Real sum_realized(const std::vector<Real>& realized) {
  return std::accumulate(realized.begin(), realized.end(), 0.0);
}

DebugStageSnapshot debug_snapshot(
    const ReplicatedDiffusionContext& context) {
  if (!debug_stage_enabled(context.chemical)) return {};
  DebugStageSnapshot snapshot{
      owned_content(context.concentration, context.domain,
                    context.cell_volume),
      sum_realized(context.sink_realized),
      context.flux.boundary_step[static_cast<size_t>(context.spec)]};
  return snapshot;
}

void print_debug_stage(
    const ChemicalSpec& chemical, const Domain& domain,
    std::string_view stage, const DebugStageSnapshot& before,
    const DebugStageSnapshot& after) {
  if (domain.rank() != 0) return;
  std::cout << "NUTRIENT_STAGE step=" << nutrient_debug_step_counter()
            << " species=" << chemical.name
            << " stage=" << stage
            << " before_content=" << debug_real(before.content)
            << " after_content=" << debug_real(after.content)
            << " content_delta=" << debug_real(
                   after.content - before.content)
            << " sink_realized_delta="
            << debug_real(after.sink_realized - before.sink_realized)
            << " boundary_delta=" << debug_real(
                   after.boundary - before.boundary)
            << '\n';
}

void emit_debug_stage(
    const ReplicatedDiffusionContext& context, std::string_view stage,
    const DebugStageSnapshot& before) {
  if (!debug_stage_enabled(context.chemical)) return;
  const DebugStageSnapshot after = debug_snapshot(context);
  print_debug_stage(
      context.chemical, context.domain, stage, before, after);
}

void prepare_replicated_diffusion(ReplicatedDiffusionContext& context) {
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    const DebugStageSnapshot before = debug_snapshot(context);
    const Real boundary = set_epithelial_boundary(
        context.concentration, context.domain,
        context.chemical.boundary_conc, context.cell_volume);
    context.flux.add_boundary(
        context.spec, boundary);
    emit_debug_stage(context, "epithelial_dirichlet_pin_prepare", before);
  }
  if (!context.preserve_gradient) return;
  const Real before_luminal = owned_content(
      context.concentration, context.domain, context.cell_volume);
  const DebugStageSnapshot luminal_snapshot = debug_snapshot(context);
  set_luminal_neumann_boundary(context.concentration, context.domain);
  context.flux.add_gradient_source(
      context.spec,
      owned_content(context.concentration, context.domain, context.cell_volume)
          - before_luminal);
  emit_debug_stage(context, "luminal_neumann_prepare", luminal_snapshot);
  const Real before_profile = owned_content(
      context.concentration, context.domain, context.cell_volume);
  const DebugStageSnapshot profile_snapshot = debug_snapshot(context);
  shift_z_gradient(
      context.concentration, context.chemical, context.domain, -1.0);
  context.flux.add_gradient_source(
      context.spec,
      owned_content(context.concentration, context.domain, context.cell_volume)
          - before_profile);
  emit_debug_stage(context, "profile_subtract", profile_snapshot);
  context.diffusion_boundary = 0.0;
}

void transport_replicated_periodic(
    const ReplicatedDiffusionContext& context) {
  if (context.delivery) {
    const DeliverySinkParameters sink_params{
        context.sink_rate, context.sink_realized, context.dt / 3.0,
        context.cell_volume,
        context.preserve_gradient ? &context.chemical : nullptr,
        context.prescribed_mass, &context.flux, context.spec};
    const DebugStageSnapshot before_x = debug_snapshot(context);
    diffuse_periodic_x_delivery(
        context.concentration, context.domain, context.alpha_x, sink_params);
    emit_debug_stage(context, "x_delivery_solve", before_x);
    const DebugStageSnapshot before_y = debug_snapshot(context);
    diffuse_periodic_y_delivery(
        context.concentration, context.domain, context.alpha_y, sink_params);
    emit_debug_stage(context, "y_delivery_solve", before_y);
  } else {
    diffuse_periodic_x(
        context.concentration, context.domain, context.alpha_x);
    diffuse_periodic_y(
        context.concentration, context.domain, context.alpha_y);
  }
}

Real transport_replicated_z(
    const ReplicatedDiffusionContext& context) {
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    if (context.delivery) {
      return diffuse_bounded_z_delivery(
          context.concentration, context.domain, context.alpha_z,
          {context.diffusion_boundary, 0.0, 0.0,
           EpithelialBoundaryMode::Dirichlet, context.cell_volume},
          {context.sink_rate, context.sink_realized, context.dt / 3.0,
           context.cell_volume,
           context.preserve_gradient ? &context.chemical : nullptr,
           context.prescribed_mass, &context.flux, context.spec});
    }
    return diffuse_bounded_z(
        context.concentration, context.domain, context.alpha_z,
        context.diffusion_boundary, context.cell_volume);
  }
  const Real beta = context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Robin
      ? context.chemical.epithelial_transfer_coeff * context.dt
          / context.domain.dx_z() : 0.0;
  const Real flux_source = context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Flux
      ? context.chemical.epithelial_flux * context.dt
          / context.domain.dx_z() : 0.0;
  if (context.delivery) {
    return diffuse_bounded_z_delivery_with_sink(
        context.concentration, context.domain, context.alpha_z,
        {context.chemical.boundary_conc, beta, flux_source,
         context.chemical.epithelial_boundary_mode, context.cell_volume},
        {context.sink_rate, context.sink_realized, context.dt / 3.0,
         context.cell_volume,
         context.preserve_gradient ? &context.chemical : nullptr,
         context.prescribed_mass, &context.flux, context.spec});
  }
  return diffuse_bounded_z_delivery(
      context.concentration, context.domain, context.alpha_z,
      {context.chemical.boundary_conc, beta, flux_source,
       context.chemical.epithelial_boundary_mode, context.cell_volume});
}

void transport_replicated_diffusion(
    ReplicatedDiffusionContext& context) {
  transport_replicated_periodic(context);
  const DebugStageSnapshot before_z = debug_snapshot(context);
  const Real boundary = transport_replicated_z(context);
  context.flux.add_boundary(
      context.spec, boundary);
  emit_debug_stage(context, "z_delivery_solve", before_z);
}

void finish_replicated_diffusion(ReplicatedDiffusionContext& context) {
  if (context.preserve_gradient) {
    const Real before_profile = owned_content(
        context.concentration, context.domain, context.cell_volume);
    const DebugStageSnapshot profile_snapshot = debug_snapshot(context);
    shift_z_gradient(
        context.concentration, context.chemical, context.domain, 1.0);
    context.flux.add_gradient_source(
        context.spec,
        owned_content(context.concentration, context.domain, context.cell_volume)
            - before_profile);
    emit_debug_stage(context, "profile_restore", profile_snapshot);
    const Real before_luminal = owned_content(
        context.concentration, context.domain, context.cell_volume);
    const DebugStageSnapshot luminal_snapshot = debug_snapshot(context);
    set_luminal_neumann_boundary(context.concentration, context.domain);
    context.flux.add_gradient_source(
        context.spec,
        owned_content(context.concentration, context.domain, context.cell_volume)
            - before_luminal);
    emit_debug_stage(context, "luminal_neumann_finish", luminal_snapshot);
  }
  if (!context.delivery) {
    clamp_nonnegative(context.concentration);
  } else if (!context.prescribed_active) {
    clamp_delivery_total(
        context.concentration, context.cell_volume,
        context.flux, context.spec);
  }
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    const DebugStageSnapshot before = debug_snapshot(context);
    const Real boundary = set_epithelial_boundary(
        context.concentration, context.domain,
        context.chemical.boundary_conc, context.cell_volume);
    context.flux.add_boundary(
        context.spec, boundary);
    emit_debug_stage(context, "epithelial_dirichlet_pin_finish", before);
  }
}

struct SlabDiffusionContext {
  std::vector<Real>& concentration;
  const Domain& domain;
  const ChemicalSpec& chemical;
  NutrientFluxAccounting& flux;
  std::vector<Real>& sink_rate;
  std::vector<Real>& sink_realized;
  Int storage_nx;
  Int halo_width;
  Int spec;
  Real dt;
  Real alpha_x;
  Real alpha_y;
  Real alpha_z;
  Real cell_volume;
  Real diffusion_boundary;
  bool preserve_gradient;
  bool delivery;
  const ChemicalSpec* gradient_spec;
  const std::vector<Real>* prescribed_mass;
  bool prescribed_active;
};

Real sum_realized(
    const std::vector<Real>& realized, const Domain& domain,
    Int storage_nx, Int halo_width) {
  Real total = 0.0;
  for (Int iz = 0; iz < domain.nz(); ++iz) {
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.local_grid_nx(); ++ix) {
        total += realized[static_cast<size_t>(slab_storage_index(
            halo_width + ix, iy, iz, storage_nx, domain.ny()))];
      }
    }
  }
  return total;
}

DebugStageSnapshot debug_snapshot(const SlabDiffusionContext& context) {
  if (!debug_stage_enabled(context.chemical)) return {};
  DebugStageSnapshot snapshot{
      owned_content_slab(
          context.concentration, context.domain, context.storage_nx,
          context.halo_width, context.cell_volume),
      sum_realized(
          context.sink_realized, context.domain, context.storage_nx,
          context.halo_width),
      context.flux.boundary_step[static_cast<size_t>(context.spec)]};
#ifdef GUTIBM_MPI
  if (context.domain.nprocs() > 1) {
    std::array<Real, 3> values{
        snapshot.content, snapshot.sink_realized, snapshot.boundary};
    MPI_Allreduce(MPI_IN_PLACE, values.data(), 3, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    snapshot = {values[0], values[1], values[2]};
  }
#endif
  return snapshot;
}

void emit_debug_stage(
    const SlabDiffusionContext& context, std::string_view stage,
    const DebugStageSnapshot& before) {
  if (!debug_stage_enabled(context.chemical)) return;
  const DebugStageSnapshot after = debug_snapshot(context);
  print_debug_stage(
      context.chemical, context.domain, stage, before, after);
}

void prepare_slab_diffusion(SlabDiffusionContext& context) {
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    const DebugStageSnapshot before = debug_snapshot(context);
    const Real boundary = set_epithelial_boundary_slab(
        context.concentration, context.domain, context.storage_nx,
        context.halo_width, context.chemical.boundary_conc,
        context.cell_volume);
    context.flux.add_boundary(
        context.spec, boundary);
    emit_debug_stage(context, "epithelial_dirichlet_pin_prepare", before);
  }
  if (!context.preserve_gradient) return;
  const Real before_luminal = owned_content_slab(
      context.concentration, context.domain, context.storage_nx,
      context.halo_width, context.cell_volume);
  const DebugStageSnapshot luminal_snapshot = debug_snapshot(context);
  set_luminal_neumann_boundary_slab(
      context.concentration, context.domain, context.storage_nx,
      context.halo_width);
  context.flux.add_gradient_source(
      context.spec,
      owned_content_slab(
          context.concentration, context.domain, context.storage_nx,
          context.halo_width, context.cell_volume) - before_luminal);
  emit_debug_stage(context, "luminal_neumann_prepare", luminal_snapshot);
  const Real before_profile = owned_content_slab(
      context.concentration, context.domain, context.storage_nx,
      context.halo_width, context.cell_volume);
  const DebugStageSnapshot profile_snapshot = debug_snapshot(context);
  shift_z_gradient_slab(
      context.concentration, context.chemical, context.domain,
      context.storage_nx, context.halo_width, -1.0);
  context.flux.add_gradient_source(
      context.spec,
      owned_content_slab(
          context.concentration, context.domain, context.storage_nx,
          context.halo_width, context.cell_volume) - before_profile);
  emit_debug_stage(context, "profile_subtract", profile_snapshot);
  context.diffusion_boundary = 0.0;
}

void transport_slab_periodic(const SlabDiffusionContext& context) {
  const SlabTransportContext transport{
      context.concentration, context.domain, context.storage_nx,
      context.halo_width, context.alpha_x,
      context.delivery ? &context.sink_rate : nullptr,
      context.delivery ? &context.sink_realized : nullptr,
      context.delivery ? context.gradient_spec : nullptr,
      context.dt / 3.0, context.cell_volume, context.prescribed_mass,
      context.delivery ? &context.flux : nullptr, context.spec};
  const DebugStageSnapshot before_x = debug_snapshot(context);
  diffuse_periodic_x_slab(transport);
  emit_debug_stage(context, "x_delivery_solve", before_x);
  if (context.delivery) {
    const DebugStageSnapshot before_y = debug_snapshot(context);
    diffuse_periodic_y_slab_delivery(
        {context.concentration, context.domain, context.storage_nx,
         context.halo_width, context.alpha_y, &context.sink_rate,
         &context.sink_realized, context.gradient_spec,
         context.dt / 3.0, context.cell_volume,
         context.prescribed_mass, &context.flux, context.spec});
    emit_debug_stage(context, "y_delivery_solve", before_y);
  } else {
    diffuse_periodic_y_slab(
        context.concentration, context.domain, context.storage_nx,
        context.halo_width, context.alpha_y);
  }
}

Real transport_slab_z(const SlabDiffusionContext& context) {
  const SlabTransportContext z_transport{
      context.concentration, context.domain, context.storage_nx,
      context.halo_width, context.alpha_z,
      context.delivery ? &context.sink_rate : nullptr,
      context.delivery ? &context.sink_realized : nullptr,
      context.delivery ? context.gradient_spec : nullptr,
      context.dt / 3.0, context.cell_volume, context.prescribed_mass,
      context.delivery ? &context.flux : nullptr, context.spec};
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    const Real boundary = context.delivery
        ? diffuse_bounded_z_slab_delivery(
              z_transport, context.diffusion_boundary)
        : diffuse_bounded_z_slab(
              context.concentration, context.domain, context.storage_nx,
              context.halo_width, context.alpha_z,
              context.diffusion_boundary, context.cell_volume);
    return boundary;
  }
  const Real beta = context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Robin
      ? context.chemical.epithelial_transfer_coeff * context.dt
          / context.domain.dx_z() : 0.0;
  const Real flux_source = context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Flux
      ? context.chemical.epithelial_flux * context.dt
          / context.domain.dx_z() : 0.0;
  const DeliveryBoundaryParameters boundary_params{
      context.chemical.boundary_conc, beta, flux_source,
      context.chemical.epithelial_boundary_mode, context.cell_volume};
  const Real boundary = context.delivery
      ? diffuse_bounded_z_delivery_with_sink_slab(
            z_transport, boundary_params)
      : diffuse_bounded_z_delivery_slab(
            context.concentration, context.domain, context.storage_nx,
            context.halo_width, context.alpha_z, boundary_params);
  return boundary;
}

void transport_slab_diffusion(SlabDiffusionContext& context) {
  transport_slab_periodic(context);
  const DebugStageSnapshot before_z = debug_snapshot(context);
  context.flux.add_boundary(context.spec, transport_slab_z(context));
  emit_debug_stage(context, "z_delivery_solve", before_z);
}

void finish_slab_diffusion(SlabDiffusionContext& context) {
  if (context.preserve_gradient) {
    const Real before_profile = owned_content_slab(
        context.concentration, context.domain, context.storage_nx,
        context.halo_width, context.cell_volume);
    const DebugStageSnapshot profile_snapshot = debug_snapshot(context);
    shift_z_gradient_slab(
        context.concentration, context.chemical, context.domain,
        context.storage_nx, context.halo_width, 1.0);
    context.flux.add_gradient_source(
        context.spec,
        owned_content_slab(
            context.concentration, context.domain, context.storage_nx,
            context.halo_width, context.cell_volume) - before_profile);
    emit_debug_stage(context, "profile_restore", profile_snapshot);
    const Real before_luminal = owned_content_slab(
        context.concentration, context.domain, context.storage_nx,
        context.halo_width, context.cell_volume);
    const DebugStageSnapshot luminal_snapshot = debug_snapshot(context);
    set_luminal_neumann_boundary_slab(
        context.concentration, context.domain, context.storage_nx,
        context.halo_width);
    context.flux.add_gradient_source(
        context.spec,
        owned_content_slab(
            context.concentration, context.domain, context.storage_nx,
            context.halo_width, context.cell_volume) - before_luminal);
    emit_debug_stage(context, "luminal_neumann_finish", luminal_snapshot);
  }
  if (!context.delivery) {
    clamp_nonnegative_slab(
        context.concentration, context.domain, context.storage_nx,
        context.halo_width);
  } else if (!context.prescribed_active) {
    clamp_delivery_total_slab(
        context.concentration, context.domain, context.storage_nx,
        context.halo_width, context.cell_volume,
        context.flux, context.spec);
  }
  if (context.chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    const DebugStageSnapshot before = debug_snapshot(context);
    const Real boundary = set_epithelial_boundary_slab(
        context.concentration, context.domain, context.storage_nx,
        context.halo_width, context.chemical.boundary_conc,
        context.cell_volume);
    context.flux.add_boundary(
        context.spec, boundary);
    emit_debug_stage(context, "epithelial_dirichlet_pin_finish", before);
  }
}

}  // namespace

void ChemicalField::apply_diffusion(const Domain& domain, Real dt) {
  if (dt <= 0.0 || domain.dx_x() <= 0.0 || domain.dx_y() <= 0.0
      || domain.dx_z() <= 0.0) return;
  if (!delivery_support_stencil_.matches(
          domain, delivery_far_field_radius_)) {
    delivery_support_stencil_ =
        make_delivery_support_stencil(domain, delivery_far_field_radius_);
    delivery_affected_mask_.assign(
        static_cast<size_t>(ncells_), 0);
    delivery_affected_cells_.clear();
    delivery_affected_cells_.reserve(static_cast<size_t>(ncells_));
  }
  sum_prescribed_sinks_across_ranks();
  if (mode_ == DecompositionMode::Slab) {
    apply_diffusion_slab(domain, dt);
    return;
  }
  for (Int s = 0; s < nspec_; ++s) {
    apply_diffusion_species(domain, dt, s);
  }
}

bool ChemicalField::apply_diffusion_gpu(
    ChemicalFieldGpu& gpu, const Domain& domain, Real dt) {
  if (!gpu.active() || !delivery_route_b_eligible(domain, *this)) {
    return false;
  }
  for (Int s = 0; s < nspec_; ++s) {
    const ChemicalSpec& chemical = specs_[static_cast<size_t>(s)];
    if (!chemical.diffuses() || !chemical.delivery_enabled) continue;
    if (!gpu_delivery_species_eligible(domain, chemical, dt)) {
      return false;
    }
  }
  bool has_delivery = false;
  for (Int s = 0; s < nspec_; ++s) {
    const auto& chemical = specs_[static_cast<size_t>(s)];
    has_delivery |= chemical.diffuses() && chemical.delivery_enabled;
  }
  if (!has_delivery) return false;
  // Non-delivery species retain the existing constant-diagonal GPU kernels.
  bool has_non_delivery_diffusion = false;
  for (Int s = 0; s < nspec_; ++s) {
    const auto& chemical = specs_[static_cast<size_t>(s)];
    has_non_delivery_diffusion |=
        chemical.diffuses() && !chemical.delivery_enabled;
  }
  if (has_non_delivery_diffusion
      && !gpu.apply_diffusion(domain, *this, dt)) {
    return false;
  }
  for (Int s = 0; s < nspec_; ++s) {
    const ChemicalSpec& chemical = specs_[static_cast<size_t>(s)];
    if (!chemical.diffuses() || !chemical.delivery_enabled) continue;

    auto& concentration = conc_[static_cast<size_t>(s)];
    auto& prescribed = prescribed_sink_[static_cast<size_t>(s)];
    const auto concentration_snapshot = concentration;
    const auto flux_snapshot = flux_accounting_;
    assert(std::ranges::all_of(
        total_sink_realized_[static_cast<size_t>(s)],
        [](const Real value) { return value == 0.0; }));
    const auto realized_snapshot =
        total_sink_realized_[static_cast<size_t>(s)];
    const auto prescribed_snapshot = prescribed;
    gpu.snapshot_delivery_species(s);
    gpu.prepare_delivery_species(
        s, sink_rate_[static_cast<size_t>(s)], prescribed);

    const auto restore = [&gpu, this, s, &concentration_snapshot,
                          &flux_snapshot, &realized_snapshot] {
      conc_[static_cast<size_t>(s)] = concentration_snapshot;
      flux_accounting_ = flux_snapshot;
      total_sink_realized_[static_cast<size_t>(s)] = realized_snapshot;
      gpu.restore_delivery_species(s);
    };
    const auto restore_original =
        [&gpu, this, s, &restore, &prescribed, &prescribed_snapshot] {
          restore();
          prescribed = prescribed_snapshot;
          gpu.upload_delivery_prescribed(prescribed);
        };
    const auto solve = [
        &gpu, this, &domain, s, dt, &chemical] {
      gpu.upload_delivery_prescribed(prescribed_sink_[static_cast<size_t>(s)]);
      gpu.reset_delivery_boundary(s);
      if (!gpu.apply_delivery_species(
              domain, chemical, s, dt,
              prescribed_active_[static_cast<size_t>(s)])) {
        throw Error(
            "CUDA delivery precondition failed for species '"
            + chemical.name + "'");
      }
      flux_accounting_.add_boundary(
          s, gpu.download_delivery_boundary(s));
      flux_accounting_.add_gradient_source(
          s, gpu.download_delivery_gradient_source());
      flux_accounting_.add_reaction_clip(
          s, gpu.download_delivery_reaction_clip());
    };
    const auto owns_cell = [this](Int cell) {
      return owns_global_cell(cell);
    };
    const auto storage_cell = [this](Int cell) {
      return global_to_storage_cell(cell);
    };
    const DeliveryRationingCallbacks callbacks{
        restore,
        restore_original,
        solve,
        [&gpu, s] {
          return gpu.delivery_has_negative(s);
        },
        [&gpu, s] {
          return gpu.delivery_negative_fraction(s);
        },
        [this, &gpu, s, &domain, &owns_cell, &storage_cell] {
          // Retry dilation intentionally uses one device-to-host transfer;
          // duplicating the established support stencil on CUDA would add a
          // second scientific implementation for a rare path.
          std::vector<Real> realized;
          gpu.download_delivery_species(
              s, conc_[static_cast<size_t>(s)], realized);
          auto& current = prescribed_sink_[static_cast<size_t>(s)];
          const Real reduced = reduce_prescribed_near_negative_cells(
              conc_[static_cast<size_t>(s)], current, domain,
              delivery_support_stencil_, delivery_affected_mask_,
              delivery_affected_cells_, global_ncells_,
              owns_cell, storage_cell);
          gpu.upload_delivery_prescribed(current);
          return reduced;
        },
        [this, s]() -> std::vector<Real>& {
          return prescribed_sink_[static_cast<size_t>(s)];
        },
        [this, &owns_cell, &storage_cell](
            const std::vector<Real>& values) {
          return owned_prescribed_sum(
              values, global_ncells_, owns_cell, storage_cell);
        },
        [this, &owns_cell, &storage_cell](
            const std::vector<Real>& original,
            const std::vector<Real>& final_values) {
          return minimum_prescribed_ratio(
              original, final_values, global_ncells_,
              owns_cell, storage_cell);
        }};
    const auto result = run_delivery_rationing(
        prescribed_snapshot, callbacks);
    std::vector<Real> realized;
    gpu.download_delivery_species(s, concentration, realized);
    auto& total_realized = total_sink_realized_[static_cast<size_t>(s)];
    assert(total_realized.size() == realized_snapshot.size());
    assert(total_realized.size() == realized.size());
    for (size_t i = 0; i < total_realized.size(); ++i) {
      total_realized[i] = realized_snapshot[i] + realized[i];
    }
    record_delivery_rationing(s, chemical, result);
  }
  return has_delivery;
}

DeliveryRetryResult ChemicalField::run_delivery_rationing_for_species(
    Int species, const Domain& domain,
    const std::vector<Real>& concentration_snapshot,
    const NutrientFluxAccounting& flux_snapshot,
    const std::vector<Real>& realized_snapshot,
    const std::vector<Real>& prescribed_snapshot,
    const std::function<void()>& solve) {
  const auto restore = [
      this, species, &concentration_snapshot, &flux_snapshot,
      &realized_snapshot] {
    conc_[static_cast<size_t>(species)] = concentration_snapshot;
    flux_accounting_ = flux_snapshot;
    total_sink_realized_[static_cast<size_t>(species)] = realized_snapshot;
  };
  const auto restore_original = [
      this, species, &restore, &prescribed_snapshot] {
    restore();
    prescribed_sink_[static_cast<size_t>(species)] = prescribed_snapshot;
  };
  const auto owns_cell = [this](Int cell) {
    return owns_global_cell(cell);
  };
  const auto storage_cell = [this](Int cell) {
    return global_to_storage_cell(cell);
  };
  const DeliveryRationingCallbacks callbacks{
      restore,
      restore_original,
      solve,
      [this, species, &owns_cell, &storage_cell] {
        return collective_negative(has_negative_owned_cell(
            conc_[static_cast<size_t>(species)], global_ncells_,
            owns_cell, storage_cell));
      },
      [this, species, &owns_cell, &storage_cell] {
        return owned_negative_fraction(
            conc_[static_cast<size_t>(species)], global_ncells_,
            owns_cell, storage_cell);
      },
      [this, species, &domain, &owns_cell, &storage_cell] {
        auto& prescribed = prescribed_sink_[static_cast<size_t>(species)];
        return reduce_prescribed_near_negative_cells(
            conc_[static_cast<size_t>(species)], prescribed, domain,
            delivery_support_stencil_, delivery_affected_mask_,
            delivery_affected_cells_, global_ncells_,
            owns_cell, storage_cell);
      },
      [this, species]() -> std::vector<Real>& {
        return prescribed_sink_[static_cast<size_t>(species)];
      },
      [this, &owns_cell, &storage_cell](
          const std::vector<Real>& values) {
        return owned_prescribed_sum(
            values, global_ncells_, owns_cell, storage_cell);
      },
      [this, &owns_cell, &storage_cell](
          const std::vector<Real>& original,
          const std::vector<Real>& final_values) {
        return minimum_prescribed_ratio(
            original, final_values, global_ncells_,
            owns_cell, storage_cell);
      }};
  return run_delivery_rationing(prescribed_snapshot, callbacks);
}

void ChemicalField::record_delivery_rationing(
    Int species, const ChemicalSpec& chemical,
    const DeliveryRetryResult& result) {
  if (result.negative_after_solve) {
    const auto [minimum, count] = negative_diagnostics(
        conc_[static_cast<size_t>(species)], global_ncells_,
        [this](Int cell) { return owns_global_cell(cell); },
        [this](Int cell) { return global_to_storage_cell(cell); },
        mode_ != DecompositionMode::Slab);
    std::cerr << "Delivery infeasible: species=" << chemical.name
              << " step=" << nutrient_debug_step_counter()
              << " minimum_concentration=" << minimum
              << " negative_cells=" << count << "\n";
    flux_accounting_.add_delivery_infeasible(species, 1.0);
  }
  add_delivery_reduction(species, result.delivery_reduction);
  flux_accounting_.add_delivery_retry_events(species, result.retry_events);
  Real rationing_factor = result.rationing_factor;
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized && !finalized) {
    MPI_Allreduce(MPI_IN_PLACE, &rationing_factor, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
  }
#endif
  flux_accounting_.add_delivery_rationing_factor(
      species, rationing_factor);
  split_delivery_sink_realized(species);
  finalize_delivery_realized(species);
}

void ChemicalField::apply_diffusion_species(
    const Domain& domain, Real dt, Int s) {
  const ChemicalSpec& chemical = specs_[s];
  if (!chemical.diffuses()) return;
  const Real effective_diffusion = chemical.diff_coeff / chemical.retardation;
  const Real alpha_x = effective_diffusion * dt
      / (domain.dx_x() * domain.dx_x());
  const Real alpha_y = effective_diffusion * dt
      / (domain.dx_y() * domain.dx_y());
  const Real alpha_z = effective_diffusion * dt
      / (domain.dx_z() * domain.dx_z());
  ReplicatedDiffusionContext context{
      conc_[s], domain, chemical, flux_accounting_,
      sink_rate_[static_cast<size_t>(s)],
      total_sink_realized_[static_cast<size_t>(s)], s, dt, alpha_x, alpha_y,
      alpha_z,
      domain.cell_volume(), chemical.boundary_conc,
      chemical.z_gradient_enabled && chemical.z_gradient_lambda > 0.0,
      chemical.delivery_enabled,
      (chemical.z_gradient_enabled && chemical.z_gradient_lambda > 0.0)
          ? &chemical : nullptr,
      &prescribed_sink_[static_cast<size_t>(s)],
      prescribed_active_[static_cast<size_t>(s)]};
  if (!context.delivery) {
    prepare_replicated_diffusion(context);
    transport_replicated_diffusion(context);
    finish_replicated_diffusion(context);
    return;
  }
  const auto concentration_snapshot = conc_[static_cast<size_t>(s)];
  const auto flux_snapshot = flux_accounting_;
  const auto realized_snapshot = total_sink_realized_[static_cast<size_t>(s)];
  const auto prescribed_snapshot = prescribed_sink_[static_cast<size_t>(s)];
  const auto solve = [&context] {
      prepare_replicated_diffusion(context);
      transport_replicated_diffusion(context);
      finish_replicated_diffusion(context);
  };
  const auto final_retry = run_delivery_rationing_for_species(
      s, domain, concentration_snapshot, flux_snapshot, realized_snapshot,
      prescribed_snapshot, solve);
  record_delivery_rationing(s, chemical, final_retry);
}

void ChemicalField::finalize_delivery_realized(Int spec) {
  Real delivery_scale = 1.0;
  if (delivery_axis_cap_enabled()) {
    const auto species_index = static_cast<size_t>(spec);
    Real demanded = 0.0;
    for (Int cell = 0; cell < global_ncells_; ++cell) {
      if (!owns_global_cell(cell)) continue;
      const Int storage = global_to_storage_cell(cell);
      if (storage >= 0) {
        demanded += prescribed_sink_[species_index][
            static_cast<size_t>(storage)];
      }
    }
    const Real capped =
        flux_accounting_.delivery_axis_prescribed_mass_step[species_index];
    if (demanded > 0.0) {
      delivery_scale = std::clamp(capped / demanded, 0.0, 1.0);
    }
  }
  for (Int cell = 0; cell < global_ncells_; ++cell) {
    if (!owns_global_cell(cell)) continue;
    const Int storage = global_to_storage_cell(cell);
    if (storage < 0) continue;
    const auto index = static_cast<size_t>(storage);
    if (prescribed_active_[static_cast<size_t>(spec)]) {
      prescribed_sink_[static_cast<size_t>(spec)][index] *= delivery_scale;
      const Real agent = prescribed_sink_[static_cast<size_t>(spec)][index];
      sink_realized_[static_cast<size_t>(spec)][index] = agent;
      total_sink_realized_[static_cast<size_t>(spec)][index] =
          vbf_sink_realized_[static_cast<size_t>(spec)][index] + agent;
    }
  }
}

void ChemicalField::apply_periodic_x_diffusion(const Domain& domain, Real dt) {
  if (dt <= 0.0 || domain.dx_x() <= 0.0) return;
  for (Int s = 0; s < nspec_; ++s) {
    apply_periodic_x_diffusion(domain, dt, s);
  }
}

void ChemicalField::apply_periodic_x_diffusion(const Domain& domain, Real dt,
                                               Int spec) {
  if (dt <= 0.0 || domain.dx_x() <= 0.0 || spec < 0 || spec >= nspec_) return;
  const ChemicalSpec& chemical = specs_[static_cast<size_t>(spec)];
  if (!chemical.diffuses()) return;
  const Real alpha = (chemical.diff_coeff / chemical.retardation) * dt
      / (domain.dx_x() * domain.dx_x());
  if (mode_ == DecompositionMode::Slab) {
    diffuse_periodic_x_slab(
        {conc_[static_cast<size_t>(spec)], domain, storage_nx_, halo_width_,
         alpha});
  } else {
    diffuse_periodic_x(conc_[static_cast<size_t>(spec)], domain, alpha);
  }
}

void ChemicalField::apply_periodic_y_diffusion(const Domain& domain, Real dt,
                                               Int spec) {
  if (dt <= 0.0 || domain.dx_y() <= 0.0 || spec < 0 || spec >= nspec_) return;
  const ChemicalSpec& chemical = specs_[static_cast<size_t>(spec)];
  if (!chemical.diffuses()) return;
  const Real alpha = (chemical.diff_coeff / chemical.retardation) * dt
      / (domain.dx_y() * domain.dx_y());
  diffuse_periodic_y(conc_[static_cast<size_t>(spec)], domain, alpha);
}

void ChemicalField::apply_bounded_z_diffusion(const Domain& domain, Real dt,
                                              Int spec) {
  if (dt <= 0.0 || domain.dx_z() <= 0.0 || spec < 0 || spec >= nspec_) return;
  const ChemicalSpec& chemical = specs_[static_cast<size_t>(spec)];
  if (!chemical.diffuses()) return;
  const Real effective_diffusion =
      chemical.diff_coeff / chemical.retardation;
  const Real alpha = effective_diffusion * dt
      / (domain.dx_z() * domain.dx_z());
  if (chemical.epithelial_boundary_mode
      == EpithelialBoundaryMode::Dirichlet) {
    flux_accounting_.add_boundary(
        spec, diffuse_bounded_z(
                  conc_[static_cast<size_t>(spec)], domain, alpha,
                  chemical.boundary_conc, domain.cell_volume()));
  } else {
    const Real beta = chemical.epithelial_boundary_mode
        == EpithelialBoundaryMode::Robin
        ? chemical.epithelial_transfer_coeff * dt / domain.dx_z() : 0.0;
    const Real flux_source = chemical.epithelial_boundary_mode
        == EpithelialBoundaryMode::Flux
        ? chemical.epithelial_flux * dt / domain.dx_z() : 0.0;
    flux_accounting_.add_boundary(
        spec, diffuse_bounded_z_delivery(
                  conc_[static_cast<size_t>(spec)], domain, alpha,
                  {chemical.boundary_conc, beta, flux_source,
                   chemical.epithelial_boundary_mode, domain.cell_volume()}));
  }
}

void ChemicalField::apply_diffusion_slab(const Domain& domain, Real dt) {
  for (Int s = 0; s < nspec_; ++s) {
    apply_diffusion_slab_species(domain, dt, s);
  }
}

void ChemicalField::apply_diffusion_slab_species(
    const Domain& domain, Real dt, Int s) {
  const ChemicalSpec& chemical = specs_[s];
  if (!chemical.diffuses()) return;
  const Real effective_diffusion = chemical.diff_coeff / chemical.retardation;
  const Real alpha_x = effective_diffusion * dt
      / (domain.dx_x() * domain.dx_x());
  const Real alpha_y = effective_diffusion * dt
      / (domain.dx_y() * domain.dx_y());
  const Real alpha_z = effective_diffusion * dt
      / (domain.dx_z() * domain.dx_z());
  SlabDiffusionContext context{
      conc_[static_cast<size_t>(s)], domain, chemical, flux_accounting_,
      sink_rate_[static_cast<size_t>(s)],
      total_sink_realized_[static_cast<size_t>(s)], storage_nx_, halo_width_, s,
      dt,
      alpha_x, alpha_y, alpha_z, domain.cell_volume(),
      chemical.boundary_conc,
      chemical.z_gradient_enabled && chemical.z_gradient_lambda > 0.0,
      chemical.delivery_enabled,
      (chemical.z_gradient_enabled && chemical.z_gradient_lambda > 0.0)
          ? &chemical : nullptr,
      &prescribed_sink_[static_cast<size_t>(s)],
      prescribed_active_[static_cast<size_t>(s)]};
  if (!context.delivery) {
    prepare_slab_diffusion(context);
    transport_slab_diffusion(context);
    finish_slab_diffusion(context);
    return;
  }
  const auto concentration_snapshot = conc_[static_cast<size_t>(s)];
  const auto flux_snapshot = flux_accounting_;
  const auto realized_snapshot = total_sink_realized_[static_cast<size_t>(s)];
  const auto prescribed_snapshot = prescribed_sink_[static_cast<size_t>(s)];
  const auto solve = [&context] {
      prepare_slab_diffusion(context);
      transport_slab_diffusion(context);
      finish_slab_diffusion(context);
  };
  const auto final_retry = run_delivery_rationing_for_species(
      s, domain, concentration_snapshot, flux_snapshot, realized_snapshot,
      prescribed_snapshot, solve);
  record_delivery_rationing(s, chemical, final_retry);
}

namespace {

void apply_epithelial_boundary_layer(
    std::vector<std::vector<Real>>& concentration,
    const Domain& domain, Int species_index, Real boundary_conc,
    NutrientFluxAccounting& flux_accounting) {
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.nx(); ++ix) {
      const Int idx = domain.cell_index(ix, iy, 0);
      const Real old = concentration[static_cast<size_t>(species_index)]
          [static_cast<size_t>(idx)];
      concentration[static_cast<size_t>(species_index)]
          [static_cast<size_t>(idx)] = boundary_conc;
      flux_accounting.add_boundary(
          species_index,
          (boundary_conc - old) * domain.cell_volume());
    }
  }
}

void mirror_non_diffusing_top_layer(std::vector<std::vector<Real>>& concentration,
                                    const Domain& domain, Int species_index) {
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.nx(); ++ix) {
      const Int top = domain.cell_index(ix, iy, domain.nz() - 1);
      const Int below = domain.cell_index(ix, iy, domain.nz() - 2);
      concentration[static_cast<size_t>(species_index)]
          [static_cast<size_t>(top)] =
              concentration[static_cast<size_t>(species_index)]
                  [static_cast<size_t>(below)];
    }
  }
}

}  // namespace

void ChemicalField::apply_boundaries(const Domain& domain) {
  if (mode_ == DecompositionMode::Slab) {
    apply_boundaries_slab(domain);
    return;
  }
  const Int nz = domain.nz();

  for (Int s = 0; s < nspec_; ++s) {
    const Real bc = specs_[s].boundary_conc;

    // z=0 (epithelial surface): Dirichlet for nutrients. When a z-gradient is
    // configured, this is the peak concentration at the epithelium.
    if (specs_[s].epithelial_boundary_mode
        == EpithelialBoundaryMode::Dirichlet) {
      apply_epithelial_boundary_layer(
          conc_, domain, s, bc, flux_accounting_);
    }

    // The implicit z solve enforces the luminal zero-flux condition directly.
    // Non-diffusing fields retain the legacy mirrored top layer.
    if (!specs_[s].diffusion_enabled && nz >= 2) {
      mirror_non_diffusing_top_layer(conc_, domain, s);
    }
  }
}

void ChemicalField::apply_boundaries_slab(const Domain& domain) {
  const Real cell_volume = domain.cell_volume();
  for (Int s = 0; s < nspec_; ++s) {
    auto& concentration = conc_[static_cast<size_t>(s)];
    if (specs_[s].epithelial_boundary_mode
        == EpithelialBoundaryMode::Dirichlet) {
      flux_accounting_.add_boundary(
          s, set_epithelial_boundary_slab(
                 concentration, domain, storage_nx_, halo_width_,
                 specs_[s].boundary_conc, cell_volume));
    }
    if (!specs_[s].diffusion_enabled && domain.nz() >= 2) {
      set_luminal_neumann_boundary_slab(
          concentration, domain, storage_nx_, halo_width_);
    }
  }
}

void ChemicalField::debug_report_step(const Domain& domain) const {
  if (!nutrient_debug_enabled() || nutrient_debug_step_counter() <= 0
      || debug_initial_content_.size() != static_cast<size_t>(nspec_)) {
    return;
  }
  const auto& flux = flux_accounting_;
  const Int oxygen = find(species::OXYGEN);
  const Int carbon = find(species::CARBON);
  for (const Int s : {oxygen, carbon}) {
    if (s < 0 || s >= nspec_) continue;
    Real initial = debug_initial_content_[static_cast<size_t>(s)];
    Real after_content = mode_ == DecompositionMode::Slab
        ? owned_content_slab(
              conc_[static_cast<size_t>(s)], domain, storage_nx_,
              halo_width_, domain.cell_volume())
        : owned_content(
              conc_[static_cast<size_t>(s)], domain, domain.cell_volume());
#ifdef GUTIBM_MPI
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized && !finalized && mode_ == DecompositionMode::Slab
        && domain.nprocs() > 1) {
      MPI_Allreduce(MPI_IN_PLACE, &initial, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &after_content, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
    }
#endif
    const auto index = static_cast<size_t>(s);
    const Real residual = initial + flux.boundary_last_step[index]
        + flux.gradient_source_last_step[index]
        + flux.vbf_source_last_step[index]
        - flux.agent_uptake_last_step[index] - flux.maintenance_last_step[index]
        - flux.vbf_sink_last_step[index]
        + flux.reaction_clip_last_step[index] - after_content;
    if (domain.rank() == 0) {
      std::cout << "NUTRIENT_STEP step=" << nutrient_debug_step_counter()
                << " species=" << specs_[index].name
                << " owned_content_before=" << debug_real(initial)
                << " owned_content_after=" << debug_real(after_content)
                << " boundary_step="
                << debug_real(flux.boundary_last_step[index])
                << " gradient_source_step="
                << debug_real(flux.gradient_source_last_step[index])
                << " vbf_source_step="
                << debug_real(flux.vbf_source_last_step[index])
                << " agent_uptake_step="
                << debug_real(flux.agent_uptake_last_step[index])
                << " maintenance_step="
                << debug_real(flux.maintenance_last_step[index])
                << " vbf_sink_step="
                << debug_real(flux.vbf_sink_last_step[index])
                << " reaction_clip_step="
                << debug_real(flux.reaction_clip_last_step[index])
                << " residual=" << debug_real(residual) << '\n';
    }
  }
}

Int ChemicalField::find(std::string_view name) const {
  for (Int i = 0; i < nspec_; ++i) {
    if (specs_[i].name == name) return i;
  }
  return -1;
}

}  // namespace gutibm
