/* -----------------------------------------------------------------------
   GutIBM – Analytical Green's function kernels for QSSA diffusion
   
   Instead of explicit FTCS grid-based PDE solvers, we compute
   steady-state concentration fields from point sources using
   analytical solutions to the advection-diffusion equation.
   
   For a point source of strength Q (mol/s) in a uniform flow U
   with effective diffusivity D, the steady-state concentration is:
   
     C(r) = (Q / 4*pi*D*r) * exp(-U*(r - x_downstream) / (2*D))
   
   where r is distance from source and x_downstream is the
   downstream projection.  This gives the characteristic
   "comet-tail" pattern without any grid timestep constraint.
   
   Method of Images enforces no-flux boundaries at z=0 (epithelium)
   and z=h (lumen), preventing artificial mass loss.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_GREENS_FUNCTION_H
#define GUTIBM_GREENS_FUNCTION_H

#include "types.h"
#include "robin_correction_table.h"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <memory>
#include <ranges>
#include <numeric>
#include <utility>
#include <vector>

namespace gutibm {

inline constexpr int kDefaultImageSeriesMaxShells = 512;
inline constexpr int kHistoricalLegacyImageSeriesShells = 3;

class Domain;
class AdvectionField;

struct GreensFunctionParams {
  Real diff_coeff = 0.0;       // effective D (m^2/s) = D_free / retardation
  Real source_rate = 0.0;      // Q (mol/s)
  Real pI = 0.0;               // isoelectric point (determines retardation)
  Real retardation = 1.0;      // mucin retardation factor
  Real decay_rate = 0.0;       // first-order toxin degradation rate (1/s)
  // Infinity retains the sealed Neumann result for direct callers. QSSA
  // populates this from toxin.lumen_transfer_length (disabled by default).
  Real lumen_transfer_length = robin::kZeroTransferLength;
  bool lumen_transfer_basis_free = false;
  Real robin_cutoff = robin::kDefaultCutoff;
  Real image_series_relative_tolerance = 1.0e-10;
  int image_series_max_shells = kDefaultImageSeriesMaxShells;
  bool image_series_max_shells_explicit = false;
  bool image_series_legacy_reflections = false;
  bool drift_correction = false;

  // NOTE: bacteriocin pI classification lives in a single source of truth,
  // `classify_by_pI()` in src/genome/plasmid.h (pI > 8.5 → CORE, pI < 7.0 →
  // HALO, else NEUTRAL). The Green's function / QSSA code consumes the
  // pre-computed `BICluster.bclass` and must never re-classify from pI here.
};

// Concurrency diagnostics updated from const kernel evaluations under OpenMP.
// Heap-owned by GreensFunction so the atomics never need to move.
struct GreensFunctionDiagnostics {
  // One slot per OpenMP thread, padded to its own cache line.
  struct alignas(64) KernelEvaluationSlot {
    std::atomic<uint64_t> count{0};
  };

  std::atomic<uint64_t> image_series_cap_hits{0};
  std::atomic<uint64_t> low_screening_evaluations{0};
  std::atomic<uint64_t> drift_envelope_evaluations{0};
  std::atomic<uint64_t> negative_field_count{0};
  std::atomic<Real> most_negative_field{0.0};
  std::atomic<uint64_t> robin_direct_evaluations{0};
  std::atomic<uint64_t> robin_host_fallback_sources{0};
  std::vector<KernelEvaluationSlot> kernel_evaluations_by_thread;

  // Lowers most_negative_field to `value` if it is smaller.
  void record_negative_field(Real value) {
    Real current = most_negative_field.load();
    while (value < current
           && !most_negative_field.compare_exchange_weak(current, value)) {
      // `current` now holds the competing value; retry the comparison.
    }
  }
};

class GreensFunction {
 public:
  GreensFunction() = default;
  GreensFunction(const GreensFunction&) = delete;
  GreensFunction& operator=(const GreensFunction&) = delete;
  GreensFunction(GreensFunction&&) noexcept = default;
  GreensFunction& operator=(GreensFunction&&) noexcept = default;

  void init(const Domain& domain, const AdvectionField& adv);

  // Steady-state concentration at `target` from a point source at `source`
  // using the advection-diffusion Green's function.
  Real concentration(const Vec3& source, const Vec3& target,
                     const GreensFunctionParams& params) const;

  // Same but with Method of Images for bounded z-domain
  Real concentration_bounded(const Vec3& source, const Vec3& target,
                              const GreensFunctionParams& params) const;

  // Superpose contributions from multiple sources onto grid cells.
  // Uses spatial hashing cutoff to limit O(N*M) to O(N*k).
  void superpose_to_grid(const std::vector<Vec3>& sources,
                          const std::vector<GreensFunctionParams>& params,
                          std::vector<Real>& grid_conc,
                          Real cutoff_radius) const;

  void superpose_to_grid(const std::vector<Vec3>& sources,
                          const std::vector<GreensFunctionParams>& params,
                          const std::vector<Real>& strength_factors,
                          std::vector<Real>& grid_conc,
                          Real cutoff_radius) const;

  void superpose_to_local_grid(
      const std::vector<Vec3>& sources,
      const std::vector<GreensFunctionParams>& params,
      const std::vector<Real>& strength_factors,
      std::vector<Real>& grid_conc,
      Real cutoff_radius, Int x_begin, Int x_end, Int storage_nx,
      Int halo_width) const;

  // Peclet number at position: Pe = U*L/D
  Real peclet(const Vec3& pos, Real D_eff, Real length_scale) const;

  uint64_t image_series_cap_hits() const {
    return diag_->image_series_cap_hits.load();
  }
  uint64_t low_screening_evaluations() const {
    return diag_->low_screening_evaluations.load();
  }
  uint64_t drift_envelope_evaluations() const {
    return diag_->drift_envelope_evaluations.load();
  }
  uint64_t negative_field_count() const {
    return diag_->negative_field_count.load();
  }
  Real most_negative_field() const {
    return diag_->most_negative_field.load();
  }
  uint64_t robin_direct_evaluations() const {
    return diag_->robin_direct_evaluations.load();
  }
  uint64_t robin_host_fallback_sources() const {
    return diag_->robin_host_fallback_sources.load();
  }
  uint64_t kernel_evaluations() const {
    uint64_t total = 0;
    for (const auto& slot : diag_->kernel_evaluations_by_thread) {
      total += slot.count.load();
    }
    return total;
  }
  void set_kernel_evaluation_counting(bool enabled);
  bool kernel_evaluation_counting_enabled() const {
    return kernel_evaluation_counting_enabled_;
  }
  void add_image_series_cap_hits(uint64_t count) const {
    diag_->image_series_cap_hits.fetch_add(count);
  }
  void add_low_screening_evaluations(uint64_t count) const {
    diag_->low_screening_evaluations.fetch_add(count);
  }
  void add_drift_envelope_evaluations(uint64_t count) const {
    diag_->drift_envelope_evaluations.fetch_add(count);
  }
  void add_negative_field_diagnostics(uint64_t count,
                                      Real most_negative) const;
  void add_kernel_evaluations(uint64_t count) const {
    if (kernel_evaluation_counting_enabled_
        && !diag_->kernel_evaluations_by_thread.empty()) {
      diag_->kernel_evaluations_by_thread[0].count.fetch_add(count);
    }
  }
  void reset_image_series_cap_hits() {
    diag_->image_series_cap_hits.store(0);
  }
  void reset_low_screening_diagnostics() {
    diag_->low_screening_evaluations.store(0);
    diag_->drift_envelope_evaluations.store(0);
    diag_->negative_field_count.store(0);
    diag_->most_negative_field.store(0.0);
  }
  void add_robin_host_fallback_sources(uint64_t count) const {
    diag_->robin_host_fallback_sources.fetch_add(count);
  }
  void reset_kernel_evaluations() {
    for (auto& slot : diag_->kernel_evaluations_by_thread) {
      slot.count.store(0);
    }
  }

 private:
  void require_init() const;
  Real concentration_sealed(const Vec3& source, const Vec3& target,
                            const GreensFunctionParams& params) const;
  std::shared_ptr<const robin::Table> robin_table(
      const GreensFunctionParams& params) const;

  // Single image contribution
  Real single_kernel(const Vec3& src, const Vec3& tgt,
                      Real D_eff, Real Q, Real decay_rate,
                      const Vec3& flow_vel) const;

  const Domain* domain_    = nullptr;
  const AdvectionField* adv_ = nullptr;

  Real z_lo_ = 0.0;
  Real z_hi_ = 100.0e-6;
  bool kernel_evaluation_counting_enabled_ = false;
  std::unique_ptr<GreensFunctionDiagnostics> diag_ =
      std::make_unique<GreensFunctionDiagnostics>();
};

}  // namespace gutibm

#endif  // GUTIBM_GREENS_FUNCTION_H
