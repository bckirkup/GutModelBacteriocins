/* -----------------------------------------------------------------------
   GutIBM – QSSA solver implementation
   ----------------------------------------------------------------------- */

#include "qssa_solver.h"
#include "metabolic_mode.h"
#include "species_names.h"
#include "fmm.h"
#include "fmm_gpu.h"
#include "domain.h"
#include "advection.h"
#include "chemical_field.h"
#include "chemical_field_gpu.h"
#include "greens_function_gpu.h"
#include "neumann_image_series.h"
#include "dispatch.h"
#include "agent.h"
#include "error.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <format>
#include <limits>
#include <numbers>
#include <numeric>
#include <ranges>
#include <iostream>
#include <type_traits>
#ifdef GUTIBM_MPI
#include <mpi.h>
#endif
#ifdef GUTIBM_OPENMP
#include <omp.h>
#endif

namespace gutibm {

namespace {

constexpr Real k_ln2 = std::numbers::ln2;

Int toxin_sample_index(ReceptorType target);

bool is_bacteriocin_spec(const ChemicalSpec& spec) {
  return spec.name == species::BACTERIOCIN_BTUB
      || spec.name == species::BACTERIOCIN_FEPA
      || spec.name == species::BACTERIOCIN_CIRA
      || spec.name == species::BACTERIOCIN_FHUA
      || spec.name == species::BACTERIOCIN_LUMPED;
}

void enforce_low_screening_policy(
    const QSSAConfig& cfg, const Domain& domain, const AdvectionField& adv,
    const std::vector<ChemicalSpec>& chemicals) {
  if (cfg.low_screening_policy == "allow") return;
  const Vec3 probe = domain.cell_center(
      domain.nx() / 2, domain.ny() / 2, domain.nz() / 2);
  const Vec3 flow = adv.velocity(probe);
  for (const auto& spec : chemicals) {
    if (!is_bacteriocin_spec(spec)) continue;
    const Real d_eff = spec.diff_coeff / spec.retardation;
    if (!(d_eff > 0.0)) continue;
    const Real flow_magnitude = std::sqrt(
        flow[0] * flow[0] + flow[1] * flow[1] + flow[2] * flow[2]);
    const Real k_h = neumann::series_screening(
        d_eff, spec.decay_rate, flow_magnitude, std::abs(flow[2]))
        * (domain.hi()[2] - domain.lo()[2]);
    if (k_h < neumann::kLowScreeningFloorThreshold) {
      const std::string message = std::format(
          "low-screening sealed Neumann image series: kH={:f}; sealed "
          "truncation error is at least approximately 13% at the kH=0.0225 "
          "threshold and grows without bound as kH approaches zero; see "
          "docs/NEUMANN_LOW_SCREENING_ENVELOPE.md",
          k_h);
      if (cfg.low_screening_policy == "error") {
        throw SimulationError(message);
      }
      std::cerr << "WARNING: " << message
                << "; use qssa.low_screening_policy=allow for deliberate "
                   "unscreened diagnostics\n";
      return;
    }
  }
}

void enforce_drift_envelope_policy(
    const QSSAConfig& cfg, const Domain& domain, const AdvectionField& adv,
    const std::vector<ChemicalSpec>* chemicals,
    const RuntimeDriftEnvelopeBasis* runtime_basis) {
  if (cfg.drift_correction) {
    // The corrected sealed field satisfies the physical zero-total-flux wall
    // law, so the first-order Pe_z envelope no longer describes its error; the
    // remaining error is table interpolation only.
    std::cerr << "NOTE: wall-normal drift correction enabled: the sealed "
                 "Neumann field uses the physical zero-total-flux wall law; "
                 "the residual error is correction-table interpolation, "
                 "measured at 1.22e-4 median, 3.88e-3 p90, and 1.55e-1 "
                 "maximum at the ColE1 basis; see "
                 "docs/NEUMANN_WALL_NORMAL_DRIFT.md\n";
    return;
  }
  if (cfg.drift_envelope_policy == "allow") return;
  const Vec3 probe = domain.cell_center(
      domain.nx() / 2, domain.ny() / 2, domain.nz() / 2);
  const Vec3 flow = adv.velocity(probe);
  const Real height = domain.hi()[2] - domain.lo()[2];
  Real worst_pe_z = 0.0;
  Real worst_d_eff = 0.0;
  std::string worst_basis;
  if (chemicals != nullptr) {
    for (const auto& spec : *chemicals) {
      if (!is_bacteriocin_spec(spec)) continue;
      const Real d_eff = spec.diff_coeff / spec.retardation;
      if (!(d_eff > 0.0)) continue;
      const Real pe_z = neumann::wall_normal_peclet(flow[2], height, d_eff);
      if (std::abs(pe_z) > std::abs(worst_pe_z)) {
        worst_pe_z = pe_z;
        worst_d_eff = d_eff;
        worst_basis = "configured-species basis: " + spec.name;
      }
    }
  }
  if (runtime_basis != nullptr && runtime_basis->d_eff > 0.0) {
    const Real pe_z = neumann::wall_normal_peclet(
        flow[2], height, runtime_basis->d_eff);
    if (std::abs(pe_z) > std::abs(worst_pe_z)) {
      worst_pe_z = pe_z;
      worst_d_eff = runtime_basis->d_eff;
      worst_basis = "configured-plasmid basis: " + runtime_basis->label;
    }
  }
  if (!neumann::drift_envelope_exceeded(
          flow[2], height, worst_d_eff)) {
    return;
  }
  const std::string error_description =
      std::abs(worst_pe_z) > 0.3
      ? "worst-case relative field error is approximately 0.44*Pe_z only in "
        "the small-Pe fit; measured worst-case field error is 20-40% for "
        "strongly retarded shipped toxins"
      : "worst-case relative field error is approximately 0.44*Pe_z";
  const std::string message = std::format(
      "wall-normal drift in sealed Neumann image series: Pe_z={:.6g} "
      "({}; D_eff={:.6g}); envelope is |Pe_z| <= 0.05; {}; see "
      "docs/NEUMANN_WALL_NORMAL_DRIFT.md",
      worst_pe_z, worst_basis, worst_d_eff, error_description);
  if (cfg.drift_envelope_policy == "error") {
    throw SimulationError(message);
  }
  std::cerr << "WARNING: " << message
            << "; use qssa.drift_envelope_policy=allow for deliberate "
               "high-wall-normal-flow diagnostics\n";
}

int image_series_shells(const QSSAConfig& cfg) {
  if (cfg.image_series_mode == "pre_fix_duplicated_reflection"
      && !cfg.image_series_max_shells_explicit) {
    return kHistoricalLegacyImageSeriesShells;
  }
  return cfg.image_series_max_shells;
}

struct FarFieldGridContext {
  const Domain& domain;
  Real fmm_theta;
  Int nx;
  Int ny;
  Int nz;
  Int x_begin = 0;
  Int x_end = 0;
  Int storage_nx = 0;
  Int halo_width = 0;
};

struct MicrocinSourceBuffers {
  std::vector<Vec3>& sources;
  std::vector<GreensFunctionParams>& params;
  std::vector<Real>& strength_factors;
  std::vector<bool>& is_nuclease;
  std::vector<ReceptorType>& targets;
};

bool source_owned_by_rank(const Domain& domain, const Vec3& position) {
  return domain.nprocs() <= 1 || domain.owner_rank(position) == domain.rank();
}

struct ToxinSourceRecord {
  Vec3 position;
  GreensFunctionParams params;
  Real strength_factor;
  int target;
  int is_nuclease;
};

static_assert(std::is_trivially_copyable_v<ToxinSourceRecord>);

bool mpi_source_exchange_available(int& ranks) {
#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return false;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  return ranks > 1;
#else
  (void)ranks;
  return false;
#endif
}

void exchange_toxin_sources(
    std::vector<Vec3>& sources,
    std::vector<GreensFunctionParams>& params,
    std::vector<Real>& strength_factors,
    std::vector<bool>& is_nuclease,
    std::vector<ReceptorType>& targets) {
  int ranks = 1;
  if (!mpi_source_exchange_available(ranks)) return;

#ifdef GUTIBM_MPI
  const size_t local_count = sources.size();
  if (const size_t max_count = static_cast<size_t>(
          std::numeric_limits<int>::max()) / sizeof(ToxinSourceRecord);
      local_count > max_count) {
    throw SimulationError("toxin source count exceeds MPI exchange capacity");
  }

  std::vector local(local_count, ToxinSourceRecord{});
  for (size_t i = 0; i < local_count; ++i) {
    local[i].position = sources[i];
    local[i].params = params[i];
    local[i].strength_factor = strength_factors[i];
    local[i].target = to_underlying(targets[i]);
    local[i].is_nuclease = is_nuclease[i] ? 1 : 0;
  }

  const auto local_bytes =
      static_cast<int>(local.size() * sizeof(ToxinSourceRecord));
  std::vector byte_counts(ranks, 0);
  MPI_Allgather(&local_bytes, 1, MPI_INT, byte_counts.data(), 1, MPI_INT,
                MPI_COMM_WORLD);

  std::vector displacements(ranks, 0);
  int total_bytes = 0;
  for (int rank = 0; rank < ranks; ++rank) {
    displacements[rank] = total_bytes;
    if (byte_counts[rank] > std::numeric_limits<int>::max() - total_bytes) {
      throw SimulationError("toxin source MPI exchange is too large");
    }
    total_bytes += byte_counts[rank];
  }

  std::vector<std::byte> gathered(static_cast<size_t>(total_bytes));
  MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE, gathered.data(),
                 byte_counts.data(), displacements.data(), MPI_BYTE,
                 MPI_COMM_WORLD);

  const size_t total_count =
      static_cast<size_t>(total_bytes) / sizeof(ToxinSourceRecord);
  std::vector<ToxinSourceRecord> combined(total_count);
  if (total_count > 0) {
    std::memcpy(combined.data(), gathered.data(),
                total_count * sizeof(ToxinSourceRecord));
  }

  sources.clear();
  params.clear();
  strength_factors.clear();
  is_nuclease.clear();
  targets.clear();
  sources.reserve(total_count);
  params.reserve(total_count);
  strength_factors.reserve(total_count);
  is_nuclease.reserve(total_count);
  targets.reserve(total_count);
  for (const ToxinSourceRecord& source : combined) {
    sources.push_back(source.position);
    params.push_back(source.params);
    strength_factors.push_back(source.strength_factor);
    is_nuclease.push_back(source.is_nuclease != 0);
    targets.push_back(static_cast<ReceptorType>(source.target));
  }
#endif
}

void accumulate_far_field_cell(
    const FMM& fmm, const GreensFunction& gf,
    const GreensFunctionParams& avg_params, const FarFieldGridContext& grid,
    std::vector<Real>& toxin_conc, Real toxin_cutoff, bool near_field_on_device,
    Int ix, Int iy, Int iz) {
  const Vec3 tgt = grid.domain.cell_center(ix, iy, iz);
  const Int idx = iz * grid.storage_nx * grid.ny
      + iy * grid.storage_nx + ix - grid.x_begin + grid.halo_width;
  Real contribution = 0.0;
  if (near_field_on_device) {
    contribution = fmm.evaluate_far_field(
        tgt, grid.fmm_theta, toxin_cutoff, gf, avg_params);
  } else {
    const Real total = fmm.evaluate_total_field(
        tgt, grid.fmm_theta, gf, avg_params);
    contribution = std::max(0.0, total - toxin_conc[idx]);
  }
  toxin_conc[idx] += contribution;
}

GreensFunctionParams weighted_avg_params(
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    const QSSAConfig& cfg,
    std::vector<Real>& strengths) {
  strengths.resize(params.size());
  GreensFunctionParams avg_params{};
  Real total_s = 0.0;
  for (size_t i = 0; i < params.size(); ++i) {
    const Real s = params[i].source_rate * strength_factors[i];
    strengths[i] = s;
    avg_params.diff_coeff  += s * params[i].diff_coeff;
    avg_params.pI          += s * params[i].pI;
    avg_params.retardation += s * params[i].retardation;
    // FMM approximates mixed-toxin screening with a source-weighted k;
    // sqrt(k) is nonlinear, so heterogeneous clusters are approximate.
    avg_params.decay_rate  += s * params[i].decay_rate;
    total_s += s;
  }
  if (total_s > 0.0) {
    avg_params.diff_coeff  /= total_s;
    avg_params.pI          /= total_s;
    avg_params.retardation /= total_s;
    avg_params.decay_rate  /= total_s;
  } else {
    avg_params.diff_coeff  = cfg.fallback_diff_coeff;
    avg_params.pI          = cfg.fallback_pI;
    avg_params.retardation = cfg.fallback_retardation;
    avg_params.decay_rate = 0.0;
  }
  avg_params.source_rate = 0.0;
  avg_params.lumen_transfer_length = cfg.lumen_transfer_length;
  avg_params.lumen_transfer_basis_free = cfg.lumen_transfer_basis == "free";
  avg_params.robin_cutoff = cfg.toxin_cutoff;
  avg_params.image_series_relative_tolerance =
      cfg.image_series_relative_tolerance;
  avg_params.image_series_max_shells = image_series_shells(cfg);
  avg_params.image_series_max_shells_explicit =
      cfg.image_series_max_shells_explicit;
  avg_params.image_series_legacy_reflections =
      cfg.image_series_mode == "pre_fix_duplicated_reflection";
  avg_params.drift_correction = cfg.drift_correction;
  return avg_params;
}

void accumulate_far_field(const FMM& fmm,
                          const GreensFunction& gf,
                          const GreensFunctionParams& avg_params,
                          const FarFieldGridContext& grid,
                          std::vector<Real>& toxin_conc,
                          Real toxin_cutoff,
                          bool near_field_on_device) {
  for (Int iz = 0; iz < grid.nz; ++iz) {
    for (Int iy = 0; iy < grid.ny; ++iy) {
      for (Int ix = grid.x_begin; ix < grid.x_end; ++ix) {
        accumulate_far_field_cell(fmm, gf, avg_params, grid, toxin_conc,
                                  toxin_cutoff, near_field_on_device, ix, iy,
                                  iz);
      }
    }
  }
}

FarFieldGridContext make_far_field_grid(const Domain& domain, Real fmm_theta) {
  return {domain, fmm_theta, domain.nx(), domain.ny(), domain.nz(),
          0, domain.nx(), domain.nx(), 0};
}

void deposit_to_chemical_field(ChemicalField& chem,
                               Int toxin_species_idx,
                               const std::vector<Real>& concentrations) {
  if (chem.slab_mode()) {
    assert(static_cast<Int>(concentrations.size()) == chem.ncells());
    #ifdef GUTIBM_OPENMP
    #pragma omp parallel for collapse(3) schedule(static)
    #endif
    for (Int iz = 0; iz < chem.global_nz(); ++iz) {
      for (Int iy = 0; iy < chem.global_ny(); ++iy) {
        for (Int ix = chem.owned_storage_x_begin();
             ix < chem.owned_storage_x_end(); ++ix) {
          const Int c = iz * chem.storage_nx() * chem.global_ny()
              + iy * chem.storage_nx() + ix;
          chem.conc(toxin_species_idx, c) = concentrations[
              static_cast<size_t>(c)];
        }
      }
    }
    return;
  }
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (size_t c = 0; c < concentrations.size(); ++c) {
    const auto global_cell = static_cast<Int>(c);
    if (!chem.owns_global_cell(global_cell)) continue;
    chem.conc_global(toxin_species_idx, global_cell) = concentrations[c];
  }
}

bool try_gpu_near_field(const Domain& domain,
                        const AdvectionField& adv,
                        const std::vector<Vec3>& sources,
                        const std::vector<GreensFunctionParams>& params,
                        const std::vector<Real>& strength_factors,
                        Real cutoff_radius,
                        ChemicalField& chem,
                        Int toxin_species_idx,
                        ChemicalFieldGpu* chem_gpu,
                        bool defer_host_sync,
                        uint64_t* cap_hits,
                        uint64_t* kernel_evaluations,
                        uint64_t* low_screening_evaluations,
                        uint64_t* negative_field_count,
                        Real* most_negative_field) {
  if (chem_gpu == nullptr || !chem_gpu->active()) return false;
  double* d_conc = chem_gpu->conc_device(toxin_species_idx);
  if (d_conc == nullptr) return false;
  if (!gpu_superpose_to_device(domain, adv, sources, params, strength_factors,
                               d_conc, cutoff_radius, cap_hits,
                               kernel_evaluations, low_screening_evaluations,
                               negative_field_count, most_negative_field)) {
    return false;
  }
  if (!defer_host_sync) {
    chem_gpu->sync_species_concentrations_to_host(chem, toxin_species_idx);
  }
  return true;
}

bool accumulate_near_field_gpu_or_cpu(const Domain& domain,
                                      const GreensFunction& gf,
                                      const AdvectionField& adv,
                                      const std::vector<Vec3>& sources,
                                      const std::vector<GreensFunctionParams>& params,
                                      const std::vector<Real>& strength_factors,
                                      Real cutoff_radius,
                                      std::vector<Real>& toxin_conc,
                                      ChemicalField& chem,
                                      Int toxin_species_idx,
                                      ChemicalFieldGpu* chem_gpu,
                                      bool defer_host_sync) {
  uint64_t gpu_cap_hits = 0;
  uint64_t gpu_kernel_evaluations = 0;
  uint64_t gpu_low_screening = 0;
  uint64_t gpu_negative_count = 0;
  Real gpu_most_negative = 0.0;
  uint64_t* kernel_evaluations = gf.kernel_evaluation_counting_enabled()
      ? &gpu_kernel_evaluations
      : nullptr;
  if (const std::vector<size_t> fallback =
          ::gutibm::robin_host_fallback_sources(domain, sources, params);
      fallback.empty()) {
    if (try_gpu_near_field(
            domain, adv, sources, params, strength_factors, cutoff_radius,
            chem, toxin_species_idx, chem_gpu, defer_host_sync,
            &gpu_cap_hits, kernel_evaluations, &gpu_low_screening,
            &gpu_negative_count, &gpu_most_negative)) {
      gf.add_image_series_cap_hits(gpu_cap_hits);
      gf.add_low_screening_evaluations(gpu_low_screening);
      gf.add_negative_field_diagnostics(gpu_negative_count,
                                        gpu_most_negative);
      if (kernel_evaluations != nullptr) {
        gf.add_kernel_evaluations(gpu_kernel_evaluations);
      }
      return true;
    }
  } else if (gpu_runtime_enabled()) {
    std::vector<Vec3> gpu_sources;
    std::vector<GreensFunctionParams> gpu_params;
    std::vector<Real> gpu_strengths;
    std::vector<Vec3> host_sources;
    std::vector<GreensFunctionParams> host_params;
    std::vector<Real> host_strengths;
    std::vector<char> is_fallback(sources.size(), 0);
    for (const size_t index : fallback) is_fallback[index] = 1;
    for (size_t index = 0; index < sources.size(); ++index) {
      if (is_fallback[index] != 0) {
        host_sources.push_back(sources[index]);
        host_params.push_back(params[index]);
        host_strengths.push_back(strength_factors[index]);
      } else {
        gpu_sources.push_back(sources[index]);
        gpu_params.push_back(params[index]);
        gpu_strengths.push_back(strength_factors[index]);
      }
    }
    if (try_gpu_near_field(
            domain, adv, gpu_sources, gpu_params, gpu_strengths,
            cutoff_radius, chem, toxin_species_idx, chem_gpu,
            true, &gpu_cap_hits, kernel_evaluations, &gpu_low_screening,
            &gpu_negative_count, &gpu_most_negative)) {
      std::vector<Real> host_grid;
      gf.superpose_to_grid(
          host_sources, host_params, host_strengths, host_grid,
          cutoff_radius);
      chem_gpu->sync_species_concentrations_to_host(chem, toxin_species_idx);
      std::vector<Real>& device_grid =
          chem.mutable_species_concentration(toxin_species_idx);
      for (size_t cell = 0; cell < device_grid.size(); ++cell) {
        device_grid[cell] += host_grid[cell];
      }
      chem_gpu->sync_species_concentrations_to_device(chem, toxin_species_idx);
      gf.add_image_series_cap_hits(gpu_cap_hits);
      gf.add_low_screening_evaluations(gpu_low_screening);
      gf.add_negative_field_diagnostics(gpu_negative_count,
                                        gpu_most_negative);
      if (kernel_evaluations != nullptr) {
        gf.add_kernel_evaluations(gpu_kernel_evaluations);
      }
      gf.add_robin_host_fallback_sources(fallback.size());
      return true;
    }
  }
  gf.superpose_to_grid(sources, params, strength_factors, toxin_conc,
                       cutoff_radius);
  return false;
}

void collect_microcin_sources(const AgentPool& agents,
                              const QSSAConfig& cfg,
                              const ProteaseConfig& protease,
                              const AdvectionField& adv,
                              const Domain& domain,
                              MicrocinSourceBuffers& out) {
  for (const Agent& a : agents) {
    if (a.state == PhenoState::DEAD || a.state == PhenoState::SOS_INDUCED) continue;
    if (!source_owned_by_rank(domain, a.x)) continue;

    for (const auto& bi : a.genome.bi_loci) {
      if (bi.release_mode != ReleaseMode::CONTINUOUS) continue;

      GreensFunctionParams gfp;
      gfp.diff_coeff   = bi.diff_coeff;
      gfp.retardation  = bi.retardation;
      gfp.pI           = bi.pI;
      gfp.source_rate  = cfg.microcin_secretion;
      gfp.lumen_transfer_length = cfg.lumen_transfer_length;
      gfp.lumen_transfer_basis_free = cfg.lumen_transfer_basis == "free";
      gfp.robin_cutoff = cfg.toxin_cutoff;
      gfp.image_series_relative_tolerance =
          cfg.image_series_relative_tolerance;
      gfp.image_series_max_shells = image_series_shells(cfg);
      gfp.image_series_max_shells_explicit =
          cfg.image_series_max_shells_explicit;
      gfp.image_series_legacy_reflections =
          cfg.image_series_mode == "pre_fix_duplicated_reflection";
      gfp.drift_correction = cfg.drift_correction;
      const Real protease_decay = (protease.enabled
                                   && bi.protease_half_life > 0.0)
          ? k_ln2 / bi.protease_half_life : 0.0;
      const Real dilution_decay = std::max(
          adv.washout_rate(a.x[2]), protease.dilution_rate);
      gfp.decay_rate = protease_decay + dilution_decay;

      out.sources.push_back(a.x);
      out.params.push_back(gfp);
      out.strength_factors.push_back(1.0);
      out.is_nuclease.push_back(bi.is_nuclease);
      out.targets.push_back(bi.target);
    }
  }
}

void append_burst_sources(const std::vector<ToxinBurstSource>& bursts,
                          Real current_time,
                          const Domain& domain,
                          MicrocinSourceBuffers& out) {
  for (const ToxinBurstSource& burst : bursts) {
    if (!source_owned_by_rank(domain, burst.pos)) continue;
    if (burst.release_tau <= 0.0) continue;
    const Real age = std::max(0.0, current_time - burst.creation_time);
    const Real factor = std::exp(-age / burst.release_tau);
    if (factor < 1.0e-12) continue;

    out.sources.push_back(burst.pos);
    out.params.push_back(burst.params);
    out.strength_factors.push_back(factor);
    out.is_nuclease.push_back(burst.is_nuclease);
    out.targets.push_back(burst.target);
  }
}

void filter_sources_by_target(const std::vector<Vec3>& sources,
                              const std::vector<GreensFunctionParams>& params,
                              const std::vector<Real>& strength_factors,
                              const std::vector<ReceptorType>& targets,
                              ReceptorType target,
                              std::vector<Vec3>& out_sources,
                              std::vector<GreensFunctionParams>& out_params,
                              std::vector<Real>& out_strengths) {
  for (size_t i = 0; i < sources.size(); ++i) {
    if (targets[i] != target) continue;
    out_sources.push_back(sources[i]);
    out_params.push_back(params[i]);
    out_strengths.push_back(strength_factors[i]);
  }
}

void filter_nuclease_sources_by_target(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    const std::vector<bool>& is_nuclease,
    const std::vector<ReceptorType>& targets,
    ReceptorType target,
    std::vector<Vec3>& out_sources,
    std::vector<GreensFunctionParams>& out_params,
    std::vector<Real>& out_strengths) {
  for (size_t i = 0; i < sources.size(); ++i) {
    if (!is_nuclease[i] || targets[i] != target) continue;
    out_sources.push_back(sources[i]);
    out_params.push_back(params[i]);
    out_strengths.push_back(strength_factors[i]);
  }
}

void zero_species_field(ChemicalField& chem,
                        Int species_idx,
                        ChemicalFieldGpu* chem_gpu) {
  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for collapse(3) schedule(static)
  #endif
  for (Int iz = 0; iz < chem.global_nz(); ++iz) {
    for (Int iy = 0; iy < chem.global_ny(); ++iy) {
      for (Int ix = chem.owned_storage_x_begin();
           ix < chem.owned_storage_x_end(); ++ix) {
        const Int c = iz * chem.storage_nx() * chem.global_ny()
            + iy * chem.storage_nx() + ix;
        chem.conc(species_idx, c) = 0.0;
      }
    }
  }
  if (chem_gpu != nullptr && chem_gpu->active()) {
    chem_gpu->sync_species_concentrations_to_device(chem, species_idx);
  }
}

}  // namespace

void QSSASolver::init(const QSSAConfig& cfg, const Domain& domain,
                       const AdvectionField& adv, bool profile_steps,
                       const std::vector<ChemicalSpec>* chemicals,
                       const RuntimeDriftEnvelopeBasis* runtime_basis) {
  cfg_    = cfg;
  domain_ = &domain;
  adv_    = &adv;
  gf_.init(domain, adv);
  gf_.set_kernel_evaluation_counting(profile_steps);
  gf_.reset_image_series_cap_hits();
  gf_.reset_low_screening_diagnostics();
  gf_.reset_kernel_evaluations();
  if (chemicals != nullptr) {
    enforce_low_screening_policy(cfg, domain, adv, *chemicals);
  }
  if (chemicals != nullptr || runtime_basis != nullptr) {
    enforce_drift_envelope_policy(cfg, domain, adv, chemicals, runtime_basis);
  }
  if (cfg.image_series_mode == "pre_fix_duplicated_reflection") {
    std::cerr
        << "WARNING: image_series_mode=pre_fix_duplicated_reflection is a "
           "benchmark-only cost reference and produces a physically wrong "
           "field; do not use it for science runs.\n";
  }
}

void QSSASolver::solve_lumped_bacteriocin_fields(
    const AgentPool& agents,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    const AdvectionField& adv,
    ChemicalField& chem,
    ChemicalFieldGpu* chem_gpu,
    bool materialize_grid) {
  const Int lumped_idx =
      chem.find(species::BACTERIOCIN_LUMPED);
  if (lumped_idx < 0) {
    throw SimulationError(
        "lumped toxin field is missing from the chemical field");
  }

  if (cfg_.toxin_evaluation == "agents") {
    sampled_species_indices_.fill(lumped_idx);
    sampled_agents_ = &agents;
    sample_bacteriocin_field(sources, params, strength_factors, agents,
                             ReceptorType::BtuB);
  }
  if (!materialize_grid || sources.empty()) {
    zero_species_field(chem, lumped_idx, chem_gpu);
  } else {
    solve_bacteriocin_field_from_sources(
        sources, params, strength_factors, adv, chem, lumped_idx, chem_gpu);
  }

  if (chem_gpu != nullptr && chem_gpu->active()) {
    chem_gpu->sync_species_concentrations_to_host(chem, lumped_idx);
  }
}

void QSSASolver::solve_bacteriocin_field(
    const AgentPool& agents,
    const std::vector<ToxinBurstSource>& bursts,
    Real current_time,
    const ProteaseConfig& protease,
    const AdvectionField& adv,
    ChemicalField& chem,
    Int toxin_species_idx,
    ReceptorType target,
    ChemicalFieldGpu* chem_gpu,
    bool materialize_grid) {

  std::vector<Vec3> all_sources;
  std::vector<GreensFunctionParams> all_params;
  std::vector<Real> all_strengths;
  std::vector<bool> is_nuclease;
  std::vector<ReceptorType> all_targets;
  MicrocinSourceBuffers buffers{all_sources, all_params, all_strengths, is_nuclease, all_targets};

  collect_microcin_sources(agents, cfg_, protease, adv, *domain_, buffers);
  append_burst_sources(bursts, current_time, *domain_, buffers);
  for (auto& param : all_params) {
    param.lumen_transfer_length = cfg_.lumen_transfer_length;
    param.lumen_transfer_basis_free = cfg_.lumen_transfer_basis == "free";
    param.robin_cutoff = cfg_.toxin_cutoff;
    param.image_series_relative_tolerance =
        cfg_.image_series_relative_tolerance;
    param.image_series_max_shells = image_series_shells(cfg_);
    param.image_series_max_shells_explicit =
        cfg_.image_series_max_shells_explicit;
    param.image_series_legacy_reflections =
        cfg_.image_series_mode == "pre_fix_duplicated_reflection";
    param.drift_correction = cfg_.drift_correction;
  }
  exchange_toxin_sources(all_sources, all_params, all_strengths, is_nuclease,
                         all_targets);
  sample_nuclease_sources(all_sources, all_params, all_strengths, is_nuclease,
                          all_targets, agents);

  if (cfg_.toxin_lumping == "lumped") {
    solve_lumped_bacteriocin_fields(
        agents, all_sources, all_params, all_strengths, adv, chem, chem_gpu,
        materialize_grid);
    return;
  }

  std::vector<Vec3> sources;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strength_factors;
  filter_sources_by_target(all_sources, all_params, all_strengths, all_targets, target,
                           sources, params, strength_factors);

  if (cfg_.toxin_evaluation == "agents") {
    sampled_species_indices_.fill(-1);
    if (const Int sample_idx = toxin_sample_index(target); sample_idx >= 0) {
      sampled_species_indices_[static_cast<size_t>(sample_idx)] =
          toxin_species_idx;
    }
    sample_bacteriocin_field(sources, params, strength_factors, agents,
                             target);
    if (materialize_grid) {
      if (sources.empty()) {
        zero_species_field(chem, toxin_species_idx, chem_gpu);
      } else {
        solve_bacteriocin_field_from_sources(
            sources, params, strength_factors, adv, chem, toxin_species_idx,
            chem_gpu);
      }
    } else {
      zero_species_field(chem, toxin_species_idx, chem_gpu);
    }
    return;
  }
  if (sources.empty()) {
    zero_species_field(chem, toxin_species_idx, chem_gpu);
    return;
  }
  solve_bacteriocin_field_from_sources(sources, params, strength_factors, adv,
                                       chem, toxin_species_idx, chem_gpu);
}

void QSSASolver::solve_bacteriocin_field_from_sources(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    const AdvectionField& adv,
    ChemicalField& chem,
    Int toxin_species_idx,
    ChemicalFieldGpu* chem_gpu) const {
  if (cfg_.use_fmm) {
    solve_bacteriocin_field_fmm(sources, params, strength_factors, adv,
                                chem, toxin_species_idx, chem_gpu);
    return;
  }

  const Int toxin_cells = chem.slab_mode() ? chem.ncells()
                                           : domain_->ncells();
  std::vector toxin_conc(static_cast<size_t>(toxin_cells), 0.0);
  const bool defer_sync = chem_gpu != nullptr && chem_gpu->active();
  const bool use_local_grid = chem.slab_mode();
  bool near_on_gpu = false;
  if (use_local_grid) {
    gf_.superpose_to_local_grid(
        sources, params, strength_factors, toxin_conc, cfg_.toxin_cutoff,
        domain_->local_grid_x_begin(), domain_->local_grid_x_end(),
        chem.storage_nx(), chem.owned_storage_x_begin());
  } else {
    near_on_gpu = accumulate_near_field_gpu_or_cpu(
        *domain_, gf_, adv, sources, params, strength_factors,
        cfg_.toxin_cutoff, toxin_conc, chem, toxin_species_idx,
        chem_gpu, defer_sync);
  }
  if (!near_on_gpu) {
    deposit_to_chemical_field(chem, toxin_species_idx, toxin_conc);
  }
}

void QSSASolver::solve_all_bacteriocin_fields(
    const AgentPool& agents,
    const std::vector<ToxinBurstSource>& bursts,
    Real current_time,
    const ProteaseConfig& protease,
    const AdvectionField& adv,
    ChemicalField& chem,
    ChemicalFieldGpu* chem_gpu,
    bool materialize_grid) {
  std::vector<Int> solved_indices;
  solved_indices.reserve(species::BACTERIOCIN_RECEPTOR_TARGETS.size());

  std::vector<Vec3> all_sources;
  std::vector<GreensFunctionParams> all_params;
  std::vector<Real> all_strengths;
  std::vector<bool> is_nuclease;
  std::vector<ReceptorType> all_targets;
  MicrocinSourceBuffers buffers{all_sources, all_params, all_strengths,
                                is_nuclease, all_targets};
  collect_microcin_sources(agents, cfg_, protease, adv, *domain_, buffers);
  append_burst_sources(bursts, current_time, *domain_, buffers);
  for (auto& param : all_params) {
    param.lumen_transfer_length = cfg_.lumen_transfer_length;
    param.lumen_transfer_basis_free = cfg_.lumen_transfer_basis == "free";
    param.robin_cutoff = cfg_.toxin_cutoff;
    param.image_series_relative_tolerance =
        cfg_.image_series_relative_tolerance;
    param.image_series_max_shells = image_series_shells(cfg_);
    param.image_series_max_shells_explicit =
        cfg_.image_series_max_shells_explicit;
    param.image_series_legacy_reflections =
        cfg_.image_series_mode == "pre_fix_duplicated_reflection";
    param.drift_correction = cfg_.drift_correction;
  }
  exchange_toxin_sources(all_sources, all_params, all_strengths, is_nuclease,
                         all_targets);
  sample_nuclease_sources(all_sources, all_params, all_strengths, is_nuclease,
                          all_targets, agents);

  if (cfg_.toxin_lumping == "lumped") {
    solve_lumped_bacteriocin_fields(
        agents, all_sources, all_params, all_strengths, adv, chem, chem_gpu,
        materialize_grid);
    return;
  }

  std::vector<Vec3> sources;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strength_factors;
  if (cfg_.toxin_evaluation == "agents") {
    sampled_species_indices_.fill(-1);
    sampled_agents_ = &agents;
    for (ReceptorType target : species::BACTERIOCIN_RECEPTOR_TARGETS) {
      const char* name = species::bacteriocin_species_for(target);
      if (name == nullptr) continue;
      Int idx = chem.find(name);
      if (idx < 0) continue;
      if (const Int sample_idx = toxin_sample_index(target);
          sample_idx >= 0) {
        sampled_species_indices_[static_cast<size_t>(sample_idx)] = idx;
      }
      filter_sources_by_target(all_sources, all_params, all_strengths,
                               all_targets, target, sources, params,
                               strength_factors);
      sample_bacteriocin_field(sources, params, strength_factors, agents,
                               target);
      if (sources.empty() || !materialize_grid) {
        zero_species_field(chem, idx, chem_gpu);
      } else {
        solve_bacteriocin_field_from_sources(
            sources, params, strength_factors, adv, chem, idx, chem_gpu);
      }
      sources.clear();
      params.clear();
      strength_factors.clear();
    }
    return;
  }

  for (ReceptorType target : species::BACTERIOCIN_RECEPTOR_TARGETS) {
    const char* name = species::bacteriocin_species_for(target);
    if (name == nullptr) continue;
    Int idx = chem.find(name);
    if (idx < 0) continue;
    filter_sources_by_target(all_sources, all_params, all_strengths, all_targets,
                             target, sources, params, strength_factors);
    if (sources.empty()) {
      zero_species_field(chem, idx, chem_gpu);
    } else {
      solve_bacteriocin_field_from_sources(
          sources, params, strength_factors, adv, chem, idx, chem_gpu);
    }
    sources.clear();
    params.clear();
    strength_factors.clear();
    if (chem_gpu != nullptr && chem_gpu->active()) {
      solved_indices.push_back(idx);
    }
  }

  if (chem_gpu != nullptr && chem_gpu->active()) {
    for (Int spec : solved_indices) {
      chem_gpu->sync_species_concentrations_to_host(chem, spec);
    }
  }
}

namespace {

Int toxin_sample_index(ReceptorType target) {
  using enum ReceptorType;
  switch (target) {
    case BtuB: return 0;
    case FepA: return 1;
    case CirA: return 2;
    case FhuA: return 3;
    default: return -1;
  }
}

}  // namespace

void QSSASolver::sample_bacteriocin_field(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    const AgentPool& agents,
    ReceptorType target) {
  const Int sample_idx = toxin_sample_index(target);
  if (sample_idx < 0) return;
  auto& field = sampled_fields_[static_cast<size_t>(sample_idx)];
  field.sources = sources;
  field.params = params;
  field.strength_factors = strength_factors;
  field.samples.assign(static_cast<size_t>(agents.size()), 0.0);
  field.fmm_ready = false;
  sampled_agents_ = &agents;
  if (sources.empty()) return;

  std::vector<Real> strengths;
  GreensFunctionParams avg_params = weighted_avg_params(
      params, strength_factors, cfg_, strengths);
  if (cfg_.use_fmm) {
    field.fmm.build(sources, strengths, *domain_, cfg_.fmm_expansion_order);
    field.fmm.compute_local_expansions(cfg_.fmm_theta, gf_, avg_params);
    field.fmm_ready = true;
  }
  for (Int i = 0; i < agents.size(); ++i) {
    field.samples[static_cast<size_t>(i)] =
        std::max(evaluate_sample(field, agents[i].x), 0.0);
  }
}

void QSSASolver::sample_nuclease_field(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    const AgentPool& agents,
    ReceptorType target) {
  auto& field = sampled_nuclease_fields_[
      static_cast<size_t>(toxin_sample_index(target))];
  field.sources = sources;
  field.params = params;
  field.strength_factors = strength_factors;
  field.samples.assign(static_cast<size_t>(agents.size()), 0.0);
  field.fmm_ready = false;
  sampled_agents_ = &agents;
  if (sources.empty()) return;

  std::vector<Real> strengths;
  const GreensFunctionParams avg_params = weighted_avg_params(
      params, strength_factors, cfg_, strengths);
  if (cfg_.use_fmm) {
    field.fmm.build(
        sources, strengths, *domain_, cfg_.fmm_expansion_order);
    field.fmm.compute_local_expansions(
        cfg_.fmm_theta, gf_, avg_params);
    field.fmm_ready = true;
  }
  for (Int i = 0; i < agents.size(); ++i) {
    field.samples[static_cast<size_t>(i)] =
        std::max(evaluate_sample(field, agents[i].x), 0.0);
  }
}

void QSSASolver::sample_nuclease_sources(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    const std::vector<bool>& is_nuclease,
    const std::vector<ReceptorType>& targets,
    const AgentPool& agents) {
  sampled_nuclease_sources_ = std::ranges::any_of(
      is_nuclease, [](bool value) {
        return value;
      });
  if (!sampled_nuclease_sources_) return;

  std::vector<Vec3> nuclease_sources;
  std::vector<GreensFunctionParams> nuclease_params;
  std::vector<Real> nuclease_strengths;
  for (ReceptorType target : species::BACTERIOCIN_RECEPTOR_TARGETS) {
    filter_nuclease_sources_by_target(
        sources, params, strength_factors, is_nuclease, targets, target,
        nuclease_sources, nuclease_params, nuclease_strengths);
    sample_nuclease_field(nuclease_sources, nuclease_params,
                          nuclease_strengths, agents, target);
    nuclease_sources.clear();
    nuclease_params.clear();
    nuclease_strengths.clear();
  }
}

Real QSSASolver::evaluate_sample(const SampledToxinField& field,
                                 const Vec3& position) const {
  if (field.sources.empty()) return 0.0;
  if (field.fmm_ready) {
    std::vector<Real> strengths;
    const GreensFunctionParams avg_params = weighted_avg_params(
        field.params, field.strength_factors, cfg_, strengths);
    return field.fmm.evaluate_field(
        position, cfg_.fmm_theta, cfg_.toxin_cutoff, gf_, field.sources,
        field.params, avg_params);
  }
  return point_concentration(position, field.sources, field.params,
                             field.strength_factors);
}

Int QSSASolver::sampled_slot_for_species(Int species_idx) const {
  const auto it = std::ranges::find(sampled_species_indices_, species_idx);
  return it == sampled_species_indices_.end()
      ? -1
      : static_cast<Int>(std::distance(sampled_species_indices_.begin(), it));
}

Real QSSASolver::sampled_toxin_conc(Int agent_index, Int species_idx) const {
  if (agent_index < 0 || species_idx < 0 || sampled_agents_ == nullptr) {
    return 0.0;
  }
  const Int sample_idx = sampled_slot_for_species(species_idx);
  if (sample_idx < 0 || agent_index >= sampled_agents_->size()) return 0.0;
  const auto& field = sampled_fields_[static_cast<size_t>(sample_idx)];
  if (agent_index >= static_cast<Int>(field.samples.size())) {
    return std::max(
        evaluate_sample(field, (*sampled_agents_)[agent_index].x), 0.0);
  }
  return field.samples[static_cast<size_t>(agent_index)];
}

Int QSSASolver::sampled_toxin_sample_count(Int species_idx) const {
  const Int sample_idx = sampled_slot_for_species(species_idx);
  if (sample_idx < 0) return 0;
  return static_cast<Int>(
      sampled_fields_[static_cast<size_t>(sample_idx)].samples.size());
}

Real QSSASolver::sampled_nuclease_conc(Int agent_index) const {
  Real total = 0.0;
  for (ReceptorType target : species::BACTERIOCIN_RECEPTOR_TARGETS) {
    total += sampled_nuclease_conc(agent_index, target);
  }
  return total;
}

Real QSSASolver::sampled_nuclease_conc(
    Int agent_index, ReceptorType target) const {
  if (!sampled_nuclease_sources_) return 0.0;
  if (agent_index < 0 || sampled_agents_ == nullptr
      || agent_index >= sampled_agents_->size()) {
    return 0.0;
  }
  const Int sample_idx = toxin_sample_index(target);
  if (sample_idx < 0) return 0.0;
  const auto& samples =
      sampled_nuclease_fields_[static_cast<size_t>(sample_idx)].samples;
  if (agent_index >= static_cast<Int>(samples.size())) return 0.0;
  return samples[static_cast<size_t>(agent_index)];
}

Real QSSASolver::sampled_toxin_max(Int species_idx) const {
  if (const Int sample_idx = sampled_slot_for_species(species_idx);
      sample_idx < 0) return 0.0;
  if (sampled_agents_ == nullptr) return 0.0;
  Real maximum = 0.0;
  for (Int i = 0; i < sampled_agents_->size(); ++i) {
    maximum = std::max(maximum, sampled_toxin_conc(i, species_idx));
  }
  return maximum;
}

void QSSASolver::solve_bacteriocin_field_fmm(
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    const AdvectionField& adv,
    ChemicalField& chem,
    Int toxin_species_idx,
    ChemicalFieldGpu* chem_gpu) const {

  std::vector<Real> strengths;
  const GreensFunctionParams avg_params =
      weighted_avg_params(params, strength_factors, cfg_, strengths);

  FMM fmm;
  fmm.build(sources, strengths, *domain_, cfg_.fmm_expansion_order);
  fmm.compute_local_expansions(cfg_.fmm_theta, gf_, avg_params);

  const Int ncells = chem.slab_mode() ? chem.ncells() : domain_->ncells();
  const FarFieldGridContext far_grid = make_far_field_grid(*domain_, cfg_.fmm_theta);

  std::vector toxin_conc(static_cast<size_t>(ncells), 0.0);
  const bool defer_sync = chem_gpu != nullptr && chem_gpu->active();
  bool near_on_gpu = false;
  if (chem.slab_mode()) {
    gf_.superpose_to_local_grid(
        sources, params, strength_factors, toxin_conc, cfg_.toxin_cutoff,
        domain_->local_grid_x_begin(), domain_->local_grid_x_end(),
        chem.storage_nx(), chem.owned_storage_x_begin());
  } else {
    near_on_gpu = accumulate_near_field_gpu_or_cpu(
        *domain_, gf_, adv, sources, params, strength_factors,
        cfg_.toxin_cutoff, toxin_conc, chem, toxin_species_idx, chem_gpu,
        defer_sync);
  }
  if (near_on_gpu) {
    if (chem_gpu == nullptr) {
      throw SimulationError(
          "GPU near-field evaluation succeeded without a GPU chemical field");
    }
    chem_gpu->sync_species_concentrations_to_host(chem, toxin_species_idx);
    for (Int c = 0; c < ncells; ++c) {
      toxin_conc[static_cast<size_t>(c)] = chem.conc(toxin_species_idx, c);
    }
  }
  if (fmm.locals_ready()
      && gpu_accumulate_far_field_local(
          fmm, *domain_, cfg_.fmm_expansion_order, toxin_conc, toxin_conc)) {
    // GPU far-field deposit complete.
  } else {
    if (chem.slab_mode()) {
      FarFieldGridContext local_grid = far_grid;
      local_grid.x_begin = domain_->local_grid_x_begin();
      local_grid.x_end = domain_->local_grid_x_end();
      local_grid.storage_nx = chem.storage_nx();
      local_grid.halo_width = chem.owned_storage_x_begin();
      accumulate_far_field(fmm, gf_, avg_params, local_grid, toxin_conc,
                           cfg_.toxin_cutoff, false);
    } else {
      accumulate_far_field(fmm, gf_, avg_params, far_grid, toxin_conc,
                           cfg_.toxin_cutoff, false);
    }
  }
  deposit_to_chemical_field(chem, toxin_species_idx, toxin_conc);
  if (chem_gpu != nullptr && chem_gpu->active()) {
    chem_gpu->sync_species_concentrations_to_device(chem, toxin_species_idx);
  }
}

void QSSASolver::solve_nutrient_depletion(
    const AgentPool& agents,
    ChemicalField& chem,
    const OxygenConfig& oxygen) const {

  Int i_oxygen = chem.find(species::OXYGEN);

  const Real cell_vol = domain_->cell_volume();

  #ifdef GUTIBM_OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (const Agent& a : agents) {
    if (a.state == PhenoState::DEAD || a.flags.is_ghost) continue;

    Int cell = a.grid_cell;
    if (cell < 0) continue;

    // Spec 6 — per-agent carbon/iron/B12 depletion is handled solely by the
    // metabolism Fix (yield-based uptake in FixMetabolism::grow_agent). The
    // carbon/iron/B12 terms formerly applied here duplicated that uptake
    // (double-counting); they have been removed. Corrinoid (B12) is no longer
    // depleted at all. This function now applies only agent O2 respiration,
    // which has no counterpart in the metabolism Fix.
    if (oxygen.enabled && !oxygen.delivery_uptake_enabled
        && i_oxygen >= 0 && cell_vol > 0.0) {
      // Pirt respiration: growth-associated + basal maintenance. The
      // maintenance term is applied per living cell regardless of growth, so
      // the O2 field tracks agent density (a non-growing cell still respires).
      const Real respiratory_fraction = oxygen.metabolic_switch_enabled
          ? 1.0 - a.realized_fermentation_fraction : 1.0;
      const Real o2_use = oxygen.q_consumption
                            * std::max(a.mu_realized, 0.0)
                            * std::clamp(respiratory_fraction, 0.0, 1.0)
                        + oxygen.q_maintenance;
      #ifdef GUTIBM_OPENMP
      #pragma omp atomic
      #endif
      chem.reac_global(i_oxygen, cell) -= o2_use / cell_vol;
    }
  }
}

Real QSSASolver::point_concentration(
    const Vec3& target,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors) const {

  Real total = 0.0;
  for (size_t s = 0; s < sources.size(); ++s) {
    GreensFunctionParams p = params[s];
    p.source_rate *= strength_factors[s];
    Real d2 = domain_->min_image_dist_sq(sources[s], target);
    if (Real cutoff2 = cfg_.toxin_cutoff * cfg_.toxin_cutoff; d2 > cutoff2) continue;
    total += gf_.concentration_bounded(sources[s], target, p);
  }
  return total;
}

}  // namespace gutibm
