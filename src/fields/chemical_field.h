/* -----------------------------------------------------------------------
   GutIBM – Grid-based chemical concentration fields
   Stores nutrient and toxin concentrations on an Eulerian grid.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_CHEMICAL_microcin_penalty_applied_H
#define GUTIBM_CHEMICAL_microcin_penalty_applied_H

#include "types.h"
#include <cassert>
#include <string>
#include <string_view>
#include <vector>

namespace gutibm {

class Domain;

struct NutrientFluxAccounting {
  std::vector<Real> boundary_interval;
  std::vector<Real> boundary_cumulative;
  std::vector<Real> vbf_source_interval;
  std::vector<Real> vbf_source_cumulative;
  std::vector<Real> vbf_sink_interval;
  std::vector<Real> vbf_sink_cumulative;
  std::vector<Real> agent_uptake_interval;
  std::vector<Real> agent_uptake_step;
  std::vector<Real> agent_uptake_cumulative;
  std::vector<Real> reaction_clip_interval;
  std::vector<Real> reaction_clip_cumulative;

  void init(size_t species_count) {
    boundary_interval.assign(species_count, 0.0);
    boundary_cumulative.assign(species_count, 0.0);
    vbf_source_interval.assign(species_count, 0.0);
    vbf_source_cumulative.assign(species_count, 0.0);
    vbf_sink_interval.assign(species_count, 0.0);
    vbf_sink_cumulative.assign(species_count, 0.0);
    agent_uptake_interval.assign(species_count, 0.0);
    agent_uptake_step.assign(species_count, 0.0);
    agent_uptake_cumulative.assign(species_count, 0.0);
    reaction_clip_interval.assign(species_count, 0.0);
    reaction_clip_cumulative.assign(species_count, 0.0);
  }

  void add_interval(Int species, Real boundary, Real source, Real sink,
                    Real uptake) {
    const auto index = static_cast<size_t>(species);
    boundary_interval[index] += boundary;
    vbf_source_interval[index] += source;
    vbf_sink_interval[index] += sink;
    agent_uptake_interval[index] += uptake;
  }

  void add_boundary(Int species, Real amount) {
    boundary_interval[static_cast<size_t>(species)] += amount;
  }

  void add_agent_uptake(Int species, Real amount) {
    #ifdef GUTIBM_OPENMP
    #pragma omp atomic
    #endif
    agent_uptake_step[static_cast<size_t>(species)] += amount;
  }

  void add_reaction_clip(Int species, Real amount) {
    #ifdef GUTIBM_OPENMP
    #pragma omp atomic
    #endif
    reaction_clip_interval[static_cast<size_t>(species)] += amount;
  }

  void commit_agent_uptake_step() {
    for (size_t i = 0; i < agent_uptake_step.size(); ++i) {
      agent_uptake_interval[i] += agent_uptake_step[i];
      agent_uptake_step[i] = 0.0;
    }
  }

  void close_interval() {
    for (size_t i = 0; i < boundary_interval.size(); ++i) {
      boundary_cumulative[i] += boundary_interval[i];
      vbf_source_cumulative[i] += vbf_source_interval[i];
      vbf_sink_cumulative[i] += vbf_sink_interval[i];
      agent_uptake_cumulative[i] += agent_uptake_interval[i];
      reaction_clip_cumulative[i] += reaction_clip_interval[i];
      boundary_interval[i] = 0.0;
      vbf_source_interval[i] = 0.0;
      vbf_sink_interval[i] = 0.0;
      agent_uptake_interval[i] = 0.0;
      reaction_clip_interval[i] = 0.0;
    }
  }
};

struct ChemicalSpec {
  std::string name;
  Real diff_coeff = 0.0;       // diffusion coefficient in free water (m^2/s)
  Real retardation = 1.0;      // mucin retardation factor (effective D = D_free / retardation)
  Real initial_conc = 0.0;     // initial bulk concentration (mol/m^3)
  Real boundary_conc = 0.0;    // Dirichlet boundary (epithelial supply or luminal)
  Real decay_rate = 0.0;       // first-order decay (1/s)

  // z-dependent gradient: C(z) = C_max * exp(-z_rel / z_gradient_lambda)
  bool z_gradient_enabled = false;
  Real z_gradient_lambda  = 25.0e-6;  // characteristic decay length (m)

  // Stable backward-Euler directional diffusion for nutrients/small molecules.
  // Toxin fields remain false because QSSA Green's functions handle them.
  bool diffusion_enabled = false;
};

class ChemicalField {
 public:
  enum class DecompositionMode { Replicated, Slab };

  ChemicalField() = default;

  void init(const Domain& domain, const std::vector<ChemicalSpec>& specs,
            std::string_view decomposition = "replicated");

  Int num_species() const { return nspec_; }
  Int ncells() const { return ncells_; }
  Int global_ncells() const { return global_ncells_; }
  Int owned_ncells() const {
    return (owned_x_end_ - owned_x_begin_) * global_ny_ * global_nz_;
  }
  Int owned_storage_x_begin() const { return halo_width_; }
  Int owned_storage_x_end() const {
    return halo_width_ + owned_x_end_ - owned_x_begin_;
  }
  Int global_ny() const { return global_ny_; }
  Int global_nz() const { return global_nz_; }
  Int storage_nx() const { return storage_nx_; }

  // Concentration accessors [species][storage cell].
  Real conc(Int spec, Int cell) const {
    assert(cell >= 0 && cell < ncells_);
    return conc_[spec][cell];
  }
  Real& conc(Int spec, Int cell) {
    assert(cell >= 0 && cell < ncells_);
    return conc_[spec][cell];
  }
  Real conc_global(Int spec, Int cell) const {
    const Int storage_cell = global_to_storage_cell(cell);
    assert(storage_cell >= 0);
    return conc_[spec][storage_cell];
  }
  Real& conc_global(Int spec, Int cell) {
    const Int storage_cell = global_to_storage_cell(cell);
    assert(storage_cell >= 0);
    return conc_[spec][storage_cell];
  }

  // Reaction rate [species][storage cell] (mol/m^3/s, negative = consumption)
  Real reac(Int spec, Int cell) const {
    assert(cell >= 0 && cell < ncells_);
    return reac_[spec][cell];
  }
  Real& reac(Int spec, Int cell) {
    assert(cell >= 0 && cell < ncells_);
    return reac_[spec][cell];
  }
  Real reac_global(Int spec, Int cell) const {
    const Int storage_cell = global_to_storage_cell(cell);
    assert(storage_cell >= 0);
    return reac_[spec][storage_cell];
  }
  Real& reac_global(Int spec, Int cell) {
    const Int storage_cell = global_to_storage_cell(cell);
    assert(storage_cell >= 0);
    return reac_[spec][storage_cell];
  }

  // Global cell mapping. In slab mode, only owned cells and configured halo
  // cells map to storage; an out-of-range global cell is a programming error.
  Int owned_global_x_begin() const { return owned_x_begin_; }
  Int owned_global_x_end() const { return owned_x_end_; }
  Int grid_halo_width() const { return halo_width_; }
  Int global_to_storage_cell(Int global_cell) const;
  Int storage_to_global_cell(Int storage_cell) const;
  bool owns_global_cell(Int global_cell) const;
  bool global_cell_in_halo(Int global_cell) const;
  bool slab_mode() const { return mode_ == DecompositionMode::Slab; }

  // Make slab concentration halos current after chemistry and before the next
  // biology pass. Replicated mode is intentionally a no-op.
  void exchange_concentration_halos();

  // Reset reaction rates to zero each timestep
  void zero_reactions();

  // Apply stable implicit diffusion for enabled nutrient species.
  void apply_diffusion(const Domain& domain, Real dt);

  // Apply boundary conditions
  void apply_boundaries(const Domain& domain);

  // Sum rank-local agent reaction fields before spatial diffusion.
  void sum_reactions_across_ranks();
  void sum_agent_uptake_across_ranks();
  void sum_accounting_across_ranks();

  // Get species index by name
  Int find(std::string_view name) const;

  const ChemicalSpec& spec(Int i) const { return specs_[i]; }
  const std::vector<ChemicalSpec>& specs() const { return specs_; }
  NutrientFluxAccounting& flux_accounting() { return flux_accounting_; }
  const NutrientFluxAccounting& flux_accounting() const {
    return flux_accounting_;
  }

  // Raw data for HDF5 output
  const std::vector<std::vector<Real>>& conc_data() const { return conc_; }

 private:
  Int nspec_  = 0;
  Int ncells_ = 0;
  Int global_ncells_ = 0;
  Int global_nx_ = 0;
  Int global_ny_ = 0;
  Int global_nz_ = 0;
  Int owned_x_begin_ = 0;
  Int owned_x_end_ = 0;
  Int halo_width_ = 0;
  Int storage_nx_ = 0;
  DecompositionMode mode_ = DecompositionMode::Replicated;
  const Domain* domain_ = nullptr;
  std::vector<ChemicalSpec> specs_;
  std::vector<std::vector<Real>> conc_;   // [nspec][ncells]
  std::vector<std::vector<Real>> reac_;   // [nspec][ncells]
  NutrientFluxAccounting flux_accounting_;

  void apply_diffusion_slab(const Domain& domain, Real dt);
  void apply_boundaries_slab(const Domain& domain);
};

}  // namespace gutibm

#endif  // GUTIBM_CHEMICAL_microcin_penalty_applied_H
