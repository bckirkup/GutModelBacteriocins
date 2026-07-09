/* -----------------------------------------------------------------------
   GutIBM – Grid-based chemical concentration fields
   Stores nutrient and toxin concentrations on an Eulerian grid.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_CHEMICAL_microcin_penalty_applied_H
#define GUTIBM_CHEMICAL_microcin_penalty_applied_H

#include "types.h"
#include <string>
#include <string_view>
#include <vector>

namespace gutibm {

class Domain;

struct ChemicalSpec {
  std::string name;
  Real diff_coeff;       // diffusion coefficient in free water (m^2/s)
  Real retardation;      // mucin retardation factor (effective D = D_free / retardation)
  Real initial_conc;     // initial bulk concentration (mol/m^3)
  Real boundary_conc;    // Dirichlet boundary (epithelial supply or luminal)
  Real decay_rate;       // first-order decay (1/s)

  // z-dependent gradient: C(z) = C_max * exp(-z_rel / z_gradient_lambda)
  bool z_gradient_enabled = false;
  Real z_gradient_lambda  = 25.0e-6;  // characteristic decay length (m)

  // Stable backward-Euler directional diffusion for nutrients/small molecules.
  // Toxin fields remain false because QSSA Green's functions handle them.
  bool diffusion_enabled = false;
};

class ChemicalField {
 public:
  ChemicalField() = default;

  void init(const Domain& domain, const std::vector<ChemicalSpec>& specs);

  Int num_species() const { return nspec_; }
  Int ncells() const { return ncells_; }

  // Concentration accessors  [species][cell]
  Real  conc(Int spec, Int cell) const { return conc_[spec][cell]; }
  Real& conc(Int spec, Int cell)       { return conc_[spec][cell]; }

  // Reaction rate  [species][cell]  (mol/m^3/s, negative = consumption)
  Real  reac(Int spec, Int cell) const { return reac_[spec][cell]; }
  Real& reac(Int spec, Int cell)       { return reac_[spec][cell]; }

  // Reset reaction rates to zero each timestep
  void zero_reactions();

  // Apply stable implicit diffusion for enabled nutrient species.
  void apply_diffusion(const Domain& domain, Real dt);

  // Apply boundary conditions
  void apply_boundaries(const Domain& domain);

  // Sum rank-local agent reaction fields before spatial diffusion.
  void sum_reactions_across_ranks();

  // Get species index by name
  Int find(std::string_view name) const;

  const ChemicalSpec& spec(Int i) const { return specs_[i]; }
  const std::vector<ChemicalSpec>& specs() const { return specs_; }

  // Raw data for HDF5 output
  const std::vector<std::vector<Real>>& conc_data() const { return conc_; }

 private:
  Int nspec_  = 0;
  Int ncells_ = 0;
  std::vector<ChemicalSpec> specs_;
  std::vector<std::vector<Real>> conc_;   // [nspec][ncells]
  std::vector<std::vector<Real>> reac_;   // [nspec][ncells]
};

}  // namespace gutibm

#endif  // GUTIBM_CHEMICAL_microcin_penalty_applied_H
