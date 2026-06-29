/* -----------------------------------------------------------------------
   GutIBM – Chemical field implementation
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "domain.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace gutibm {

void ChemicalField::init(const Domain& domain,
                          const std::vector<ChemicalSpec>& specs) {
  specs_  = specs;
  nspec_  = static_cast<Int>(specs.size());
  ncells_ = domain.ncells();

  conc_.resize(nspec_);
  reac_.resize(nspec_);
  for (Int s = 0; s < nspec_; ++s) {
    conc_[s].assign(ncells_, specs_[s].initial_conc);
    reac_[s].assign(ncells_, 0.0);

    if (specs_[s].z_gradient_enabled) {
      const Int nx = domain.nx();
      const Int ny = domain.ny();
      const Int nz = domain.nz();
      for (Int iz = 0; iz < nz; ++iz) {
        Real z_rel = (iz + 0.5) * domain.dx();
        Real factor = std::exp(-z_rel / specs_[s].z_gradient_lambda);
        for (Int iy = 0; iy < ny; ++iy) {
          for (Int ix = 0; ix < nx; ++ix) {
            Int cell = domain.cell_index(ix, iy, iz);
            conc_[s][cell] = specs_[s].initial_conc * factor;
          }
        }
      }
    }
  }
}

void ChemicalField::zero_reactions() {
  for (Int s = 0; s < nspec_; ++s) {
    std::ranges::fill(reac_[s], 0.0);
  }
}

void ChemicalField::apply_boundaries(const Domain& domain) {
  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();

  for (Int s = 0; s < nspec_; ++s) {
    Real bc = specs_[s].boundary_conc;

    // z=0 (epithelial surface): Dirichlet for nutrients
    // When z_gradient_enabled, bc is the peak concentration at the epithelium
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        Int idx = domain.cell_index(ix, iy, 0);
        conc_[s][idx] = bc;
      }
    }

    // z=nz-1 (luminal surface): open boundary (zero-gradient Neumann)
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        Int top = domain.cell_index(ix, iy, nz - 1);
        Int below = domain.cell_index(ix, iy, nz - 2);
        conc_[s][top] = conc_[s][below];
      }
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
