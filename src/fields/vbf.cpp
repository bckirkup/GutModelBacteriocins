/* -----------------------------------------------------------------------
   GutIBM – VBF implementation
   ----------------------------------------------------------------------- */

#include "vbf.h"
#include "domain.h"
#include "chemical_field.h"
#include <cmath>

namespace gutibm {

void VBF::init(const VBFConfig& cfg, const Domain& domain) {
  cfg_          = cfg;
  ncells_       = domain.ncells();
  carrying_cap_ = cfg.carrying_cap;
}

Real VBF::mucin_rate(Real z_rel) const {
  if (!cfg_.mucin_z_gradient_enabled) {
    return cfg_.mucin_liberation;
  }
  return cfg_.mucin_liberation * std::exp(-z_rel / cfg_.mucin_z_gradient_lambda);
}

void VBF::apply_nutrient_coupling(ChemicalField& chem, const Domain& domain,
                                   Real /*dt*/) const {
  Int i_carbon = chem.find("carbon");
  Int i_iron   = chem.find("iron");

  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();

  for (Int iz = 0; iz < nz; ++iz) {
    Real z_rel = (iz + 0.5) * domain.dx();
    Real liberation = mucin_rate(z_rel);

    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        Int c = domain.cell_index(ix, iy, iz);

        if (i_carbon >= 0) {
          chem.reac(i_carbon, c) += liberation;
        }
        if (i_iron >= 0) {
          chem.reac(i_iron, c) -= cfg_.nutrient_sink * chem.conc(i_iron, c);
        }
      }
    }
  }
}

Vec3 VBF::drag_force(const Vec3& agent_vel) const {
  return {
    -cfg_.drag_coeff * agent_vel[0],
    -cfg_.drag_coeff * agent_vel[1],
    -cfg_.drag_coeff * agent_vel[2]
  };
}

}  // namespace gutibm
