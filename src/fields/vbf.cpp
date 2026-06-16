/* -----------------------------------------------------------------------
   GutIBM – VBF implementation
   ----------------------------------------------------------------------- */

#include "vbf.h"
#include "domain.h"
#include "chemical_field.h"

namespace gutibm {

void VBF::init(const VBFConfig& cfg, const Domain& domain) {
  cfg_          = cfg;
  ncells_       = domain.ncells();
  carrying_cap_ = cfg.carrying_cap;
}

void VBF::apply_nutrient_coupling(ChemicalField& chem, Real dt) const {
  // Carbon source index (monosaccharides from mucin degradation)
  Int i_carbon = chem.find("carbon");
  // Iron sink (background consumption)
  Int i_iron   = chem.find("iron");

  for (Int c = 0; c < chem.ncells(); ++c) {
    // Monosaccharide liberation from mucin by anaerobic background
    if (i_carbon >= 0) {
      chem.reac(i_carbon, c) += cfg_.mucin_liberation;
    }
    // Background iron consumption (competitive depletion)
    if (i_iron >= 0) {
      chem.reac(i_iron, c) -= cfg_.nutrient_sink * chem.conc(i_iron, c);
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
