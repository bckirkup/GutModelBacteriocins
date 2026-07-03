/* -----------------------------------------------------------------------
   GutIBM – VBF implementation
   ----------------------------------------------------------------------- */

#include "vbf.h"
#include "domain.h"
#include "chemical_field.h"
#include "chem_environment_config.h"
#include <cmath>

namespace gutibm {

namespace {

Real dynamic_mucin_liberation(Real mucin_conc,
                              const VBFConfig& vbf_cfg,
                              const MucinConfig& mucin_cfg) {
  const Real substrate = mucin_conc / (mucin_cfg.Km_degradation + mucin_conc);
  return mucin_cfg.k_liberation * vbf_cfg.density * substrate;
}

}  // namespace

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
                                   Real /*dt*/,
                                   const OxygenConfig& oxygen,
                                   const AcetateConfig& acetate,
                                   const MucinConfig& mucin) const {
  Int i_carbon = chem.find("carbon");
  Int i_iron   = chem.find("iron");
  Int i_oxygen = chem.find("oxygen");
  Int i_acetate = chem.find("acetate");
  Int i_mucin  = chem.find("mucin");

  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();

  for (Int iz = 0; iz < nz; ++iz) {
    Real z_rel = (iz + 0.5) * domain.dx();
    Real z_weight = cfg_.mucin_z_gradient_enabled
        ? std::exp(-z_rel / cfg_.mucin_z_gradient_lambda)
        : 1.0;

    const Real static_liberation = cfg_.use_dynamic_mucin ? 0.0 : mucin_rate(z_rel);

    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        Int c = domain.cell_index(ix, iy, iz);

        if (cfg_.use_dynamic_mucin && mucin.enabled && i_mucin >= 0) {
          const Real liberation =
              dynamic_mucin_liberation(chem.conc(i_mucin, c), cfg_, mucin);
          chem.reac(i_mucin, c) -= liberation;
          if (i_carbon >= 0) {
            chem.reac(i_carbon, c) += liberation;
          }
        } else if (i_carbon >= 0) {
          chem.reac(i_carbon, c) += static_liberation;
        }

        if (i_iron >= 0) {
          chem.reac(i_iron, c) -= cfg_.nutrient_sink * chem.conc(i_iron, c);
        }

        if (oxygen.enabled && i_oxygen >= 0) {
          chem.reac(i_oxygen, c) -= oxygen.vbf_sink;
        }

        if (acetate.enabled && i_acetate >= 0) {
          chem.reac(i_acetate, c) += acetate.vbf_production * z_weight;
          chem.reac(i_acetate, c) -= acetate.vbf_consumption;
          if (iz == 0) {
            chem.reac(i_acetate, c) -= acetate.epithelial_uptake;
          }
        }

        if (mucin.enabled && i_mucin >= 0 && iz == 0) {
          chem.reac(i_mucin, c) += mucin.secretion_rate;
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
