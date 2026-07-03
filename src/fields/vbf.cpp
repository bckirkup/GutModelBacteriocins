/* -----------------------------------------------------------------------
   GutIBM – VBF implementation
   ----------------------------------------------------------------------- */

#include "vbf.h"
#include "species_names.h"
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

struct VbfSpeciesIndices {
  Int carbon = -1;
  Int iron = -1;
  Int oxygen = -1;
  Int acetate = -1;
  Int mucin = -1;
};

struct VbfCellContext {
  const VBFConfig& cfg;
  const OxygenConfig& oxygen;
  const AcetateConfig& acetate;
  const MucinConfig& mucin;
  VbfSpeciesIndices idx;
  Real static_liberation;
  Real z_weight;
  Int iz;
};

void apply_carbon_source(ChemicalField& chem, Int cell, const VbfCellContext& ctx) {
  if (ctx.cfg.use_dynamic_mucin && ctx.mucin.enabled && ctx.idx.mucin >= 0) {
    const Real liberation =
        dynamic_mucin_liberation(chem.conc(ctx.idx.mucin, cell), ctx.cfg, ctx.mucin);
    chem.reac(ctx.idx.mucin, cell) -= liberation;
    if (ctx.idx.carbon >= 0) {
      chem.reac(ctx.idx.carbon, cell) += liberation;
    }
    return;
  }
  if (ctx.idx.carbon >= 0) {
    chem.reac(ctx.idx.carbon, cell) += ctx.static_liberation;
  }
}

void apply_iron_sink(ChemicalField& chem, Int cell, const VbfCellContext& ctx) {
  if (ctx.idx.iron < 0) return;
  chem.reac(ctx.idx.iron, cell) -= ctx.cfg.nutrient_sink * chem.conc(ctx.idx.iron, cell);
}

void apply_oxygen_sink(ChemicalField& chem, Int cell, const VbfCellContext& ctx) {
  if (!ctx.oxygen.enabled || ctx.idx.oxygen < 0) return;
  chem.reac(ctx.idx.oxygen, cell) -= ctx.oxygen.vbf_sink;
}

void apply_acetate_coupling(ChemicalField& chem, Int cell, const VbfCellContext& ctx) {
  if (!ctx.acetate.enabled || ctx.idx.acetate < 0) return;
  chem.reac(ctx.idx.acetate, cell) += ctx.acetate.vbf_production * ctx.z_weight;
  chem.reac(ctx.idx.acetate, cell) -= ctx.acetate.vbf_consumption;
  if (ctx.iz == 0) {
    chem.reac(ctx.idx.acetate, cell) -= ctx.acetate.epithelial_uptake;
  }
}

void apply_mucin_secretion(ChemicalField& chem, Int cell, const VbfCellContext& ctx) {
  if (!ctx.mucin.enabled || ctx.idx.mucin < 0 || ctx.iz != 0) return;
  chem.reac(ctx.idx.mucin, cell) += ctx.mucin.secretion_rate;
}

void apply_vbf_at_cell(ChemicalField& chem, Int cell, const VbfCellContext& ctx) {
  apply_carbon_source(chem, cell, ctx);
  apply_iron_sink(chem, cell, ctx);
  apply_oxygen_sink(chem, cell, ctx);
  apply_acetate_coupling(chem, cell, ctx);
  apply_mucin_secretion(chem, cell, ctx);
}

VbfSpeciesIndices find_vbf_species(const ChemicalField& chem) {
  return {
    chem.find(species::CARBON),
    chem.find(species::IRON),
    chem.find(species::OXYGEN),
    chem.find(species::ACETATE),
    chem.find(species::MUCIN),
  };
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
  const VbfSpeciesIndices idx = find_vbf_species(chem);

  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();

  for (Int iz = 0; iz < nz; ++iz) {
    const Real z_rel = (iz + 0.5) * domain.dx();
    const Real z_weight = cfg_.mucin_z_gradient_enabled
        ? std::exp(-z_rel / cfg_.mucin_z_gradient_lambda)
        : 1.0;
    const Real static_liberation = cfg_.use_dynamic_mucin ? 0.0 : mucin_rate(z_rel);

    const VbfCellContext ctx{cfg_, oxygen, acetate, mucin, idx,
                             static_liberation, z_weight, iz};

    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        apply_vbf_at_cell(chem, domain.cell_index(ix, iy, iz), ctx);
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
