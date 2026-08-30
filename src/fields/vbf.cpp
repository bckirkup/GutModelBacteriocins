/* -----------------------------------------------------------------------
   GutIBM – VBF implementation
   ----------------------------------------------------------------------- */

#include "vbf.h"
#include "species_names.h"
#include "domain.h"
#include "chemical_field.h"
#include "chem_environment_config.h"
#include "vbf_carbon_sink.h"
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

Int agent_count_for_cell(const std::vector<Int>& agent_counts,
                         Int global_cell) {
  if (global_cell < 0
      || global_cell >= static_cast<Int>(agent_counts.size())) {
    return 0;
  }
  return agent_counts[static_cast<size_t>(global_cell)];
}

void apply_carbon_source(ChemicalField& chem, Int cell,
                         const VbfCellContext& ctx,
                         VbfFluxTotals* totals,
                         Real cell_volume,
                         Real dt) {
  if (ctx.cfg.use_dynamic_mucin && ctx.mucin.enabled && ctx.idx.mucin >= 0) {
    const Real liberation =
        dynamic_mucin_liberation(chem.conc(ctx.idx.mucin, cell), ctx.cfg, ctx.mucin);
    chem.reac(ctx.idx.mucin, cell) -= liberation;
    if (ctx.idx.carbon >= 0) {
      chem.reac(ctx.idx.carbon, cell) += liberation;
      if (totals != nullptr) {
        totals->carbon_source += liberation * cell_volume * dt;
      }
    }
    return;
  }
  if (ctx.idx.carbon >= 0) {
    chem.reac(ctx.idx.carbon, cell) += ctx.static_liberation;
    if (totals != nullptr) {
      totals->carbon_source += ctx.static_liberation * cell_volume * dt;
    }
  }
}

// Spec 5 §1 — VBF carbon consumption (Restaurant Hypothesis). Monod-saturating
// sink: near-complete uptake at high [C], leaving a thin residual at low [C].
void apply_carbon_sink(ChemicalField& chem, Int cell,
                       const VbfCellContext& ctx,
                       VbfFluxTotals* totals,
                       Real cell_volume,
                       Real dt,
                       Int agent_count) {
  if (ctx.idx.carbon < 0) return;
  const Real c = chem.conc(ctx.idx.carbon, cell);
  const Real vmax = ctx.cfg.carbon_sink_vmax
      + ctx.cfg.agent_carbon_coupling * static_cast<Real>(agent_count)
          / cell_volume;
  if (vmax <= 0.0) return;
  const Real sink = vbf::implicit_carbon_sink(
      c, vmax, ctx.cfg.carbon_sink_km, dt);
  chem.reac(ctx.idx.carbon, cell) -= sink;
  if (totals != nullptr) totals->carbon_sink += sink * cell_volume * dt;
}

void apply_iron_sink(ChemicalField& chem, Int cell,
                     const VbfCellContext& ctx,
                     VbfFluxTotals* totals,
                     Real cell_volume,
                     Real dt) {
  if (ctx.idx.iron < 0) return;
  // First-order (concentration-dependent) uptake: nutrient_sink is a rate
  // constant (1/s). See VBFConfig::nutrient_sink for why this is not a
  // zero-order mol/m^3/s removal.
  const Real sink = ctx.cfg.nutrient_sink * chem.conc(ctx.idx.iron, cell);
  chem.reac(ctx.idx.iron, cell) -= sink;
  if (totals != nullptr) totals->iron_sink += sink * cell_volume * dt;
}

void apply_oxygen_sink(ChemicalField& chem, Int cell,
                       const VbfCellContext& ctx,
                       VbfFluxTotals* totals,
                       Real cell_volume,
                       Real dt) {
  if (!ctx.oxygen.enabled || ctx.idx.oxygen < 0) return;
  if (const ChemicalSpec& oxygen_spec = chem.spec(ctx.idx.oxygen);
      oxygen_spec.delivery_enabled) {
    chem.add_vbf_sink_rate_global(
        ctx.idx.oxygen, cell, ctx.oxygen.vbf_sink);
    return;
  }
  // First-order background O2 consumption by the anaerobic majority:
  // reac -= vbf_sink * [O2] (1/s rate constant), mirroring the iron sink. A
  // zero-order (constant mol/m^3/s) removal removes O2 that may not be present,
  // driving the interior to a hard zero in a single bio step and masking the
  // per-agent respiration signal entirely — the density-tracking bug Edison
  // reported. A first-order sink is self-limiting: it scales with the local O2
  // it can actually consume, so a smooth gradient survives and agent
  // respiration remains visible on top of it.
  const Real sink = ctx.oxygen.vbf_sink * chem.conc(ctx.idx.oxygen, cell);
  chem.reac(ctx.idx.oxygen, cell) -= sink;
  if (totals != nullptr) totals->oxygen_sink += sink * cell_volume * dt;
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

void apply_vbf_at_cell(ChemicalField& chem, Int cell,
                       const VbfCellContext& ctx,
                       VbfFluxTotals* totals,
                       Real cell_volume,
                       Real dt,
                       Int agent_count) {
  apply_carbon_source(chem, cell, ctx, totals, cell_volume, dt);
  apply_carbon_sink(chem, cell, ctx, totals, cell_volume, dt, agent_count);
  apply_iron_sink(chem, cell, ctx, totals, cell_volume, dt);
  apply_oxygen_sink(chem, cell, ctx, totals, cell_volume, dt);
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
}

Real VBF::mucin_rate(Real z_rel) const {
  if (!cfg_.mucin_z_gradient_enabled) {
    return cfg_.mucin_liberation;
  }
  return cfg_.mucin_liberation * std::exp(-z_rel / cfg_.mucin_z_gradient_lambda);
}

void VBF::apply_nutrient_coupling(ChemicalField& chem, const Domain& domain,
                                   Real dt,
                                   const OxygenConfig& oxygen,
                                   const AcetateConfig& acetate,
                                   const MucinConfig& mucin,
                                   VbfFluxTotals* totals,
                                   const std::vector<Int>& agent_counts) const {
  const VbfSpeciesIndices idx = find_vbf_species(chem);

  const Int nx = domain.nx();
  const Int ny = domain.ny();
  const Int nz = domain.nz();
  const Real cell_volume = domain.cell_volume();

  if (chem.slab_mode()) {
    for (Int iz = 0; iz < nz; ++iz) {
      const Real z_rel = (iz + 0.5) * domain.dx_z();
      const Real z_weight = cfg_.mucin_z_gradient_enabled
          ? std::exp(-z_rel / cfg_.mucin_z_gradient_lambda) : 1.0;
      const Real static_liberation =
          cfg_.use_dynamic_mucin ? 0.0 : mucin_rate(z_rel);
      const VbfCellContext ctx{cfg_, oxygen, acetate, mucin, idx,
                               static_liberation, z_weight, iz};
      for (Int iy = 0; iy < ny; ++iy) {
        for (Int ix = domain.local_grid_x_begin();
             ix < domain.local_grid_x_end(); ++ix) {
          const Int global_cell = domain.cell_index(ix, iy, iz);
          const Int storage_cell = chem.global_to_storage_cell(global_cell);
          const Int agent_count =
              agent_count_for_cell(agent_counts, global_cell);
          apply_vbf_at_cell(chem, storage_cell, ctx, totals,
                            cell_volume, dt, agent_count);
        }
      }
    }
    return;
  }

  for (Int iz = 0; iz < nz; ++iz) {
    const Real z_rel = (iz + 0.5) * domain.dx_z();
    const Real z_weight = cfg_.mucin_z_gradient_enabled
        ? std::exp(-z_rel / cfg_.mucin_z_gradient_lambda)
        : 1.0;
    const Real static_liberation = cfg_.use_dynamic_mucin ? 0.0 : mucin_rate(z_rel);

    const VbfCellContext ctx{cfg_, oxygen, acetate, mucin, idx,
                             static_liberation, z_weight, iz};

    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        const Int global_cell = domain.cell_index(ix, iy, iz);
        apply_vbf_at_cell(
            chem, global_cell, ctx, totals, cell_volume, dt,
            agent_count_for_cell(agent_counts, global_cell));
      }
    }
  }
}

}  // namespace gutibm
