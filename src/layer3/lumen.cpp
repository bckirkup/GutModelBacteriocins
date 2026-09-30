/* -----------------------------------------------------------------------
   GutIBM – Spec 13 Phase 3: per-region luminal population compartment.
   ----------------------------------------------------------------------- */

#include "lumen.h"

#include "input_parser.h"
#include "metabolic_mode.h"

#include <algorithm>
#include <cmath>

namespace gutibm {

namespace {

uint64_t splitmix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

uint64_t hash_combine(uint64_t h, uint64_t v) {
  return splitmix64(h ^ v);
}

// Order-stable draw seed: identical for a given (step, region, lineage)
// regardless of processing order (G3).
uint64_t draw_seed(uint64_t seed, Int step_index, int region_index,
                   TagID lineage_id, uint64_t salt) {
  uint64_t h = seed;
  h = hash_combine(h, static_cast<uint64_t>(step_index));
  h = hash_combine(h, static_cast<uint64_t>(region_index));
  h = hash_combine(h, static_cast<uint64_t>(lineage_id));
  return hash_combine(h, salt);
}

}  // namespace

void LumenCompartment::init(const Layer3RegionConfig& rc,
                            const std::vector<Agent>& seed_exemplars,
                            Real seed_cells) {
  lineages_.clear();
  scalars_.carbon_mol_m3 = rc.lumen_carbon_mol_m3;
  scalars_.acetate_mol_m3 = 0.0;
  if (seed_cells <= 0.0 || seed_exemplars.empty()) return;
  const Real per_strain = seed_cells /
                          static_cast<Real>(seed_exemplars.size());
  for (const Agent& a : seed_exemplars) {
    LuminalLineage lin;
    lin.exemplar = a;
    lin.count = per_strain;
    lineages_[a.genome.lineage_id] = std::move(lin);
  }
}

void LumenCompartment::add_cells(const std::vector<Agent>& cells) {
  for (const Agent& a : cells) {
    auto& lin = lineages_[a.genome.lineage_id];
    if (lin.count <= 0.0) {
      lin.exemplar = a;  // new lineage: carry the genotype forward
    }
    lin.count += 1.0;
  }
}

LumenStepPlan LumenCompartment::compute(
    Real dt, const Layer3RegionConfig& rc, const Layer3Config& global,
    const LumenPhysio& physio, const LuminalScalars& upstream_scalars,
    Real occupancy_fraction, Int step_index, int region_index,
    uint64_t seed) const {
  LumenStepPlan plan;
  const Real tau = std::max(rc.transit_h, 1e-6) * 3600.0;
  const Real washout = 1.0 - std::exp(-dt / tau);
  const Real volume_m3 = std::max(rc.lumen_volume_L, 1e-9) * 1.0e-3;
  const Real carbon_cost = physio.yield_carbon * physio.carbon_cost_factor;
  const Real acid_block =
      metabolic_mode::acid_inhibition(scalars_.acetate_mol_m3, rc.ph,
                                      physio.acetate_pka,
                                      physio.acid_inhibition_ki,
                                      physio.acid_inhibition_max);
  const Real p_attach_cell =
      global.reattach_p0 * std::max(0.0, 1.0 - occupancy_fraction) *
      (dt / tau);

  Real carbon_consumed_mol = 0.0;
  Real mu_weighted = 0.0;
  for (const auto& [id, lin] : lineages_) {
    if (lin.count <= 0.0) continue;
    const Real mu = lin.exemplar.mu_max * physio.anaerobic_mu_factor *
                    metabolic_mode::monod(scalars_.carbon_mol_m3,
                                          lin.exemplar.km.km_carbon) *
                    (1.0 - acid_block);
    const Real grown = lin.count * std::exp(mu * dt);
    const Real births = grown - lin.count;
    const Real out = grown * washout;
    const Real resident = grown - out;
    plan.resident[id] = resident;
    if (out > 0.0) {
      plan.outflow.emplace(id, std::make_pair(out, lin.exemplar));
    }
    plan.births += births;
    plan.outflow_cells += out;
    mu_weighted += mu * grown;
    carbon_consumed_mol += births * lin.exemplar.biomass * carbon_cost;

    if (p_attach_cell > 0.0 && resident >= 1.0) {
      RNG draw;
      draw.seed(draw_seed(seed, step_index, region_index, id,
                          0xA77AC4ULL));
      const Int attempts = std::min<Int>(
          draw.poisson(p_attach_cell * resident),
          static_cast<Int>(resident));
      if (attempts > 0) {
        plan.reattach_events.emplace_back(id, attempts);
      }
    }
  }
  Real grown_total = plan.outflow_cells;
  for (const auto& [id, c] : plan.resident) {
    grown_total += c;
  }
  plan.mean_mu_per_s = grown_total > 0.0 ? mu_weighted / grown_total : 0.0;

  // Scalar transport: each region mixes toward its upstream inflow at
  // the transit rate, then growth consumption/acid production apply.
  const Real mix = washout;
  Real c_next = upstream_scalars.carbon_mol_m3 +
                (scalars_.carbon_mol_m3 - upstream_scalars.carbon_mol_m3) *
                    (1.0 - mix);
  c_next = std::max(0.0, c_next - carbon_consumed_mol / volume_m3);
  plan.scalars.carbon_mol_m3 = c_next;
  const Real a_next =
      upstream_scalars.acetate_mol_m3 +
      (scalars_.acetate_mol_m3 - upstream_scalars.acetate_mol_m3) *
          (1.0 - mix) +
      carbon_consumed_mol * physio.ferm_acid_yield / volume_m3;
  plan.scalars.acetate_mol_m3 = a_next;
  return plan;
}

void LumenCompartment::apply(const LumenStepPlan& plan,
                             const LumenStepPlan* upstream,
                             Real upstream_survival,
                             const std::map<TagID, Int>& deposited) {
  // Resident exemplars come from the pre-step state; upstream arrivals
  // carry their own. Map iteration is key-sorted in both passes, so
  // the result is permutation-invariant by construction (G3).
  std::map<TagID, LuminalLineage> next;
  for (const auto& [id, resident] : plan.resident) {
    Real count = resident;
    if (auto it = deposited.find(id); it != deposited.end()) {
      count -= static_cast<Real>(it->second);
    }
    if (count <= 0.0) continue;
    auto it = lineages_.find(id);
    LuminalLineage lin;
    lin.exemplar = it->second.exemplar;
    lin.count = count;
    next[id] = std::move(lin);
  }
  if (upstream != nullptr) {
    for (const auto& [id, pair] : upstream->outflow) {
      const Real arriving = pair.first * upstream_survival;
      if (arriving <= 0.0) continue;
      auto& lin = next[id];
      if (lin.count <= 0.0) {
        lin.exemplar = pair.second;
      }
      lin.count += arriving;
    }
  }
  lineages_ = std::move(next);
  scalars_ = plan.scalars;
}

Real LumenCompartment::purge(Real fraction) {
  Real removed = 0.0;
  for (auto& [id, lin] : lineages_) {
    removed += lin.count * fraction;
    lin.count *= (1.0 - fraction);
  }
  return removed;
}

Real LumenCompartment::total_cells() const {
  Real total = 0.0;
  for (const auto& [id, lin] : lineages_) {
    total += lin.count;
  }
  return total;
}

uint64_t LumenCompartment::fingerprint() const {
  uint64_t h = 0x9e3779b97f4a7c15ULL;
  for (const auto& [id, lin] : lineages_) {
    h = hash_combine(h, static_cast<uint64_t>(id));
    h = hash_combine(h, static_cast<uint64_t>(lin.count * 1e6));
    h = hash_combine(h, static_cast<uint64_t>(
                            lin.exemplar.genome.lineage_id));
  }
  h = hash_combine(h, static_cast<uint64_t>(scalars_.carbon_mol_m3 * 1e12));
  h = hash_combine(h,
                   static_cast<uint64_t>(scalars_.acetate_mol_m3 * 1e12));
  return h;
}

}  // namespace gutibm
