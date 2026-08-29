/* -----------------------------------------------------------------------
   GutIBM – Quasi-Steady-State Approximation (QSSA) diffusion solver
   
   Because small-molecule diffusion (microseconds) is orders of
   magnitude faster than cell division and advection (minutes–hours),
   we decouple these timescales:
   
   1. At each biological timestep, collect all active toxin sources
   2. Compute steady-state concentration field via superposition
      of Green's function kernels
   3. Apply the resulting field to agent receptor-binding kinetics
   
   This avoids the CFL stability limitation of explicit FTCS solvers
   that would require sub-millisecond timesteps at 1 um resolution.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_QSSA_SOLVER_H
#define GUTIBM_QSSA_SOLVER_H

#include "types.h"
#include "greens_function.h"
#include "robin_correction_table.h"
#include "fmm.h"
#include "chem_environment_config.h"
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace gutibm {

class Domain;
class AdvectionField;
class ChemicalField;
class ChemicalFieldGpu;
class AgentPool;

struct QSSAConfig {
  // Cutoff radius for Green's function evaluation (m)
  // Beyond this, contribution is negligible
  Real toxin_cutoff     = 200.0e-6;   // 200 um for lethal halos
  // Lumen unstirred-layer transfer length; infinity disables the Robin
  // correction and preserves the sealed boundary result.
  Real lumen_transfer_length = robin::kZeroTransferLength;
  // effective uses k_c = D_eff / length; free uses D_free / length.
  std::string lumen_transfer_basis = "effective";
  Real nutrient_cutoff  = 50.0e-6;    // 50 um for nutrient depletion zones

  // Bacteriocin source parameters
  Real colicin_release_rate = 1.0e-18; // mol/s per lysed cell (burst)
  Real microcin_secretion   = 1.0e-20; // mol/s continuous secretion

  // FMM / Barnes-Hut acceleration (far-field aggregation)
  bool use_fmm  = false;   // enable octree far-field acceleration
  std::string toxin_evaluation = "grid";
  std::string toxin_lumping = "per_receptor";
  Real fmm_theta = 0.5;    // opening angle (0→exact, 1→fast/approximate)
  int  fmm_expansion_order = 2;  // 1=monopole, 2=dipole+quadrupole, 3=octupole

  // Spec 6 — per-agent carbon/iron/B12 depletion is now owned solely by the
  // metabolism Fix (yield-based uptake). The former iron/b12/carbon
  // stoichiometry knobs, which drove the removed QSSA nutrient terms, no longer
  // exist. solve_nutrient_depletion applies only agent O2 respiration.

  // Fallback averaged Green's function params used when there are toxin
  // sources but their total weighted strength is zero (degenerate case).
  Real fallback_diff_coeff  = 4.0e-11;
  Real fallback_pI          = 7.0;
  Real fallback_retardation = 5.0;
};

// Persistent burst from SOS lysis with finite inventory and release timescale.
struct ToxinBurstSource {
  Vec3 pos;
  GreensFunctionParams params;
  Real creation_time = 0.0;
  Real release_tau = 300.0;  // exponential release timescale (s)
  bool is_nuclease = false;
  ReceptorType target = ReceptorType::BtuB;
};

class QSSASolver {
 public:
  QSSASolver() = default;

  void init(const QSSAConfig& cfg, const Domain& domain,
            const AdvectionField& adv);

  // Compute steady-state bacteriocin field from current toxin sources
  // and deposit onto chemical field grid
  void solve_bacteriocin_field(
      const AgentPool& agents,
      const std::vector<ToxinBurstSource>& bursts,
      Real current_time,
      const ProteaseConfig& protease,
      const AdvectionField& adv,
      ChemicalField& chem,
      Int toxin_species_idx,
      ReceptorType target,
      ChemicalFieldGpu* chem_gpu = nullptr,
      bool materialize_grid = true);

  // Solve all four per-receptor bacteriocin fields.
  void solve_all_bacteriocin_fields(
      const AgentPool& agents,
      const std::vector<ToxinBurstSource>& bursts,
      Real current_time,
      const ProteaseConfig& protease,
      const AdvectionField& adv,
      ChemicalField& chem,
      ChemicalFieldGpu* chem_gpu = nullptr,
      bool materialize_grid = true);

  // Compute nutrient depletion zones around colonies
  void solve_nutrient_depletion(
      const AgentPool& agents,
      ChemicalField& chem,
      const OxygenConfig& oxygen) const;

  // Concentration at a specific point from all nearby sources
  Real point_concentration(
      const Vec3& target,
      const std::vector<Vec3>& sources,
      const std::vector<GreensFunctionParams>& params,
      const std::vector<Real>& strength_factors) const;

  Real sampled_toxin_conc(Int agent_index, Int species_idx) const;
  Real sampled_nuclease_conc(Int agent_index) const;
  Real sampled_nuclease_conc(Int agent_index, ReceptorType target) const;
  Real sampled_toxin_max(Int species_idx) const;
  bool agent_sampling() const { return cfg_.toxin_evaluation == "agents"; }
  bool toxin_lumping() const { return cfg_.toxin_lumping == "lumped"; }

  const GreensFunction& gf() const { return gf_; }

 private:
  struct SampledToxinField {
    std::vector<Vec3> sources;
    std::vector<GreensFunctionParams> params;
    std::vector<Real> strength_factors;
    std::vector<Real> samples;
    FMM fmm;
    bool fmm_ready = false;
  };

  void solve_bacteriocin_field_fmm(
      const std::vector<Vec3>& sources,
      const std::vector<GreensFunctionParams>& params,
      const std::vector<Real>& strength_factors,
      const AdvectionField& adv,
      ChemicalField& chem,
      Int toxin_species_idx,
      ChemicalFieldGpu* chem_gpu = nullptr) const;

  void sample_bacteriocin_field(
      const std::vector<Vec3>& sources,
      const std::vector<GreensFunctionParams>& params,
      const std::vector<Real>& strength_factors,
      const AgentPool& agents,
      ReceptorType target);
  void sample_nuclease_field(
      const std::vector<Vec3>& sources,
      const std::vector<GreensFunctionParams>& params,
      const std::vector<Real>& strength_factors,
      const AgentPool& agents,
      ReceptorType target);
  void sample_nuclease_sources(
      const std::vector<Vec3>& sources,
      const std::vector<GreensFunctionParams>& params,
      const std::vector<Real>& strength_factors,
      const std::vector<bool>& is_nuclease,
      const std::vector<ReceptorType>& targets,
      const AgentPool& agents);

  // Lumped mode has one field over all sources, so both entry points do the
  // same work regardless of which receptor target was requested.
  void solve_lumped_bacteriocin_fields(
      const AgentPool& agents,
      const std::vector<Vec3>& sources,
      const std::vector<GreensFunctionParams>& params,
      const std::vector<Real>& strength_factors,
      const AdvectionField& adv,
      ChemicalField& chem,
      ChemicalFieldGpu* chem_gpu,
      bool materialize_grid);

  Int sampled_slot_for_species(Int species_idx) const;
  Real evaluate_sample(const SampledToxinField& field,
                       const Vec3& position) const;

  void solve_bacteriocin_field_from_sources(
      const std::vector<Vec3>& sources,
      const std::vector<GreensFunctionParams>& params,
      const std::vector<Real>& strength_factors,
      const AdvectionField& adv,
      ChemicalField& chem,
      Int toxin_species_idx,
      ChemicalFieldGpu* chem_gpu = nullptr) const;

  QSSAConfig cfg_;
  GreensFunction gf_;
  const Domain* domain_ = nullptr;
  const AdvectionField* adv_ = nullptr;
  const AgentPool* sampled_agents_ = nullptr;
  mutable std::array<SampledToxinField, 4> sampled_fields_;
  mutable std::array<SampledToxinField, 4> sampled_nuclease_fields_;
  bool sampled_nuclease_sources_ = false;
  std::array<Int, 4> sampled_species_indices_{{-1, -1, -1, -1}};
};

}  // namespace gutibm

#endif  // GUTIBM_QSSA_SOLVER_H
