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
#include "fmm.h"
#include "chem_environment_config.h"
#include <vector>

namespace gutibm {

class Domain;
class AdvectionField;
class ChemicalField;
class AgentPool;

struct QSSAConfig {
  // Cutoff radius for Green's function evaluation (m)
  // Beyond this, contribution is negligible
  Real toxin_cutoff     = 200.0e-6;   // 200 um for lethal halos
  Real nutrient_cutoff  = 50.0e-6;    // 50 um for nutrient depletion zones

  // Bacteriocin source parameters
  Real colicin_release_rate = 1.0e-18; // mol/s per lysed cell (burst)
  Real microcin_secretion   = 1.0e-20; // mol/s continuous secretion

  // FMM / Barnes-Hut acceleration (far-field aggregation)
  bool use_fmm  = false;   // enable octree far-field acceleration
  Real fmm_theta = 0.5;    // opening angle (0→exact, 1→fast/approximate)
  int  fmm_expansion_order = 2;  // 1=monopole, 2=dipole+quadrupole, 3=octupole
};

// Persistent burst from SOS lysis (Spec 1 protease decay)
struct ToxinBurstSource {
  Vec3 pos;
  GreensFunctionParams params;
  Real creation_time = 0.0;
  Real decay_rate = 0.0;   // ln(2) / protease_half_life (1/s)
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
      Int toxin_species_idx) const;

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

  const GreensFunction& gf() const { return gf_; }

 private:
  void solve_bacteriocin_field_fmm(
      const std::vector<Vec3>& sources,
      const std::vector<GreensFunctionParams>& params,
      const std::vector<Real>& strength_factors,
      ChemicalField& chem,
      Int toxin_species_idx) const;

  QSSAConfig cfg_;
  GreensFunction gf_;
  const Domain* domain_ = nullptr;
  const AdvectionField* adv_ = nullptr;
};

}  // namespace gutibm

#endif  // GUTIBM_QSSA_SOLVER_H
