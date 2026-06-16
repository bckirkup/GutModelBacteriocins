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
      ChemicalField& chem,
      Int toxin_species_idx) const;

  // Compute nutrient depletion zones around colonies
  void solve_nutrient_depletion(
      const AgentPool& agents,
      ChemicalField& chem) const;

  // Concentration at a specific point from all nearby sources
  Real point_concentration(
      const Vec3& target,
      const std::vector<Vec3>& sources,
      const std::vector<GreensFunctionParams>& params) const;

  const GreensFunction& gf() const { return gf_; }

 private:
  QSSAConfig cfg_;
  GreensFunction gf_;
  const Domain* domain_ = nullptr;
};

}  // namespace gutibm

#endif  // GUTIBM_QSSA_SOLVER_H
