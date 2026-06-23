/* -----------------------------------------------------------------------
   GutIBM – Soft-sphere mechanical repulsion (Hertzian contact model)

   Implements proper Hertzian contact mechanics (F = k * overlap^(3/2))
   between overlapping cells. Optionally includes EPS-mediated adhesion
   for biofilm-like cell clustering.

   Replaces the simplistic linear overlap force previously inline in
   simulation.cpp.

   References:
   - Hertz (1882) contact theory: F = (4/3)*E*sqrt(R)*delta^(3/2)
   - Typical bacterial elastic modulus: ~0.1-1 MPa (AFM measurements)
   - NUFEB heritage: nufeb/fix_eps_extract.cpp

   Issue: #16
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_FIX_MECHANICS_H
#define GUTIBM_FIX_MECHANICS_H

#include "fix.h"
#include "agent.h"

namespace gutibm {

struct MechanicsConfig {
  Real hertz_k = 1.0e-6;            // Hertzian spring constant (N/m^1.5)
  bool hertzian_enabled = true;      // Use Hertzian (true) vs linear (false)
  bool adhesion_enabled = false;     // EPS-mediated adhesion
  Real adhesion_strength = 1.0e-12;  // EPS adhesion force magnitude (N)
  Real adhesion_range = 0.5e-6;      // Adhesion range beyond contact (m)
};

class FixMechanics : public Fix {
 public:
  FixMechanics(Simulation& sim, const MechanicsConfig& cfg)
      : Fix("mechanics", sim), cfg_(cfg) {}

  void compute(Real dt) override;

 private:
  MechanicsConfig cfg_;
};

}  // namespace gutibm

#endif  // GUTIBM_FIX_MECHANICS_H
