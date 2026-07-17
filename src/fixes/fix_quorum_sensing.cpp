/* -----------------------------------------------------------------------
   GutIBM – AI-2 quorum sensing: LuxS production + Lsr import (Spec 11)
   ----------------------------------------------------------------------- */

#include "fix_quorum_sensing.h"
#include "species_names.h"
#include "simulation.h"
#include <algorithm>

namespace gutibm {

FixQuorumSensing::FixQuorumSensing(Simulation& sim, const QuorumSensingConfig& cfg)
    : Fix("quorum_sensing", sim), cfg_(cfg) {}

void FixQuorumSensing::compute(Real /*dt*/) {
  if (!cfg_.enabled) return;

  auto& chem = sim_.chemical_field();
  const Int i_ai2 = chem.find(species::AI2);
  if (i_ai2 < 0) return;

  const Real dx = sim_.domain().dx();
  const Real cell_vol = dx * dx * dx;
  if (cell_vol <= 0.0) return;

  // Agent production (LuxS) and Lsr import
  for (const Agent& agent : sim_.agents()) {
    if (agent.state == PhenoState::DEAD) continue;
    if (agent.grid_cell < 0) continue;

    const Real production = cfg_.ai2_basal_rate
        + cfg_.ai2_growth_coupled * std::max(agent.mu_realized, 0.0);
    chem.reac(i_ai2, agent.grid_cell) += production / cell_vol;

    const Real ai2_conc = chem.conc(i_ai2, agent.grid_cell);
    const Real import_rate =
        cfg_.lsr_vmax * ai2_conc / (cfg_.lsr_km + ai2_conc);
    chem.reac(i_ai2, agent.grid_cell) -= import_rate / cell_vol;
  }

  // First-order background decay (uses ChemicalSpec.decay_rate when set)
  const Real decay = std::max(chem.spec(i_ai2).decay_rate, 0.0);
  if (decay > 0.0) {
    for (Int c = 0; c < chem.ncells(); ++c) {
      chem.reac(i_ai2, c) -= decay * chem.conc(i_ai2, c);
    }
  }
}

}  // namespace gutibm
