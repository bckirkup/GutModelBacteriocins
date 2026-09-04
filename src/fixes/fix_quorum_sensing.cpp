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

  const Real cell_vol = sim_.domain().cell_volume();
  if (cell_vol <= 0.0) return;

  // Agent production (LuxS) and Lsr import
  for (const Agent& agent : sim_.agents()) {
    if (agent.state == PhenoState::DEAD || agent.flags.is_ghost) continue;
    if (agent.grid_cell < 0) continue;

    const Real production = cfg_.ai2_basal_rate
        + cfg_.ai2_growth_coupled * std::max(agent.mu_realized, 0.0);
    chem.reac_global(i_ai2, agent.grid_cell) += production / cell_vol;

    const Real ai2_conc = chem.conc_global(i_ai2, agent.grid_cell);
    const Real import_rate =
        cfg_.lsr_vmax * ai2_conc / (cfg_.lsr_km + ai2_conc);
    chem.reac_global(i_ai2, agent.grid_cell) -= import_rate / cell_vol;
  }

  // First-order background decay (uses ChemicalSpec.decay_rate when set)
  const Real decay = std::max(chem.spec(i_ai2).decay_rate, 0.0);
  if (decay > 0.0) {
    const Int x_begin = chem.slab_mode()
        ? chem.owned_storage_x_begin()
        : sim_.domain().local_grid_x_begin();
    const Int x_end = chem.slab_mode()
        ? chem.owned_storage_x_end()
        : sim_.domain().local_grid_x_end();
    const Int storage_nx = chem.slab_mode()
        ? chem.storage_nx() : sim_.domain().nx();
    for (Int iz = 0; iz < chem.global_nz(); ++iz) {
      for (Int iy = 0; iy < chem.global_ny(); ++iy) {
        for (Int ix = x_begin; ix < x_end; ++ix) {
          const Int c = iz * storage_nx * chem.global_ny()
              + iy * storage_nx + ix;
          chem.reac(i_ai2, c) -= decay * chem.conc(i_ai2, c);
        }
      }
    }
  }
  chem.mark_host_reac_dirty(i_ai2);
}

}  // namespace gutibm
