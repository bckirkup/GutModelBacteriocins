#include "chemistry_pipeline.h"
#include "agent.h"
#include "chemical_field.h"
#include "chemical_field_gpu.h"
#include "dispatch.h"
#include "qssa_gpu.h"
#include "qssa_solver.h"
#include "vbf.h"
#include "vbf_gpu.h"

namespace gutibm {

ChemistryPipelineResult run_chemistry_pipeline(ChemistryPipelineInput& in, Real dt) {
  ChemistryPipelineResult result;

  bool applied_o2_on_gpu = false;
  if (in.gpu_active) {
    applied_o2_on_gpu = gpu_solve_nutrient_depletion(
        in.agents_gpu, in.num_agents, in.chem_gpu, in.chem,
        in.oxygen, in.domain);
  }

  if (!applied_o2_on_gpu) {
    if (in.gpu_active) {
      in.chem_gpu.sync_reactions_to_host(in.chem);
    }
    in.qssa.solve_nutrient_depletion(in.agents, in.chem, in.oxygen);
  } else if (in.gpu_active) {
    in.chem_gpu.sync_reactions_to_host(in.chem);
  }

  // Every rank holds the full chemical grid but only its local agents. Sum the
  // rank-local agent reaction fields before adding the identical global VBF.
  in.chem.sum_reactions_across_ranks();

  bool reactions_on_device = false;
  bool applied_vbf_on_gpu = false;
  if (in.gpu_active && applied_o2_on_gpu) {
    in.chem_gpu.sync_reactions_to_device(in.chem);
    reactions_on_device = true;
    applied_vbf_on_gpu = gpu_apply_vbf_coupling(
        in.chem_gpu, in.chem, in.domain, in.vbf,
        in.oxygen, in.acetate, in.mucin);
  }
  if (!applied_vbf_on_gpu) {
    in.vbf.apply_nutrient_coupling(in.chem, in.domain, dt,
                                   in.oxygen, in.acetate, in.mucin);
    if (in.gpu_active) {
      in.chem_gpu.sync_reactions_to_device(in.chem);
      reactions_on_device = true;
    }
  }

  if (in.gpu_active) {
    if (!reactions_on_device) {
      in.chem_gpu.sync_reactions_to_device(in.chem);
    }
    result.reactions_on_gpu = in.chem_gpu.apply_reactions(dt, in.domain);
  }

  if (!result.reactions_on_gpu) {
    Int s = 0;
    for (const auto& conc_row : in.chem.conc_data()) {
      (void)conc_row;
      #ifdef GUTIBM_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (Int c = 0; c < in.chem.ncells(); ++c) {
        in.chem.conc(s, c) += in.chem.reac(s, c) * dt;
        in.chem.conc(s, c) = std::max(in.chem.conc(s, c), 0.0);
      }
      ++s;
    }
  }

  if (in.gpu_active && result.reactions_on_gpu) {
    result.diffusion_on_gpu =
        in.chem_gpu.apply_diffusion(in.domain, in.chem, dt);
    if (result.diffusion_on_gpu) {
      in.chem_gpu.apply_boundaries(in.domain, in.chem);
      in.chem_gpu.sync_concentrations_to_host(in.chem);
    }
  }

  if (!result.diffusion_on_gpu) {
    if (in.gpu_active && result.reactions_on_gpu) {
      in.chem_gpu.sync_concentrations_to_host(in.chem);
    }
    in.chem.apply_diffusion(in.domain, dt);
    in.chem.apply_boundaries(in.domain);
    if (in.gpu_active) {
      in.chem_gpu.sync_concentrations_to_device(in.chem);
    }
  }

  if (in.gpu_active) {
    gpu_sync_compute();
  }

  return result;
}

}  // namespace gutibm
