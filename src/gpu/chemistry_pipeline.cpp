#include "chemistry_pipeline.h"
#include "agent.h"
#include "chemical_field.h"
#include "chemical_field_gpu.h"
#include "domain.h"
#include "dispatch.h"
#include "step_profiler.h"
#include "qssa_gpu.h"
#include "qssa_solver.h"
#include "species_names.h"
#include "vbf.h"
#include "vbf_gpu.h"

#include <chrono>
#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {

namespace {

bool sum_reactions_with_optional_device(ChemistryPipelineInput& in,
                                        bool reactions_on_device) {
  const auto t0 = std::chrono::steady_clock::now();
  if (reactions_on_device
      && in.chem_gpu.try_sum_reactions_on_device(in.chem)) {
    if (in.step_profile != nullptr) {
      const auto t1 = std::chrono::steady_clock::now();
      in.step_profile->mpi_reaction_reduce_s +=
          std::chrono::duration<double>(t1 - t0).count();
    }
    return true;
  }
  if (reactions_on_device) {
    in.chem_gpu.sync_reactions_to_host(in.chem);
  }
  in.chem.sum_reactions_across_ranks();
  if (in.step_profile != nullptr) {
    const auto t1 = std::chrono::steady_clock::now();
    in.step_profile->mpi_reaction_reduce_s +=
        std::chrono::duration<double>(t1 - t0).count();
  }
  return false;
}

}  // namespace

ChemistryPipelineResult run_chemistry_pipeline(ChemistryPipelineInput& in, Real dt) {
  ChemistryPipelineResult result;
  bool reactions_on_device = in.gpu_active && in.metabolism_on_gpu;
  if (in.gpu_active && !reactions_on_device) {
    in.chem_gpu.sync_reactions_to_device(in.chem);
    reactions_on_device = true;
  }

  bool oxygen_on_gpu = false;
  if (in.gpu_active) {
    oxygen_on_gpu = gpu_solve_nutrient_depletion(
        in.agents_gpu, in.num_agents, in.chem_gpu, in.chem,
        in.oxygen, in.domain);
  }

  if (!oxygen_on_gpu) {
    if (reactions_on_device) {
      in.chem_gpu.sync_reactions_to_host(in.chem);
      reactions_on_device = false;
    }
    in.qssa.solve_nutrient_depletion(in.agents, in.chem, in.oxygen);
  }

  // Sum rank-local agent reaction fields before adding the VBF contribution.
  const bool reactions_reduced_on_device =
      sum_reactions_with_optional_device(in, reactions_on_device);
  reactions_on_device = reactions_reduced_on_device;
  const Int carbon = in.chem.find(species::CARBON);
  const Int iron = in.chem.find(species::IRON);
  const Int oxygen = in.chem.find(species::OXYGEN);
  VbfFluxTotals vbf_totals;
  in.chem.sum_agent_uptake_across_ranks();
  in.chem.flux_accounting().commit_agent_uptake_step();
  bool applied_vbf_on_gpu = false;
  if (in.gpu_active && !reactions_on_device) {
    in.chem_gpu.sync_reactions_to_device(in.chem);
  }
  if (in.gpu_active) {
    in.chem_gpu.reset_vbf_totals();
    applied_vbf_on_gpu = gpu_apply_vbf_coupling(
        in.chem_gpu, in.chem, in.domain, in.vbf,
        in.oxygen, in.acetate, in.mucin, vbf_totals, dt);
  }
  if (!applied_vbf_on_gpu) {
    if (reactions_on_device) {
      in.chem_gpu.sync_reactions_to_host(in.chem);
    }
    in.vbf.apply_nutrient_coupling(in.chem, in.domain, dt,
                                   in.oxygen, in.acetate, in.mucin,
                                   &vbf_totals);
  }
  if (carbon >= 0) {
#ifdef GUTIBM_MPI
    if (in.chem.slab_mode() && in.domain.nprocs() > 1) {
      MPI_Allreduce(MPI_IN_PLACE, &vbf_totals.carbon_source, 1, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &vbf_totals.carbon_sink, 1, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    in.flux_accounting.add_interval(
        carbon, 0.0, vbf_totals.carbon_source,
        vbf_totals.carbon_sink, 0.0);
  }
  if (iron >= 0) {
#ifdef GUTIBM_MPI
    if (in.chem.slab_mode() && in.domain.nprocs() > 1) {
      MPI_Allreduce(MPI_IN_PLACE, &vbf_totals.iron_sink, 1, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    in.flux_accounting.add_interval(
        iron, 0.0, 0.0, vbf_totals.iron_sink, 0.0);
  }
  if (oxygen >= 0) {
#ifdef GUTIBM_MPI
    if (in.chem.slab_mode() && in.domain.nprocs() > 1) {
      MPI_Allreduce(MPI_IN_PLACE, &vbf_totals.oxygen_sink, 1, MPI_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    in.flux_accounting.add_interval(
        oxygen, 0.0, 0.0, vbf_totals.oxygen_sink, 0.0);
  }
  if (in.gpu_active && applied_vbf_on_gpu) {
    result.reactions_on_gpu = in.chem_gpu.apply_reactions(dt, in.domain);
    if (result.reactions_on_gpu) {
      in.chem_gpu.download_reaction_clip(in.chem);
    }
  }

  if (!result.reactions_on_gpu) {
    const Real cell_volume = in.domain.cell_volume();
    Int s = 0;
    for (const auto& conc_row : in.chem.conc_data()) {
      (void)conc_row;
      #ifdef GUTIBM_OPENMP
      #pragma omp parallel for collapse(3) schedule(static)
      #endif
      for (Int iz = 0; iz < in.chem.global_nz(); ++iz) {
        for (Int iy = 0; iy < in.chem.global_ny(); ++iy) {
          for (Int ix = in.chem.owned_storage_x_begin();
               ix < in.chem.owned_storage_x_end(); ++ix) {
            const Int c = iz * in.chem.storage_nx() * in.chem.global_ny()
                + iy * in.chem.storage_nx() + ix;
            const Real updated = in.chem.conc(s, c)
                + in.chem.reac(s, c) * dt;
            if (updated < 0.0) {
              in.flux_accounting.add_reaction_clip(
                  s, -updated * cell_volume);
            }
            in.chem.conc(s, c) = std::max(updated, 0.0);
          }
        }
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

  in.chem.sum_accounting_across_ranks();
  in.chem.flux_accounting().commit_boundary_and_reaction_step();

  return result;
}

}  // namespace gutibm
