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

namespace gutibm {

namespace {

void sum_reactions_with_optional_device(ChemistryPipelineInput& in) {
  const auto t0 = std::chrono::steady_clock::now();
  if (in.gpu_active && in.chem_gpu.try_sum_reactions_on_device(in.chem)) {
    if (in.step_profile != nullptr) {
      const auto t1 = std::chrono::steady_clock::now();
      in.step_profile->mpi_reaction_reduce_s +=
          std::chrono::duration<double>(t1 - t0).count();
    }
    return;
  }
  in.chem.sum_reactions_across_ranks();
  if (in.step_profile != nullptr) {
    const auto t1 = std::chrono::steady_clock::now();
    in.step_profile->mpi_reaction_reduce_s +=
        std::chrono::duration<double>(t1 - t0).count();
  }
}

}  // namespace

ChemistryPipelineResult run_chemistry_pipeline(ChemistryPipelineInput& in, Real dt) {
  ChemistryPipelineResult result;
  std::vector<Real> agent_reaction_amount(in.chem.num_species(), 0.0);

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
  sum_reactions_with_optional_device(in);
  for (Int s = 0; s < in.chem.num_species(); ++s) {
    for (Int cell = 0; cell < in.chem.ncells(); ++cell) {
      agent_reaction_amount[static_cast<size_t>(s)] +=
          in.chem.reac(s, cell) * in.domain.dx() * in.domain.dx()
          * in.domain.dx() * dt;
    }
  }

  const Int carbon = in.chem.find(species::CARBON);
  const Int iron = in.chem.find(species::IRON);
  const Int oxygen = in.chem.find(species::OXYGEN);
  VbfFluxTotals vbf_totals;
  const auto record_uptake = [&](Int species_index, Real amount) {
    if (species_index < 0) return;
    const Real uptake = std::max(-amount, 0.0);
    in.flux_accounting.add_interval(species_index, 0.0, 0.0, 0.0, uptake);
  };
  if (carbon >= 0) {
    record_uptake(carbon, agent_reaction_amount[static_cast<size_t>(carbon)]);
  }
  if (iron >= 0) {
    record_uptake(iron, agent_reaction_amount[static_cast<size_t>(iron)]);
  }
  bool reactions_on_device = false;
  // Keep VBF coupling on the host so the accounting uses the same expression
  // as the reaction update. The resulting reaction field is uploaded before
  // the existing device concentration integration.
  in.vbf.apply_nutrient_coupling(in.chem, in.domain, dt,
                                 in.oxygen, in.acetate, in.mucin,
                                 &vbf_totals);
  if (carbon >= 0) {
    in.flux_accounting.add_interval(
        carbon, 0.0, vbf_totals.carbon_source,
        vbf_totals.carbon_sink, 0.0);
  }
  if (iron >= 0) {
    in.flux_accounting.add_interval(
        iron, 0.0, 0.0, vbf_totals.iron_sink, 0.0);
  }
  if (oxygen >= 0) {
    in.flux_accounting.add_interval(
        oxygen, 0.0, 0.0, vbf_totals.oxygen_sink, 0.0);
  }
  if (in.gpu_active) {
    in.chem_gpu.sync_reactions_to_device(in.chem);
    reactions_on_device = true;
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
