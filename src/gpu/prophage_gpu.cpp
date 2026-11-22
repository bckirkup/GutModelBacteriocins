/* -----------------------------------------------------------------------
   GutIBM – GPU dispatch implementation for Spec 14 temperate prophage
   ----------------------------------------------------------------------- */

#include "prophage_gpu.h"
#include "advection.h"
#include "agent.h"
#include "agent_pool_gpu.h"
#include "chemical_field.h"
#include "chemical_field_gpu.h"
#include "device_memory.h"
#include "dispatch.h"
#include "domain.h"
#include "fix_prophage.h"
#include "gpu_kernels.h"

#include <vector>

namespace gutibm {

bool gpu_compute_prophage_induction_probs(
    const AgentPoolGpu& agents,
    const std::vector<double>& damage_signals,
    const ProphageConfig& cfg,
    double dt,
    std::vector<double>& induction_probs_out) {
#ifndef GUTIBM_CUDA
  (void)agents;
  (void)damage_signals;
  (void)cfg;
  (void)dt;
  (void)induction_probs_out;
  return false;
#else
  const Int n = agents.size();
  if (!gpu_runtime_enabled() || n <= 0 || dt <= 0.0
      || static_cast<Int>(damage_signals.size()) != n) {
    return false;
  }

  GpuTransferSite site("prophage_induction");
  DeviceBuffer<double> d_damage;
  DeviceBuffer<double> d_probs;
  d_damage.upload(damage_signals);
  d_probs.allocate(static_cast<size_t>(n));

  gpu::launch_prophage_induction_prob_kernel(
      agents.state(), agents.is_ghost(), agents.has_inducible_prophage(),
      agents.prophage_induced(), d_damage.data(), d_probs.data(),
      static_cast<int>(n), cfg.k_spont, cfg.k_max_damage, cfg.K_damage,
      cfg.hill_n, dt, gpu_compute_stream());
  gpu_sync_compute();
  gpu_check_error("prophage_induction_prob_kernel");

  induction_probs_out.resize(static_cast<size_t>(n));
  d_probs.download(induction_probs_out);
  return true;
#endif
}

bool gpu_compute_prophage_adsorption(
    const AgentPoolGpu& agents,
    const AgentPool& pool,
    const ChemicalFieldGpu& chem_gpu,
    const ChemicalField& chem,
    const Domain& domain,
    const AdvectionField& advection,
    const ProphageConfig& cfg,
    Int virion_species_idx,
    double dt,
    std::vector<double>& infect_probs_out,
    std::vector<double>& reac_delta_out,
    double& total_decayed_out,
    double& total_adsorbed_out) {
#ifndef GUTIBM_CUDA
  (void)agents;
  (void)pool;
  (void)chem_gpu;
  (void)chem;
  (void)domain;
  (void)advection;
  (void)cfg;
  (void)virion_species_idx;
  (void)dt;
  (void)infect_probs_out;
  (void)reac_delta_out;
  (void)total_decayed_out;
  (void)total_adsorbed_out;
  return false;
#else
  const Int n = pool.size();
  const Int ncells = chem.ncells();
  const double cell_vol = domain.cell_volume();
  if (!gpu_runtime_enabled() || !chem_gpu.active()
      || virion_species_idx < 0 || ncells <= 0 || cell_vol <= 0.0
      || dt <= 0.0 || agents.size() != n) {
    return false;
  }

  GpuTransferSite site("prophage_adsorption");
  const Int nz = chem.global_nz();
  std::vector<double> washout_by_z(static_cast<size_t>(nz), 0.0);
  if (cfg.virion_washout_scale > 0.0) {
    for (Int iz = 0; iz < nz; ++iz) {
      const double z_center = domain.lo()[2]
          + (static_cast<double>(iz) + 0.5) * domain.dx_z();
      washout_by_z[static_cast<size_t>(iz)] =
          cfg.virion_washout_scale * advection.washout_rate(z_center);
    }
  }

  const auto cell_count = static_cast<size_t>(ncells);
  DeviceBuffer<double> d_cell_adsorption_rate(cell_count);
  DeviceBuffer<double> d_reac_delta(cell_count);
  DeviceBuffer<double> d_cell_adsorb_factor(cell_count);
  DeviceBuffer<double> d_washout_by_z;
  DeviceBuffer<double> d_totals(2);

  d_cell_adsorption_rate.upload(std::vector(cell_count, 0.0));
  d_reac_delta.upload(std::vector(cell_count, 0.0));
  d_cell_adsorb_factor.upload(std::vector(cell_count, 0.0));
  d_washout_by_z.upload(washout_by_z);
  d_totals.upload(std::vector(2, 0.0));

  const Int storage_nx = chem.slab_mode() ? chem.storage_nx() : domain.nx();
  const Int owned_global_x_begin = chem.slab_mode()
      ? chem.owned_global_x_begin() : domain.local_grid_x_begin();
  const Int owned_global_x_end = chem.slab_mode()
      ? chem.owned_global_x_end() : domain.local_grid_x_end();
  const Int owned_storage_x_begin = chem.slab_mode()
      ? chem.owned_storage_x_begin() : domain.local_grid_x_begin();
  const Int owned_storage_x_end = chem.slab_mode()
      ? chem.owned_storage_x_end() : domain.local_grid_x_end();
  const Int receptor_idx = prophage_adsorption_receptor_index(cfg);

  if (n > 0 && cfg.k_adsorb > 0.0) {
    gpu::launch_prophage_adsorption_accum_kernel(
        agents.grid_cell(), agents.state(), agents.is_ghost(),
        agents.receptor_expr(), d_cell_adsorption_rate.data(),
        static_cast<int>(n), static_cast<int>(receptor_idx), cfg.k_adsorb,
        cell_vol, domain.nx(), domain.ny(), storage_nx,
        owned_global_x_begin, owned_global_x_end, owned_storage_x_begin,
        gpu_compute_stream());
  }

  gpu::launch_prophage_cell_removal_kernel(
      chem_gpu.conc_device(virion_species_idx), d_cell_adsorption_rate.data(),
      d_washout_by_z.data(), d_reac_delta.data(), d_cell_adsorb_factor.data(),
      d_totals.data(), storage_nx, chem.global_ny(), nz,
      owned_storage_x_begin, owned_storage_x_end, cfg.k_virion_decay,
      cell_vol, dt, gpu_compute_stream());

  infect_probs_out.assign(static_cast<size_t>(n), 0.0);
  if (n > 0) {
    DeviceBuffer<double> d_infect_probs(static_cast<size_t>(n));
    gpu::launch_prophage_infect_prob_kernel(
        agents.grid_cell(), agents.state(), agents.is_ghost(),
        agents.has_prophage_immunity(), agents.receptor_expr(),
        d_cell_adsorb_factor.data(), d_infect_probs.data(),
        static_cast<int>(n), static_cast<int>(receptor_idx), cfg.k_adsorb,
        cell_vol, domain.nx(), domain.ny(), storage_nx,
        owned_global_x_begin, owned_global_x_end, owned_storage_x_begin,
        gpu_compute_stream());
    gpu_sync_compute();
    gpu_check_error("prophage_adsorption_kernels");
    d_infect_probs.download(infect_probs_out);
  } else {
    gpu_sync_compute();
    gpu_check_error("prophage_cell_removal_kernel");
  }

  reac_delta_out.resize(cell_count);
  d_reac_delta.download(reac_delta_out);

  std::vector<double> totals_host(2, 0.0);
  d_totals.download(totals_host);
  total_decayed_out = totals_host[0];
  total_adsorbed_out = totals_host[1];
  return true;
#endif
}

}  // namespace gutibm
