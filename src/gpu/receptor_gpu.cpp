#include "receptor_gpu.h"
#include "agent_pool_gpu.h"
#include "agent.h"
#include "chemical_field_gpu.h"
#include "chemical_field.h"
#include "fix_receptor.h"
#include "species_names.h"
#include "dispatch.h"
#include "device_memory.h"
#include "gpu_kernels.h"

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

#include <vector>

namespace gutibm {

namespace {

constexpr int kImmunityReceptors = 4;

void pack_immunity_eff(const AgentPool& pool, const ReceptorConfig& cfg,
                       std::vector<double>& out) {
  const Int n = pool.size();
  out.assign(static_cast<size_t>(n) * kImmunityReceptors, 1.0);
  for (Int i = 0; i < n; ++i) {
    const Agent& a = pool[i];
    auto set_eff = [&](ReceptorType target, int slot) {
      double eff = 1.0;
      for (const auto& bi : a.genome.bi_loci) {
        if (bi.target == target) {
          const double candidate = cfg.immunity_factor * bi.immunity_binding_affinity;
          if (candidate < eff) eff = candidate;
        }
      }
      out[static_cast<size_t>(i) * kImmunityReceptors + slot] = eff;
    };
    set_eff(ReceptorType::BtuB, 0);
    set_eff(ReceptorType::FepA, 1);
    set_eff(ReceptorType::CirA, 2);
    set_eff(ReceptorType::FhuA, 3);
  }
}

void pack_toxin_affinity(const AgentPool& pool, std::vector<double>& out) {
  const Int n = pool.size();
  out.assign(static_cast<size_t>(n) * NUM_RECEPTORS, 1.0);
  for (Int i = 0; i < n; ++i) {
    for (int r = 0; r < NUM_RECEPTORS; ++r) {
      out[static_cast<size_t>(i) * NUM_RECEPTORS + r] = pool[i].genome.toxin_affinity[r];
    }
  }
}

}  // namespace

bool gpu_compute_receptor_kill_probs_host_packed(
    const AgentPoolGpu& agents,
    const AgentPool& pool,
    const ChemicalFieldGpu& chem_gpu,
    const ChemicalField& chem,
    const ReceptorConfig& cfg,
    double dt,
    std::vector<double>& kill_probs_out) {

#ifndef GUTIBM_CUDA
  (void)agents;
  (void)pool;
  (void)chem_gpu;
  (void)chem;
  (void)cfg;
  (void)dt;
  (void)kill_probs_out;
  return false;
#else
  const Int n = pool.size();
  if (!gpu_runtime_enabled() || !chem_gpu.active() || n <= 0) return false;

  std::vector<double> immunity;
  std::vector<double> toxin_aff;
  pack_immunity_eff(pool, cfg, immunity);
  pack_toxin_affinity(pool, toxin_aff);

  DeviceBuffer<double> d_immunity;
  DeviceBuffer<double> d_toxin_aff;
  DeviceBuffer<double> d_kill;
  d_immunity.upload(immunity);
  d_toxin_aff.upload(toxin_aff);
  d_kill.allocate(static_cast<size_t>(n));

  Int i_btuB = chem.find(species::BACTERIOCIN_BTUB);
  Int i_fepA = chem.find(species::BACTERIOCIN_FEPA);
  Int i_cirA = chem.find(species::BACTERIOCIN_CIRA);
  Int i_fhuA = chem.find(species::BACTERIOCIN_FHUA);
  Int i_b12 = chem.find(species::B12);
  Int i_iron = chem.find(species::IRON);
  Int i_sid = chem.find(species::SIDEROPHORE);

  gpu::launch_receptor_kill_prob_kernel(
      agents.grid_cell(), agents.state(),
      agents.receptor_expr(), agents.ligand_affinity(),
      d_toxin_aff.data(), d_immunity.data(),
      i_btuB >= 0 ? chem_gpu.conc_device(i_btuB) : nullptr,
      i_fepA >= 0 ? chem_gpu.conc_device(i_fepA) : nullptr,
      i_cirA >= 0 ? chem_gpu.conc_device(i_cirA) : nullptr,
      i_fhuA >= 0 ? chem_gpu.conc_device(i_fhuA) : nullptr,
      i_b12 >= 0 ? chem_gpu.conc_device(i_b12) : nullptr,
      i_iron >= 0 ? chem_gpu.conc_device(i_iron) : nullptr,
      i_sid >= 0 ? chem_gpu.conc_device(i_sid) : nullptr,
      d_kill.data(), static_cast<int>(n), dt,
      cfg.kd_b12_btuB, cfg.kd_colicinE_btuB, cfg.kd_enterobactin,
      cfg.kd_colicinB_fepA, cfg.kd_lin_enterobactin, cfg.kd_colicinIa_cirA,
      cfg.kd_colicinM_fhuA, cfg.kd_ferrichrome, cfg.cirA_linearized_fraction,
      cfg.kill_rate_colicin, cfg.kill_rate_microcin, gpu_compute_stream());
  gpu_sync_compute();
  gpu_check_error("receptor_kill_prob_kernel");

  kill_probs_out.resize(static_cast<size_t>(n));
  d_kill.download(kill_probs_out);
  return true;
#endif
}

}  // namespace gutibm
