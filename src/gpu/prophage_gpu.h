/* -----------------------------------------------------------------------
   GutIBM – GPU dispatch interface for Spec 14 temperate prophage kernels
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_PROPHAGE_GPU_H
#define GUTIBM_PROPHAGE_GPU_H

#include "types.h"
#include <vector>

namespace gutibm {

class AdvectionField;
class AgentPool;
class AgentPoolGpu;
class ChemicalField;
class ChemicalFieldGpu;
class Domain;
struct ProphageConfig;

// Compute per-agent Hill induction probabilities on GPU.
// Returns false when CUDA is unavailable or inactive (caller falls back to CPU).
bool gpu_compute_prophage_induction_probs(
    const AgentPoolGpu& agents,
    const std::vector<double>& damage_signals,
    const ProphageConfig& cfg,
    double dt,
    std::vector<double>& induction_probs_out);

// Compute per-cell virion decay/washout/adsorption split and per-agent
// susceptible infection probabilities on GPU.
// Returns false when CUDA is unavailable or inactive (caller falls back to CPU).
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
    double& total_adsorbed_out);

}  // namespace gutibm

#endif  // GUTIBM_PROPHAGE_GPU_H
