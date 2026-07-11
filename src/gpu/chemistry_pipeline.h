#ifndef GUTIBM_CHEMISTRY_PIPELINE_H
#define GUTIBM_CHEMISTRY_PIPELINE_H

#include "types.h"
#include "agent.h"
#include "chem_environment_config.h"

namespace gutibm {

class AgentPool;

class AgentPoolGpu;
class ChemicalFieldGpu;
class ChemicalField;
class Domain;
class VBF;
class QSSASolver;

struct ChemistryPipelineInput {
  bool gpu_active = false;
  AgentPoolGpu& agents_gpu;
  ChemicalFieldGpu& chem_gpu;
  ChemicalField& chem;
  const Domain& domain;
  const VBF& vbf;
  QSSASolver& qssa;
  const AgentPool& agents;
  const OxygenConfig& oxygen;
  const AcetateConfig& acetate;
  const MucinConfig& mucin;
  Int num_agents = 0;
};

struct ChemistryPipelineResult {
  bool reactions_on_gpu = false;
  bool diffusion_on_gpu = false;
};

// Runs nutrient depletion, MPI reaction sum, VBF coupling, reaction integration,
// and implicit diffusion with GPU-first dispatch and CPU fallbacks.
ChemistryPipelineResult run_chemistry_pipeline(ChemistryPipelineInput& in, Real dt);

}  // namespace gutibm

#endif  // GUTIBM_CHEMISTRY_PIPELINE_H
