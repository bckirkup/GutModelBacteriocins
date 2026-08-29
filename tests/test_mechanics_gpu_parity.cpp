/* -----------------------------------------------------------------------
   GutIBM – GPU mechanics parity (issue #157)
   ----------------------------------------------------------------------- */

#include "fix_mechanics.h"
#include "simulation.h"
#include "input_parser.h"
#include "dispatch.h"
#include "device.h"
#include "mechanics_test_helpers.h"
#include "gpu_test_support.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace gutibm;
using gutibm::test::make_two_agent_sim;

namespace {

void test_gpu_overlapping_agents_match_cpu() {
#ifdef GUTIBM_CUDA
  GpuConfig gcfg;
  gcfg.enabled = true;
  gpu_set_config(gcfg);
  assert(gpu_init_for_rank(0, 1));

  Real r = CELL_RADIUS_DEFAULT;
  Real separation = 2 * r - 0.2e-6;
  Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b = {50e-6 + separation, 50e-6, 50e-6};

  MechanicsConfig mcfg;
  mcfg.hertzian_enabled = true;
  mcfg.hertz_k = 1.0e-6;
  mcfg.adhesion_enabled = false;

  Simulation sim_cpu = make_two_agent_sim(pos_a, pos_b, mcfg, false);
  FixMechanics fix_cpu(sim_cpu, mcfg);
  fix_cpu.compute(1.0);

  Simulation sim_gpu = make_two_agent_sim(pos_a, pos_b, mcfg, true);
  assert(sim_gpu.gpu_active());
  FixMechanics fix_gpu(sim_gpu, mcfg);
  fix_gpu.compute(1.0);

  for (Int i = 0; i < 2; ++i) {
    for (int d = 0; d < 3; ++d) {
      const Real diff = std::abs(sim_cpu.agents()[i].x[d] - sim_gpu.agents()[i].x[d]);
      assert(diff < 1e-12);
    }
  }

  std::cout << "  test_gpu_overlapping_agents_match_cpu: PASSED\n";
#endif
}

void test_gpu_adhesion_match_cpu() {
#ifdef GUTIBM_CUDA
  GpuConfig gcfg;
  gcfg.enabled = true;
  gpu_set_config(gcfg);
  assert(gpu_init_for_rank(0, 1));

  Real r = CELL_RADIUS_DEFAULT;
  Real gap = 0.2e-6;
  Real separation = 2 * r + gap;
  Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  Vec3 pos_b = {50e-6 + separation, 50e-6, 50e-6};

  MechanicsConfig mcfg;
  mcfg.hertzian_enabled = true;
  mcfg.adhesion_enabled = true;
  mcfg.adhesion_strength = 1.0e-12;
  mcfg.adhesion_range = 0.5e-6;

  Simulation sim_cpu = make_two_agent_sim(pos_a, pos_b, mcfg, false);
  FixMechanics fix_cpu(sim_cpu, mcfg);
  fix_cpu.compute(1.0);

  Simulation sim_gpu = make_two_agent_sim(pos_a, pos_b, mcfg, true);
  assert(sim_gpu.gpu_active());
  FixMechanics fix_gpu(sim_gpu, mcfg);
  fix_gpu.compute(1.0);

  const Real dist_cpu = sim_cpu.agents()[1].x[0] - sim_cpu.agents()[0].x[0];
  const Real dist_gpu = sim_gpu.agents()[1].x[0] - sim_gpu.agents()[0].x[0];
  assert(dist_cpu < separation);
  assert(dist_gpu < separation);
  assert(std::abs(dist_cpu - dist_gpu) < 1e-12);

  std::cout << "  test_gpu_adhesion_match_cpu: PASSED\n";
#endif
}

#ifdef GUTIBM_CUDA
void mark_corpse(Simulation& sim, Real death_time) {
  sim.agents()[1].state = PhenoState::DEAD;
  sim.agents()[1].timers.death_time = death_time;
}

CdiConfig corpse_cdi_config() {
  CdiConfig cdi;
  cdi.enabled = true;
  cdi.corpse_persistence = 300.0;
  return cdi;
}
#endif

void test_gpu_fresh_corpse_match_cpu() {
#ifdef GUTIBM_CUDA
  GpuConfig gcfg;
  gcfg.enabled = true;
  gpu_set_config(gcfg);
  assert(gpu_init_for_rank(0, 1));

  const Real r = CELL_RADIUS_DEFAULT;
  const Real separation = 2 * r - 0.2e-6;
  const Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  const Vec3 pos_b = {50e-6 + separation, 50e-6, 50e-6};

  MechanicsConfig mcfg;
  mcfg.hertzian_enabled = true;
  mcfg.hertz_k = 1.0e-6;
  mcfg.adhesion_enabled = false;
  const CdiConfig cdi = corpse_cdi_config();

  Simulation sim_cpu =
      make_two_agent_sim(pos_a, pos_b, mcfg, false, 0.01, cdi);
  mark_corpse(sim_cpu, 0.0);
  FixMechanics fix_cpu(sim_cpu, mcfg);
  fix_cpu.compute(1.0);

  Simulation sim_gpu =
      make_two_agent_sim(pos_a, pos_b, mcfg, true, 0.01, cdi);
  mark_corpse(sim_gpu, 0.0);
  assert(sim_gpu.gpu_active());
  FixMechanics fix_gpu(sim_gpu, mcfg);
  fix_gpu.compute(1.0);

  const Real host_displacement =
      std::abs(sim_cpu.agents()[0].x[0] - pos_a[0])
      + std::abs(sim_cpu.agents()[1].x[0] - pos_b[0]);
  const Real device_displacement =
      std::abs(sim_gpu.agents()[0].x[0] - pos_a[0])
      + std::abs(sim_gpu.agents()[1].x[0] - pos_b[0]);
  assert(host_displacement > 1.0e-12);
  assert(device_displacement > 1.0e-12);
  for (Int i = 0; i < 2; ++i) {
    for (int d = 0; d < 3; ++d) {
      assert(std::abs(sim_cpu.agents()[i].x[d]
                      - sim_gpu.agents()[i].x[d]) < 1e-12);
    }
  }

  std::cout << "  test_gpu_fresh_corpse_match_cpu: PASSED\n";
#endif
}

void test_gpu_expired_corpse_match_cpu() {
#ifdef GUTIBM_CUDA
  GpuConfig gcfg;
  gcfg.enabled = true;
  gpu_set_config(gcfg);
  assert(gpu_init_for_rank(0, 1));

  const Real r = CELL_RADIUS_DEFAULT;
  const Real separation = 2 * r - 0.2e-6;
  const Vec3 pos_a = {50e-6, 50e-6, 50e-6};
  const Vec3 pos_b = {50e-6 + separation, 50e-6, 50e-6};
  MechanicsConfig mcfg;
  mcfg.hertzian_enabled = true;
  mcfg.hertz_k = 1.0e-6;
  mcfg.adhesion_enabled = false;
  const CdiConfig cdi = corpse_cdi_config();

  Simulation sim_cpu =
      make_two_agent_sim(pos_a, pos_b, mcfg, false, 0.01, cdi);
  mark_corpse(sim_cpu, -(cdi.corpse_persistence + 1.0));
  FixMechanics fix_cpu(sim_cpu, mcfg);
  fix_cpu.compute(1.0);

  Simulation sim_gpu =
      make_two_agent_sim(pos_a, pos_b, mcfg, true, 0.01, cdi);
  mark_corpse(sim_gpu, -(cdi.corpse_persistence + 1.0));
  assert(sim_gpu.gpu_active());
  FixMechanics fix_gpu(sim_gpu, mcfg);
  fix_gpu.compute(1.0);

  for (Int i = 0; i < 2; ++i) {
    for (int d = 0; d < 3; ++d) {
      assert(std::abs(sim_cpu.agents()[i].x[d]
                      - (i == 0 ? pos_a[d] : pos_b[d])) < 1e-12);
      assert(std::abs(sim_cpu.agents()[i].x[d]
                      - sim_gpu.agents()[i].x[d]) < 1e-12);
    }
  }
  std::cout << "  test_gpu_expired_corpse_match_cpu: PASSED\n";
#endif
}

}  // namespace

int main() {
  std::cout << "=== Mechanics GPU Parity Tests ===\n";
  const int gpu_status = test::require_gpu("mechanics_gpu_parity");
  if (gpu_status != 0) return gpu_status;
  test_gpu_overlapping_agents_match_cpu();
  test_gpu_adhesion_match_cpu();
  test_gpu_fresh_corpse_match_cpu();
  test_gpu_expired_corpse_match_cpu();
  std::cout << "All mechanics GPU parity tests passed.\n";
  return 0;
}
