/* -----------------------------------------------------------------------
   GutIBM – GPU FMM far-field parity & device-residency test (Issue #159)

   Verifies:
   1. GPU vs CPU FMM parity (<= 1e-10 relative L2 vs CPU FMM) and exact
      Green's function accuracy (<= 2e-2 relative L2 for orders 2, 3) for
      BOTH locals_ready() == true (Mode A) and locals_ready() == false
      (Mode B), across expansion orders 1, 2, 3, periodic cross-seam
      sources, and BOTH replicated and slab chemical storage layouts.
   2. Zero full-grid H2D or D2H concentration transfers during
      QSSASolver::solve_bacteriocin_field with use_fmm = true and
      defer_sync = true (device residency enabled), audited via GpuProfile.
   ----------------------------------------------------------------------- */

#include "advection.h"
#include "agent.h"
#include "chemical_field.h"
#include "chemical_field_gpu.h"
#include "device_memory.h"
#include "dispatch.h"
#include "domain.h"
#include "fmm.h"
#include "fmm_gpu.h"
#include "gpu_profile.h"
#include "gpu_test_support.h"
#include "greens_function.h"
#include "input_parser.h"
#include "plasmid.h"
#include "qssa_solver.h"
#include "random.h"
#include "species_names.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace gutibm;

#ifdef GUTIBM_CUDA
namespace {

struct CrossSeamFixture {
  Domain domain;
  AdvectionField adv;
  GreensFunction gf;
  std::vector<Vec3> positions;
  std::vector<GreensFunctionParams> params;
  std::vector<Real> strengths;
  GreensFunctionParams avg_params;
};

CrossSeamFixture make_cross_seam_fixture() {
  CrossSeamFixture fx;
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {200.0e-6, 200.0e-6, 100.0e-6};
  dcfg.periodic = {true, true, false};
  dcfg.grid_dx = 10.0e-6;  // 20 x 20 x 10 = 4,000 voxels
  fx.domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 1e20;
  acfg.distal_transit_time = 1e20;
  acfg.mucus_thickness = 100.0e-6;
  acfg.distal_length = 200.0e-6;
  fx.adv.init(acfg, fx.domain);
  fx.gf.init(fx.domain, fx.adv);

  RNG rng(15977);
  auto add_cluster = [&](int count, Real x0, Real x1, Real y0, Real y1,
                         Real z0, Real z1) {
    for (int i = 0; i < count; ++i) {
      fx.positions.push_back({rng.uniform(x0, x1), rng.uniform(y0, y1),
                              rng.uniform(z0, z1)});
      GreensFunctionParams p;
      p.diff_coeff = 4.0e-11;
      p.source_rate = 1.0e-18;
      p.pI = 7.0;
      p.retardation = 5.0;
      p.decay_rate = 1.0e-3;
      fx.params.push_back(p);
      fx.strengths.push_back(p.source_rate);
    }
  };

  // Clusters near x = lo[0], x = hi[0], y = lo[1], and interior.
  add_cluster(18, 2.0e-6, 9.0e-6, 85.0e-6, 115.0e-6, 35.0e-6, 65.0e-6);
  add_cluster(12, 2.0e-6, 9.0e-6, 2.0e-6, 9.0e-6, 40.0e-6, 60.0e-6);
  add_cluster(18, 92.0e-6, 108.0e-6, 92.0e-6, 108.0e-6, 35.0e-6, 65.0e-6);
  add_cluster(12, 191.0e-6, 198.0e-6, 140.0e-6, 160.0e-6, 40.0e-6, 60.0e-6);

  fx.avg_params = fx.params[0];
  fx.avg_params.source_rate = 0.0;
  return fx;
}

Real relative_l2(const std::vector<Real>& approx,
                 const std::vector<Real>& ref) {
  assert(approx.size() == ref.size());
  Real num = 0.0;
  Real den = 0.0;
  for (size_t i = 0; i < ref.size(); ++i) {
    const Real d = approx[i] - ref[i];
    num += d * d;
    den += ref[i] * ref[i];
  }
  return den > 0.0 ? std::sqrt(num / den) : std::sqrt(num);
}

void compute_exact_and_near_grid(const CrossSeamFixture& fx,
                                 Real near_cutoff,
                                 const FmmDeviceLayout& layout,
                                 std::vector<Real>& exact_owned,
                                 std::vector<Real>& near_storage,
                                 Real halo_sentinel) {
  const Int ny = fx.domain.ny();
  const Int nz = fx.domain.nz();
  const Int owned_nx = layout.owned_storage_x_end - layout.owned_storage_x_begin;
  const Int storage_cells = layout.storage_nx * ny * nz;
  near_storage.assign(static_cast<size_t>(storage_cells), halo_sentinel);
  exact_owned.assign(static_cast<size_t>(owned_nx * ny * nz), 0.0);

  const Real cutoff2 = near_cutoff * near_cutoff;
  size_t owned_idx = 0;
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int local = 0; local < owned_nx; ++local) {
        const Int ix_storage = layout.owned_storage_x_begin + local;
        const Int ix_global = layout.local_global_x_begin + local;
        const Vec3 tgt = fx.domain.cell_center(ix_global, iy, iz);
        Real exact_sum = 0.0;
        Real near_sum = 0.0;
        for (size_t s = 0; s < fx.positions.size(); ++s) {
          const Real c = fx.gf.concentration_bounded(
              fx.positions[s], tgt, fx.params[s]);
          exact_sum += c;
          if (fx.domain.min_image_dist_sq(fx.positions[s], tgt) <= cutoff2) {
            near_sum += c;
          }
        }
        const Int storage_cell = iz * layout.storage_nx * ny
                               + iy * layout.storage_nx + ix_storage;
        near_storage[static_cast<size_t>(storage_cell)] = near_sum;
        exact_owned[owned_idx++] = exact_sum;
      }
    }
  }
}

void verify_parity_for_layout(const CrossSeamFixture& fx,
                              const std::string& layout_name,
                              const FmmDeviceLayout& layout,
                              bool check_halo_sentinel) {
  constexpr Real kCutoff = 35.0e-6;
  constexpr Real kTheta = 0.5;
  constexpr Real kHaloSentinel = -123.5;

  std::vector<Real> exact_owned;
  std::vector<Real> near_storage;
  compute_exact_and_near_grid(
      fx, kCutoff, layout, exact_owned, near_storage, kHaloSentinel);

  const Int ny = fx.domain.ny();
  const Int nz = fx.domain.nz();
  const Int owned_nx = layout.owned_storage_x_end - layout.owned_storage_x_begin;

  for (int order : {1, 2, 3}) {
    for (bool use_locals : {false, true}) {
      FMM fmm;
      fmm.build(fx.positions, fx.strengths, fx.domain, order);
      if (use_locals) {
        fmm.compute_local_expansions(kTheta, fx.gf, fx.avg_params, kCutoff);
        assert(fmm.locals_ready());
      } else {
        assert(!fmm.locals_ready());
      }

      // CPU reference: near + FMM far-field over owned cells.
      std::vector<Real> cpu_owned(static_cast<size_t>(owned_nx * ny * nz), 0.0);
      size_t owned_idx = 0;
      for (Int iz = 0; iz < nz; ++iz) {
        for (Int iy = 0; iy < ny; ++iy) {
          for (Int local = 0; local < owned_nx; ++local) {
            const Int ix_storage = layout.owned_storage_x_begin + local;
            const Int ix_global = layout.local_global_x_begin + local;
            const Vec3 tgt = fx.domain.cell_center(ix_global, iy, iz);
            const Int storage_cell = iz * layout.storage_nx * ny
                                   + iy * layout.storage_nx + ix_storage;
            const Real far = fmm.evaluate_far_field(
                tgt, kTheta, kCutoff, fx.gf, fx.avg_params);
            cpu_owned[owned_idx++] =
                near_storage[static_cast<size_t>(storage_cell)] + far;
          }
        }
      }

      // GPU evaluation directly into device buffer holding near_storage.
      DeviceBuffer<double> d_conc;
      d_conc.upload(near_storage);
      const bool ok = gpu_accumulate_far_field_device(
          fmm, fx.domain, fx.adv, fx.avg_params, kTheta, kCutoff,
          d_conc.data(), layout);
      assert(ok);

      std::vector<double> gpu_storage;
      d_conc.download(gpu_storage);
      assert(gpu_storage.size() == near_storage.size());

      std::vector<Real> gpu_owned(static_cast<size_t>(owned_nx * ny * nz), 0.0);
      owned_idx = 0;
      for (Int iz = 0; iz < nz; ++iz) {
        for (Int iy = 0; iy < ny; ++iy) {
          for (Int ix_s = 0; ix_s < layout.storage_nx; ++ix_s) {
            const Int storage_cell = iz * layout.storage_nx * ny
                                   + iy * layout.storage_nx + ix_s;
            const double val = gpu_storage[static_cast<size_t>(storage_cell)];
            assert(std::isfinite(val));
            if (ix_s < layout.owned_storage_x_begin
                || ix_s >= layout.owned_storage_x_end) {
              if (check_halo_sentinel) {
                assert(std::abs(val - kHaloSentinel) <= 1.0e-12);
              }
            } else {
              assert(val >= 0.0);
              gpu_owned[owned_idx++] = val;
            }
          }
        }
      }

      const Real rel_vs_cpu = relative_l2(gpu_owned, cpu_owned);
      const Real rel_vs_exact = relative_l2(gpu_owned, exact_owned);
      std::cout << "    layout=" << layout_name << " order=" << order
                << " locals_ready=" << (use_locals ? 1 : 0)
                << " rel_l2_vs_cpu=" << rel_vs_cpu
                << " rel_l2_vs_exact=" << rel_vs_exact << std::endl;
      assert(rel_vs_cpu <= 1.0e-10);
      if (order >= 2) {
        assert(rel_vs_exact <= 2.0e-2);
      }
    }
  }
}

void test_gpu_fmm_replicated_and_slab_parity() {
  const CrossSeamFixture fx = make_cross_seam_fixture();

  FmmDeviceLayout replicated{};
  replicated.storage_nx = fx.domain.nx();
  replicated.owned_storage_x_begin = 0;
  replicated.owned_storage_x_end = fx.domain.nx();
  replicated.local_global_x_begin = 0;
  verify_parity_for_layout(fx, "replicated", replicated, false);

  // Slab layout with 2-cell halos on both sides covering an interior+seam slice.
  FmmDeviceLayout slab{};
  slab.owned_storage_x_begin = 2;
  slab.owned_storage_x_end = 12;
  slab.storage_nx = 14;
  slab.local_global_x_begin = 10;  // global x in [10, 20), touching x = hi[0] seam
  verify_parity_for_layout(fx, "slab", slab, true);

  std::cout << "  test_gpu_fmm_replicated_and_slab_parity: PASSED" << std::endl;
}

void test_gpu_fmm_zero_full_grid_transfers() {
  DomainConfig dcfg;
  dcfg.lo = {0.0, 0.0, 0.0};
  dcfg.hi = {200.0e-6, 200.0e-6, 100.0e-6};
  dcfg.periodic = {true, true, false};
  dcfg.grid_dx = 5.0e-6;  // 40 x 40 x 20 = 32,000 cells (256 KB per species)
  Domain domain;
  domain.init(dcfg);

  AdvectionConfig acfg;
  acfg.radial_turnover = 1e20;
  acfg.distal_transit_time = 1e20;
  acfg.mucus_thickness = 100.0e-6;
  acfg.distal_length = 200.0e-6;
  AdvectionField adv;
  adv.init(acfg, domain);

  QSSAConfig qcfg;
  qcfg.use_fmm = true;
  qcfg.fmm_theta = 0.5;
  qcfg.fmm_expansion_order = 2;
  qcfg.toxin_cutoff = 40.0e-6;
  qcfg.colicin_release_rate = 1.0e-18;

  QSSASolver qssa_cpu;
  qssa_cpu.init(qcfg, domain, adv);
  QSSASolver qssa_gpu;
  qssa_gpu.init(qcfg, domain, adv);

  const SimulationConfig default_cfg = InputParser::default_config();
  ChemicalField chem_cpu;
  chem_cpu.init(domain, default_cfg.chemicals);
  ChemicalField chem_gpu_host;
  chem_gpu_host.init(domain, default_cfg.chemicals);

  const Int btub_idx = chem_cpu.find(species::BACTERIOCIN_BTUB);
  assert(btub_idx >= 0);

  ChemicalFieldGpu chem_gpu;
  chem_gpu.init(chem_gpu_host);
  assert(chem_gpu.active());
  assert(ChemicalFieldGpu::device_conc_residency_enabled());

  // Populate ColE1-carrying agents and active ColE1 lysis bursts near x=lo,
  // interior, and x=hi to exercise periodic cross-seam FMM on GPU.
  AgentPool agents;
  std::vector<ToxinBurstSource> bursts;
  const PlasmidEntry* cole1 = PlasmidLibrary::find("ColE1");
  assert(cole1 != nullptr);
  RNG rng(15999);
  for (int i = 0; i < 48; ++i) {
    Agent a{};
    a.state = PhenoState::NORMAL;
    const Real x = (i < 16) ? rng.uniform(2.0e-6, 9.0e-6)
                 : (i < 32) ? rng.uniform(90.0e-6, 110.0e-6)
                            : rng.uniform(191.0e-6, 198.0e-6);
    a.x = {x, rng.uniform(10.0e-6, 190.0e-6), rng.uniform(20.0e-6, 80.0e-6)};
    a.genome.bi_loci.push_back(cole1->cluster);
    assert(!a.genome.bi_loci.empty());
    agents.push_back(a);

    ToxinBurstSource burst{};
    burst.pos = a.x;
    burst.params.diff_coeff = cole1->cluster.diff_coeff;
    burst.params.source_rate = qcfg.colicin_release_rate;
    burst.params.pI = cole1->cluster.pI;
    burst.params.retardation = cole1->cluster.retardation;
    burst.params.decay_rate = 1.0e-3;
    burst.creation_time = 0.0;
    burst.release_tau = 300.0;
    burst.is_nuclease = false;
    burst.target = ReceptorType::BtuB;
    bursts.push_back(burst);
  }

  const ProteaseConfig protease{};

  // CPU reference solve.
  qssa_cpu.solve_bacteriocin_field(
      agents, bursts, 0.0, protease, adv, chem_cpu, btub_idx,
      ReceptorType::BtuB, nullptr, true);

  // GPU device-resident solve with transfer profiling.
  gpu_transfer_profile_set_enabled(true);
  gpu_transfer_profile_reset();

  qssa_gpu.solve_bacteriocin_field(
      agents, bursts, 0.0, protease, adv, chem_gpu_host, btub_idx,
      ReceptorType::BtuB, &chem_gpu, true);

  const GpuProfile prof = gpu_transfer_profile_snapshot();
  const auto sites = gpu_transfer_site_profiles();
  gpu_transfer_profile_set_enabled(false);

  const unsigned long long full_grid_bytes =
      static_cast<unsigned long long>(domain.ncells()) * sizeof(double);
  unsigned long long fmm_tree_h2d = 0;
  unsigned long long fmm_tree_d2h = 0;
  unsigned long long host_sync_bytes = 0;
  unsigned long long fallback_bytes = 0;
  for (const auto& s : sites) {
    if (s.label == "qssa_fmm_tree") {
      fmm_tree_h2d += s.h2d_bytes;
      fmm_tree_d2h += s.d2h_bytes;
    } else if (s.label == "qssa_toxin_host_sync") {
      host_sync_bytes += s.h2d_bytes + s.d2h_bytes;
    } else if (s.label == "qssa_fmm_fallback") {
      fallback_bytes += s.h2d_bytes + s.d2h_bytes;
    }
  }

  std::cout << "    solve_bacteriocin_field (use_fmm=true, defer_sync=true): "
            << "full_grid_bytes=" << full_grid_bytes
            << " total_h2d=" << prof.h2d_bytes
            << " total_d2h=" << prof.d2h_bytes
            << " fmm_tree_h2d=" << fmm_tree_h2d
            << " fmm_tree_d2h=" << fmm_tree_d2h << std::endl;

  // Zero full-grid concentration H2D or D2H transfers during the solve.
  assert(fmm_tree_h2d > 0 && fmm_tree_h2d < full_grid_bytes);
  assert(fmm_tree_d2h == 0);
  assert(host_sync_bytes == 0);
  assert(fallback_bytes == 0);
  assert(prof.d2h_bytes < full_grid_bytes);
  assert(prof.d2h_bytes <= 128ULL);
  assert(prof.h2d_bytes < full_grid_bytes);

  // Now materialize GPU result lazily via host read and verify CPU/GPU parity.
  const Int ncells = domain.ncells();
  std::vector<Real> cpu_field(static_cast<size_t>(ncells), 0.0);
  std::vector<Real> gpu_field(static_cast<size_t>(ncells), 0.0);
  for (Int c = 0; c < ncells; ++c) {
    cpu_field[static_cast<size_t>(c)] = chem_cpu.conc(btub_idx, c);
    gpu_field[static_cast<size_t>(c)] = chem_gpu_host.conc(btub_idx, c);
  }
  const Real rel_l2 = relative_l2(gpu_field, cpu_field);
  std::cout << "    solve_bacteriocin_field CPU vs GPU rel_l2=" << rel_l2
            << std::endl;
  assert(rel_l2 <= 1.0e-10);

  std::cout << "  test_gpu_fmm_zero_full_grid_transfers: PASSED" << std::endl;
}

}  // namespace
#endif

int main() {
  std::cout << "=== GPU FMM Far-Field Parity & Residency Tests ===\n";
  if (const int gpu_status = test::require_gpu("gpu_fmm_parity");
      gpu_status != 0) {
    return gpu_status;
  }

#ifndef GUTIBM_CUDA
  std::cout << "  SKIPPED (CUDA not compiled in)\n";
  return 0;
#else
  GpuConfig gcfg;
  gcfg.enabled = true;
  gcfg.device_id = 0;
  gpu_set_config(gcfg);
  if (!gpu_init_for_rank(0, 1)) return 1;

  test_gpu_fmm_replicated_and_slab_parity();
  test_gpu_fmm_zero_full_grid_transfers();

  std::cout << "All GPU FMM parity tests passed.\n";
  return 0;
#endif
}
