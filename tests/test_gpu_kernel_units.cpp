#include "gpu_kernels.h"
#include "gpu_test_support.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#ifdef GUTIBM_CUDA
#include "device_memory.h"
#include <cuda_runtime.h>
#endif

namespace {

#ifdef GUTIBM_CUDA

using gutibm::DeviceBuffer;
using gutibm::gpu::AdvectionParams;
using gutibm::gpu::DomainParams;
using gutibm::gpu::GfSourceParams;
using gutibm::gpu::MechanicsLaunchParams;
using gutibm::gpu::VbfLaunchParams;

constexpr int kNx = 4;
constexpr int kNy = 4;
constexpr int kNz = 4;
constexpr int kCells = kNx * kNy * kNz;

template <typename T>
std::vector<T> download(const DeviceBuffer<T>& device, int count) {
  std::vector<T> host(static_cast<size_t>(count));
  device.download(host.data(), static_cast<size_t>(count));
  return host;
}

void synchronize() {
  assert(cudaDeviceSynchronize() == cudaSuccess);
}

DomainParams domain_params() {
  DomainParams dom{};
  dom.nx = kNx;
  dom.ny = kNy;
  dom.nz = kNz;
  dom.dx_x = 1.0;
  dom.dx_y = 1.0;
  dom.dx_z = 1.0;
  dom.lo = {0.0, 0.0, 0.0};
  dom.extent = {4.0, 4.0, 4.0};
  dom.periodic = {true, true, false};
  dom.z_lo = 0.0;
  dom.z_hi = 4.0;
  return dom;
}

AdvectionParams zero_advection() {
  AdvectionParams adv{};
  adv.h = 4.0;
  adv.lo_z = 0.0;
  adv.profile_alpha = 1.0;
  return adv;
}

void test_field_and_boundaries() {
  DeviceBuffer<double> conc(kCells);
  DeviceBuffer<double> reac(kCells);
  DeviceBuffer<double> clip(1);
  std::vector<double> host_conc(kCells, 1.0);
  std::vector<double> host_reac(kCells, 0.0);
  host_reac[0] = -2.0;
  conc.upload(host_conc);
  reac.upload(host_reac);
  clip.upload(std::vector<double>{0.0});
  gutibm::gpu::launch_field_update_kernel(
      conc.data(), reac.data(), kCells, 1, 1.0, clip.data(), 1.0,
      kNx, kNy, kNz, 0, kNx, nullptr);
  synchronize();
  host_conc = download(conc, kCells);
  const auto clip_host = download(clip, 1);
  assert(host_conc[0] == 0.0);
  assert(std::abs(clip_host[0] - 1.0) < 1.0e-12);
  for (int cell = 1; cell < kCells; ++cell) assert(host_conc[cell] == 1.0);
  gutibm::gpu::launch_field_update_kernel(
      conc.data(), reac.data(), kCells, 1, 1.0, clip.data(), 1.0,
      kNx, kNy, kNz, 0, kNx, nullptr);
  synchronize();
  const auto doubled_clip = download(clip, 1);
  assert(doubled_clip[0] > clip_host[0]);

  DeviceBuffer<double> graded_conc(kCells);
  DeviceBuffer<double> graded_reac(kCells);
  graded_conc.upload(std::vector<double>(kCells, 1.0));
  graded_reac.upload(std::vector<double>(kCells, 0.5));
  const double dts[] = {0.5, 1.0, 2.0};
  double previous = 1.0;
  for (const double dt : dts) {
    graded_conc.upload(std::vector<double>(kCells, 1.0));
    gutibm::gpu::launch_field_update_kernel(
        graded_conc.data(), graded_reac.data(), kCells, 1, dt, nullptr, 1.0,
        kNx, kNy, kNz, 0, kNx, nullptr);
    synchronize();
    const auto graded = download(graded_conc, kCells);
    assert(graded[0] > previous);
    previous = graded[0];
  }
  graded_conc.upload(std::vector<double>(kCells, 1.0));
  gutibm::gpu::launch_field_update_kernel(
      graded_conc.data(), graded_reac.data(), kCells, 1, 0.0, nullptr, 1.0,
      kNx, kNy, kNz, 0, kNx, nullptr);
  synchronize();
  assert(download(graded_conc, kCells)[0] == 1.0);

  DeviceBuffer<double> boundary(1);
  boundary.upload(std::vector<double>{0.5});
  gutibm::gpu::launch_apply_boundaries_kernel(
      conc.data(), kNx, kNy, kNz, 1, boundary.data(), nullptr);
  synchronize();
  host_conc = download(conc, kCells);
  for (int iy = 0; iy < kNy; ++iy) {
    for (int ix = 0; ix < kNx; ++ix) {
      const int bottom = iy * kNx + ix;
      const int top = (kNz - 1) * kNx * kNy + iy * kNx + ix;
      assert(host_conc[bottom] == 0.5);
      assert(host_conc[top] == host_conc[(kNz - 2) * kNx * kNy + iy * kNx + ix]);
    }
  }
}

void test_grid_and_metabolism() {
  constexpr int agents = 2;
  DeviceBuffer<double> x(agents);
  DeviceBuffer<double> y(agents);
  DeviceBuffer<double> z(agents);
  DeviceBuffer<int> cells(agents);
  DeviceBuffer<int> state(agents);
  x.upload(std::vector<double>{0.1, 2.1});
  y.upload(std::vector<double>{1.1, 2.1});
  z.upload(std::vector<double>{0.1, 3.1});
  state.upload(std::vector<int>{0, 3});
  cells.upload(std::vector<int>{-1, -1});
  gutibm::gpu::launch_grid_coupling_kernel(
      x.data(), y.data(), z.data(), cells.data(), state.data(), 0.0, 0.0, 0.0,
      1.0, 1.0, 1.0, kNx, kNy, kNz, agents, nullptr);
  synchronize();
  const auto mapped = download(cells, agents);
  assert(mapped[0] == 1 * kNx + 0);
  assert(mapped[1] == -1);

  DeviceBuffer<double> c(kCells);
  DeviceBuffer<double> iron(kCells);
  DeviceBuffer<double> b12(kCells);
  DeviceBuffer<double> acetate(kCells);
  DeviceBuffer<double> eut(kCells);
  DeviceBuffer<double> oxygen(kCells);
  c.upload(std::vector<double>(kCells, 1.0));
  iron.upload(std::vector<double>(kCells, 1.0));
  b12.upload(std::vector<double>(kCells, 1.0));
  acetate.upload(std::vector<double>(kCells, 0.0));
  eut.upload(std::vector<double>(kCells, 0.0));
  oxygen.upload(std::vector<double>(kCells, 0.0));
  DeviceBuffer<double> reac_c(kCells);
  DeviceBuffer<double> reac_i(kCells);
  DeviceBuffer<double> reac_b(kCells);
  reac_c.upload(std::vector<double>(kCells, 0.25));
  reac_i.upload(std::vector<double>(kCells, 0.5));
  reac_b.upload(std::vector<double>(kCells, 0.0));
  DeviceBuffer<double> mu(agents);
  DeviceBuffer<double> biomass(agents);
  DeviceBuffer<double> radius(agents);
  DeviceBuffer<double> mass(agents);
  DeviceBuffer<double> age(agents);
  DeviceBuffer<double> mu_max(agents);
  DeviceBuffer<double> km_b12(agents);
  DeviceBuffer<double> km_carbon(agents);
  DeviceBuffer<double> receptor(8 * agents);
  DeviceBuffer<double> ligand(8 * agents);
  DeviceBuffer<int> loci(agents);
  DeviceBuffer<double> amelioration(agents);
  DeviceBuffer<double> uptake(2);
  mu_max.upload(std::vector<double>(agents, 1.0e-3));
  km_b12.upload(std::vector<double>(agents, 0.1));
  km_carbon.upload(std::vector<double>(agents, 0.1));
  receptor.upload(std::vector<double>(8 * agents, 1.0));
  ligand.upload(std::vector<double>(8 * agents, 1.0));
  loci.upload(std::vector<int>(agents, 0));
  amelioration.upload(std::vector<double>(agents, 0.0));
  biomass.upload(std::vector<double>(agents, 1.0));
  radius.upload(std::vector<double>(agents, 1.0));
  mass.upload(std::vector<double>(agents, 1.0));
  age.upload(std::vector<double>(agents, 0.0));
  uptake.upload(std::vector<double>(2, 0.0));
  gutibm::gpu::launch_metabolism_kernel(
      c.data(), iron.data(), b12.data(), acetate.data(), eut.data(),
      reac_c.data(), reac_i.data(), reac_b.data(), mu.data(), biomass.data(),
      radius.data(), mass.data(), age.data(), cells.data(), state.data(),
      mu_max.data(), km_b12.data(), km_carbon.data(), receptor.data(),
      ligand.data(), loci.data(), amelioration.data(), agents, 1.0, 1.0, 1.0,
      0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2,
      0.1, 0.1, 0, 0.0, 1.0, oxygen.data(), uptake.data(), kNx, kNy, kNx,
      0, kNx, 0, nullptr);
  synchronize();
  const auto reaction = download(reac_c, kCells);
  const auto uptake_host = download(uptake, 2);
  assert(reaction[1 * kNx] > 0.0);
  assert(uptake_host[0] > 0.0);
  assert(uptake_host[1] > 0.0);
}

void test_diffusion_boundaries_and_gradient() {
  DeviceBuffer<double> field(kCells);
  const std::vector<double> uniform(kCells, 2.0);
  field.upload(uniform);
  DeviceBuffer<double> correction(kNx);
  correction.upload(std::vector<double>(kNx, 0.0));
  gutibm::gpu::launch_diffuse_x_periodic(
      field.data(), kNx, kNy, kNz, 0.1, 0.1, 0.1, 1.0, correction.data(),
      nullptr);
  gutibm::gpu::launch_diffuse_y_periodic(
      field.data(), kNx, kNy, kNz, 0, kNx, 0.1, 0.1, 0.1, 1.0,
      correction.data(), nullptr);
  synchronize();
  for (const double value : download(field, kCells)) {
    assert(std::isfinite(value));
    assert(value >= 0.0);
  }

  DeviceBuffer<double> exchange(1);
  exchange.upload(std::vector<double>{0.0});
  gutibm::gpu::launch_diffuse_z_bounded(
      field.data(), kNx, kNy, kNz, 0, kNx, 0.1, 1.0, 1.0, exchange.data(),
      nullptr);
  gutibm::gpu::launch_set_epithelial_boundary(
      field.data(), kNx, kNy, 0, kNx, 1.0, 1.0, exchange.data(), nullptr);
  gutibm::gpu::launch_set_luminal_neumann(
      field.data(), kNx, kNy, kNz, 0, kNx, nullptr);
  gutibm::gpu::launch_shift_z_gradient(
      field.data(), kNx, kNy, kNz, 0, kNx, 1.0, 1.0, 2.0, 1.0, 1.0,
      nullptr);
  gutibm::gpu::launch_clamp_nonneg(
      field.data(), kNx, kNy, kNz, 0, kNx, nullptr);
  synchronize();
  for (const double value : download(field, kCells)) {
    assert(std::isfinite(value));
    assert(value >= 0.0);
  }
}

void test_sources_and_reactions() {
  const DomainParams dom = domain_params();
  const AdvectionParams adv = zero_advection();
  DeviceBuffer<double> sx(1);
  DeviceBuffer<double> sy(1);
  DeviceBuffer<double> sz(1);
  DeviceBuffer<GfSourceParams> params(1);
  DeviceBuffer<double> grid(kCells);
  sx.upload(std::vector<double>{1.5});
  sy.upload(std::vector<double>{1.5});
  sz.upload(std::vector<double>{1.5});
  params.upload(std::vector<GfSourceParams>{{1.0, 1.0, 1.0, 0.0}});
  grid.upload(std::vector<double>(kCells, 0.0));
  gutibm::gpu::launch_superpose_kernel(
      sx.data(), sy.data(), sz.data(), params.data(), grid.data(), dom, adv, 1,
      1, 1, 1, nullptr);
  synchronize();
  const auto source_field = download(grid, kCells);
  assert(std::isfinite(*std::max_element(source_field.begin(), source_field.end())));
  assert(*std::max_element(source_field.begin(), source_field.end()) > 0.0);

  DeviceBuffer<double> leaf_local(1);
  DeviceBuffer<double> leaf_center(3);
  DeviceBuffer<int> cell_leaf(kCells);
  DeviceBuffer<double> out(kCells);
  leaf_local.upload(std::vector<double>{1.0});
  leaf_center.upload(std::vector<double>{1.5, 1.5, 1.5});
  cell_leaf.upload(std::vector<int>(kCells, 0));
  out.upload(std::vector<double>(kCells, 0.0));
  gutibm::gpu::launch_fmm_far_local_kernel(
      leaf_local.data(), leaf_center.data(), cell_leaf.data(), nullptr,
      out.data(), kCells, 1, 1, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, kNx, kNy,
      nullptr);
  synchronize();
  for (const double value : download(out, kCells)) assert(value == 1.0);
  gutibm::gpu::launch_fmm_far_local_kernel(
      leaf_local.data(), leaf_center.data(), cell_leaf.data(), nullptr,
      out.data(), kCells, 1, 1, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, kNx, kNy,
      nullptr);
  synchronize();
  for (const double value : download(out, kCells)) assert(value == 2.0);

  constexpr int storage_nx = 6;
  const int storage_cells = storage_nx * kNy * kNz;
  DeviceBuffer<double> reac_c(storage_cells);
  DeviceBuffer<double> reac_i(storage_cells);
  DeviceBuffer<double> reac_o(storage_cells);
  DeviceBuffer<double> reac_a(storage_cells);
  DeviceBuffer<double> reac_m(storage_cells);
  DeviceBuffer<double> conc(storage_cells);
  reac_c.upload(std::vector<double>(storage_cells, 0.0));
  reac_i.upload(std::vector<double>(storage_cells, 0.0));
  reac_o.upload(std::vector<double>(storage_cells, 0.0));
  reac_a.upload(std::vector<double>(storage_cells, 0.0));
  reac_m.upload(std::vector<double>(storage_cells, 0.0));
  conc.upload(std::vector<double>(storage_cells, 1.0));
  VbfLaunchParams vbf{};
  vbf.storage_nx = storage_nx;
  vbf.owned_x_begin = 1;
  vbf.owned_x_end = 3;
  vbf.global_nx = kNx;
  vbf.ny = kNy;
  vbf.nz = kNz;
  vbf.dx_x = 1.0;
  vbf.dx_y = 1.0;
  vbf.dx_z = 1.0;
  vbf.mucin_liberation = 0.1;
  vbf.vbf_density = 1.0;
  gutibm::gpu::launch_vbf_coupling_kernel(
      2 * kNy * kNz, vbf, reac_c.data(), conc.data(), reac_i.data(), conc.data(),
      reac_o.data(), conc.data(), reac_a.data(), reac_m.data(), conc.data(),
      nullptr, 1.0, nullptr);
  synchronize();
  const auto vbf_reaction = download(reac_c, storage_cells);
  for (int iz = 0; iz < kNz; ++iz) {
    for (int iy = 0; iy < kNy; ++iy) {
      for (int ix = 1; ix < 3; ++ix) {
        const int cell = iz * (storage_nx * kNy) + iy * storage_nx + ix;
        assert(vbf_reaction[cell] > 0.0);
      }
    }
  }

  DeviceBuffer<double> mu(2);
  DeviceBuffer<int> agent_cell(2);
  DeviceBuffer<int> agent_state(2);
  DeviceBuffer<double> o2_reaction(storage_cells);
  mu.upload(std::vector<double>{1.0, 1.0});
  agent_cell.upload(std::vector<int>{1, 1});
  agent_state.upload(std::vector<int>{0, 0});
  o2_reaction.upload(std::vector<double>(storage_cells, 0.0));
  gutibm::gpu::launch_o2_depletion_kernel(
      o2_reaction.data(), mu.data(), agent_cell.data(), agent_state.data(), 2,
      0.5, 0.1, 1.0, kNx, kNy, storage_nx, 1, 3, 1, nullptr);
  synchronize();
  const auto o2 = download(o2_reaction, storage_cells);
  assert(o2[1] < 0.0);
  assert(std::abs(o2[1] + 1.2) < 1.0e-12);
}

void test_hash_mechanics_and_receptor() {
  constexpr int agents = 2;
  DeviceBuffer<double> x(agents);
  DeviceBuffer<double> y(agents);
  DeviceBuffer<double> z(agents);
  DeviceBuffer<int> state(agents);
  DeviceBuffer<int> keys(agents);
  DeviceBuffer<int> sorted(agents);
  x.upload(std::vector<double>{0.2, 1.2});
  y.upload(std::vector<double>{0.2, 1.2});
  z.upload(std::vector<double>{0.2, 1.2});
  state.upload(std::vector<int>{0, 0});
  gutibm::gpu::launch_spatial_hash_build_kernel(
      x.data(), y.data(), z.data(), state.data(), keys.data(), sorted.data(),
      agents, 0.0, 0.0, 0.0, 1.0, kNx, kNy, kNz, nullptr);
  synchronize();
  const auto hash_keys = download(keys, agents);
  const auto hash_sorted = download(sorted, agents);
  assert(hash_keys[0] == 0);
  assert(hash_keys[1] == 1 + kNx + kNx * kNy);
  assert(hash_sorted[0] == 0);
  assert(hash_sorted[1] == 1);

  DeviceBuffer<double> radius(agents);
  DeviceBuffer<double> dx(agents);
  DeviceBuffer<double> dy(agents);
  DeviceBuffer<double> dz(agents);
  DeviceBuffer<int> offsets(2);
  DeviceBuffer<int> clamp(1);
  radius.upload(std::vector<double>{0.5, 0.5});
  dx.upload(std::vector<double>{3.0, 3.0});
  dy.upload(std::vector<double>{0.0, 0.0});
  dz.upload(std::vector<double>{0.0, 0.0});
  offsets.upload(std::vector<int>{0, 2});
  clamp.upload(std::vector<int>{0});
  gutibm::gpu::launch_mechanics_clear_kernel(
      dx.data(), dy.data(), dz.data(), clamp.data(), agents, nullptr);
  synchronize();
  const auto cleared_dx = download(dx, agents);
  assert(std::all_of(cleared_dx.begin(), cleared_dx.end(),
                     [](double value) { return value == 0.0; }));

  dx.upload(std::vector<double>{0.0, 0.0});
  dy.upload(std::vector<double>{0.0, 0.0});
  dz.upload(std::vector<double>{0.0, 0.0});
  x.upload(std::vector<double>{0.5, 0.8});
  y.upload(std::vector<double>{0.5, 0.5});
  z.upload(std::vector<double>{0.5, 0.5});
  MechanicsLaunchParams mech{};
  mech.hertzian_enabled = 1;
  mech.hertz_k = 1.0;
  mech.viscosity = 1.0;
  mech.dt = 1.0;
  mech.max_displacement_fraction = 0.1;
  mech.cell_size = 1.0;
  mech.nx_cells = 1;
  mech.ny_cells = 1;
  mech.nz_cells = 1;
  mech.lo0 = 0.0;
  mech.lo1 = 0.0;
  mech.lo2 = 0.0;
  mech.hi0 = 4.0;
  mech.hi1 = 4.0;
  mech.hi2 = 4.0;
  gutibm::gpu::launch_mechanics_forces_kernel(
      x.data(), y.data(), z.data(), radius.data(), state.data(), offsets.data(),
      sorted.data(), dx.data(), dy.data(), dz.data(), agents, mech, nullptr);
  synchronize();
  const auto force = download(dx, agents);
  assert(force[0] < 0.0);
  assert(force[1] > 0.0);
  gutibm::gpu::launch_mechanics_apply_kernel(
      x.data(), y.data(), z.data(), radius.data(), dx.data(), dy.data(),
      dz.data(), clamp.data(), agents, mech, nullptr);
  synchronize();
  const auto moved = download(x, agents);
  assert(moved[0] < 0.5);
  assert(moved[1] > 0.8);
  assert(std::abs(moved[0] - 0.5) <= 0.05);
  assert(std::abs(moved[1] - 0.8) <= 0.05);

  DeviceBuffer<int> grid_cell(1);
  DeviceBuffer<int> receptor_state(1);
  DeviceBuffer<double> receptors(8);
  DeviceBuffer<double> ligands(8);
  DeviceBuffer<double> toxin_affinity(8);
  DeviceBuffer<double> immunity(4);
  DeviceBuffer<double> toxin(1);
  DeviceBuffer<double> no_toxin(1);
  DeviceBuffer<double> kill(1);
  grid_cell.upload(std::vector<int>{0});
  receptor_state.upload(std::vector<int>{0});
  receptors.upload(std::vector<double>(8, 1.0));
  ligands.upload(std::vector<double>(8, 1.0));
  toxin_affinity.upload(std::vector<double>(8, 1.0));
  immunity.upload(std::vector<double>(4, 1.0));
  toxin.upload(std::vector<double>{1.0});
  no_toxin.upload(std::vector<double>{0.0});
  gutibm::gpu::launch_receptor_kill_prob_kernel(
      grid_cell.data(), receptor_state.data(), receptors.data(), ligands.data(),
      toxin_affinity.data(), immunity.data(), toxin.data(), no_toxin.data(),
      no_toxin.data(), no_toxin.data(), no_toxin.data(), no_toxin.data(),
      no_toxin.data(), kill.data(), 1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, kNx, kNy, kNx, 0, kNx, 0, nullptr);
  synchronize();
  const auto kill_host = download(kill, 1);
  assert(kill_host[0] > 0.0 && kill_host[0] < 1.0);
  DeviceBuffer<double> competitor(1);
  competitor.upload(std::vector<double>{10.0});
  toxin.upload(std::vector<double>{1.0});
  gutibm::gpu::launch_receptor_kill_prob_kernel(
      grid_cell.data(), receptor_state.data(), receptors.data(), ligands.data(),
      toxin_affinity.data(), immunity.data(), toxin.data(), no_toxin.data(),
      no_toxin.data(), no_toxin.data(), competitor.data(), no_toxin.data(),
      no_toxin.data(), kill.data(), 1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, kNx, kNy, kNx, 0, kNx, 0, nullptr);
  synchronize();
  assert(download(kill, 1)[0] < kill_host[0]);
  toxin.upload(std::vector<double>{0.5});
  gutibm::gpu::launch_receptor_kill_prob_kernel(
      grid_cell.data(), receptor_state.data(), receptors.data(), ligands.data(),
      toxin_affinity.data(), immunity.data(), toxin.data(), no_toxin.data(),
      no_toxin.data(), no_toxin.data(), no_toxin.data(), no_toxin.data(),
      no_toxin.data(), kill.data(), 1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, kNx, kNy, kNx, 0, kNx, 0, nullptr);
  synchronize();
  assert(download(kill, 1)[0] < kill_host[0]);
}

#endif

}  // namespace

int main() {
  std::cout << "=== Direct GPU Kernel Unit Tests ===\n";
  const int gpu_status = gutibm::test::require_gpu("gpu_kernel_units");
  if (gpu_status != 0) return gpu_status;
#ifdef GUTIBM_CUDA
  test_field_and_boundaries();
  test_grid_and_metabolism();
  test_diffusion_boundaries_and_gradient();
  test_sources_and_reactions();
  test_hash_mechanics_and_receptor();
  std::cout << "All 20 GPU launch entry points were invoked directly.\n";
  return 0;
#else
  return 77;
#endif
}
