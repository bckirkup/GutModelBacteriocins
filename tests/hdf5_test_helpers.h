/* -----------------------------------------------------------------------
   GutIBM – Shared HDF5 test helpers
   ----------------------------------------------------------------------- */

#pragma once

#include "agent.h"
#include "hdf5_reader.h"
#include "simulation.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

namespace gutibm {
namespace test {

constexpr double kAgentSnapshotTol = 1e-12;

struct AgentSnapshot {
  int64_t id;
  int32_t type;
  int32_t state;
  double x;
  double y;
  double z;
  double radius;
  double biomass;
  double mu;
  int64_t lineage;
};

inline std::vector<AgentSnapshot> collect_agent_snapshots(const Simulation& sim) {
  std::vector<AgentSnapshot> out;
  out.reserve(sim.agents().size());
  for (const Agent& a : sim.agents()) {
    out.emplace_back(
        a.identity.tag,
        a.identity.type,
        static_cast<int32_t>(to_underlying(a.state)),
        a.x[0], a.x[1], a.x[2],
        a.radius,
        a.biomass,
        a.mu_realized,
        a.genome.lineage_id);
  }
  std::ranges::sort(out, [](const AgentSnapshot& lhs, const AgentSnapshot& rhs) {
    return lhs.id < rhs.id;
  });
  return out;
}

inline std::vector<AgentSnapshot> agent_snapshots_from_checkpoint(
    const HDF5CheckpointAgents& atoms) {
  const size_t n = atoms.id.size();
  std::vector<AgentSnapshot> out(n);
  size_t i = 0;
  for (int64_t id : atoms.id) {
    (void)id;
    out[i] = AgentSnapshot{
        atoms.id[i], atoms.type[i], atoms.state[i],
        atoms.x[i], atoms.y[i], atoms.z[i],
        atoms.radius[i], atoms.biomass[i], atoms.mu[i], atoms.lineage[i],
    };
    ++i;
  }
  std::ranges::sort(out, [](const AgentSnapshot& lhs, const AgentSnapshot& rhs) {
    return lhs.id < rhs.id;
  });
  return out;
}

inline void compare_agent_snapshots(const std::vector<AgentSnapshot>& expected,
                                    const std::vector<AgentSnapshot>& actual,
                                    double tol = kAgentSnapshotTol) {
  assert(expected.size() == actual.size());
  auto it_actual = actual.begin();
  for (const AgentSnapshot& exp : expected) {
    const AgentSnapshot& act = *it_actual++;
    assert(exp.id == act.id);
    assert(exp.type == act.type);
    assert(exp.state == act.state);
    assert(std::abs(exp.x - act.x) < tol);
    assert(std::abs(exp.y - act.y) < tol);
    assert(std::abs(exp.z - act.z) < tol);
    assert(std::abs(exp.radius - act.radius) < tol);
    assert(std::abs(exp.biomass - act.biomass) < tol);
    assert(std::abs(exp.mu - act.mu) < tol);
    assert(exp.lineage == act.lineage);
  }
}

#ifdef GUTIBM_HDF5

inline bool hdf5_dataset_exists(hid_t file, const std::string& path) {
  return H5Lexists(file, path.c_str(), H5P_DEFAULT) > 0;
}

template <typename T>
inline std::vector<T> hdf5_read_dataset_1d(hid_t file, const std::string& path,
                                           hid_t h5_type) {
  hid_t dset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  assert(dset >= 0);

  hid_t space = H5Dget_space(dset);
  int ndims = H5Sget_simple_extent_ndims(space);
  assert(ndims == 1);
  hsize_t len = 0;
  H5Sget_simple_extent_dims(space, &len, nullptr);

  std::vector<T> data(static_cast<size_t>(len));
  H5Dread(dset, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

  H5Sclose(space);
  H5Dclose(dset);
  return data;
}

template <typename T>
inline T hdf5_read_scalar(hid_t file, const std::string& path, hid_t h5_type) {
  auto data = hdf5_read_dataset_1d<T>(file, path, h5_type);
  assert(data.size() == 1);
  return data[0];
}

inline std::vector<AgentSnapshot> read_agent_snapshots(hid_t file,
                                                       const std::string& step) {
  const std::string prefix = "agents/" + step + "/";
  auto ids = hdf5_read_dataset_1d<int64_t>(file, prefix + "id", H5T_NATIVE_INT64);
  auto types = hdf5_read_dataset_1d<int32_t>(file, prefix + "type", H5T_NATIVE_INT32);
  auto states = hdf5_read_dataset_1d<int32_t>(file, prefix + "state", H5T_NATIVE_INT32);
  auto xs = hdf5_read_dataset_1d<double>(file, prefix + "x", H5T_NATIVE_DOUBLE);
  auto ys = hdf5_read_dataset_1d<double>(file, prefix + "y", H5T_NATIVE_DOUBLE);
  auto zs = hdf5_read_dataset_1d<double>(file, prefix + "z", H5T_NATIVE_DOUBLE);
  auto radii = hdf5_read_dataset_1d<double>(file, prefix + "radius", H5T_NATIVE_DOUBLE);
  auto biomass = hdf5_read_dataset_1d<double>(file, prefix + "biomass", H5T_NATIVE_DOUBLE);
  auto mus = hdf5_read_dataset_1d<double>(file, prefix + "mu_realized", H5T_NATIVE_DOUBLE);
  auto lineages =
      hdf5_read_dataset_1d<int64_t>(file, prefix + "lineage_id", H5T_NATIVE_INT64);

  const size_t n = ids.size();
  assert(types.size() == n);
  assert(states.size() == n);
  assert(xs.size() == n);
  assert(ys.size() == n);
  assert(zs.size() == n);
  assert(radii.size() == n);
  assert(biomass.size() == n);
  assert(mus.size() == n);
  assert(lineages.size() == n);

  std::vector<AgentSnapshot> out(n);
  size_t i = 0;
  for (int64_t id : ids) {
    (void)id;
    out[i] = AgentSnapshot{
        ids[i], types[i], states[i],
        xs[i], ys[i], zs[i],
        radii[i], biomass[i], mus[i], lineages[i],
    };
    ++i;
  }
  std::ranges::sort(out, [](const AgentSnapshot& lhs, const AgentSnapshot& rhs) {
    return lhs.id < rhs.id;
  });
  return out;
}

#endif  // GUTIBM_HDF5

}  // namespace test
}  // namespace gutibm
