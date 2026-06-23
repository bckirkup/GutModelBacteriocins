/* -----------------------------------------------------------------------
   GutIBM – HDF5 writer implementation
   ----------------------------------------------------------------------- */

#include "hdf5_writer.h"
#include "simulation.h"

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

namespace gutibm {

void HDF5Writer::init(const HDF5Config& cfg) {
  cfg_ = cfg;

  if (!cfg_.enabled) {
    enabled_ = false;
    return;
  }

#ifdef GUTIBM_HDF5
  enabled_ = true;

  hid_t plist = H5P_DEFAULT;
#ifdef GUTIBM_MPI
  if (cfg_.parallel) {
    plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);
  }
#endif

  file_id_ = static_cast<int64_t>(
      H5Fcreate(cfg_.filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist));

#ifdef GUTIBM_MPI
  if (cfg_.parallel && plist != H5P_DEFAULT) {
    H5Pclose(plist);
  }
#endif

  if (file_id_ < 0) {
    enabled_ = false;
  }
#endif
}

void HDF5Writer::write_step(const Simulation& sim, Int step, Real time) {
#ifdef GUTIBM_HDF5
  if (!enabled_) return;
  if (step % cfg_.dump_every != 0) return;

  hid_t fid = static_cast<hid_t>(file_id_);

  // Create step group: /step_NNNNNN
  std::ostringstream oss;
  oss << "step_" << std::setw(6) << std::setfill('0') << step;
  std::string gname = oss.str();

  hid_t step_group = H5Gcreate2(fid, gname.c_str(), H5P_DEFAULT,
                                  H5P_DEFAULT, H5P_DEFAULT);

  write_agents(sim, gname);
  write_grid(sim, gname);
  write_metadata(sim, gname, step, time);
  write_lineage(sim, gname);

  H5Gclose(step_group);
#else
  (void)sim; (void)step; (void)time;
#endif
}

void HDF5Writer::write_agents(const Simulation& sim, const std::string& group) {
#ifdef GUTIBM_HDF5
  hid_t fid = static_cast<hid_t>(file_id_);
  std::string agroup = group + "/atoms";
  hid_t ag = H5Gcreate2(fid, agroup.c_str(), H5P_DEFAULT,
                          H5P_DEFAULT, H5P_DEFAULT);

  const auto& agents = sim.agents();
  Int n = agents.size();

  // Collect local arrays
  std::vector<int64_t> ids(n);
  std::vector<int32_t> types(n), states(n);
  std::vector<double>  x(n), y(n), z(n), radius(n), biomass(n), mu(n);
  std::vector<int64_t> lineage(n);

  for (Int i = 0; i < n; ++i) {
    const Agent& a = agents[i];
    ids[i]     = a.tag;
    types[i]   = a.type;
    states[i]  = static_cast<int32_t>(a.state);
    x[i]       = a.x[0];
    y[i]       = a.x[1];
    z[i]       = a.x[2];
    radius[i]  = a.radius;
    biomass[i] = a.biomass;
    mu[i]      = a.mu_realized;
    lineage[i] = a.genome.lineage_id;
  }

#ifdef GUTIBM_MPI
  if (cfg_.parallel) {
    // Parallel HDF5: each rank writes its own slice
    // Gather counts from all ranks to compute offsets
    int local_n = static_cast<int>(n);
    int nprocs = 1;
    int my_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    std::vector<int> counts(nprocs);
    MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);

    hsize_t total = 0;
    hsize_t offset = 0;
    for (int r = 0; r < nprocs; ++r) {
      if (r < my_rank) offset += counts[r];
      total += counts[r];
    }

    hsize_t local_count = static_cast<hsize_t>(local_n);
    hid_t filespace = H5Screate_simple(1, &total, nullptr);
    hid_t memspace  = H5Screate_simple(1, &local_count, nullptr);

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, nullptr,
                        &local_count, nullptr);

    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);

    auto write_ds_par = [&](const char* name, hid_t type, const void* data) {
      std::string dsname = agroup + "/" + name;
      hid_t ds = H5Dcreate2(fid, dsname.c_str(), type, filespace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(ds, type, memspace, filespace, xfer_plist, data);
      H5Dclose(ds);
    };

    write_ds_par("id",       H5T_NATIVE_INT64,  ids.data());
    write_ds_par("type",     H5T_NATIVE_INT32,  types.data());
    write_ds_par("state",    H5T_NATIVE_INT32,  states.data());
    write_ds_par("x",        H5T_NATIVE_DOUBLE, x.data());
    write_ds_par("y",        H5T_NATIVE_DOUBLE, y.data());
    write_ds_par("z",        H5T_NATIVE_DOUBLE, z.data());
    write_ds_par("radius",   H5T_NATIVE_DOUBLE, radius.data());
    write_ds_par("biomass",  H5T_NATIVE_DOUBLE, biomass.data());
    write_ds_par("mu",       H5T_NATIVE_DOUBLE, mu.data());
    write_ds_par("lineage",  H5T_NATIVE_INT64,  lineage.data());

    H5Pclose(xfer_plist);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Gclose(ag);
    return;
  }
#endif

  // Serial path: single-rank writes
  hsize_t dims[1] = {static_cast<hsize_t>(n)};
  hid_t space = H5Screate_simple(1, dims, nullptr);

  auto write_ds = [&](const char* name, hid_t type, const void* data) {
    std::string dsname = agroup + "/" + name;
    hid_t ds = H5Dcreate2(fid, dsname.c_str(), type, space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(ds);
  };

  write_ds("id",       H5T_NATIVE_INT64,  ids.data());
  write_ds("type",     H5T_NATIVE_INT32,  types.data());
  write_ds("state",    H5T_NATIVE_INT32,  states.data());
  write_ds("x",        H5T_NATIVE_DOUBLE, x.data());
  write_ds("y",        H5T_NATIVE_DOUBLE, y.data());
  write_ds("z",        H5T_NATIVE_DOUBLE, z.data());
  write_ds("radius",   H5T_NATIVE_DOUBLE, radius.data());
  write_ds("biomass",  H5T_NATIVE_DOUBLE, biomass.data());
  write_ds("mu",       H5T_NATIVE_DOUBLE, mu.data());
  write_ds("lineage",  H5T_NATIVE_INT64,  lineage.data());

  H5Sclose(space);
  H5Gclose(ag);
#else
  (void)sim; (void)group;
#endif
}

void HDF5Writer::write_grid(const Simulation& sim, const std::string& group) {
#ifdef GUTIBM_HDF5
  hid_t fid = static_cast<hid_t>(file_id_);
  std::string ggroup = group + "/grid";
  hid_t gg = H5Gcreate2(fid, ggroup.c_str(), H5P_DEFAULT,
                          H5P_DEFAULT, H5P_DEFAULT);

  // In parallel mode, only rank 0 writes full grid data
#ifdef GUTIBM_MPI
  if (cfg_.parallel) {
    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank != 0) {
      H5Gclose(gg);
      return;
    }
  }
#endif

  const auto& chem = sim.chemical_field();
  Int ncells = chem.ncells();
  hsize_t dims[1] = {static_cast<hsize_t>(ncells)};
  hid_t space = H5Screate_simple(1, dims, nullptr);

  for (Int s = 0; s < chem.num_species(); ++s) {
    std::string dsname = ggroup + "/" + chem.spec(s).name;
    hid_t ds = H5Dcreate2(fid, dsname.c_str(), H5T_NATIVE_DOUBLE, space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
              chem.conc_data()[s].data());
    H5Dclose(ds);
  }

  H5Sclose(space);
  H5Gclose(gg);
#else
  (void)sim; (void)group;
#endif
}

void HDF5Writer::write_metadata(const Simulation& sim, const std::string& group,
                                  Int step, Real time) {
#ifdef GUTIBM_HDF5
  hid_t fid = static_cast<hid_t>(file_id_);
  std::string mgroup = group + "/metadata";
  hid_t mg = H5Gcreate2(fid, mgroup.c_str(), H5P_DEFAULT,
                          H5P_DEFAULT, H5P_DEFAULT);

  // In parallel mode, only rank 0 writes metadata
#ifdef GUTIBM_MPI
  if (cfg_.parallel) {
    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank != 0) {
      H5Gclose(mg);
      return;
    }
  }
#endif

  hsize_t one = 1;
  hid_t scalar = H5Screate_simple(1, &one, nullptr);

  auto write_attr = [&](const char* name, hid_t type, const void* val) {
    std::string dsname = mgroup + "/" + name;
    hid_t ds = H5Dcreate2(fid, dsname.c_str(), type, scalar,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
    H5Dclose(ds);
  };

  double t = time;
  int32_t n_agents = static_cast<int32_t>(sim.global_agent_count());
  int32_t n_lin    = static_cast<int32_t>(
      sim.lineage_tracker().snapshots().empty() ? 0
      : sim.lineage_tracker().snapshots().back().num_lineages);
  int32_t s = step;

  write_attr("time",         H5T_NATIVE_DOUBLE, &t);
  write_attr("step",         H5T_NATIVE_INT32,  &s);
  write_attr("num_agents",   H5T_NATIVE_INT32,  &n_agents);
  write_attr("num_lineages", H5T_NATIVE_INT32,  &n_lin);

  H5Sclose(scalar);
  H5Gclose(mg);
#else
  (void)sim; (void)group; (void)step; (void)time;
#endif
}

void HDF5Writer::write_lineage(const Simulation& sim, const std::string& group) {
#ifdef GUTIBM_HDF5
  hid_t fid = static_cast<hid_t>(file_id_);
  std::string lgroup = group + "/lineage";
  hid_t lg = H5Gcreate2(fid, lgroup.c_str(), H5P_DEFAULT,
                          H5P_DEFAULT, H5P_DEFAULT);

  // In parallel mode, only rank 0 writes lineage data
#ifdef GUTIBM_MPI
  if (cfg_.parallel) {
    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank != 0) {
      H5Gclose(lg);
      return;
    }
  }
#endif

  // Write per-agent receptor expression
  const auto& agents = sim.agents();
  Int n = agents.size();

  std::vector<double> btuB_expr(n), fepA_expr(n);
  std::vector<int32_t> n_bi(n), generation(n);

  for (Int i = 0; i < n; ++i) {
    const Agent& a = agents[i];
    btuB_expr[i]  = a.receptor_expr[static_cast<int>(ReceptorType::BtuB)];
    fepA_expr[i]  = a.receptor_expr[static_cast<int>(ReceptorType::FepA)];
    n_bi[i]       = static_cast<int32_t>(a.genome.bi_loci.size());
    generation[i] = static_cast<int32_t>(a.genome.generation);
  }

  hsize_t dims[1] = {static_cast<hsize_t>(n)};
  hid_t space = H5Screate_simple(1, dims, nullptr);

  auto write_ds = [&](const char* name, hid_t type, const void* data) {
    std::string dsname = lgroup + "/" + name;
    hid_t ds = H5Dcreate2(fid, dsname.c_str(), type, space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(ds);
  };

  write_ds("btuB_expression", H5T_NATIVE_DOUBLE, btuB_expr.data());
  write_ds("fepA_expression", H5T_NATIVE_DOUBLE, fepA_expr.data());
  write_ds("num_bi_loci",     H5T_NATIVE_INT32,  n_bi.data());
  write_ds("generation",      H5T_NATIVE_INT32,  generation.data());

  H5Sclose(space);
  H5Gclose(lg);
#else
  (void)sim; (void)group;
#endif
}

void HDF5Writer::finalize() {
#ifdef GUTIBM_HDF5
  if (enabled_ && file_id_ >= 0) {
    H5Fclose(static_cast<hid_t>(file_id_));
    file_id_ = -1;
  }
#endif
}

}  // namespace gutibm
