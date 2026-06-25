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
#include <numeric>

namespace gutibm {

namespace {

#ifdef GUTIBM_HDF5

#ifdef GUTIBM_MPI
bool mpi_is_active() {
  int initialized = 0;
  MPI_Initialized(&initialized);
  return initialized != 0;
}

int mpi_rank_world() {
  if (!mpi_is_active()) return 0;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int mpi_nprocs_world() {
  if (!mpi_is_active()) return 1;
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  return nprocs;
}

bool mpi_multi_rank() {
  return mpi_is_active() && mpi_nprocs_world() > 1;
}

void mpi_barrier_world() {
  if (mpi_is_active()) {
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
#else
int mpi_rank_world() { return 0; }
int mpi_nprocs_world() { return 1; }
bool mpi_multi_rank() { return false; }
void mpi_barrier_world() {}
#endif

int io_rank(const HDF5Config& cfg) {
  (void)cfg;
  return mpi_rank_world();
}

void mpi_barrier(const HDF5Config& cfg) {
  (void)cfg;
  if (mpi_multi_rank()) mpi_barrier_world();
}

void ensure_group(hid_t fid, const std::string& path, const HDF5Config& cfg) {
  if (io_rank(cfg) == 0 && fid >= 0) {
    hid_t g = H5Gcreate2(fid, path.c_str(), H5P_DEFAULT,
                          H5P_DEFAULT, H5P_DEFAULT);
    if (g < 0) {
      H5Eclear2(H5E_DEFAULT);
      g = H5Gopen2(fid, path.c_str(), H5P_DEFAULT);
    }
    if (g >= 0) H5Gclose(g);
  }
  mpi_barrier(cfg);
}

struct ParSlice {
  hsize_t offset = 0;
  hsize_t count  = 0;
  hsize_t total  = 0;
};

ParSlice compute_par_slice(hsize_t local_count, const HDF5Config& cfg) {
  (void)cfg;
  ParSlice slice;
  slice.count = local_count;
#ifdef GUTIBM_MPI
  if (mpi_multi_rank()) {
    int local_n = static_cast<int>(local_count);
    int nprocs = mpi_nprocs_world();
    int my_rank = mpi_rank_world();
    std::vector<int> counts(static_cast<size_t>(nprocs));
    MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    for (int r = 0; r < nprocs; ++r) {
      if (r < my_rank) slice.offset += static_cast<hsize_t>(counts[r]);
      slice.total += static_cast<hsize_t>(counts[r]);
    }
    return slice;
  }
#endif
  slice.offset = 0;
  slice.total = local_count;
  return slice;
}

template <typename T>
void write_dataset_1d(hid_t fid, const std::string& path, hid_t h5_type,
                      const T* local_data, hsize_t local_len,
                      const HDF5Config& cfg) {
  ParSlice slice = compute_par_slice(local_len, cfg);

#ifdef GUTIBM_MPI
  if (mpi_multi_rank()) {
    int rank = mpi_rank_world();
    int nprocs = mpi_nprocs_world();

    const int local_n = static_cast<int>(local_len);
    std::vector<int> counts(static_cast<size_t>(nprocs));
    std::vector<int> displs(static_cast<size_t>(nprocs));
    MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    displs[0] = 0;
    for (int r = 1; r < nprocs; ++r) {
      displs[r] = displs[r - 1] + counts[r - 1];
    }

    std::vector<T> global(static_cast<size_t>(slice.total));
    std::vector<int> byte_counts(static_cast<size_t>(nprocs));
    std::vector<int> byte_displs(static_cast<size_t>(nprocs));
    for (int r = 0; r < nprocs; ++r) {
      byte_counts[r] = counts[r] * static_cast<int>(sizeof(T));
      byte_displs[r] = displs[r] * static_cast<int>(sizeof(T));
    }
    T dummy{};
    const T* send_ptr = (local_n > 0 && local_data != nullptr) ? local_data : &dummy;
    MPI_Gatherv(send_ptr, local_n * static_cast<int>(sizeof(T)), MPI_BYTE,
                global.data(), byte_counts.data(), byte_displs.data(),
                MPI_BYTE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      hsize_t dims[1] = {slice.total};
      hid_t space = H5Screate_simple(1, dims, nullptr);
      hid_t ds = H5Dcreate2(fid, path.c_str(), h5_type, space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (ds < 0) {
        H5Eclear2(H5E_DEFAULT);
        ds = H5Dopen2(fid, path.c_str(), H5P_DEFAULT);
      }
      if (slice.total > 0) {
        H5Dwrite(ds, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, global.data());
      }
      H5Dclose(ds);
      H5Sclose(space);
    }
    mpi_barrier(cfg);
    return;
  }
#endif

  hsize_t dims[1] = {slice.total};
  hid_t space = H5Screate_simple(1, dims, nullptr);
  hid_t ds = H5Dcreate2(fid, path.c_str(), h5_type, space,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, local_data);
  H5Dclose(ds);
  H5Sclose(space);
}

#endif  // GUTIBM_HDF5

}  // namespace

void HDF5Writer::init(const HDF5Config& cfg) {
  cfg_ = cfg;

  if (!cfg_.enabled) {
    enabled_ = false;
    return;
  }

#ifdef GUTIBM_HDF5
  enabled_ = true;

#ifdef GUTIBM_MPI
  if (mpi_multi_rank() && mpi_rank_world() != 0) {
    file_id_ = -1;
    return;
  }
#endif

  file_id_ = static_cast<int64_t>(
      H5Fcreate(cfg_.filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));

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

  std::ostringstream oss;
  oss << "step_" << std::setw(6) << std::setfill('0') << step;
  std::string gname = oss.str();

  int skip = 0;
  if (io_rank(cfg_) == 0 && fid >= 0) {
    skip = (H5Lexists(fid, gname.c_str(), H5P_DEFAULT) > 0) ? 1 : 0;
  }
#ifdef GUTIBM_MPI
  if (mpi_multi_rank()) {
    MPI_Bcast(&skip, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif
  if (skip) return;

  ensure_group(fid, gname, cfg_);

  write_agents(sim, gname);
  write_lineage(sim, gname);
  write_genome(sim, gname);
  write_grid(sim, gname);
  write_metadata(sim, gname, step, time);
#else
  (void)sim; (void)step; (void)time;
#endif
}

void HDF5Writer::write_agents(const Simulation& sim, const std::string& group) {
#ifdef GUTIBM_HDF5
  hid_t fid = static_cast<hid_t>(file_id_);
  std::string agroup = group + "/atoms";
  ensure_group(fid, agroup, cfg_);

  const auto& agents = sim.agents();
  Int n = agents.size();

  std::vector<int64_t> ids(static_cast<size_t>(n));
  std::vector<int32_t> types(static_cast<size_t>(n)), states(static_cast<size_t>(n));
  std::vector<double>  x(static_cast<size_t>(n)), y(static_cast<size_t>(n)),
                       z(static_cast<size_t>(n)), radius(static_cast<size_t>(n)),
                       biomass(static_cast<size_t>(n)), mu(static_cast<size_t>(n));
  std::vector<int64_t> lineage(static_cast<size_t>(n));

  for (Int i = 0; i < n; ++i) {
    const Agent& a = agents[i];
    ids[static_cast<size_t>(i)]     = a.tag;
    types[static_cast<size_t>(i)]   = a.type;
    states[static_cast<size_t>(i)]  = static_cast<int32_t>(a.state);
    x[static_cast<size_t>(i)]       = a.x[0];
    y[static_cast<size_t>(i)]       = a.x[1];
    z[static_cast<size_t>(i)]       = a.x[2];
    radius[static_cast<size_t>(i)]  = a.radius;
    biomass[static_cast<size_t>(i)] = a.biomass;
    mu[static_cast<size_t>(i)]      = a.mu_realized;
    lineage[static_cast<size_t>(i)] = a.genome.lineage_id;
  }

  const hsize_t local_n = static_cast<hsize_t>(n);
  write_dataset_1d(fid, agroup + "/id",       H5T_NATIVE_INT64,  ids.data(), local_n, cfg_);
  write_dataset_1d(fid, agroup + "/type",     H5T_NATIVE_INT32,  types.data(), local_n, cfg_);
  write_dataset_1d(fid, agroup + "/state",    H5T_NATIVE_INT32,  states.data(), local_n, cfg_);
  write_dataset_1d(fid, agroup + "/x",        H5T_NATIVE_DOUBLE, x.data(), local_n, cfg_);
  write_dataset_1d(fid, agroup + "/y",        H5T_NATIVE_DOUBLE, y.data(), local_n, cfg_);
  write_dataset_1d(fid, agroup + "/z",        H5T_NATIVE_DOUBLE, z.data(), local_n, cfg_);
  write_dataset_1d(fid, agroup + "/radius",   H5T_NATIVE_DOUBLE, radius.data(), local_n, cfg_);
  write_dataset_1d(fid, agroup + "/biomass",  H5T_NATIVE_DOUBLE, biomass.data(), local_n, cfg_);
  write_dataset_1d(fid, agroup + "/mu",       H5T_NATIVE_DOUBLE, mu.data(), local_n, cfg_);
  write_dataset_1d(fid, agroup + "/lineage",  H5T_NATIVE_INT64,  lineage.data(), local_n, cfg_);
#else
  (void)sim; (void)group;
#endif
}

void HDF5Writer::write_grid(const Simulation& sim, const std::string& group) {
#ifdef GUTIBM_HDF5
  hid_t fid = static_cast<hid_t>(file_id_);
  std::string ggroup = group + "/grid";
  ensure_group(fid, ggroup, cfg_);

  if (io_rank(cfg_) == 0 && fid >= 0) {
    const auto& chem = sim.chemical_field();
    Int ncells = chem.ncells();
    hsize_t dims[1] = {static_cast<hsize_t>(ncells)};
    hid_t space = H5Screate_simple(1, dims, nullptr);

    for (Int s = 0; s < chem.num_species(); ++s) {
      std::string dsname = ggroup + "/" + chem.spec(s).name;
      hid_t ds = H5Dcreate2(fid, dsname.c_str(), H5T_NATIVE_DOUBLE, space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (ds < 0) {
        H5Eclear2(H5E_DEFAULT);
        ds = H5Dopen2(fid, dsname.c_str(), H5P_DEFAULT);
      }
      H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               chem.conc_data()[s].data());
      H5Dclose(ds);
    }

    H5Sclose(space);
  }
  mpi_barrier(cfg_);
#else
  (void)sim; (void)group;
#endif
}

void HDF5Writer::write_metadata(const Simulation& sim, const std::string& group,
                                  Int step, Real time) {
#ifdef GUTIBM_HDF5
  hid_t fid = static_cast<hid_t>(file_id_);
  std::string mgroup = group + "/metadata";
  ensure_group(fid, mgroup, cfg_);

  if (io_rank(cfg_) == 0 && fid >= 0) {
    hsize_t one = 1;
    hid_t scalar = H5Screate_simple(1, &one, nullptr);

    double t = time;
    int32_t n_agents = static_cast<int32_t>(sim.global_agent_count());
    int32_t n_lin    = static_cast<int32_t>(
        sim.lineage_tracker().snapshots().empty() ? 0
        : sim.lineage_tracker().snapshots().back().num_lineages);
    int32_t s = step;

    auto write_scalar = [&](const char* name, hid_t type, const void* val) {
      std::string dsname = mgroup + "/" + name;
      hid_t ds = H5Dcreate2(fid, dsname.c_str(), type, scalar,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (ds < 0) {
        H5Eclear2(H5E_DEFAULT);
        ds = H5Dopen2(fid, dsname.c_str(), H5P_DEFAULT);
      }
      H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
      H5Dclose(ds);
    };

    write_scalar("time",         H5T_NATIVE_DOUBLE, &t);
    write_scalar("step",         H5T_NATIVE_INT32,  &s);
    write_scalar("num_agents",   H5T_NATIVE_INT32,  &n_agents);
    write_scalar("num_lineages", H5T_NATIVE_INT32,  &n_lin);

    H5Sclose(scalar);
  }
  mpi_barrier(cfg_);
#else
  (void)sim; (void)group; (void)step; (void)time;
#endif
}

void HDF5Writer::write_lineage(const Simulation& sim, const std::string& group) {
#ifdef GUTIBM_HDF5
  hid_t fid = static_cast<hid_t>(file_id_);
  std::string lgroup = group + "/lineage";
  ensure_group(fid, lgroup, cfg_);

  const auto& agents = sim.agents();
  Int n = agents.size();

  std::vector<double> btuB_expr(static_cast<size_t>(n)), fepA_expr(static_cast<size_t>(n));
  std::vector<int32_t> n_bi(static_cast<size_t>(n)), generation(static_cast<size_t>(n));

  for (Int i = 0; i < n; ++i) {
    const Agent& a = agents[i];
    btuB_expr[static_cast<size_t>(i)]  = a.receptor_expr[static_cast<int>(ReceptorType::BtuB)];
    fepA_expr[static_cast<size_t>(i)]  = a.receptor_expr[static_cast<int>(ReceptorType::FepA)];
    n_bi[static_cast<size_t>(i)]       = static_cast<int32_t>(a.genome.bi_loci.size());
    generation[static_cast<size_t>(i)] = static_cast<int32_t>(a.genome.generation);
  }

  const hsize_t local_n = static_cast<hsize_t>(n);
  write_dataset_1d(fid, lgroup + "/btuB_expression", H5T_NATIVE_DOUBLE,
                    btuB_expr.data(), local_n, cfg_);
  write_dataset_1d(fid, lgroup + "/fepA_expression", H5T_NATIVE_DOUBLE,
                    fepA_expr.data(), local_n, cfg_);
  write_dataset_1d(fid, lgroup + "/num_bi_loci", H5T_NATIVE_INT32,
                    n_bi.data(), local_n, cfg_);
  write_dataset_1d(fid, lgroup + "/generation", H5T_NATIVE_INT32,
                    generation.data(), local_n, cfg_);
#else
  (void)sim; (void)group;
#endif
}

void HDF5Writer::write_genome(const Simulation& sim, const std::string& group) {
#ifdef GUTIBM_HDF5
  hid_t fid = static_cast<hid_t>(file_id_);
  std::string ggroup = group + "/genome";
  ensure_group(fid, ggroup, cfg_);

  const auto& agents = sim.agents();
  Int n = agents.size();
  const hsize_t local_n = static_cast<hsize_t>(n);

  std::vector<int64_t> parent_id(static_cast<size_t>(n));
  std::vector<int32_t> mutations(static_cast<size_t>(n)),
                       has_conjugative(static_cast<size_t>(n));
  std::vector<double> plasmid_amel(static_cast<size_t>(n));
  std::vector<double> receptor_expr(static_cast<size_t>(n) * NUM_RECEPTORS),
                      toxin_aff(static_cast<size_t>(n) * NUM_RECEPTORS),
                      ligand_aff(static_cast<size_t>(n) * NUM_RECEPTORS);

  std::vector<int32_t> bi_toxin_id, bi_immunity_id, bi_target, bi_bclass;
  std::vector<double> bi_pI, bi_diff, bi_ret, bi_mw, bi_imm_aff;

  for (Int i = 0; i < n; ++i) {
    const Agent& a = agents[i];
    const size_t idx = static_cast<size_t>(i);
    parent_id[idx] = a.genome.parent_id;
    mutations[idx] = static_cast<int32_t>(a.genome.mutations);
    has_conjugative[idx] = a.genome.has_conjugative_plasmid ? 1 : 0;
    plasmid_amel[idx] = a.genome.plasmid_cost_amelioration;
    for (Int r = 0; r < NUM_RECEPTORS; ++r) {
      receptor_expr[idx * NUM_RECEPTORS + static_cast<size_t>(r)] =
          a.genome.receptor_expression[r];
      toxin_aff[idx * NUM_RECEPTORS + static_cast<size_t>(r)] =
          a.genome.toxin_affinity[r];
      ligand_aff[idx * NUM_RECEPTORS + static_cast<size_t>(r)] =
          a.genome.ligand_affinity[r];
    }
    for (const auto& bi : a.genome.bi_loci) {
      bi_toxin_id.push_back(static_cast<int32_t>(bi.toxin_id));
      bi_immunity_id.push_back(static_cast<int32_t>(bi.immunity_id));
      bi_target.push_back(static_cast<int32_t>(bi.target));
      bi_bclass.push_back(static_cast<int32_t>(bi.bclass));
      bi_pI.push_back(bi.pI);
      bi_diff.push_back(bi.diff_coeff);
      bi_ret.push_back(bi.retardation);
      bi_mw.push_back(bi.molecular_weight);
      bi_imm_aff.push_back(bi.immunity_binding_affinity);
    }
  }

  write_dataset_1d(fid, ggroup + "/parent_id", H5T_NATIVE_INT64,
                   parent_id.data(), local_n, cfg_);
  write_dataset_1d(fid, ggroup + "/mutations", H5T_NATIVE_INT32,
                   mutations.data(), local_n, cfg_);
  write_dataset_1d(fid, ggroup + "/has_conjugative_plasmid", H5T_NATIVE_INT32,
                   has_conjugative.data(), local_n, cfg_);
  write_dataset_1d(fid, ggroup + "/plasmid_cost_amelioration", H5T_NATIVE_DOUBLE,
                   plasmid_amel.data(), local_n, cfg_);
  write_dataset_1d(fid, ggroup + "/receptor_expression", H5T_NATIVE_DOUBLE,
                   receptor_expr.data(), local_n * NUM_RECEPTORS, cfg_);
  write_dataset_1d(fid, ggroup + "/toxin_affinity", H5T_NATIVE_DOUBLE,
                   toxin_aff.data(), local_n * NUM_RECEPTORS, cfg_);
  write_dataset_1d(fid, ggroup + "/ligand_affinity", H5T_NATIVE_DOUBLE,
                   ligand_aff.data(), local_n * NUM_RECEPTORS, cfg_);

  const hsize_t local_bi = static_cast<hsize_t>(bi_toxin_id.size());
  write_dataset_1d(fid, ggroup + "/bi_toxin_id", H5T_NATIVE_INT32,
                   bi_toxin_id.data(), local_bi, cfg_);
  write_dataset_1d(fid, ggroup + "/bi_immunity_id", H5T_NATIVE_INT32,
                   bi_immunity_id.data(), local_bi, cfg_);
  write_dataset_1d(fid, ggroup + "/bi_target", H5T_NATIVE_INT32,
                   bi_target.data(), local_bi, cfg_);
  write_dataset_1d(fid, ggroup + "/bi_bclass", H5T_NATIVE_INT32,
                   bi_bclass.data(), local_bi, cfg_);
  write_dataset_1d(fid, ggroup + "/bi_pI", H5T_NATIVE_DOUBLE,
                   bi_pI.data(), local_bi, cfg_);
  write_dataset_1d(fid, ggroup + "/bi_diff_coeff", H5T_NATIVE_DOUBLE,
                   bi_diff.data(), local_bi, cfg_);
  write_dataset_1d(fid, ggroup + "/bi_retardation", H5T_NATIVE_DOUBLE,
                   bi_ret.data(), local_bi, cfg_);
  write_dataset_1d(fid, ggroup + "/bi_molecular_weight", H5T_NATIVE_DOUBLE,
                   bi_mw.data(), local_bi, cfg_);
  write_dataset_1d(fid, ggroup + "/bi_immunity_binding_affinity", H5T_NATIVE_DOUBLE,
                   bi_imm_aff.data(), local_bi, cfg_);
#else
  (void)sim; (void)group;
#endif
}

void HDF5Writer::finalize() {
#ifdef GUTIBM_HDF5
  if (mpi_multi_rank()) {
    mpi_barrier_world();
  }
  if (enabled_ && file_id_ >= 0) {
    H5Fclose(static_cast<hid_t>(file_id_));
    file_id_ = -1;
  }
  if (mpi_multi_rank()) {
    mpi_barrier_world();
  }
#endif
}

}  // namespace gutibm
