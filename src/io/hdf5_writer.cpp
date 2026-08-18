/* -----------------------------------------------------------------------
   GutIBM – HDF5 writer implementation (Spec 4 layered schema)
   ----------------------------------------------------------------------- */

#include "hdf5_writer.h"
#include "path_utils.h"
#include "simulation.h"
#include "species_names.h"
#include "step_events.h"
#include "config_json.h"
#include "error.h"

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <format>
#include <iostream>
#include <limits>
#include <numeric>
#include <ranges>
#include <set>
#include <system_error>
#include <string>
#include <vector>
#include <cstdlib>
#include "error.h"

namespace gutibm {

namespace {

constexpr int k_max_types = 8;
constexpr int k_num_pheno_states = 4;

bool schedule_has_output(const HDF5Schedule& sched) {
  return sched.summary > 0 || sched.agents > 0 || sched.grid > 0 ||
         sched.lineage > 0 || sched.genome > 0 || sched.provenance > 0;
}

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
    hid_t g = -1;
    if (H5Lexists(fid, path.c_str(), H5P_DEFAULT) > 0) {
      g = H5Gopen2(fid, path.c_str(), H5P_DEFAULT);
    } else {
      g = H5Gcreate2(fid, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (g >= 0) H5Gclose(g);
  }
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
    auto local_n = static_cast<int>(local_count);
    auto nprocs = mpi_nprocs_world();
    auto my_rank = mpi_rank_world();
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
void write_dataset_1d_serial(hid_t fid, const std::string& path, hid_t h5_type,
                             const T* data, hsize_t len) {
  std::array<hsize_t, 1> dims = {len};
  hid_t space = H5Screate_simple(1, dims.data(), nullptr);
  hid_t ds = H5Dcreate2(fid, path.c_str(), h5_type, space,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (ds < 0) {
    H5Eclear2(H5E_DEFAULT);
    ds = H5Dopen2(fid, path.c_str(), H5P_DEFAULT);
  }
  if (len > 0) {
    H5Dwrite(ds, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  }
  H5Dclose(ds);
  H5Sclose(space);
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

    const auto local_n = static_cast<int>(local_len);
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
      std::array<hsize_t, 1> dims = {slice.total};
      hid_t space = H5Screate_simple(1, dims.data(), nullptr);
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

  std::array<hsize_t, 1> dims = {slice.total};
  hid_t space = H5Screate_simple(1, dims.data(), nullptr);
  hid_t ds = H5Dcreate2(fid, path.c_str(), h5_type, space,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (slice.total > 0) {
    H5Dwrite(ds, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, local_data);
  }
  H5Dclose(ds);
  H5Sclose(space);
}

void write_scalar_dataset(hid_t fid, const std::string& path, hid_t h5_type,
                          const void* value) {
  hsize_t one = 1;
  hid_t space = H5Screate_simple(1, &one, nullptr);
  hid_t ds = H5Dcreate2(fid, path.c_str(), h5_type, space,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (ds < 0) {
    H5Eclear2(H5E_DEFAULT);
    ds = H5Dopen2(fid, path.c_str(), H5P_DEFAULT);
  }
  H5Dwrite(ds, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
  H5Dclose(ds);
  H5Sclose(space);
}

void write_file_attr(hid_t fid, const char* name, hid_t type, const void* value) {
  hid_t space = H5Screate(H5S_SCALAR);
  hid_t attr = H5Acreate2(fid, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, type, value);
  H5Aclose(attr);
  H5Sclose(space);
}

void write_string_dataset(hid_t fid, const std::string& path,
                          const std::string& value) {
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, H5T_VARIABLE);
  hid_t space = H5Screate(H5S_SCALAR);
  hid_t ds = H5Dcreate2(fid, path.c_str(), type, space,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (ds < 0) {
    H5Eclear2(H5E_DEFAULT);
    ds = H5Dopen2(fid, path.c_str(), H5P_DEFAULT);
  }
  const char* text = value.c_str();
  if (ds >= 0) H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &text);
  if (ds >= 0) H5Dclose(ds);
  H5Sclose(space);
  H5Tclose(type);
}

hid_t make_dataset_plist(const HDF5Config& cfg, const hsize_t* chunk_dims,
                         int rank) {
  hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
  if (cfg.compression == "gzip") {
    H5Pset_chunk(plist, rank, chunk_dims);
    const auto level = static_cast<unsigned>(
        std::clamp(cfg.compression_level, 0, 9));
    H5Pset_deflate(plist, level);
  }
  return plist;
}

void pack_grid_species(const ChemicalField& chem, const Domain& domain,
                       Int species_index, Int nx, Int ny, Int nz,
                       std::vector<double>& grid3d) {
  for (Int iz = 0; iz < nz; ++iz) {
    for (Int iy = 0; iy < ny; ++iy) {
      for (Int ix = 0; ix < nx; ++ix) {
        const Int flat = domain.cell_index(ix, iy, iz);
        const size_t idx = static_cast<size_t>(iz) * static_cast<size_t>(nx * ny)
            + static_cast<size_t>(iy) * static_cast<size_t>(nx)
            + static_cast<size_t>(ix);
        grid3d[idx] = chem.conc(species_index, flat);
      }
    }
  }
}

Real field_mean(const ChemicalField& chem, Int species_idx) {
  if (species_idx < 0) return 0.0;
  Real sum = 0.0;
  Int n = 0;
  for (Int iz = 0; iz < chem.global_nz(); ++iz) {
    for (Int iy = 0; iy < chem.global_ny(); ++iy) {
      for (Int ix = chem.owned_storage_x_begin();
           ix < chem.owned_storage_x_end(); ++ix) {
        const Int c = iz * chem.storage_nx() * chem.global_ny()
            + iy * chem.storage_nx() + ix;
        sum += chem.conc(species_idx, c);
        ++n;
      }
    }
  }
#ifdef GUTIBM_MPI
  if (chem.slab_mode()) {
    Real global_sum = 0.0;
    Int global_n = 0;
    MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&n, &global_n, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    sum = global_sum;
    n = global_n;
  }
#endif
  return n > 0 ? sum / static_cast<Real>(n) : 0.0;
}

Real field_max(const ChemicalField& chem, Int species_idx) {
  if (species_idx < 0) return 0.0;
  Real mx = 0.0;
  for (Int iz = 0; iz < chem.global_nz(); ++iz) {
    for (Int iy = 0; iy < chem.global_ny(); ++iy) {
      for (Int ix = chem.owned_storage_x_begin();
           ix < chem.owned_storage_x_end(); ++ix) {
        const Int c = iz * chem.storage_nx() * chem.global_ny()
            + iy * chem.storage_nx() + ix;
        mx = std::max(mx, chem.conc(species_idx, c));
      }
    }
  }
#ifdef GUTIBM_MPI
  if (chem.slab_mode()) {
    Real global_max = 0.0;
    MPI_Allreduce(&mx, &global_max, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    mx = global_max;
  }
#endif
  return mx;
}

int32_t count_live_lineages(const AgentPool& agents) {
  std::set<int64_t> unique_lineages;
  for (const Agent& agent : agents) {
    if (agent.state != PhenoState::DEAD) {
      unique_lineages.insert(static_cast<int64_t>(agent.genome.lineage_id));
    }
  }
  return static_cast<int32_t>(unique_lineages.size());
}

std::vector<char> species_name_table(const ChemicalField& chem) {
  constexpr size_t kSpeciesNameWidth = 48;
  std::vector names(chem.specs().size() * kSpeciesNameWidth, '\0');
  for (size_t i = 0; i < chem.specs().size(); ++i) {
    const std::string& name = chem.specs()[i].name;
    const size_t count = std::min(name.size(), kSpeciesNameWidth - 1);
    std::copy_n(name.data(), count, names.begin() + i * kSpeciesNameWidth);
  }
  return names;
}

std::vector<Real> boundary_flux_per_area(
    const std::vector<Real>& interval, Real area, Real interval_time) {
  std::vector result(interval.size(), 0.0);
  if (area <= 0.0 || interval_time <= 0.0) return result;
  for (size_t i = 0; i < result.size(); ++i) {
    result[i] = interval[i] / (area * interval_time);
  }
  return result;
}

void summarize_agents(
    const AgentPool& agents, std::array<int32_t, k_max_types>& n_by_type,
    std::array<int32_t, k_max_types>& n_in_crypt,
    std::array<int32_t, k_num_pheno_states>& n_by_state,
    std::array<double, k_max_types>& mean_z,
    std::array<double, k_max_types>& mean_mu) {
  std::array<int32_t, k_max_types> count_by_type{};
  for (const Agent& agent : agents) {
    if (agent.state == PhenoState::DEAD) continue;
    const Int tidx = std::clamp(agent.identity.type, 0, k_max_types - 1);
    n_by_type[static_cast<size_t>(tidx)]++;
    if (agent.flags.in_crypt) n_in_crypt[static_cast<size_t>(tidx)]++;
    mean_z[static_cast<size_t>(tidx)] += agent.x[2];
    mean_mu[static_cast<size_t>(tidx)] += agent.mu_realized;
    count_by_type[static_cast<size_t>(tidx)]++;
    const Int sidx = std::clamp(
        static_cast<Int>(to_underlying(agent.state)), 0,
        k_num_pheno_states - 1);
    n_by_state[static_cast<size_t>(sidx)]++;
  }
  for (Int type_idx = 0; type_idx < k_max_types; ++type_idx) {
    const auto index = static_cast<size_t>(type_idx);
    if (count_by_type[index] > 0) {
      const Real inv = 1.0 / static_cast<Real>(count_by_type[index]);
      mean_z[index] *= inv;
      mean_mu[index] *= inv;
    }
  }
}

std::vector<double> mean_receptor_expression(const AgentPool& agents) {
  std::vector result(NUM_RECEPTORS, 0.0);
  Int live = 0;
  for (const Agent& agent : agents) {
    if (agent.state == PhenoState::DEAD) continue;
    ++live;
    for (Int receptor = 0; receptor < NUM_RECEPTORS; ++receptor) {
      result[static_cast<size_t>(receptor)] += agent.receptor_expr[receptor];
    }
  }
  if (live > 0) {
    for (double& value : result) value /= static_cast<double>(live);
  }
  return result;
}

#endif  // GUTIBM_HDF5

}  // namespace

bool HDF5Writer::layer_due(Int interval, Int step) const {
  return interval > 0 && step % interval == 0;
}

bool HDF5Writer::should_write_species(const std::string& name) const {
  const auto& species = cfg_.schedule.grid_species;
  if (species.empty()) return false;
  return std::ranges::find(species, std::string("all")) != species.end()
      || std::ranges::find(species, name) != species.end();
}

#ifdef GUTIBM_HDF5
void HDF5Writer::initialize_file() {
  if (enabled_ && io_rank(cfg_) == 0) {
    // Parallel HDF5 builds can return invalid FILE_CREATE property lists from
    // H5Pcreate on some platforms; H5P_DEFAULT is reliable for rank-0 serial I/O.
    hid_t fcpl = H5P_DEFAULT;
    hid_t created_fcpl = H5Pcreate(H5P_FILE_CREATE);
    if (created_fcpl >= 0 && H5Iis_valid(created_fcpl) &&
        H5Pget_class(created_fcpl) == H5P_FILE_CREATE) {
      if (H5Pset_libver_bounds(created_fcpl, H5F_LIBVER_V18, H5F_LIBVER_LATEST) >= 0) {
        fcpl = created_fcpl;
        created_fcpl = H5I_INVALID_HID;
      } else {
        H5Eclear2(H5E_DEFAULT);
      }
    }
    if (created_fcpl >= 0 && H5Iis_valid(created_fcpl)) {
      H5Pclose(created_fcpl);
    }

    file_id_ = static_cast<int64_t>(
        H5Fcreate(cfg_.filename.c_str(), H5F_ACC_TRUNC, fcpl, H5P_DEFAULT));
    if (fcpl != H5P_DEFAULT && fcpl >= 0 && H5Iis_valid(fcpl)) {
      H5Pclose(fcpl);
    }

    if (file_id_ < 0 || !H5Iis_valid(static_cast<hid_t>(file_id_)) ||
        H5Fis_hdf5(cfg_.filename.c_str()) <= 0) {
      H5Eclear2(H5E_DEFAULT);
      if (file_id_ >= 0 && H5Iis_valid(static_cast<hid_t>(file_id_))) {
        H5Fclose(static_cast<hid_t>(file_id_));
      }
      file_id_ = -1;
      enabled_ = false;
      std::remove(cfg_.filename.c_str());
    } else {
      auto fid = static_cast<hid_t>(file_id_);
      const int32_t nx_attr = nx_;
      const int32_t ny_attr = ny_;
      const int32_t nz_attr = nz_;
      const double dx_attr = grid_dx_[0];
      const double dx_x_attr = grid_dx_[0];
      const double dx_y_attr = grid_dx_[1];
      const double dx_z_attr = grid_dx_[2];
      const int32_t version = 4;
      write_file_attr(fid, "gutibm_version", H5T_NATIVE_INT32, &version);
      write_file_attr(fid, "nx", H5T_NATIVE_INT32, &nx_attr);
      write_file_attr(fid, "ny", H5T_NATIVE_INT32, &ny_attr);
      write_file_attr(fid, "nz", H5T_NATIVE_INT32, &nz_attr);
      write_file_attr(fid, "grid_dx", H5T_NATIVE_DOUBLE, &dx_attr);
      write_file_attr(fid, "grid_dx_x", H5T_NATIVE_DOUBLE, &dx_x_attr);
      write_file_attr(fid, "grid_dx_y", H5T_NATIVE_DOUBLE, &dx_y_attr);
      write_file_attr(fid, "grid_dx_z", H5T_NATIVE_DOUBLE, &dx_z_attr);
      write_file_attr(fid, "domain_lo_x", H5T_NATIVE_DOUBLE, &domain_lo_[0]);
      write_file_attr(fid, "domain_lo_y", H5T_NATIVE_DOUBLE, &domain_lo_[1]);
      write_file_attr(fid, "domain_lo_z", H5T_NATIVE_DOUBLE, &domain_lo_[2]);
      write_file_attr(fid, "domain_hi_x", H5T_NATIVE_DOUBLE, &domain_hi_[0]);
      write_file_attr(fid, "domain_hi_y", H5T_NATIVE_DOUBLE, &domain_hi_[1]);
      write_file_attr(fid, "domain_hi_z", H5T_NATIVE_DOUBLE, &domain_hi_[2]);
    }
  }

}
#endif

void HDF5Writer::init(const HDF5Config& cfg, const Domain& domain) {
  cfg_ = cfg;
  nx_ = domain.nx();
  ny_ = domain.ny();
  nz_ = domain.nz();
  grid_dx_ = {domain.dx_x(), domain.dx_y(), domain.dx_z()};
  domain_lo_ = domain.lo();
  domain_hi_ = domain.hi();

  if (!cfg_.enabled || !schedule_has_output(cfg_.schedule)) {
    enabled_ = false;
    return;
  }

#ifdef GUTIBM_HDF5
  enabled_ = true;
  file_id_ = -1;

  if (io_rank(cfg_) == 0) {
    try {
      validate_output_file_path(cfg_.filename);
    } catch (const IOError& ex) {
      std::cerr << "Warning: invalid HDF5 output path '" << cfg_.filename
                << "': " << ex.what() << " — HDF5 output disabled\n";
      enabled_ = false;
    }
  }

  initialize_file();
  mpi_barrier(cfg_);
#ifdef GUTIBM_MPI
  if (mpi_multi_rank()) {
    int enabled_flag = enabled_ ? 1 : 0;
    MPI_Bcast(&enabled_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    enabled_ = enabled_flag != 0;
  }
#endif
#endif
}

void HDF5Writer::write_run_provenance(const Simulation& sim) const {
#ifdef GUTIBM_HDF5
  if (run_provenance_written_ || !enabled_ || io_rank(cfg_) != 0
      || file_id_ < 0) {
    return;
  }
  const auto fid = static_cast<hid_t>(file_id_);
  ensure_group(fid, "run_provenance", cfg_);
  const auto config = ConfigJson::serialize_document(sim.config());
  write_string_dataset(fid, "run_provenance/resolved_config", config);
  write_string_dataset(fid, "run_provenance/git_sha", GUTIBM_GIT_SHA);
  write_string_dataset(fid, "run_provenance/version", GUTIBM_VERSION);
  const int32_t mpi_compiled =
#ifdef GUTIBM_MPI
      1;
#else
      0;
#endif
  const int32_t hdf5_compiled = 1;
  const int32_t gpu_compiled =
#ifdef GUTIBM_CUDA
      1;
#else
      0;
#endif
  const int32_t openmp_compiled =
#ifdef GUTIBM_OPENMP
      1;
#else
      0;
#endif
  write_scalar_dataset(fid, "run_provenance/mpi_compiled",
                       H5T_NATIVE_INT32, &mpi_compiled);
  write_scalar_dataset(fid, "run_provenance/hdf5_compiled",
                       H5T_NATIVE_INT32, &hdf5_compiled);
  write_scalar_dataset(fid, "run_provenance/gpu_compiled",
                       H5T_NATIVE_INT32, &gpu_compiled);
  write_scalar_dataset(fid, "run_provenance/openmp_compiled",
                       H5T_NATIVE_INT32, &openmp_compiled);
  const int32_t rank_count = static_cast<int32_t>(mpi_nprocs_world());
  write_scalar_dataset(fid, "run_provenance/mpi_rank_count",
                       H5T_NATIVE_INT32, &rank_count);
  const auto optional_env = [&fid](const char* name, const char* env_name) {
    if (const char* value = std::getenv(env_name);
        value != nullptr && value[0] != '\0') {
      write_string_dataset(fid, std::string("run_provenance/") + name, value);
    }
  };
  optional_env("container_image_digest", "GUTIBM_IMAGE_DIGEST");
  optional_env("job_id", "AWS_BATCH_JOB_ID");
  run_provenance_written_ = true;
#else
  (void)sim;
#endif
}

void HDF5Writer::write_halt_metadata(const Simulation& sim, Int step) const {
#ifdef GUTIBM_HDF5
  if (enabled_ && io_rank(cfg_) == 0 && file_id_ >= 0) {
    write_run_provenance(sim);
    const auto fid = static_cast<hid_t>(file_id_);
    const std::string group = "summary/" + std::format("step_{:06}", step);
    ensure_group(fid, "summary", cfg_);
    ensure_group(fid, group, cfg_);
    const int32_t halt_reason = sim.halted_for_dysbiosis() ? 1 : 0;
    const double halt_density = sim.halt_density_cells_per_mL();
    write_scalar_dataset(fid, group + "/halt_reason_code",
                         H5T_NATIVE_INT32, &halt_reason);
    write_scalar_dataset(fid, group + "/halt_density_cells_per_mL",
                         H5T_NATIVE_DOUBLE, &halt_density);
    H5Fflush(fid, H5F_SCOPE_LOCAL);
  }
  mpi_barrier(cfg_);
#else
  (void)sim;
  (void)step;
#endif
}

void HDF5Writer::write_run_termination(const Simulation& sim, Int step,
                                       Real time) const {
#ifdef GUTIBM_HDF5
  if (enabled_ && io_rank(cfg_) == 0 && file_id_ >= 0) {
    write_run_provenance(sim);
    const auto fid = static_cast<hid_t>(file_id_);
    ensure_group(fid, "run_provenance", cfg_);
    const int32_t halt_reason = sim.halted_for_dysbiosis() ? 1 : 0;
    const double halt_density = sim.halt_density_cells_per_mL();
    const int32_t halt_step = sim.halted_for_dysbiosis() ? step : 0;
    const double halt_time = sim.halted_for_dysbiosis() ? time : 0.0;
    const int32_t completed_total_time =
        time >= sim.config().time.total_time ? 1 : 0;
    const int32_t termination_reason =
        sim.halted_for_dysbiosis() ? 1 : (completed_total_time != 0 ? 0 : 2);
    const int32_t termination_step = step;
    const double termination_time = time;
    write_scalar_dataset(fid, "run_provenance/halt_reason_code",
                         H5T_NATIVE_INT32, &halt_reason);
    write_scalar_dataset(fid, "run_provenance/halt_density_cells_per_mL",
                         H5T_NATIVE_DOUBLE, &halt_density);
    write_scalar_dataset(fid, "run_provenance/halt_step",
                         H5T_NATIVE_INT32, &halt_step);
    write_scalar_dataset(fid, "run_provenance/halt_time",
                         H5T_NATIVE_DOUBLE, &halt_time);
    write_scalar_dataset(fid, "run_provenance/completed_total_time",
                         H5T_NATIVE_INT32, &completed_total_time);
    write_scalar_dataset(fid, "run_provenance/termination_reason_code",
                         H5T_NATIVE_INT32, &termination_reason);
    write_scalar_dataset(fid, "run_provenance/termination_step",
                         H5T_NATIVE_INT32, &termination_step);
    write_scalar_dataset(fid, "run_provenance/termination_time",
                         H5T_NATIVE_DOUBLE, &termination_time);
    H5Fflush(fid, H5F_SCOPE_LOCAL);
  }
  mpi_barrier(cfg_);
#else
  (void)sim;
  (void)step;
  (void)time;
#endif
}

void HDF5Writer::write_step(Simulation& sim, Int step, Real time, Real dt) const {
#ifdef GUTIBM_HDF5
  if (!enabled_) return;

  write_run_provenance(sim);

  const std::string step_name = std::format("step_{:06}", step);
  auto fid = static_cast<hid_t>(file_id_);
  const bool summary_due = layer_due(cfg_.schedule.summary, step);

  if (summary_due) {
    sim.prepare_step_events_for_summary();
    sim.prepare_mechanics_stats_for_summary();
    sim.prepare_population_stocks_for_summary();
    const std::string path = "summary/" + step_name;
    ensure_group(fid, "summary", cfg_);
    ensure_group(fid, path, cfg_);
    write_summary(sim, path, step, time, dt);
  }
  if (layer_due(cfg_.schedule.agents, step)) {
    const std::string path = "agents/" + step_name;
    ensure_group(fid, "agents", cfg_);
    ensure_group(fid, path, cfg_);
    write_agents_layer(sim, path);
  }
  if (layer_due(cfg_.schedule.grid, step)) {
    const std::string path = "grid/" + step_name;
    ensure_group(fid, "grid", cfg_);
    ensure_group(fid, path, cfg_);
    write_grid_layer(sim, path);
  }
  if (layer_due(cfg_.schedule.lineage, step)) {
    const std::string path = "lineage/" + step_name;
    ensure_group(fid, "lineage", cfg_);
    ensure_group(fid, path, cfg_);
    write_lineage_layer(sim, path);
  }
  if (layer_due(cfg_.schedule.genome, step)) {
    const std::string path = "genome/" + step_name;
    ensure_group(fid, "genome", cfg_);
    ensure_group(fid, path, cfg_);
    write_genome_layer(sim, path);
  }
  if (layer_due(cfg_.schedule.provenance, step)) {
    const std::string path = "provenance/" + step_name;
    ensure_group(fid, "provenance", cfg_);
    ensure_group(fid, path, cfg_);
    write_provenance_layer(sim, path);
  }
  if (summary_due) {
    sim.commit_step_events_after_summary(step, time);
    sim.chemical_field().flux_accounting().close_interval();
  }

  // Flush so a concurrent checkpoint copy (entry.sh SIGSTOP+cp) sees consistent
  // metadata; without this, mid-write S3 uploads can be unreadable on resume.
  if (io_rank(cfg_) == 0 && file_id_ >= 0) {
    H5Fflush(fid, H5F_SCOPE_LOCAL);
  }

  mpi_barrier(cfg_);
#else
  (void)sim; (void)step; (void)time; (void)dt;
#endif
}

void HDF5Writer::write_summary(Simulation& sim, const std::string& group,
                                Int step, Real time, Real dt) const {
#ifdef GUTIBM_HDF5
  const auto& chem = sim.chemical_field();
  const double mean_carbon = field_mean(chem, chem.find(species::CARBON));
  const double mean_iron = field_mean(chem, chem.find(species::IRON));
  const double mean_oxygen = field_mean(chem, chem.find(species::OXYGEN));
  const auto toxin_max = [&sim, &chem](const char* name) {
    const Int idx = chem.find(name);
    Real value = sim.qssa().agent_sampling()
        ? sim.qssa().sampled_toxin_max(idx)
        : field_max(chem, idx);
#ifdef GUTIBM_MPI
    if (sim.qssa().agent_sampling()) {
      Real global_value = 0.0;
      MPI_Allreduce(&value, &global_value, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD);
      value = global_value;
    }
#endif
    return value;
  };
  const bool lumped = sim.qssa().toxin_lumping();
  const double max_btuB = toxin_max(
      species::bacteriocin_species_for(ReceptorType::BtuB, lumped));
  const double max_fepA = toxin_max(
      species::bacteriocin_species_for(ReceptorType::FepA, lumped));
  const double max_cirA = toxin_max(
      species::bacteriocin_species_for(ReceptorType::CirA, lumped));
  const double max_fhuA = toxin_max(
      species::bacteriocin_species_for(ReceptorType::FhuA, lumped));

  if (io_rank(cfg_) == 0 && file_id_ >= 0) {
  auto fid = static_cast<hid_t>(file_id_);
  const auto& agents = sim.agents();
  const auto& events = sim.summary_events();

  const double t = time;
  const double dt_val = dt;
  const int32_t step_val = step;
  const auto n_total = static_cast<int32_t>(sim.global_agent_count());

  const auto num_lineages = count_live_lineages(agents);
  const auto& stocks = sim.population_stocks();
  const auto& mechanics = sim.mechanics_summary_stats();

  write_scalar_dataset(fid, group + "/time", H5T_NATIVE_DOUBLE, &t);
  write_scalar_dataset(fid, group + "/dt", H5T_NATIVE_DOUBLE, &dt_val);
  write_scalar_dataset(fid, group + "/step", H5T_NATIVE_INT32, &step_val);
  write_scalar_dataset(fid, group + "/n_total", H5T_NATIVE_INT32, &n_total);
  write_scalar_dataset(fid, group + "/num_lineages", H5T_NATIVE_INT32, &num_lineages);
  write_scalar_dataset(fid, group + "/num_agents", H5T_NATIVE_INT32, &n_total);
  ensure_group(fid, group + "/stocks", cfg_);
  const int32_t bacteriostatic_live = stocks.bacteriostatic_live;
  const int32_t washout_trapped_live = stocks.washout_trapped_live;
  write_scalar_dataset(fid, group + "/stocks/bacteriostatic_live_agents",
                       H5T_NATIVE_INT32, &bacteriostatic_live);
  write_scalar_dataset(fid, group + "/stocks/washout_trapped_live_agents",
                       H5T_NATIVE_INT32, &washout_trapped_live);
  ensure_group(fid, group + "/mechanics", cfg_);
  const int32_t displacement_clamps = mechanics.displacement_clamps;
  const int32_t cumulative_displacement_clamps =
      sim.mechanics_cumulative_stats().displacement_clamps
      + displacement_clamps;
  write_scalar_dataset(fid, group + "/mechanics/displacement_clamps",
                       H5T_NATIVE_INT32, &displacement_clamps);
  write_scalar_dataset(fid,
                       group + "/mechanics/cumulative_displacement_clamps",
                       H5T_NATIVE_INT32, &cumulative_displacement_clamps);
  const int32_t halt_reason = sim.halted_for_dysbiosis() ? 1 : 0;
  const double halt_density = sim.halt_density_cells_per_mL();
  write_scalar_dataset(fid, group + "/halt_reason_code",
                       H5T_NATIVE_INT32, &halt_reason);
  write_scalar_dataset(fid, group + "/halt_density_cells_per_mL",
                       H5T_NATIVE_DOUBLE, &halt_density);
  ensure_group(fid, group + "/nutrient_flux", cfg_);
  const auto& flux = sim.chemical_field().flux_accounting();
  const auto write_flux = [&](const char* name,
                              const std::vector<Real>& values) {
    write_dataset_1d_serial(fid, group + "/nutrient_flux/" + name,
                             H5T_NATIVE_DOUBLE, values.data(), values.size());
  };
  constexpr size_t kSpeciesNameWidth = 48;
  const auto species_names = species_name_table(chem);
  std::array<hsize_t, 2> name_dims = {
      static_cast<hsize_t>(chem.specs().size()),
      static_cast<hsize_t>(kSpeciesNameWidth)};
  hid_t name_space = H5Screate_simple(2, name_dims.data(), nullptr);
  if (hid_t name_ds = H5Dcreate2(
          fid, (group + "/nutrient_flux/species_names").c_str(),
          H5T_NATIVE_CHAR, name_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      name_ds >= 0) {
    if (!species_names.empty()) {
      H5Dwrite(name_ds, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               species_names.data());
    }
    H5Dclose(name_ds);
  }
  H5Sclose(name_space);
  write_flux("boundary_interval", flux.boundary_interval);
  const auto add_cumulative = [](const std::vector<Real>& prior,
                                 const std::vector<Real>& interval) {
    std::vector values(prior.size(), 0.0);
    for (size_t i = 0; i < values.size(); ++i) {
      values[i] = prior[i] + interval[i];
    }
    return values;
  };
  write_flux("boundary_cumulative",
             add_cumulative(flux.boundary_cumulative, flux.boundary_interval));
  write_flux("vbf_source_interval", flux.vbf_source_interval);
  write_flux("vbf_source_cumulative",
             add_cumulative(flux.vbf_source_cumulative, flux.vbf_source_interval));
  write_flux("vbf_sink_interval", flux.vbf_sink_interval);
  write_flux("vbf_sink_cumulative",
             add_cumulative(flux.vbf_sink_cumulative, flux.vbf_sink_interval));
  write_flux("agent_uptake_interval", flux.agent_uptake_interval);
  write_flux("agent_uptake_cumulative",
             add_cumulative(flux.agent_uptake_cumulative,
                           flux.agent_uptake_interval));
  write_flux("reaction_clip_interval", flux.reaction_clip_interval);
  write_flux("reaction_clip_cumulative",
             add_cumulative(flux.reaction_clip_cumulative,
                           flux.reaction_clip_interval));
  const Real area = (sim.domain().hi()[0] - sim.domain().lo()[0])
      * (sim.domain().hi()[1] - sim.domain().lo()[1]);
  const Real interval_time = std::max(
      sim.time() - sim.event_window_start_time(), 0.0);
  const auto boundary_area_flux = boundary_flux_per_area(
      flux.boundary_interval, area, interval_time);
  write_flux("boundary_area_flux_interval", boundary_area_flux);
  const auto write_flux_bound = [&](const char* name, Real value) {
    write_scalar_dataset(fid, group + "/nutrient_flux/" + name,
                         H5T_NATIVE_DOUBLE, &value);
  };
  write_flux_bound("interval_start_step",
                   static_cast<Real>(sim.event_window_start_step()));
  write_flux_bound("interval_end_step", static_cast<Real>(step));
  write_flux_bound("interval_start_time", sim.event_window_start_time());
  write_flux_bound("interval_end_time", time);
  std::array<int32_t, k_max_types> n_by_type{};
  std::array<int32_t, k_max_types> n_in_crypt{};
  std::array<int32_t, k_num_pheno_states> n_by_state{};
  std::array<double, k_max_types> mean_z{};
  std::array<double, k_max_types> mean_mu{};

  summarize_agents(agents, n_by_type, n_in_crypt, n_by_state, mean_z, mean_mu);

  write_dataset_1d_serial(fid, group + "/n_by_type", H5T_NATIVE_INT32,
                          n_by_type.data(), k_max_types);
  write_dataset_1d_serial(fid, group + "/n_in_crypt", H5T_NATIVE_INT32,
                          n_in_crypt.data(), k_max_types);
  write_dataset_1d_serial(fid, group + "/n_by_state", H5T_NATIVE_INT32,
                          n_by_state.data(), 4);
  write_dataset_1d_serial(fid, group + "/mean_z_by_type", H5T_NATIVE_DOUBLE,
                          mean_z.data(), k_max_types);
  write_dataset_1d_serial(fid, group + "/mean_mu_by_type", H5T_NATIVE_DOUBLE,
                          mean_mu.data(), k_max_types);

  const auto mean_receptor = mean_receptor_expression(agents);
  write_dataset_1d_serial(fid, group + "/mean_receptor_expr", H5T_NATIVE_DOUBLE,
                          mean_receptor.data(), NUM_RECEPTORS);

  ensure_group(fid, group + "/events", cfg_);
  const auto write_event = [&](const char* name, Int val) {
    const int32_t v = val;
    write_scalar_dataset(fid, group + std::string("/events/") + name,
                         H5T_NATIVE_INT32, &v);
  };
  write_event("sos_inductions", events.sos_inductions);
  write_event("phage_inductions", events.phage_inductions);
  write_event("mortality_colicin", events.mortality_colicin);
  write_event("mortality_cdi", events.mortality_cdi);
  write_event("outflow_washout", events.outflow_washout);
  write_event("outflow_boundary", events.outflow_boundary);
  write_event("mortality_lysis", events.mortality_lysis);
  write_event("divisions", events.divisions);
  write_event("conjugation_transfers", events.conjugation_transfers);
  write_event("mutations", events.mutations);
  write_event("immigrations", events.immigrations);

  StepEvents cumulative = sim.cumulative_events();
  cumulative.add(events);
  const auto write_cumulative_event = [&](const char* name, Int val) {
    const int32_t v = val;
    write_scalar_dataset(fid, group + std::string("/events/cumulative_") + name,
                         H5T_NATIVE_INT32, &v);
  };
  write_cumulative_event("sos_inductions", cumulative.sos_inductions);
  write_cumulative_event("phage_inductions", cumulative.phage_inductions);
  write_cumulative_event("mortality_colicin", cumulative.mortality_colicin);
  write_cumulative_event("mortality_cdi", cumulative.mortality_cdi);
  write_cumulative_event("outflow_washout", cumulative.outflow_washout);
  write_cumulative_event("outflow_boundary", cumulative.outflow_boundary);
  write_cumulative_event("mortality_lysis", cumulative.mortality_lysis);
  write_cumulative_event("divisions", cumulative.divisions);
  write_cumulative_event("conjugation_transfers", cumulative.conjugation_transfers);
  write_cumulative_event("mutations", cumulative.mutations);
  write_cumulative_event("immigrations", cumulative.immigrations);

  const int32_t interval_start_step = sim.event_window_start_step();
  const int32_t interval_end_step = step;
  const double interval_start_time = sim.event_window_start_time();
  const double interval_end_time = time;
  write_scalar_dataset(fid, group + "/events/interval_start_step",
                       H5T_NATIVE_INT32, &interval_start_step);
  write_scalar_dataset(fid, group + "/events/interval_end_step",
                       H5T_NATIVE_INT32, &interval_end_step);
  write_scalar_dataset(fid, group + "/events/interval_start_time",
                       H5T_NATIVE_DOUBLE, &interval_start_time);
  write_scalar_dataset(fid, group + "/events/interval_end_time",
                       H5T_NATIVE_DOUBLE, &interval_end_time);

  ensure_group(fid, group + "/chem", cfg_);
  write_scalar_dataset(fid, group + "/chem/mean_carbon", H5T_NATIVE_DOUBLE, &mean_carbon);
  write_scalar_dataset(fid, group + "/chem/mean_iron", H5T_NATIVE_DOUBLE, &mean_iron);
  write_scalar_dataset(fid, group + "/chem/mean_oxygen", H5T_NATIVE_DOUBLE, &mean_oxygen);
  write_scalar_dataset(fid, group + "/chem/max_toxin_BtuB", H5T_NATIVE_DOUBLE, &max_btuB);
  write_scalar_dataset(fid, group + "/chem/max_toxin_FepA", H5T_NATIVE_DOUBLE, &max_fepA);
  write_scalar_dataset(fid, group + "/chem/max_toxin_CirA", H5T_NATIVE_DOUBLE, &max_cirA);
  write_scalar_dataset(fid, group + "/chem/max_toxin_FhuA", H5T_NATIVE_DOUBLE, &max_fhuA);

  }

  mpi_barrier(cfg_);
#else
  (void)sim; (void)group; (void)step; (void)time; (void)dt;
#endif
}

void HDF5Writer::write_provenance_layer(Simulation& sim,
                                         const std::string& group) const {
#ifdef GUTIBM_HDF5
  auto fid = static_cast<hid_t>(file_id_);
  const auto& events = sim.kill_provenance();
  const auto local_n = static_cast<hsize_t>(events.size());

  std::vector<int64_t> victim_id(events.size());
  std::vector<double> x(events.size());
  std::vector<double> y(events.size());
  std::vector<double> z(events.size());
  std::vector<int32_t> strain(events.size());
  std::vector<int32_t> cause(events.size());
  std::vector<int64_t> cdi_attacker_id(events.size());
  std::vector<int32_t> cdi_attacker_known(events.size());
  std::vector<double> toxin_concentration(events.size() * 4);
  std::vector<double> toxin_occupancy(events.size() * 4);
  std::vector<double> toxin_hazard(events.size() * 4);

  for (size_t i = 0; i < events.size(); ++i) {
    const auto& event = events[i];
    victim_id[i] = event.victim_id;
    x[i] = event.position[0];
    y[i] = event.position[1];
    z[i] = event.position[2];
    strain[i] = event.strain;
    cause[i] = to_underlying(event.cause);
    cdi_attacker_id[i] = event.cdi_attacker_id;
    cdi_attacker_known[i] = event.cdi_attacker_known ? 1 : 0;
    for (size_t toxin = 0; toxin < 4; ++toxin) {
      const size_t flat = i * 4 + toxin;
      toxin_concentration[flat] = event.toxin_concentration[toxin];
      toxin_occupancy[flat] = event.toxin_occupancy[toxin];
      toxin_hazard[flat] = event.toxin_hazard[toxin];
    }
  }

  write_dataset_1d(fid, group + "/victim_id", H5T_NATIVE_INT64,
                   victim_id.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/x", H5T_NATIVE_DOUBLE, x.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/y", H5T_NATIVE_DOUBLE, y.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/z", H5T_NATIVE_DOUBLE, z.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/strain", H5T_NATIVE_INT32,
                   strain.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/cause", H5T_NATIVE_INT32,
                   cause.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/cdi_attacker_id", H5T_NATIVE_INT64,
                   cdi_attacker_id.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/cdi_attacker_known", H5T_NATIVE_INT32,
                   cdi_attacker_known.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/toxin_concentration", H5T_NATIVE_DOUBLE,
                   toxin_concentration.data(), local_n * 4, cfg_);
  write_dataset_1d(fid, group + "/toxin_occupancy", H5T_NATIVE_DOUBLE,
                   toxin_occupancy.data(), local_n * 4, cfg_);
  write_dataset_1d(fid, group + "/toxin_hazard", H5T_NATIVE_DOUBLE,
                   toxin_hazard.data(), local_n * 4, cfg_);
  sim.clear_kill_provenance();
  mpi_barrier(cfg_);
#else
  (void)sim;
  (void)group;
#endif
}

void HDF5Writer::write_agents_layer(const Simulation& sim,
                                     const std::string& group) const {
#ifdef GUTIBM_HDF5
  auto fid = static_cast<hid_t>(file_id_);
  const auto agents = output_agents(sim);
  const auto n = static_cast<Int>(agents.size());

  std::vector<int64_t> ids(static_cast<size_t>(n));
  std::vector<int32_t> types(static_cast<size_t>(n));
  std::vector<int32_t> states(static_cast<size_t>(n));
  std::vector<double> x(static_cast<size_t>(n));
  std::vector<double> y(static_cast<size_t>(n));
  std::vector<double> z(static_cast<size_t>(n));
  std::vector<double> mu(static_cast<size_t>(n));
  std::vector<double> mu_max(static_cast<size_t>(n));
  std::vector<double> biomass(static_cast<size_t>(n));
  std::vector<int32_t> in_crypt(static_cast<size_t>(n));
  std::vector<int32_t> n_bi(static_cast<size_t>(n));
  std::vector<double> radius(static_cast<size_t>(n));
  std::vector<int64_t> lineage_id(static_cast<size_t>(n));
  std::vector<double> receptor_expr(static_cast<size_t>(n) * NUM_RECEPTORS);

  for (Int i = 0; i < n; ++i) {
    const Agent& a = *agents[static_cast<size_t>(i)];
    const auto idx = static_cast<size_t>(i);
    ids[idx] = a.identity.tag;
    types[idx] = a.identity.type;
    states[idx] = static_cast<int32_t>(to_underlying(a.state));
    x[idx] = a.x[0];
    y[idx] = a.x[1];
    z[idx] = a.x[2];
    mu[idx] = a.mu_realized;
    mu_max[idx] = a.mu_max;
    biomass[idx] = a.biomass;
    in_crypt[idx] = a.flags.in_crypt ? 1 : 0;
    n_bi[idx] = static_cast<int32_t>(a.genome.bi_loci.size());
    radius[idx] = a.radius;
    lineage_id[idx] = static_cast<int64_t>(a.genome.lineage_id);
    for (Int r = 0; r < NUM_RECEPTORS; ++r) {
      receptor_expr[idx * NUM_RECEPTORS + static_cast<size_t>(r)] = a.receptor_expr[r];
    }
  }

  const auto local_n = static_cast<hsize_t>(n);
  write_dataset_1d(fid, group + "/id", H5T_NATIVE_INT64, ids.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/type", H5T_NATIVE_INT32, types.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/state", H5T_NATIVE_INT32, states.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/x", H5T_NATIVE_DOUBLE, x.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/y", H5T_NATIVE_DOUBLE, y.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/z", H5T_NATIVE_DOUBLE, z.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/mu_realized", H5T_NATIVE_DOUBLE, mu.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/mu_max", H5T_NATIVE_DOUBLE, mu_max.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/biomass", H5T_NATIVE_DOUBLE, biomass.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/in_crypt", H5T_NATIVE_INT32, in_crypt.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/n_bi_loci", H5T_NATIVE_INT32, n_bi.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/radius", H5T_NATIVE_DOUBLE, radius.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/lineage_id", H5T_NATIVE_INT64, lineage_id.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/receptor_expr", H5T_NATIVE_DOUBLE,
                   receptor_expr.data(), local_n * NUM_RECEPTORS, cfg_);
  mpi_barrier(cfg_);
#else
  (void)sim; (void)group;
#endif
}

void HDF5Writer::write_grid_layer(const Simulation& sim,
                                   const std::string& group) const {
#ifdef GUTIBM_HDF5
  auto fid = static_cast<hid_t>(file_id_);
  const auto& chem = sim.chemical_field();
  const auto& domain = sim.domain();

  std::array<hsize_t, 3> dims = {static_cast<hsize_t>(nz_),
                                 static_cast<hsize_t>(ny_),
                                 static_cast<hsize_t>(nx_)};
  std::array<hsize_t, 3> chunk = {
      static_cast<hsize_t>(std::min(nz_, 32)),
      static_cast<hsize_t>(std::min(ny_, 32)),
      static_cast<hsize_t>(std::min(nx_, 32))};
  hid_t plist = -1;
  hid_t space = -1;
  std::vector<double> grid3d;
  if (io_rank(cfg_) == 0 && file_id_ >= 0) {
    plist = make_dataset_plist(cfg_, chunk.data(), 3);
    space = H5Screate_simple(3, dims.data(), nullptr);
    grid3d.resize(static_cast<size_t>(nx_ * ny_ * nz_));
  }

  Int local_begin = 0;
  Int local_nx = nx_;
#ifdef GUTIBM_MPI
  int rank = 0;
  int nprocs = 1;
  std::vector<int> counts;
  std::vector<int> displacements;
  std::vector<int> begins;
  std::vector<int> widths;
  int total_count = 0;
  const bool gather_slab = chem.slab_mode() && mpi_multi_rank();
  if (chem.slab_mode()) {
    local_begin = chem.owned_storage_x_begin();
    local_nx = chem.owned_storage_x_end() - local_begin;
  }
  if (gather_slab) {
    rank = mpi_rank_world();
    nprocs = mpi_nprocs_world();
    const int local_count = static_cast<int>(local_nx * ny_ * nz_);
    counts.resize(static_cast<size_t>(nprocs));
    displacements.resize(static_cast<size_t>(nprocs));
    begins.resize(static_cast<size_t>(nprocs));
    widths.resize(static_cast<size_t>(nprocs));
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    const int global_begin = domain.local_grid_x_begin();
    const int global_width = domain.local_grid_nx();
    MPI_Allgather(&global_begin, 1, MPI_INT, begins.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    MPI_Allgather(&global_width, 1, MPI_INT, widths.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    Int gathered_width = 0;
    for (int r = 0; r < nprocs; ++r) {
      displacements[static_cast<size_t>(r)] = total_count;
      total_count += counts[static_cast<size_t>(r)];
      gathered_width += widths[static_cast<size_t>(r)];
    }
    if (gathered_width != nx_) {
      throw HDF5Error("slab grid gather width does not cover global nx");
    }
    if (total_count != nx_ * ny_ * nz_) {
      throw HDF5Error("slab grid gather count does not cover global grid");
    }
  }
#else
  if (chem.slab_mode()) {
    local_begin = chem.owned_storage_x_begin();
    local_nx = chem.owned_storage_x_end() - local_begin;
  }
#endif

  for (Int s = 0; s < chem.num_species(); ++s) {
    const std::string name = chem.spec(s).name;
    if (!should_write_species(name)) continue;

    if (chem.slab_mode()) {
      std::vector<double> local(
          static_cast<size_t>(local_nx * ny_ * nz_));
      for (Int iz = 0; iz < nz_; ++iz) {
        for (Int iy = 0; iy < ny_; ++iy) {
          for (Int ix = 0; ix < local_nx; ++ix) {
            const Int storage_cell =
                (iz * chem.storage_nx() * ny_) + iy * chem.storage_nx()
                + local_begin + ix;
            const size_t packed =
                static_cast<size_t>(iz) * static_cast<size_t>(ny_ * local_nx)
                + static_cast<size_t>(iy * local_nx + ix);
            local[packed] = chem.conc(s, storage_cell);
          }
        }
      }
#ifdef GUTIBM_MPI
      if (gather_slab) {
        const int local_count = static_cast<int>(local.size());
        std::vector<double> gathered;
        if (rank == 0) gathered.resize(static_cast<size_t>(total_count));
        double dummy = 0.0;
        const double* send = local.empty() ? &dummy : local.data();
        MPI_Gatherv(send, local_count, MPI_DOUBLE,
                    rank == 0 ? gathered.data() : nullptr, counts.data(),
                    displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank == 0) {
          std::fill(grid3d.begin(), grid3d.end(), 0.0);
          for (int r = 0; r < nprocs; ++r) {
            const int width = widths[static_cast<size_t>(r)];
            const int source_offset = displacements[static_cast<size_t>(r)];
            for (Int iz = 0; iz < nz_; ++iz) {
              for (Int iy = 0; iy < ny_; ++iy) {
                for (Int ix = 0; ix < width; ++ix) {
                  const size_t source = static_cast<size_t>(source_offset)
                      + static_cast<size_t>(iz * ny_ * width + iy * width + ix);
                  const size_t target = static_cast<size_t>(iz * ny_ * nx_
                      + iy * nx_ + begins[static_cast<size_t>(r)] + ix);
                  grid3d[target] = gathered[source];
                }
              }
            }
          }
        }
      } else {
        for (Int iz = 0; iz < nz_; ++iz) {
          for (Int iy = 0; iy < ny_; ++iy) {
            for (Int ix = 0; ix < local_nx; ++ix) {
              const size_t source =
                  static_cast<size_t>(iz * ny_ * local_nx + iy * local_nx + ix);
              const size_t target = static_cast<size_t>(
                  iz * ny_ * nx_ + iy * nx_ + domain.local_grid_x_begin() + ix);
              grid3d[target] = local[source];
            }
          }
        }
      }
#else
      grid3d = std::move(local);
#endif
    } else if (io_rank(cfg_) == 0) {
      pack_grid_species(chem, domain, s, nx_, ny_, nz_, grid3d);
    }

    if (io_rank(cfg_) == 0 && file_id_ >= 0) {
      const std::string dsname = group + "/" + name;
      hid_t ds = H5Dcreate2(fid, dsname.c_str(), H5T_NATIVE_DOUBLE, space,
                            H5P_DEFAULT, plist, H5P_DEFAULT);
      if (ds < 0) {
        H5Eclear2(H5E_DEFAULT);
        ds = H5Dopen2(fid, dsname.c_str(), H5P_DEFAULT);
      }
      if (ds < 0) {
        std::cerr << "Warning: cannot create/open grid dataset '" << dsname
                  << "'\n";
        continue;
      }
      if (H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   grid3d.data()) < 0) {
        std::cerr << "Warning: H5Dwrite failed for grid dataset '" << dsname
                  << "'\n";
        H5Eclear2(H5E_DEFAULT);
      }
      H5Dclose(ds);
    }
  }

  if (io_rank(cfg_) == 0 && file_id_ >= 0) {
    H5Pclose(plist);
    H5Sclose(space);
  }

  mpi_barrier(cfg_);
#else
  (void)sim; (void)group;
#endif
}

void HDF5Writer::write_lineage_layer(const Simulation& sim,
                                      const std::string& group) const {
#ifdef GUTIBM_HDF5
  auto fid = static_cast<hid_t>(file_id_);
  const auto agents = output_agents(sim);
  const auto n = static_cast<Int>(agents.size());

  std::vector<double> btuB_expr(static_cast<size_t>(n));
  std::vector<double> fepA_expr(static_cast<size_t>(n));
  std::vector<int32_t> n_bi(static_cast<size_t>(n));
  std::vector<int32_t> generation(static_cast<size_t>(n));

  for (Int i = 0; i < n; ++i) {
    const Agent& a = *agents[static_cast<size_t>(i)];
    btuB_expr[static_cast<size_t>(i)] = a.receptor_expr[to_underlying(ReceptorType::BtuB)];
    fepA_expr[static_cast<size_t>(i)] = a.receptor_expr[to_underlying(ReceptorType::FepA)];
    n_bi[static_cast<size_t>(i)] = static_cast<int32_t>(a.genome.bi_loci.size());
    generation[static_cast<size_t>(i)] = static_cast<int32_t>(a.genome.generation);
  }

  const auto local_n = static_cast<hsize_t>(n);
  write_dataset_1d(fid, group + "/btuB_expression", H5T_NATIVE_DOUBLE,
                    btuB_expr.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/fepA_expression", H5T_NATIVE_DOUBLE,
                    fepA_expr.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/num_bi_loci", H5T_NATIVE_INT32,
                    n_bi.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/generation", H5T_NATIVE_INT32,
                    generation.data(), local_n, cfg_);
  mpi_barrier(cfg_);
#else
  (void)sim; (void)group;
#endif
}

void HDF5Writer::write_genome_layer(const Simulation& sim,
                                     const std::string& group) const {
#ifdef GUTIBM_HDF5
  auto fid = static_cast<hid_t>(file_id_);
  const auto agents = output_agents(sim);
  const auto n = static_cast<Int>(agents.size());
  const auto local_n = static_cast<hsize_t>(n);

  std::vector<int64_t> ids(static_cast<size_t>(n));
  std::vector<int64_t> bi_offsets(static_cast<size_t>(n));
  std::vector<int32_t> bi_counts(static_cast<size_t>(n));
  std::vector<int64_t> parent_id(static_cast<size_t>(n));
  std::vector<int32_t> mutations(static_cast<size_t>(n));
  std::vector<int32_t> has_conjugative(static_cast<size_t>(n));
  std::vector<double> plasmid_amel(static_cast<size_t>(n));
  std::vector<int32_t> cdi_type(static_cast<size_t>(n));
  std::vector<int32_t> cdi_immunity(static_cast<size_t>(n));
  std::vector<double> receptor_expr(static_cast<size_t>(n) * NUM_RECEPTORS);
  std::vector<double> toxin_aff(static_cast<size_t>(n) * NUM_RECEPTORS);
  std::vector<double> ligand_aff(static_cast<size_t>(n) * NUM_RECEPTORS);

  std::vector<int32_t> bi_toxin_id;
  std::vector<int32_t> bi_immunity_id;
  std::vector<int32_t> bi_target;
  std::vector<int32_t> bi_bclass;
  std::vector<double> bi_pI;
  std::vector<double> bi_diff;
  std::vector<double> bi_ret;
  std::vector<double> bi_mw;
  std::vector<double> bi_imm_aff;

  for (Int i = 0; i < n; ++i) {
    const Agent& a = *agents[static_cast<size_t>(i)];
    const auto idx = static_cast<size_t>(i);
    ids[idx] = a.identity.tag;
    bi_counts[idx] = static_cast<int32_t>(a.genome.bi_loci.size());
    parent_id[idx] = a.genome.parent_id;
    mutations[idx] = static_cast<int32_t>(a.genome.mutations);
    has_conjugative[idx] = a.genome.has_conjugative_plasmid ? 1 : 0;
    plasmid_amel[idx] = a.genome.plasmid_cost_amelioration;
    cdi_type[idx] = static_cast<int32_t>(a.genome.cdi_type);
    cdi_immunity[idx] = static_cast<int32_t>(a.genome.cdi_immunity);
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
      bi_target.push_back(static_cast<int32_t>(to_underlying(bi.target)));
      bi_bclass.push_back(static_cast<int32_t>(to_underlying(bi.bclass)));
      bi_pI.push_back(bi.pI);
      bi_diff.push_back(bi.diff_coeff);
      bi_ret.push_back(bi.retardation);
      bi_mw.push_back(bi.molecular_weight);
      bi_imm_aff.push_back(bi.immunity_binding_affinity);
    }
  }

  const auto local_bi = static_cast<hsize_t>(bi_toxin_id.size());
  const ParSlice bi_slice = compute_par_slice(local_bi, cfg_);
  hsize_t local_offset = bi_slice.offset;
  for (size_t i = 0; i < bi_offsets.size(); ++i) {
    bi_offsets[i] = static_cast<int64_t>(local_offset);
    local_offset += static_cast<hsize_t>(bi_counts[i]);
  }

  write_dataset_1d(fid, group + "/id", H5T_NATIVE_INT64,
                   ids.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/bi_offset", H5T_NATIVE_INT64,
                   bi_offsets.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/bi_count", H5T_NATIVE_INT32,
                   bi_counts.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/parent_id", H5T_NATIVE_INT64,
                   parent_id.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/mutations", H5T_NATIVE_INT32,
                   mutations.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/has_conjugative_plasmid", H5T_NATIVE_INT32,
                   has_conjugative.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/plasmid_cost_amelioration", H5T_NATIVE_DOUBLE,
                   plasmid_amel.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/cdi_type", H5T_NATIVE_INT32,
                   cdi_type.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/cdi_immunity", H5T_NATIVE_INT32,
                   cdi_immunity.data(), local_n, cfg_);
  write_dataset_1d(fid, group + "/receptor_expression", H5T_NATIVE_DOUBLE,
                   receptor_expr.data(), local_n * NUM_RECEPTORS, cfg_);
  write_dataset_1d(fid, group + "/toxin_affinity", H5T_NATIVE_DOUBLE,
                   toxin_aff.data(), local_n * NUM_RECEPTORS, cfg_);
  write_dataset_1d(fid, group + "/ligand_affinity", H5T_NATIVE_DOUBLE,
                   ligand_aff.data(), local_n * NUM_RECEPTORS, cfg_);

  write_dataset_1d(fid, group + "/bi_toxin_id", H5T_NATIVE_INT32,
                   bi_toxin_id.data(), local_bi, cfg_);
  write_dataset_1d(fid, group + "/bi_immunity_id", H5T_NATIVE_INT32,
                   bi_immunity_id.data(), local_bi, cfg_);
  write_dataset_1d(fid, group + "/bi_target", H5T_NATIVE_INT32,
                   bi_target.data(), local_bi, cfg_);
  write_dataset_1d(fid, group + "/bi_bclass", H5T_NATIVE_INT32,
                   bi_bclass.data(), local_bi, cfg_);
  write_dataset_1d(fid, group + "/bi_pI", H5T_NATIVE_DOUBLE,
                   bi_pI.data(), local_bi, cfg_);
  write_dataset_1d(fid, group + "/bi_diff_coeff", H5T_NATIVE_DOUBLE,
                   bi_diff.data(), local_bi, cfg_);
  write_dataset_1d(fid, group + "/bi_retardation", H5T_NATIVE_DOUBLE,
                   bi_ret.data(), local_bi, cfg_);
  write_dataset_1d(fid, group + "/bi_molecular_weight", H5T_NATIVE_DOUBLE,
                   bi_mw.data(), local_bi, cfg_);
  write_dataset_1d(fid, group + "/bi_immunity_binding_affinity", H5T_NATIVE_DOUBLE,
                   bi_imm_aff.data(), local_bi, cfg_);
  mpi_barrier(cfg_);
#else
  (void)sim; (void)group;
#endif
}

std::vector<const Agent*> HDF5Writer::output_agents(const Simulation& sim) const {
  std::vector<const Agent*> output;
  output.reserve(sim.agents().size());
  for (const Agent& agent : sim.agents()) {
    if (!include_dead_agents_ && agent.state == PhenoState::DEAD) continue;
    output.push_back(&agent);
  }
  return output;
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

bool HDF5Writer::write_closed_restart(Simulation& sim, const std::string& path,
                                      Int step, Real time, Real dt,
                                      bool preserve_event_counters) {
#ifndef GUTIBM_HDF5
  (void)sim;
  (void)path;
  (void)step;
  (void)time;
  (void)dt;
  (void)preserve_event_counters;
  return false;
#else
  namespace fs = std::filesystem;
  try {
    validate_path_syntax(path);
  } catch (const IOError& ex) {
    std::cerr << "Warning: invalid restart path '" << path << "': " << ex.what()
              << "\n";
    return false;
  }

  const fs::path out(path);
  const fs::path tmp = out.string() + ".tmp";
  const bool writer_rank = io_rank(HDF5Config{}) == 0;
  bool preparation_ok = true;
  if (writer_rank && out.has_parent_path()) {
    std::error_code ec;
    fs::create_directories(out.parent_path(), ec);
    if (ec) {
      std::cerr << "Warning: cannot create restart directory '"
                << out.parent_path().string() << "': " << ec.message() << "\n";
      preparation_ok = false;
    }
  }

  // Stage-3 grids are ~4.4 GB uncompressed (50M cells × 11 species). Gzip keeps
  // Spot artifacts small enough for container scratch disks and S3 upload.
  const auto nx = static_cast<std::uint64_t>(sim.domain().nx());
  const auto ny = static_cast<std::uint64_t>(sim.domain().ny());
  const auto nz = static_cast<std::uint64_t>(sim.domain().nz());
  const auto nspecies =
      static_cast<std::uint64_t>(std::max(sim.chemical_field().num_species(), 1));
  const std::uint64_t uncompressed_est =
      nx * ny * nz * nspecies * sizeof(double) + (64ULL << 20);
  {
    std::error_code space_ec;
    const fs::path space_root =
        out.has_parent_path() ? out.parent_path() : fs::current_path();
    const fs::space_info si = fs::space(space_root, space_ec);
    // Gzip shrinks sparse toxin fields a lot; still keep multi-GiB headroom for
    // Stage-3 chemistry dumps on small container scratch volumes.
    constexpr std::uint64_t kHeadroom = 1ULL << 30;
    const std::uint64_t need =
        std::max<std::uint64_t>((uncompressed_est / 4) + kHeadroom, 2ULL << 30);
    if (writer_rank && !space_ec && si.available < need) {
      std::cerr << "Warning: refusing restart write to '" << path
                << "': only " << (si.available >> 20)
                << " MiB free, need ~" << (need >> 20) << " MiB\n";
      preparation_ok = false;
    }
  }
#ifdef GUTIBM_MPI
  if (mpi_is_active()) {
    MPI_Bcast(&preparation_ok, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  }
  if (!preparation_ok) return false;
#else
  if (!preparation_ok) return false;
#endif

  HDF5Config cfg;
  cfg.filename = tmp.string();
  cfg.enabled = true;
  cfg.compression = "gzip";
  cfg.compression_level = 4;
  cfg.schedule.summary = 1;
  cfg.schedule.agents = 1;
  cfg.schedule.grid = 1;
  cfg.schedule.lineage = 1;
  cfg.schedule.genome = 1;
  cfg.schedule.provenance = 1;
  cfg.schedule.grid_species = {"all"};

  HDF5Writer writer;
  writer.init(cfg, sim.domain());
  writer.include_dead_agents_ = true;
  if (!writer.is_enabled()) {
    std::error_code ec;
    fs::remove(tmp, ec);
    return false;
  }
  const StepEvents saved_step_events = sim.step_events();
  const StepEvents saved_summary_events = sim.summary_events();
  const StepEvents saved_cumulative_events = sim.cumulative_events();
  const Int saved_window_start_step = sim.event_window_start_step();
  const Real saved_window_start_time = sim.event_window_start_time();
  const auto saved_provenance = sim.kill_provenance();
  const auto saved_flux_accounting = sim.chemical_field().flux_accounting();
  sim.materialize_bacteriocin_fields_for_output();
  writer.write_step(sim, step, time, dt);
  if (preserve_event_counters) {
    sim.step_events() = saved_step_events;
    sim.summary_events() = saved_summary_events;
    sim.cumulative_events() = saved_cumulative_events;
    sim.set_event_window_start(saved_window_start_step, saved_window_start_time);
  }
  sim.kill_provenance() = saved_provenance;
  sim.chemical_field().flux_accounting() = saved_flux_accounting;
  writer.finalize();

  bool published = true;
  if (writer_rank) {
    std::error_code sz_ec;
    if (const auto tmp_bytes = fs::file_size(tmp, sz_ec);
        sz_ec || tmp_bytes < 4096 || H5Fis_hdf5(tmp.string().c_str()) <= 0) {
      std::cerr << "Warning: restart tmp '" << tmp.string()
                << "' is missing/unreadable after write (size="
                << (sz_ec ? 0 : tmp_bytes) << ")\n";
      fs::remove(tmp, sz_ec);
      published = false;
    } else {
      std::error_code rename_ec;
      fs::rename(tmp, out, rename_ec);
      if (rename_ec) {
        std::cerr << "Warning: failed to publish restart '" << path
                  << "': " << rename_ec.message() << "\n";
        fs::remove(tmp, rename_ec);
        published = false;
      }
    }
  }
#ifdef GUTIBM_MPI
  if (mpi_is_active()) {
    MPI_Bcast(&published, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    mpi_barrier(cfg);
  }
#endif
  return published;
#endif
}

}  // namespace gutibm
