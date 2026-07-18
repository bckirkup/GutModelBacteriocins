/* -----------------------------------------------------------------------
   GutIBM – HDF5 reader implementation (Spec 4 layered schema)
   ----------------------------------------------------------------------- */

#include "hdf5_reader.h"
#include "path_utils.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <string_view>
#include "error.h"

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}
#endif

namespace gutibm {

namespace {

#ifdef GUTIBM_HDF5

bool link_exists(hid_t loc, const std::string& path) {
  return H5Lexists(loc, path.c_str(), H5P_DEFAULT) > 0;
}

template <typename T>
std::vector<T> read_dataset_1d(hid_t file, const std::string& path, hid_t h5_type) {
  if (!link_exists(file, path)) {
    throw HDF5Error("missing HDF5 dataset: " + path);
  }

  hid_t dset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
  if (dset < 0) {
    throw HDF5Error("failed to open HDF5 dataset: " + path);
  }

  hid_t space = H5Dget_space(dset);
  if (int ndims = H5Sget_simple_extent_ndims(space); ndims != 1) {
    H5Sclose(space);
    H5Dclose(dset);
    throw HDF5Error("expected 1D dataset: " + path);
  }

  hsize_t len = 0;
  H5Sget_simple_extent_dims(space, &len, nullptr);

  std::vector<T> data(static_cast<size_t>(len));
  if (len > 0) {
    herr_t status = H5Dread(dset, h5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    if (status < 0) {
      H5Sclose(space);
      H5Dclose(dset);
      throw HDF5Error("failed to read HDF5 dataset: " + path);
    }
  }

  H5Sclose(space);
  H5Dclose(dset);
  return data;
}

template <typename T>
T read_scalar(hid_t file, const std::string& path, hid_t h5_type) {
  auto data = read_dataset_1d<T>(file, path, h5_type);
  if (data.size() != 1) {
    throw HDF5Error("expected scalar dataset: " + path);
  }
  return data[0];
}

struct StepListData {
  std::vector<std::string>* steps;
};

herr_t collect_step_link(hid_t /*group*/, const char* name, const H5L_info_t* /*info*/,
                         void* op_data) {
  auto* data = static_cast<StepListData*>(op_data);
  if (std::string link(name); link.rfind("step_", 0) == 0) {
    data->steps->push_back(std::move(link));
  }
  return 0;
}

std::vector<std::string> list_step_groups(hid_t file, const char* layer) {
  StepListData data;
  std::vector<std::string> steps;
  data.steps = &steps;
  if (!link_exists(file, layer)) {
    return steps;
  }
  hid_t grp = H5Gopen2(file, layer, H5P_DEFAULT);
  if (grp < 0) return steps;
  H5Literate(grp, H5_INDEX_NAME, H5_ITER_INC, nullptr, collect_step_link, &data);
  H5Gclose(grp);
  std::ranges::sort(steps, [](std::string_view sa, std::string_view sb) {
              const auto na = std::stoi(std::string(sa.substr(5)));
              const auto nb = std::stoi(std::string(sb.substr(5)));
              return na < nb;
            });
  return steps;
}

int32_t read_file_attr_int(hid_t file, const char* name) {
  if (H5Aexists(file, name) <= 0) return 0;
  hid_t attr = H5Aopen(file, name, H5P_DEFAULT);
  int32_t value = 0;
  H5Aread(attr, H5T_NATIVE_INT32, &value);
  H5Aclose(attr);
  return value;
}

HDF5CheckpointAgents read_agents(hid_t file, const std::string& step) {
  const std::string prefix = "agents/" + step + "/";
  HDF5CheckpointAgents out;
  out.id       = read_dataset_1d<int64_t>(file, prefix + "id",      H5T_NATIVE_INT64);
  out.type     = read_dataset_1d<int32_t>(file, prefix + "type",    H5T_NATIVE_INT32);
  out.state    = read_dataset_1d<int32_t>(file, prefix + "state",   H5T_NATIVE_INT32);
  out.x        = read_dataset_1d<double>(file, prefix + "x",       H5T_NATIVE_DOUBLE);
  out.y        = read_dataset_1d<double>(file, prefix + "y",       H5T_NATIVE_DOUBLE);
  out.z        = read_dataset_1d<double>(file, prefix + "z",       H5T_NATIVE_DOUBLE);
  out.radius   = read_dataset_1d<double>(file, prefix + "radius",   H5T_NATIVE_DOUBLE);
  out.biomass  = read_dataset_1d<double>(file, prefix + "biomass", H5T_NATIVE_DOUBLE);

  const std::string mu_path = link_exists(file, prefix + "mu_realized")
      ? prefix + "mu_realized" : prefix + "mu";
  out.mu = read_dataset_1d<double>(file, mu_path, H5T_NATIVE_DOUBLE);

  const std::string lineage_path = link_exists(file, prefix + "lineage_id")
      ? prefix + "lineage_id" : prefix + "lineage";
  out.lineage = read_dataset_1d<int64_t>(file, lineage_path, H5T_NATIVE_INT64);

  if (const auto n = out.id.size();
      out.type.size() != n || out.state.size() != n || out.x.size() != n ||
      out.y.size() != n || out.z.size() != n || out.radius.size() != n ||
      out.biomass.size() != n || out.mu.size() != n || out.lineage.size() != n) {
    throw HDF5Error("inconsistent agent dataset lengths in agents/" + step);
  }
  return out;
}

HDF5CheckpointLineage read_lineage(hid_t file, const std::string& step) {
  const std::string prefix = "lineage/" + step + "/";
  HDF5CheckpointLineage out;
  out.btuB_expression = read_dataset_1d<double>(file, prefix + "btuB_expression",
                                                H5T_NATIVE_DOUBLE);
  out.fepA_expression = read_dataset_1d<double>(file, prefix + "fepA_expression",
                                                H5T_NATIVE_DOUBLE);
  out.num_bi_loci     = read_dataset_1d<int32_t>(file, prefix + "num_bi_loci",
                                                 H5T_NATIVE_INT32);
  out.generation      = read_dataset_1d<int32_t>(file, prefix + "generation",
                                                 H5T_NATIVE_INT32);

  if (const auto n = out.btuB_expression.size();
      out.fepA_expression.size() != n || out.num_bi_loci.size() != n ||
      out.generation.size() != n) {
    throw HDF5Error("inconsistent lineage dataset lengths in lineage/" + step);
  }
  return out;
}

HDF5CheckpointGenome read_genome(hid_t file, const std::string& step) {
  HDF5CheckpointGenome out;
  const std::string prefix = "genome/" + step + "/";
  if (!link_exists(file, "genome/" + step)) {
    return out;
  }

  out.present = true;
  out.parent_id = read_dataset_1d<int64_t>(file, prefix + "parent_id", H5T_NATIVE_INT64);
  out.mutations = read_dataset_1d<int32_t>(file, prefix + "mutations", H5T_NATIVE_INT32);
  out.has_conjugative_plasmid =
      read_dataset_1d<int32_t>(file, prefix + "has_conjugative_plasmid", H5T_NATIVE_INT32);
  out.plasmid_cost_amelioration =
      read_dataset_1d<double>(file, prefix + "plasmid_cost_amelioration", H5T_NATIVE_DOUBLE);
  if (link_exists(file, prefix + "cdi_type")) {
    out.cdi_type = read_dataset_1d<int32_t>(file, prefix + "cdi_type", H5T_NATIVE_INT32);
    out.cdi_immunity =
        read_dataset_1d<int32_t>(file, prefix + "cdi_immunity", H5T_NATIVE_INT32);
  }
  out.receptor_expression =
      read_dataset_1d<double>(file, prefix + "receptor_expression", H5T_NATIVE_DOUBLE);
  out.toxin_affinity =
      read_dataset_1d<double>(file, prefix + "toxin_affinity", H5T_NATIVE_DOUBLE);
  out.ligand_affinity =
      read_dataset_1d<double>(file, prefix + "ligand_affinity", H5T_NATIVE_DOUBLE);

  out.bi_toxin_id =
      read_dataset_1d<int32_t>(file, prefix + "bi_toxin_id", H5T_NATIVE_INT32);
  out.bi_immunity_id =
      read_dataset_1d<int32_t>(file, prefix + "bi_immunity_id", H5T_NATIVE_INT32);
  out.bi_target =
      read_dataset_1d<int32_t>(file, prefix + "bi_target", H5T_NATIVE_INT32);
  out.bi_bclass =
      read_dataset_1d<int32_t>(file, prefix + "bi_bclass", H5T_NATIVE_INT32);
  out.bi_pI =
      read_dataset_1d<double>(file, prefix + "bi_pI", H5T_NATIVE_DOUBLE);
  out.bi_diff_coeff =
      read_dataset_1d<double>(file, prefix + "bi_diff_coeff", H5T_NATIVE_DOUBLE);
  out.bi_retardation =
      read_dataset_1d<double>(file, prefix + "bi_retardation", H5T_NATIVE_DOUBLE);
  out.bi_molecular_weight =
      read_dataset_1d<double>(file, prefix + "bi_molecular_weight", H5T_NATIVE_DOUBLE);
  out.bi_immunity_binding_affinity =
      read_dataset_1d<double>(file, prefix + "bi_immunity_binding_affinity",
                              H5T_NATIVE_DOUBLE);

  if (const auto n = out.parent_id.size();
      out.mutations.size() != n || out.has_conjugative_plasmid.size() != n ||
      out.plasmid_cost_amelioration.size() != n ||
      out.receptor_expression.size() != n * NUM_RECEPTORS ||
      out.toxin_affinity.size() != n * NUM_RECEPTORS ||
      out.ligand_affinity.size() != n * NUM_RECEPTORS) {
    throw HDF5Error("inconsistent genome per-agent dataset lengths in genome/" + step);
  }

  if (const auto n_bi = out.bi_toxin_id.size();
      out.bi_immunity_id.size() != n_bi || out.bi_target.size() != n_bi ||
      out.bi_bclass.size() != n_bi || out.bi_pI.size() != n_bi ||
      out.bi_diff_coeff.size() != n_bi || out.bi_retardation.size() != n_bi ||
      out.bi_molecular_weight.size() != n_bi ||
      out.bi_immunity_binding_affinity.size() != n_bi) {
    throw HDF5Error("inconsistent BI locus dataset lengths in genome/" + step);
  }

  return out;
}

HDF5CheckpointMetadata read_metadata(hid_t file, const std::string& step) {
  const std::string prefix = "summary/" + step + "/";
  HDF5CheckpointMetadata meta;
  meta.time = read_scalar<double>(file, prefix + "time", H5T_NATIVE_DOUBLE);
  meta.step = read_scalar<int32_t>(file, prefix + "step", H5T_NATIVE_INT32);
  if (link_exists(file, prefix + "num_agents")) {
    meta.num_agents = read_scalar<int32_t>(file, prefix + "num_agents", H5T_NATIVE_INT32);
  } else {
    meta.num_agents = read_scalar<int32_t>(file, prefix + "n_total", H5T_NATIVE_INT32);
  }
  if (link_exists(file, prefix + "num_lineages")) {
    meta.num_lineages = read_scalar<int32_t>(file, prefix + "num_lineages", H5T_NATIVE_INT32);
  }
  return meta;
}

struct LinkNameCollector {
  std::vector<std::string>* names;
};

herr_t collect_link_name(hid_t /*group*/, const char* name, const H5L_info_t* /*info*/,
                         void* op_data) {
  static_cast<LinkNameCollector*>(op_data)->names->emplace_back(name);
  return 0;
}

HDF5CheckpointGrid read_grid(hid_t file, const std::string& step) {
  const std::string grid_path = "grid/" + step;
  HDF5CheckpointGrid grid;
  if (!link_exists(file, grid_path)) {
    return grid;
  }

  const int32_t nx = read_file_attr_int(file, "nx");
  const int32_t ny = read_file_attr_int(file, "ny");
  const int32_t nz = read_file_attr_int(file, "nz");
  const size_t flat_cells = static_cast<size_t>(nx) * static_cast<size_t>(ny)
      * static_cast<size_t>(nz);

  hid_t grid_group = H5Gopen2(file, grid_path.c_str(), H5P_DEFAULT);
  if (grid_group < 0) {
    throw HDF5Error("failed to open grid group: " + grid_path);
  }

  std::vector<std::string> datasets;
  LinkNameCollector collector{&datasets};
  H5Literate(grid_group, H5_INDEX_NAME, H5_ITER_INC, nullptr, collect_link_name, &collector);
  H5Gclose(grid_group);

  for (const auto& ds_name : datasets) {
    const std::string path = grid_path + "/" + ds_name;
    hid_t dset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
    if (dset < 0) continue;

    hid_t space = H5Dget_space(dset);
    const int ndims = H5Sget_simple_extent_ndims(space);
    std::array<hsize_t, 3> dims{0, 0, 0};
    H5Sget_simple_extent_dims(space, dims.data(), nullptr);

    std::vector<double> values;
    if (ndims == 3) {
      const size_t total = static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1])
          * static_cast<size_t>(dims[2]);
      values.resize(total);
      H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
      if (flat_cells > 0 && values.size() != flat_cells) {
        H5Sclose(space);
        H5Dclose(dset);
        throw HDF5Error("grid dataset size mismatch for species: " + ds_name);
      }
    } else if (ndims == 1) {
      values = read_dataset_1d<double>(file, path, H5T_NATIVE_DOUBLE);
    } else {
      H5Sclose(space);
      H5Dclose(dset);
      throw HDF5Error("unsupported grid dataset rank for species: " + ds_name);
    }

    H5Sclose(space);
    H5Dclose(dset);
    grid.species[ds_name] = std::move(values);
  }
  return grid;
}

#endif  // GUTIBM_HDF5

}  // namespace

HDF5Reader::~HDF5Reader() {
  close();
}

bool HDF5Reader::open(const std::string& filename) {
  close();

  try {
    filename_ = validate_input_file_path(filename);
  } catch (const IOError& ) {
    open_ = false;
    file_id_ = -1;
    return false;
  }

#ifdef GUTIBM_HDF5
  H5Eclear2(H5E_DEFAULT);

  hid_t fid = H5Fopen(filename_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fid < 0) {
    open_ = false;
    file_id_ = -1;
    return false;
  }
  file_id_ = static_cast<int64_t>(fid);
  open_ = true;
  return true;
#else
  (void)filename;
  return false;
#endif
}

void HDF5Reader::close() {
#ifdef GUTIBM_HDF5
  if (open_ && file_id_ >= 0) {
    H5Fclose(static_cast<hid_t>(file_id_));
    file_id_ = -1;
  }
#endif
  open_ = false;
  filename_.clear();
}

std::vector<std::string> HDF5Reader::list_steps() const {
#ifdef GUTIBM_HDF5
  if (!open_) return {};
  auto steps = list_step_groups(static_cast<hid_t>(file_id_), "summary");
  if (steps.empty()) {
    steps = list_step_groups(static_cast<hid_t>(file_id_), "agents");
  }
  return steps;
#else
  return {};
#endif
}

std::string HDF5Reader::latest_step() const {
  auto steps = list_steps();
  if (steps.empty()) return {};
  return steps.back();
}

HDF5CheckpointSnapshot HDF5Reader::load_step(const std::string& step_name) const {
#ifdef GUTIBM_HDF5
  if (!open_) {
    throw HDF5Error("HDF5Reader is not open");
  }

  const bool has_agents = link_exists(static_cast<hid_t>(file_id_), "agents/" + step_name);
  const bool has_summary = link_exists(static_cast<hid_t>(file_id_), "summary/" + step_name);
  if (!has_agents && !has_summary) {
    throw HDF5Error("missing step group: " + step_name);
  }

  HDF5CheckpointSnapshot snap;
  snap.step_name = step_name;
  snap.metadata = read_metadata(static_cast<hid_t>(file_id_), step_name);
  snap.agents   = read_agents(static_cast<hid_t>(file_id_), step_name);
  snap.lineage  = read_lineage(static_cast<hid_t>(file_id_), step_name);
  snap.genome   = read_genome(static_cast<hid_t>(file_id_), step_name);
  snap.grid     = read_grid(static_cast<hid_t>(file_id_), step_name);

  if (snap.metadata.num_agents != static_cast<Int>(snap.agents.id.size())) {
    throw HDF5Error("metadata num_agents does not match agents dataset");
  }
  if (snap.metadata.num_agents != static_cast<Int>(snap.lineage.generation.size())) {
    throw HDF5Error("metadata num_agents does not match lineage dataset");
  }
  return snap;
#else
  (void)step_name;
  throw HDF5Error("HDF5 support not compiled in");
#endif
}

HDF5CheckpointSnapshot HDF5Reader::load_snapshot(const std::string& filename,
                                                 const std::string& step_name) {
  HDF5Reader reader;
  if (!reader.open(filename)) {
    throw HDF5Error("failed to open HDF5 file: " + filename);
  }

  std::string step = step_name;
  if (step.empty()) {
    step = reader.latest_step();
    if (step.empty()) {
      throw HDF5Error("no step groups found in: " + filename);
    }
  }
  return reader.load_step(step);
}

}  // namespace gutibm
