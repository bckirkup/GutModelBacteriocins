/* -----------------------------------------------------------------------
   GutIBM – HDF5 reader implementation (Spec 4 layered schema)
   ----------------------------------------------------------------------- */

#include "hdf5_reader.h"
#include "path_utils.h"
#include "step_events.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <string_view>
#include "error.h"
#include "h5_util.h"
#include "species_names.h"

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

double read_file_attr_real(hid_t file, const char* name) {
  if (H5Aexists(file, name) <= 0) return 0.0;
  hid_t attr = H5Aopen(file, name, H5P_DEFAULT);
  double value = 0.0;
  H5Aread(attr, H5T_NATIVE_DOUBLE, &value);
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

  if (link_exists(file, prefix + "mu_max")) {
    out.mu_max = read_dataset_1d<double>(file, prefix + "mu_max", H5T_NATIVE_DOUBLE);
  }
  if (link_exists(file, prefix + "realized_fermentation_fraction")) {
    out.realized_fermentation_fraction = read_dataset_1d<double>(
        file, prefix + "realized_fermentation_fraction", H5T_NATIVE_DOUBLE);
  }
  if (link_exists(file, prefix + "respired_oxygen_rate")) {
    out.respired_oxygen_rate = read_dataset_1d<double>(
        file, prefix + "respired_oxygen_rate", H5T_NATIVE_DOUBLE);
  }
  if (link_exists(file, prefix + "in_crypt")) {
    out.in_crypt =
        read_dataset_1d<int32_t>(file, prefix + "in_crypt", H5T_NATIVE_INT32);
  }
  if (link_exists(file, prefix + "microcin_penalty_applied")) {
    out.microcin_penalty_applied = read_dataset_1d<int32_t>(
        file, prefix + "microcin_penalty_applied", H5T_NATIVE_INT32);
  }
  const auto read_opt_double = [&](const std::string& name, std::vector<double>& dst) {
    if (link_exists(file, prefix + name)) {
      dst = read_dataset_1d<double>(file, prefix + name, H5T_NATIVE_DOUBLE);
    }
  };
  const auto read_opt_i32 = [&](const std::string& name, std::vector<int32_t>& dst) {
    if (link_exists(file, prefix + name)) {
      dst = read_dataset_1d<int32_t>(file, prefix + name, H5T_NATIVE_INT32);
    }
  };
  read_opt_double("mass", out.mass);
  read_opt_double("outer_radius", out.outer_radius);
  read_opt_double("age", out.age);
  read_opt_double("sos_timer", out.sos_timer);
  read_opt_double("death_time", out.death_time);
  read_opt_double("pending_growth_chemistry_biomass",
                  out.pending_growth_chemistry_biomass);
  read_opt_double("receptor_expr", out.receptor_expr);
  read_opt_double("receptor_expr_base", out.receptor_expr_base);
  read_opt_double("motility_dir_x", out.motility_dir_x);
  read_opt_double("motility_dir_y", out.motility_dir_y);
  read_opt_double("motility_dir_z", out.motility_dir_z);
  read_opt_double("motility_swim_speed", out.motility_swim_speed);
  read_opt_double("motility_run_timer", out.motility_run_timer);
  read_opt_i32("motility_is_stopped", out.motility_is_stopped);
  read_opt_double("motility_stop_timer", out.motility_stop_timer);
  read_opt_double("motility_prev_carbon", out.motility_prev_carbon);
  read_opt_double("motility_prev_oxygen", out.motility_prev_oxygen);
  read_opt_double("motility_prev_ai2", out.motility_prev_ai2);
  read_opt_double("sampled_nuclease_conc", out.sampled_nuclease_conc);

  const std::string lineage_path = link_exists(file, prefix + "lineage_id")
      ? prefix + "lineage_id" : prefix + "lineage";
  out.lineage = read_dataset_1d<int64_t>(file, lineage_path, H5T_NATIVE_INT64);

  if (const auto n = out.id.size();
      out.type.size() != n || out.state.size() != n || out.x.size() != n ||
      out.y.size() != n || out.z.size() != n || out.radius.size() != n ||
      out.biomass.size() != n || out.mu.size() != n || out.lineage.size() != n) {
    throw HDF5Error("inconsistent agent dataset lengths in agents/" + step);
  }
  const auto check_opt_len = [&](size_t sz, size_t expected, const char* label) {
    if (sz != 0 && sz != expected) {
      throw HDF5Error(std::string("inconsistent ") + label + " length in agents/" + step);
    }
  };
  const auto n = out.id.size();
  check_opt_len(out.mu_max.size(), n, "mu_max");
  check_opt_len(out.realized_fermentation_fraction.size(), n, "realized fermentation");
  check_opt_len(out.respired_oxygen_rate.size(), n, "respired oxygen rate");
  check_opt_len(out.in_crypt.size(), n, "in_crypt");
  check_opt_len(out.microcin_penalty_applied.size(), n, "microcin penalty");
  check_opt_len(out.mass.size(), n, "mass");
  check_opt_len(out.outer_radius.size(), n, "outer_radius");
  check_opt_len(out.age.size(), n, "age");
  check_opt_len(out.sos_timer.size(), n, "sos_timer");
  check_opt_len(out.death_time.size(), n, "death_time");
  check_opt_len(out.pending_growth_chemistry_biomass.size(), n,
                "pending_growth_chemistry_biomass");
  check_opt_len(out.receptor_expr.size(), n * NUM_RECEPTORS, "receptor_expr");
  check_opt_len(out.receptor_expr_base.size(), n * NUM_RECEPTORS, "receptor_expr_base");
  check_opt_len(out.motility_dir_x.size(), n, "motility_dir_x");
  check_opt_len(out.motility_dir_y.size(), n, "motility_dir_y");
  check_opt_len(out.motility_dir_z.size(), n, "motility_dir_z");
  check_opt_len(out.motility_swim_speed.size(), n, "motility_swim_speed");
  check_opt_len(out.motility_run_timer.size(), n, "motility_run_timer");
  check_opt_len(out.motility_is_stopped.size(), n, "motility_is_stopped");
  check_opt_len(out.motility_stop_timer.size(), n, "motility_stop_timer");
  check_opt_len(out.motility_prev_carbon.size(), n, "motility_prev_carbon");
  check_opt_len(out.motility_prev_oxygen.size(), n, "motility_prev_oxygen");
  check_opt_len(out.motility_prev_ai2.size(), n, "motility_prev_ai2");
  check_opt_len(out.sampled_nuclease_conc.size(),
                n * species::BACTERIOCIN_RECEPTOR_TARGETS.size(),
                "sampled_nuclease_conc");
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
  if (link_exists(file, prefix + "bi_release_mode")) {
    out.bi_release_mode =
        read_dataset_1d<int32_t>(file, prefix + "bi_release_mode",
                                 H5T_NATIVE_INT32);
  }
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
  if (link_exists(file, prefix + "bi_protease_half_life")) {
    out.bi_protease_half_life = read_dataset_1d<double>(
        file, prefix + "bi_protease_half_life", H5T_NATIVE_DOUBLE);
  }
  if (link_exists(file, prefix + "bi_is_nuclease")) {
    out.bi_is_nuclease = read_dataset_1d<int32_t>(
        file, prefix + "bi_is_nuclease", H5T_NATIVE_INT32);
  }
  if (link_exists(file, prefix + "bi_burst_size")) {
    out.bi_burst_size = read_dataset_1d<double>(
        file, prefix + "bi_burst_size", H5T_NATIVE_DOUBLE);
  }
  if (link_exists(file, prefix + "bi_phage_induction_rate")) {
    out.bi_phage_induction_rate = read_dataset_1d<double>(
        file, prefix + "bi_phage_induction_rate", H5T_NATIVE_DOUBLE);
  }
  if (link_exists(file, prefix + "bi_phage_burst_size")) {
    out.bi_phage_burst_size = read_dataset_1d<double>(
        file, prefix + "bi_phage_burst_size", H5T_NATIVE_DOUBLE);
  }
  if (link_exists(file, prefix + "bi_phage_lysogeny_rate")) {
    out.bi_phage_lysogeny_rate = read_dataset_1d<double>(
        file, prefix + "bi_phage_lysogeny_rate", H5T_NATIVE_DOUBLE);
  }

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
      out.bi_immunity_binding_affinity.size() != n_bi ||
      (!out.bi_release_mode.empty() && out.bi_release_mode.size() != n_bi) ||
      (!out.bi_protease_half_life.empty() && out.bi_protease_half_life.size() != n_bi) ||
      (!out.bi_is_nuclease.empty() && out.bi_is_nuclease.size() != n_bi) ||
      (!out.bi_burst_size.empty() && out.bi_burst_size.size() != n_bi) ||
      (!out.bi_phage_induction_rate.empty() && out.bi_phage_induction_rate.size() != n_bi) ||
      (!out.bi_phage_burst_size.empty() && out.bi_phage_burst_size.size() != n_bi) ||
      (!out.bi_phage_lysogeny_rate.empty() && out.bi_phage_lysogeny_rate.size() != n_bi)) {
    throw HDF5Error("inconsistent BI locus dataset lengths in genome/" + step);
  }

  return out;
}

HDF5CheckpointMetadata read_metadata(hid_t file, const std::string& step) {
  const std::string prefix = "summary/" + step + "/";
  HDF5CheckpointMetadata meta;
  const double legacy_dx = read_file_attr_real(file, "grid_dx");
  meta.grid_spacing = {
      read_file_attr_real(file, "grid_dx_x"),
      read_file_attr_real(file, "grid_dx_y"),
      read_file_attr_real(file, "grid_dx_z")};
  if (meta.grid_spacing[0] <= 0.0) meta.grid_spacing[0] = legacy_dx;
  if (meta.grid_spacing[1] <= 0.0) meta.grid_spacing[1] = legacy_dx;
  if (meta.grid_spacing[2] <= 0.0) meta.grid_spacing[2] = legacy_dx;
  meta.has_grid_spacing = meta.grid_spacing[0] > 0.0
      && meta.grid_spacing[1] > 0.0 && meta.grid_spacing[2] > 0.0;
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
  const std::string events = prefix + "events/";
  const auto read_event = [&](const std::string& name, Int& value) {
    if (link_exists(file, events + name)) {
      value = read_scalar<int32_t>(file, events + name, H5T_NATIVE_INT32);
    }
  };
  read_event("sos_inductions", meta.interval_events.sos_inductions);
  read_event("phage_inductions", meta.interval_events.phage_inductions);
  read_event("mortality_colicin", meta.interval_events.mortality_colicin);
  read_event("mortality_cdi", meta.interval_events.mortality_cdi);
  read_event("outflow_washout", meta.interval_events.outflow_washout);
  read_event("outflow_boundary", meta.interval_events.outflow_boundary);
  read_event("mortality_lysis", meta.interval_events.mortality_lysis);
  read_event("divisions", meta.interval_events.divisions);
  if (link_exists(file, events + "divisions_by_type")) {
    auto values = read_dataset_1d<int32_t>(file, events + "divisions_by_type",
                                           H5T_NATIVE_INT32);
    for (size_t i = 0; i < values.size() && i < MAX_AGENT_TYPES; ++i) {
      meta.interval_events.divisions_by_type[i] = values[i];
    }
  }
  read_event("conjugation_transfers", meta.interval_events.conjugation_transfers);
  read_event("mutations", meta.interval_events.mutations);
  read_event("immigrations", meta.interval_events.immigrations);
  read_event("cumulative_sos_inductions", meta.cumulative_events.sos_inductions);
  read_event("cumulative_phage_inductions", meta.cumulative_events.phage_inductions);
  read_event("cumulative_mortality_colicin", meta.cumulative_events.mortality_colicin);
  read_event("cumulative_mortality_cdi", meta.cumulative_events.mortality_cdi);
  read_event("cumulative_outflow_washout", meta.cumulative_events.outflow_washout);
  read_event("cumulative_outflow_boundary", meta.cumulative_events.outflow_boundary);
  read_event("cumulative_mortality_lysis", meta.cumulative_events.mortality_lysis);
  read_event("cumulative_divisions", meta.cumulative_events.divisions);
  if (link_exists(file, events + "cumulative_divisions_by_type")) {
    auto values = read_dataset_1d<int32_t>(
        file, events + "cumulative_divisions_by_type", H5T_NATIVE_INT32);
    for (size_t i = 0; i < values.size() && i < MAX_AGENT_TYPES; ++i) {
      meta.cumulative_events.divisions_by_type[i] = values[i];
    }
  }
  read_event("cumulative_conjugation_transfers",
             meta.cumulative_events.conjugation_transfers);
  read_event("cumulative_mutations", meta.cumulative_events.mutations);
  read_event("cumulative_immigrations", meta.cumulative_events.immigrations);
  const auto read_rate = [&](const std::string& name, Real& value) {
    if (link_exists(file, events + name)) {
      value = read_scalar<double>(file, events + name, H5T_NATIVE_DOUBLE);
    }
  };
  read_rate("sos_basal_rate", meta.interval_events.sos_basal_rate);
  read_rate("sos_post_division_rate",
            meta.interval_events.sos_post_division_rate);
  read_rate("sos_nuclease_cross_induction_rate",
            meta.interval_events.sos_nuclease_cross_induction_rate);
  read_rate("sos_ros_rate", meta.interval_events.sos_ros_rate);
  read_rate("cumulative_sos_basal_rate",
            meta.cumulative_events.sos_basal_rate);
  read_rate("cumulative_sos_post_division_rate",
            meta.cumulative_events.sos_post_division_rate);
  read_rate("cumulative_sos_nuclease_cross_induction_rate",
            meta.cumulative_events.sos_nuclease_cross_induction_rate);
  read_rate("cumulative_sos_ros_rate",
            meta.cumulative_events.sos_ros_rate);
  if (link_exists(file, events + "interval_start_step")) {
    meta.event_window_start_step =
        read_scalar<int32_t>(file, events + "interval_start_step", H5T_NATIVE_INT32);
    meta.event_window_end_step =
        read_scalar<int32_t>(file, events + "interval_end_step", H5T_NATIVE_INT32);
    meta.event_window_start_time =
        read_scalar<double>(file, events + "interval_start_time", H5T_NATIVE_DOUBLE);
    meta.event_window_end_time =
        read_scalar<double>(file, events + "interval_end_time", H5T_NATIVE_DOUBLE);
  }
  const std::string flux = prefix + "nutrient_flux/";
  const auto read_flux = [&](const std::string& name, std::vector<Real>& values) {
    if (link_exists(file, flux + name)) {
      values = read_dataset_1d<double>(file, flux + name, H5T_NATIVE_DOUBLE);
    }
  };
  read_flux("boundary_interval", meta.flux_accounting.boundary_interval);
  read_flux("boundary_cumulative", meta.flux_accounting.boundary_cumulative);
  read_flux("gradient_source_interval",
            meta.flux_accounting.gradient_source_interval);
  read_flux("gradient_source_cumulative",
            meta.flux_accounting.gradient_source_cumulative);
  read_flux("vbf_source_interval", meta.flux_accounting.vbf_source_interval);
  read_flux("vbf_source_cumulative", meta.flux_accounting.vbf_source_cumulative);
  read_flux("vbf_sink_interval", meta.flux_accounting.vbf_sink_interval);
  read_flux("vbf_sink_cumulative", meta.flux_accounting.vbf_sink_cumulative);
  read_flux("agent_uptake_interval", meta.flux_accounting.agent_uptake_interval);
  read_flux("agent_uptake_cumulative", meta.flux_accounting.agent_uptake_cumulative);
  read_flux("nutrient_blocking_fraction",
            meta.flux_accounting.nutrient_blocking_fraction);
  read_flux("maintenance_interval", meta.flux_accounting.maintenance_interval);
  read_flux("maintenance_cumulative",
            meta.flux_accounting.maintenance_cumulative);
  read_flux("maintenance_shortfall_interval",
            meta.flux_accounting.maintenance_shortfall_interval);
  read_flux("maintenance_shortfall_cumulative",
            meta.flux_accounting.maintenance_shortfall_cumulative);
  read_flux("maintenance_limited_agents_interval",
            meta.flux_accounting.maintenance_limited_agents_interval);
  read_flux("maintenance_limited_agents_cumulative",
            meta.flux_accounting.maintenance_limited_agents_cumulative);
  read_flux("uptake_demand_interval", meta.flux_accounting.uptake_demand_interval);
  read_flux("uptake_demand_cumulative",
            meta.flux_accounting.uptake_demand_cumulative);
  read_flux("uptake_shortfall_interval",
            meta.flux_accounting.uptake_shortfall_interval);
  read_flux("uptake_shortfall_cumulative",
            meta.flux_accounting.uptake_shortfall_cumulative);
  read_flux("uptake_limited_agents_interval",
            meta.flux_accounting.uptake_limited_interval);
  read_flux("uptake_limited_agents_cumulative",
            meta.flux_accounting.uptake_limited_cumulative);
  read_flux("delivery_reduction_interval",
            meta.flux_accounting.delivery_reduction_interval);
  read_flux("delivery_reduction_cumulative",
            meta.flux_accounting.delivery_reduction_cumulative);
  read_flux("delivery_retry_events_interval",
            meta.flux_accounting.delivery_retry_events_interval);
  read_flux("delivery_retry_events_cumulative",
            meta.flux_accounting.delivery_retry_events_cumulative);
  read_flux("delivery_rationing_factor_interval",
            meta.flux_accounting.delivery_rationing_factor_interval);
  read_flux("delivery_rationing_factor_cumulative",
            meta.flux_accounting.delivery_rationing_factor_cumulative);
  read_flux("delivery_infeasible_interval",
            meta.flux_accounting.delivery_infeasible_interval);
  read_flux("delivery_infeasible_cumulative",
            meta.flux_accounting.delivery_infeasible_cumulative);
  read_flux("reaction_clip_interval",
            meta.flux_accounting.reaction_clip_interval);
  read_flux("reaction_clip_cumulative",
            meta.flux_accounting.reaction_clip_cumulative);

  const std::string rng_prefix = prefix + "rng/";
  if (link_exists(file, prefix + "rng")
      && link_exists(file, rng_prefix + "simulation_state")) {
    meta.simulation_rng_states =
        h5::read_str_vec(file, (rng_prefix + "simulation_state").c_str());
    meta.has_rng_state = !meta.simulation_rng_states.empty();
    if (link_exists(file, rng_prefix + "immigration_state")) {
      const auto imm =
          h5::read_str_vec(file, (rng_prefix + "immigration_state").c_str());
      if (!imm.empty()) {
        meta.immigration_rng_state = imm.front();
      }
    }
    if (link_exists(file, rng_prefix + "seed")) {
      meta.rng_seed = static_cast<uint64_t>(read_scalar<unsigned long long>(
          file, rng_prefix + "seed", H5T_NATIVE_ULLONG));
      meta.has_rng_seed = true;
    }
    if (link_exists(file, rng_prefix + "immigration_enabled")) {
      meta.immigration_was_enabled = read_scalar<int32_t>(
          file, rng_prefix + "immigration_enabled", H5T_NATIVE_INT32) != 0;
      meta.has_immigration_enabled = true;
    }
    if (link_exists(file, rng_prefix + "immigration_start_step")) {
      meta.immigration_start_step = read_scalar<int32_t>(
          file, rng_prefix + "immigration_start_step", H5T_NATIVE_INT32);
      meta.has_immigration_start_step = true;
    }
    if (link_exists(file, rng_prefix + "next_tag")) {
      meta.next_tags =
          read_dataset_1d<int64_t>(file, rng_prefix + "next_tag", H5T_NATIVE_INT64);
    }
  }

  const std::string burst_prefix = prefix + "toxin_bursts/";
  if (link_exists(file, prefix + "toxin_bursts")
      && link_exists(file, burst_prefix + "pos_x")) {
    meta.has_toxin_bursts = true;
    const auto pos_x = read_dataset_1d<double>(file, burst_prefix + "pos_x", H5T_NATIVE_DOUBLE);
    const auto pos_y = read_dataset_1d<double>(file, burst_prefix + "pos_y", H5T_NATIVE_DOUBLE);
    const auto pos_z = read_dataset_1d<double>(file, burst_prefix + "pos_z", H5T_NATIVE_DOUBLE);
    const auto diff = read_dataset_1d<double>(file, burst_prefix + "diff_coeff", H5T_NATIVE_DOUBLE);
    const auto source = read_dataset_1d<double>(file, burst_prefix + "source_rate", H5T_NATIVE_DOUBLE);
    const auto pI = read_dataset_1d<double>(file, burst_prefix + "pI", H5T_NATIVE_DOUBLE);
    const auto ret = read_dataset_1d<double>(file, burst_prefix + "retardation", H5T_NATIVE_DOUBLE);
    const auto decay = read_dataset_1d<double>(file, burst_prefix + "decay_rate", H5T_NATIVE_DOUBLE);
    const auto lumen_len = read_dataset_1d<double>(file, burst_prefix + "lumen_transfer_length", H5T_NATIVE_DOUBLE);
    const auto lumen_free = read_dataset_1d<int32_t>(file, burst_prefix + "lumen_transfer_basis_free", H5T_NATIVE_INT32);
    const auto robin_cut = read_dataset_1d<double>(file, burst_prefix + "robin_cutoff", H5T_NATIVE_DOUBLE);
    const auto img_tol = read_dataset_1d<double>(file, burst_prefix + "image_series_relative_tolerance", H5T_NATIVE_DOUBLE);
    const auto img_shells = read_dataset_1d<int32_t>(file, burst_prefix + "image_series_max_shells", H5T_NATIVE_INT32);
    const auto img_explicit = read_dataset_1d<int32_t>(file, burst_prefix + "image_series_max_shells_explicit", H5T_NATIVE_INT32);
    const auto img_legacy = read_dataset_1d<int32_t>(file, burst_prefix + "image_series_legacy_reflections", H5T_NATIVE_INT32);
    const auto drift_corr = read_dataset_1d<int32_t>(file, burst_prefix + "drift_correction", H5T_NATIVE_INT32);
    const auto creation = read_dataset_1d<double>(file, burst_prefix + "creation_time", H5T_NATIVE_DOUBLE);
    const auto tau = read_dataset_1d<double>(file, burst_prefix + "release_tau", H5T_NATIVE_DOUBLE);
    const auto nuclease = read_dataset_1d<int32_t>(file, burst_prefix + "is_nuclease", H5T_NATIVE_INT32);
    const auto target = read_dataset_1d<int32_t>(file, burst_prefix + "target", H5T_NATIVE_INT32);
    const size_t nb = pos_x.size();
    if (pos_y.size() != nb || pos_z.size() != nb || diff.size() != nb
        || source.size() != nb || pI.size() != nb || ret.size() != nb
        || decay.size() != nb || lumen_len.size() != nb || lumen_free.size() != nb
        || robin_cut.size() != nb || img_tol.size() != nb || img_shells.size() != nb
        || img_explicit.size() != nb || img_legacy.size() != nb || drift_corr.size() != nb
        || creation.size() != nb || tau.size() != nb || nuclease.size() != nb
        || target.size() != nb) {
      throw HDF5Error("inconsistent toxin_bursts dataset lengths in summary/" + step);
    }
    meta.toxin_bursts.resize(nb);
    for (size_t i = 0; i < nb; ++i) {
      ToxinBurstSource& b = meta.toxin_bursts[i];
      b.pos = {pos_x[i], pos_y[i], pos_z[i]};
      b.params.diff_coeff = diff[i];
      b.params.source_rate = source[i];
      b.params.pI = pI[i];
      b.params.retardation = ret[i];
      b.params.decay_rate = decay[i];
      b.params.lumen_transfer_length = lumen_len[i];
      b.params.lumen_transfer_basis_free = (lumen_free[i] != 0);
      b.params.robin_cutoff = robin_cut[i];
      b.params.image_series_relative_tolerance = img_tol[i];
      b.params.image_series_max_shells = img_shells[i];
      b.params.image_series_max_shells_explicit = (img_explicit[i] != 0);
      b.params.image_series_legacy_reflections = (img_legacy[i] != 0);
      b.params.drift_correction = (drift_corr[i] != 0);
      b.creation_time = creation[i];
      b.release_tau = tau[i];
      b.is_nuclease = (nuclease[i] != 0);
      b.target = static_cast<ReceptorType>(target[i]);
    }
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
#ifdef GUTIBM_HDF5
    file_id_ = -1;
#endif
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
  if (const bool has_summary =
          link_exists(static_cast<hid_t>(file_id_), "summary/" + step_name);
      !has_agents && !has_summary) {
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
