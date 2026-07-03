/* -----------------------------------------------------------------------
   GutIBM – Input parser implementation
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "config_json.h"
#include "path_utils.h"
#include "plasmid.h"
#include <iostream>
#include <algorithm>
#include "error.h"
#include <sstream>
#include <stdexcept>
#include <cstdlib>

namespace gutibm {

namespace {

bool strict_config_enabled() {
  const char* env = std::getenv("GUTIBM_STRICT_CONFIG");
  if (!env || env[0] == '\0') return false;
  return env[0] != '0';
}

void warn_parse_failure(const char* kind,
                        const std::string& key,
                        const std::string& val) {
  std::cerr << "Warning: config key '" << key << "' has invalid " << kind
            << " value '" << val << "' — using 0\n";
  if (strict_config_enabled()) {
    throw ConfigError("invalid config value for key '" + key + "'");
  }
}

std::string trim_config(std::string_view s) {
  auto start = s.find_first_not_of(" \t\r\n");
  if (start == std::string_view::npos) return "";
  auto end = s.find_last_not_of(" \t\r\n");
  return std::string(s.substr(start, end - start + 1));
}

Real parse_config_real(const std::string& key, const std::string& val) {
  const std::string trimmed = trim_config(val);
  if (trimmed.empty()) {
    warn_parse_failure("numeric", key, val);
    return 0.0;
  }

  try {
    size_t consumed = 0;
    const Real result = std::stod(trimmed, &consumed);
    if (consumed != trimmed.size()) {
      warn_parse_failure("numeric", key, val);
      return 0.0;
    }
    return result;
  } catch (const std::invalid_argument&) {
    warn_parse_failure("numeric", key, val);
    return 0.0;
  } catch (const std::out_of_range&) {
    warn_parse_failure("numeric", key, val);
    return 0.0;
  }
}

Int parse_config_int(const std::string& key, const std::string& val) {
  const std::string trimmed = trim_config(val);
  if (trimmed.empty()) {
    warn_parse_failure("integer", key, val);
    return 0;
  }

  try {
    size_t consumed = 0;
    const Int result = std::stoi(trimmed, &consumed);
    if (consumed != trimmed.size()) {
      warn_parse_failure("integer", key, val);
      return 0;
    }
    return result;
  } catch (const std::invalid_argument&) {
    warn_parse_failure("integer", key, val);
    return 0;
  } catch (const std::out_of_range&) {
    warn_parse_failure("integer", key, val);
    return 0;
  }
}

}  // namespace

SimulationConfig InputParser::default_config() {
  SimulationConfig cfg;

  // Default chemical species
  // Carbon gets z-gradient enabled (mucin-derived monosaccharides peak at epithelium)
  cfg.chemicals = {
    {"carbon",      1.0e-9, 1.0,  5.0e-3, 5.0e-3,  0.0, true,  25.0e-6},
    {"iron",        1.0e-9, 1.0,  1.0e-4, 1.0e-4,  0.0, false, 25.0e-6},
    {"b12",         1.0e-9, 1.0,  1.0e-9, 1.0e-9,  0.0, false, 25.0e-6},
    {"bacteriocin", 4.0e-11, 10.0, 0.0,    0.0,     1.0e-4, false, 25.0e-6},
    {"acetate",     1.2e-9,  1.0,  80.0,   80.0,    0.0, false, 25.0e-6},
    {"ethanolamine", 1.0e-9, 1.0, 0.5e-3, 0.5e-3, 0.0, false, 25.0e-6},
  };

  // VBF mucin z-gradient enabled by default (consistent with carbon gradient)
  cfg.vbf.mucin_z_gradient_enabled = true;
  cfg.vbf.mucin_z_gradient_lambda  = 25.0e-6;

  // Default initial population: resident (B2 phylogroup) + immigrant
  SimulationConfig::InitialStrain resident;
  resident.type       = 1;
  resident.count      = 500;
  resident.mu_max     = 5.5e-4;  // ~30 min doubling in rich media
  resident.plasmids   = {"ColE1", "ColB"};
  resident.conjugative = true;

  SimulationConfig::InitialStrain immigrant;
  immigrant.type       = 2;
  immigrant.count      = 100;
  immigrant.mu_max     = 5.0e-4;
  immigrant.plasmids   = {};
  immigrant.conjugative = false;

  cfg.initial_strains = {resident, immigrant};

  finalize_config(cfg);
  return cfg;
}

namespace {

Int find_chemical_spec(std::vector<ChemicalSpec>& chemicals, std::string_view name) {
  for (Int i = 0; i < static_cast<Int>(chemicals.size()); ++i) {
    if (chemicals[static_cast<size_t>(i)].name == name) return i;
  }
  return -1;
}

bool parse_bool_config(const std::string& val) {
  return val == "true" || val == "1";
}

}  // namespace

void InputParser::finalize_config(SimulationConfig& cfg) {
  constexpr Real k_z_lambda = 25.0e-6;

  if (cfg.oxygen.enabled) {
    if (find_chemical_spec(cfg.chemicals, "oxygen") < 0) {
      const Real c0 = cfg.oxygen.epithelial_conc;
      cfg.chemicals.push_back({
          "oxygen", cfg.oxygen.D_free, 1.0, c0, c0, 0.0, true, k_z_lambda});
    }
  }

  if (cfg.acetate.enabled) {
    Int idx = find_chemical_spec(cfg.chemicals, "acetate");
    if (idx < 0) {
      cfg.chemicals.push_back({
          "acetate", cfg.acetate.D_free, 1.0, 0.0, 0.0, 0.0, false, k_z_lambda});
    } else {
      auto& spec = cfg.chemicals[static_cast<size_t>(idx)];
      spec.diff_coeff = cfg.acetate.D_free;
      spec.initial_conc = 0.0;
      spec.boundary_conc = 0.0;
    }
  }

  if (cfg.mucin.enabled) {
    if (find_chemical_spec(cfg.chemicals, "mucin") < 0) {
      const Real c0 = cfg.mucin.initial_conc;
      cfg.chemicals.push_back({
          "mucin", cfg.mucin.D_free, cfg.mucin.retardation,
          c0, c0, 0.0, false, k_z_lambda});
    }
    cfg.vbf.use_dynamic_mucin = true;
  }
}

namespace {

using FlatKeyHandler = bool (*)(SimulationConfig&, const std::string&, const std::string&);

bool apply_time_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "total_time")           { cfg.total_time = parse_config_real(key, val); return true; }
  if (key == "bio_dt")               { cfg.bio_dt = parse_config_real(key, val); return true; }
  if (key == "output_interval")      { cfg.output_interval = parse_config_real(key, val); return true; }
  if (key == "seed")                 { cfg.seed = static_cast<uint64_t>(parse_config_int(key, val)); return true; }
  return false;
}

bool apply_domain_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "grid_dx")              { cfg.domain.grid_dx = parse_config_real(key, val); return true; }
  if (key == "domain_x")             { cfg.domain.hi[0] = parse_config_real(key, val); return true; }
  if (key == "domain_y")             { cfg.domain.hi[1] = parse_config_real(key, val); return true; }
  if (key == "domain_z")             { cfg.domain.hi[2] = parse_config_real(key, val); return true; }
  return false;
}

bool apply_advection_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "mucus_thickness")      { cfg.advection.mucus_thickness = parse_config_real(key, val); return true; }
  if (key == "radial_turnover")      { cfg.advection.radial_turnover = parse_config_real(key, val); return true; }
  if (key == "distal_transit")       { cfg.advection.distal_transit_time = parse_config_real(key, val); return true; }
  if (key == "peristaltic_enabled")  { cfg.advection.peristaltic_enabled = (val == "true" || val == "1"); return true; }
  if (key == "peristaltic_period")   { cfg.advection.peristaltic_period = parse_config_real(key, val); return true; }
  if (key == "peristaltic_amplitude") { cfg.advection.peristaltic_amplitude = parse_config_real(key, val); return true; }
  if (key == "peristaltic_wavelength") { cfg.advection.peristaltic_wavelength = parse_config_real(key, val); return true; }
  if (key == "crypts_enabled")       { cfg.advection.crypts_enabled = (val == "true" || val == "1"); return true; }
  if (key == "crypt_depth")          { cfg.advection.crypt_depth = parse_config_real(key, val); return true; }
  if (key == "crypt_exit_rate")      { cfg.advection.crypt_exit_rate = parse_config_real(key, val); return true; }
  if (key == "crypt_entry_rate")     { cfg.advection.crypt_entry_rate = parse_config_real(key, val); return true; }
  if (key == "crypt_carrying_capacity") { cfg.advection.crypt_carrying_capacity = parse_config_int(key, val); return true; }
  return false;
}

bool apply_qssa_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "toxin_cutoff")         { cfg.qssa.toxin_cutoff = parse_config_real(key, val); return true; }
  if (key == "nutrient_cutoff")      { cfg.qssa.nutrient_cutoff = parse_config_real(key, val); return true; }
  if (key == "use_fmm")              { cfg.qssa.use_fmm = (val == "true" || val == "1"); return true; }
  if (key == "fmm_theta")            { cfg.qssa.fmm_theta = parse_config_real(key, val); return true; }
  if (key == "fmm_expansion_order")  { cfg.qssa.fmm_expansion_order = parse_config_int(key, val); return true; }
  return false;
}

bool apply_vbf_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "vbf_density")          { cfg.vbf.density = parse_config_real(key, val); return true; }
  if (key == "vbf_viscosity")        { cfg.vbf.viscosity = parse_config_real(key, val); return true; }
  if (key == "vbf_mucin_z_gradient") { cfg.vbf.mucin_z_gradient_enabled = (val == "true" || val == "1"); return true; }
  if (key == "vbf_mucin_z_lambda")   { cfg.vbf.mucin_z_gradient_lambda = parse_config_real(key, val); return true; }
  return false;
}

bool apply_chemical_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "carbon_z_gradient") {
    for (auto& c : cfg.chemicals) {
      if (c.name == "carbon") { c.z_gradient_enabled = (val == "true" || val == "1"); return true; }
    }
    return true;
  }
  if (key == "carbon_z_lambda") {
    for (auto& c : cfg.chemicals) {
      if (c.name == "carbon") { c.z_gradient_lambda = parse_config_real(key, val); return true; }
    }
    return true;
  }
  if (key == "sos_lysis_prob")       { cfg.bacteriocin.sos_lysis_prob = parse_config_real(key, val); return true; }
  return false;
}

bool apply_receptor_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "kd_b12_btuB")          { cfg.receptor.kd_b12_btuB = parse_config_real(key, val); return true; }
  if (key == "kd_colicinE_btuB")     { cfg.receptor.kd_colicinE_btuB = parse_config_real(key, val); return true; }
  if (key == "kd_enterobactin")       { cfg.receptor.kd_enterobactin = parse_config_real(key, val); return true; }
  if (key == "kd_colicinB_fepA")      { cfg.receptor.kd_colicinB_fepA = parse_config_real(key, val); return true; }
  if (key == "kd_lin_enterobactin")   { cfg.receptor.kd_lin_enterobactin = parse_config_real(key, val); return true; }
  if (key == "kd_colicinIa_cirA")    { cfg.receptor.kd_colicinIa_cirA = parse_config_real(key, val); return true; }
  if (key == "kill_rate_colicin")     { cfg.receptor.kill_rate_colicin = parse_config_real(key, val); return true; }
  if (key == "kill_rate_microcin")    { cfg.receptor.kill_rate_microcin = parse_config_real(key, val); return true; }
  if (key == "immunity_factor")       { cfg.receptor.immunity_factor = parse_config_real(key, val); return true; }
  return false;
}

bool apply_conjugation_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "pili_length")           { cfg.conjugation.pili_length = parse_config_real(key, val); return true; }
  if (key == "base_transfer_rate")    { cfg.conjugation.base_transfer_rate = parse_config_real(key, val); return true; }
  if (key == "shear_critical")          { cfg.conjugation.shear_critical = parse_config_real(key, val); return true; }
  if (key == "plasmid_copy_cost")     { cfg.conjugation.plasmid_copy_cost = parse_config_real(key, val); return true; }
  if (key == "pili_heterogeneity")    { cfg.conjugation.pili_heterogeneity = (val == "true" || val == "1"); return true; }
  if (key == "pili_length_min")       { cfg.conjugation.pili_length_min = parse_config_real(key, val); return true; }
  if (key == "pili_length_max")       { cfg.conjugation.pili_length_max = parse_config_real(key, val); return true; }
  return false;
}

bool apply_mutation_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "bi_duplication_rate")     { cfg.mutation.bi_duplication_rate = parse_config_real(key, val); return true; }
  if (key == "bi_recombination_rate")   { cfg.mutation.bi_recombination_rate = parse_config_real(key, val); return true; }
  if (key == "receptor_mutation_rate")  { cfg.mutation.receptor_mutation_rate = parse_config_real(key, val); return true; }
  if (key == "super_killer_rate")       { cfg.mutation.super_killer_rate = parse_config_real(key, val); return true; }
  if (key == "partial_resistance_rate") { cfg.mutation.partial_resistance_rate = parse_config_real(key, val); return true; }
  if (key == "receptor_reduction")      { cfg.mutation.receptor_reduction = parse_config_real(key, val); return true; }
  if (key == "max_bi_loci")             { cfg.mutation.max_bi_loci = parse_config_int(key, val); return true; }
  if (key == "immunity_escape_prob")    { cfg.mutation.immunity_escape_prob = parse_config_real(key, val); return true; }
  if (key == "escape_affinity_lo")      { cfg.mutation.escape_affinity_lo = parse_config_real(key, val); return true; }
  if (key == "escape_affinity_hi")      { cfg.mutation.escape_affinity_hi = parse_config_real(key, val); return true; }
  if (key == "compensatory_rate")       { cfg.mutation.compensatory_rate = parse_config_real(key, val); return true; }
  if (key == "compensatory_reduction")  { cfg.mutation.compensatory_reduction = parse_config_real(key, val); return true; }
  return false;
}

bool apply_io_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "hdf5_file") {
    validate_path_syntax(val);
    cfg.hdf5.filename = val;
    return true;
  }
  if (key == "hdf5_every")           { cfg.hdf5.dump_every = parse_config_int(key, val); return true; }
  if (key == "checkpoint_file") {
    validate_path_syntax(val);
    cfg.checkpoint.file = val;
    return true;
  }
  if (key == "checkpoint_step")        { cfg.checkpoint.step = val; return true; }
  return false;
}

bool apply_dt_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "adaptive_dt_enabled")  { cfg.adaptive_dt_enabled = (val == "true" || val == "1"); return true; }
  if (key == "dt_min")               { cfg.dt_min = parse_config_real(key, val); return true; }
  if (key == "dt_max")               { cfg.dt_max = parse_config_real(key, val); return true; }
  if (key == "dt_safety")            { cfg.dt_safety = parse_config_real(key, val); return true; }
  if (key == "dt_growth_limit")      { cfg.dt_growth_limit = parse_config_real(key, val); return true; }
  return false;
}

bool apply_misc_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "gpu_enabled")          { cfg.gpu.enabled = (val == "true" || val == "1"); return true; }
  if (key == "gpu_device_id")        { cfg.gpu.device_id = parse_config_int(key, val); return true; }
  if (key == "profile_steps")        { cfg.profile_steps = (val == "true" || val == "1"); return true; }
  return false;
}

bool apply_oxygen_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "oxygen.enabled" || key == "oxygen_enabled") {
    cfg.oxygen.enabled = parse_bool_config(val); return true;
  }
  if (key == "oxygen.epithelial_conc" || key == "oxygen_epithelial_conc") {
    cfg.oxygen.epithelial_conc = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.D_free" || key == "oxygen_D_free") {
    cfg.oxygen.D_free = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.Km" || key == "oxygen_Km") {
    cfg.oxygen.Km = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.boost_max" || key == "oxygen_boost_max") {
    cfg.oxygen.boost_max = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.q_consumption" || key == "oxygen_q_consumption") {
    cfg.oxygen.q_consumption = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.vbf_sink" || key == "oxygen_vbf_sink") {
    cfg.oxygen.vbf_sink = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.k_ROS" || key == "oxygen_k_ROS") {
    cfg.oxygen.k_ROS = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_acetate_dyn_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "acetate.enabled" || key == "acetate_enabled") {
    cfg.acetate.enabled = parse_bool_config(val); return true;
  }
  if (key == "acetate.D_free" || key == "acetate_D_free") {
    cfg.acetate.D_free = parse_config_real(key, val); return true;
  }
  if (key == "acetate.vbf_production" || key == "acetate_vbf_production") {
    cfg.acetate.vbf_production = parse_config_real(key, val); return true;
  }
  if (key == "acetate.vbf_consumption" || key == "acetate_vbf_consumption") {
    cfg.acetate.vbf_consumption = parse_config_real(key, val); return true;
  }
  if (key == "acetate.overflow_threshold" || key == "acetate_overflow_threshold") {
    cfg.acetate.overflow_threshold = parse_config_real(key, val); return true;
  }
  if (key == "acetate.overflow_rate" || key == "acetate_overflow_rate") {
    cfg.acetate.overflow_rate = parse_config_real(key, val); return true;
  }
  if (key == "acetate.scavenge_Km" || key == "acetate_scavenge_Km") {
    cfg.acetate.scavenge_Km = parse_config_real(key, val); return true;
  }
  if (key == "acetate.scavenge_rate" || key == "acetate_scavenge_rate") {
    cfg.acetate.scavenge_rate = parse_config_real(key, val); return true;
  }
  if (key == "acetate.epithelial_uptake" || key == "acetate_epithelial_uptake") {
    cfg.acetate.epithelial_uptake = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_mucin_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "mucin.enabled" || key == "mucin_enabled") {
    cfg.mucin.enabled = parse_bool_config(val); return true;
  }
  if (key == "mucin.secretion_rate" || key == "mucin_secretion_rate") {
    cfg.mucin.secretion_rate = parse_config_real(key, val); return true;
  }
  if (key == "mucin.Km_degradation" || key == "mucin_Km_degradation") {
    cfg.mucin.Km_degradation = parse_config_real(key, val); return true;
  }
  if (key == "mucin.k_liberation" || key == "mucin_k_liberation") {
    cfg.mucin.k_liberation = parse_config_real(key, val); return true;
  }
  if (key == "mucin.initial_conc" || key == "mucin_initial_conc") {
    cfg.mucin.initial_conc = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_protease_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "protease.enabled" || key == "protease_enabled") {
    cfg.protease.enabled = parse_bool_config(val); return true;
  }
  if (key == "protease.default_half_life" || key == "protease_default_half_life") {
    cfg.protease.default_half_life = parse_config_real(key, val); return true;
  }
  if (key == "protease.dilution_rate" || key == "protease_dilution_rate") {
    cfg.protease.dilution_rate = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_fur_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "fur.enabled" || key == "fur_enabled") {
    cfg.fur.enabled = parse_bool_config(val); return true;
  }
  if (key == "fur.Km" || key == "fur_Km") {
    cfg.fur.Km = parse_config_real(key, val); return true;
  }
  if (key == "fur.upregulation_max" || key == "fur_upregulation_max") {
    cfg.fur.upregulation_max = parse_config_real(key, val); return true;
  }
  if (key == "fur.receptor_max" || key == "fur_receptor_max") {
    cfg.fur.receptor_max = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_cdi_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "cdi.enabled" || key == "cdi_enabled") {
    cfg.cdi.enabled = parse_bool_config(val); return true;
  }
  if (key == "cdi.kill_rate" || key == "cdi_kill_rate") {
    cfg.cdi.kill_rate = parse_config_real(key, val); return true;
  }
  if (key == "cdi.contact_radius" || key == "cdi_contact_radius") {
    cfg.cdi.contact_radius = parse_config_real(key, val); return true;
  }
  if (key == "cdi.corpse_persistence" || key == "cdi_corpse_persistence") {
    cfg.cdi.corpse_persistence = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_motility_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "motility.enabled" || key == "motility_enabled") {
    cfg.motility.enabled = parse_bool_config(val); return true;
  }
  if (key == "motility.swim_speed" || key == "motility_swim_speed") {
    cfg.motility.swim_speed = parse_config_real(key, val); return true;
  }
  if (key == "motility.run_mean_duration" || key == "motility_run_mean_duration") {
    cfg.motility.run_mean_duration = parse_config_real(key, val); return true;
  }
  if (key == "motility.stop_probability" || key == "motility_stop_probability") {
    cfg.motility.stop_probability = parse_config_real(key, val); return true;
  }
  if (key == "motility.stop_duration" || key == "motility_stop_duration") {
    cfg.motility.stop_duration = parse_config_real(key, val); return true;
  }
  if (key == "motility.chemotaxis_enabled" || key == "motility_chemotaxis_enabled") {
    cfg.motility.chemotaxis_enabled = parse_bool_config(val); return true;
  }
  if (key == "motility.chi_carbon" || key == "motility_chi_carbon") {
    cfg.motility.chi_carbon = parse_config_real(key, val); return true;
  }
  if (key == "motility.chi_oxygen" || key == "motility_chi_oxygen") {
    cfg.motility.chi_oxygen = parse_config_real(key, val); return true;
  }
  if (key == "motility.cluster_suppress_radius" || key == "motility_cluster_suppress_radius") {
    cfg.motility.cluster_suppress_radius = parse_config_real(key, val); return true;
  }
  if (key == "motility.cluster_suppress_threshold" || key == "motility_cluster_suppress_threshold") {
    cfg.motility.cluster_suppress_threshold = parse_config_int(key, val); return true;
  }
  if (key == "motility.cluster_tumble_factor" || key == "motility_cluster_tumble_factor") {
    cfg.motility.cluster_tumble_factor = parse_config_real(key, val); return true;
  }
  return false;
}

constexpr FlatKeyHandler k_flat_key_handlers[] = {
  apply_time_key,
  apply_domain_key,
  apply_advection_key,
  apply_qssa_key,
  apply_vbf_key,
  apply_chemical_key,
  apply_receptor_key,
  apply_conjugation_key,
  apply_mutation_key,
  apply_io_key,
  apply_dt_key,
  apply_misc_key,
  apply_oxygen_key,
  apply_acetate_dyn_key,
  apply_mucin_key,
  apply_protease_key,
  apply_fur_key,
  apply_cdi_key,
  apply_motility_key,
};

bool parse_legacy_key_value(const std::string& line,
                            std::string& key_out,
                            std::string& val_out) {
  std::string trimmed = trim_config(line);
  if (trimmed.empty() || trimmed[0] == '#' || trimmed[0] == '/' ||
      trimmed[0] == '{' || trimmed[0] == '}') {
    return false;
  }

  auto colon = trimmed.find(':');
  if (colon == std::string::npos) return false;

  std::string key = trim_config(trimmed.substr(0, colon));
  std::string val = trim_config(trimmed.substr(colon + 1));
  if (key.empty() || key.front() == '_') return false;

  if (key.size() >= 2 && key.front() == '"' && key.back() == '"') {
    key = key.substr(1, key.size() - 2);
    if (!key.empty() && key.front() == '_') return false;
  }

  if (!val.empty() && val.back() == ',') val.pop_back();
  val = trim_config(val);
  if (val.size() >= 2 && val.front() == '"' && val.back() == '"') {
    val = val.substr(1, val.size() - 2);
  }

  key_out = std::move(key);
  val_out = std::move(val);
  return true;
}

void parse_legacy_flat_keys(const std::string& content, SimulationConfig& cfg) {
  std::istringstream lines(content);
  std::string line;
  while (std::getline(lines, line)) {
    std::string key;
    std::string val;
    if (!parse_legacy_key_value(line, key, val)) continue;
    InputParser::apply_flat_key(cfg, key, val);
  }
}

}  // namespace

void InputParser::apply_flat_key(SimulationConfig& cfg,
                                 const std::string& key,
                                 const std::string& val) {
  for (FlatKeyHandler handler : k_flat_key_handlers) {
    if (handler(cfg, key, val)) return;
  }
}

SimulationConfig InputParser::parse(const std::string& filename) {
  SimulationConfig cfg = default_config();

  std::string config_path;
  try {
    config_path = validate_input_file_path(filename);
  } catch (const IOError& ex) {
    std::cerr << "Warning: cannot open config file '" << filename
              << "' (" << ex.what() << "), using defaults\n";
    return cfg;
  }

  std::ifstream ifs(config_path);
  if (!ifs.is_open()) {
    std::cerr << "Warning: cannot open config file '" << config_path
              << "', using defaults\n";
    return cfg;
  }

  std::ostringstream content_stream;
  content_stream << ifs.rdbuf();
  const std::string content = content_stream.str();

  if (ConfigJson::parse_document(cfg, content)) {
    finalize_config(cfg);
    return cfg;
  }

  parse_legacy_flat_keys(content, cfg);

  if (auto strains = ConfigJson::parse_initial_strains(content); strains.found) {
    cfg.initial_strains = std::move(strains.strains);
  }

  if (auto fixes = ConfigJson::parse_enabled_fixes(content); fixes.found) {
    cfg.enabled_fixes = std::move(fixes.names);
  }

  finalize_config(cfg);
  return cfg;
}

std::string InputParser::trim(std::string_view s) {
  return trim_config(s);
}

Real InputParser::parse_real(const std::string& key, const std::string& val) {
  return parse_config_real(key, val);
}

Int InputParser::parse_int(const std::string& key, const std::string& val) {
  return parse_config_int(key, val);
}

}  // namespace gutibm
