/* -----------------------------------------------------------------------
   GutIBM – Input parser implementation
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "species_names.h"
#include "config_json.h"
#include "path_utils.h"
#include "plasmid.h"
#include <iostream>
#include <algorithm>
#include "error.h"
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <cstdlib>
#include <array>
#include <cctype>

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
    {species::CARBON,      1.0e-9, 1.0,  5.0e-3, 5.0e-3,  0.0, true,  25.0e-6},
    {species::IRON,        1.0e-9, 1.0,  1.0e-4, 1.0e-4,  0.0, false, 25.0e-6},
    {species::B12,         1.0e-9, 1.0,  1.0e-9, 1.0e-9,  0.0, false, 25.0e-6},
    {species::BACTERIOCIN_BTUB, 4.0e-11, 10.0, 0.0, 0.0, 1.0e-4, false, 25.0e-6},
    {species::BACTERIOCIN_FEPA, 4.0e-11, 10.0, 0.0, 0.0, 1.0e-4, false, 25.0e-6},
    {species::BACTERIOCIN_CIRA, 4.0e-11, 10.0, 0.0, 0.0, 1.0e-4, false, 25.0e-6},
    {species::BACTERIOCIN_FHUA, 4.0e-11, 10.0, 0.0, 0.0, 1.0e-4, false, 25.0e-6},
    {species::ACETATE,     1.2e-9,  1.0,  80.0,   80.0,    0.0, false, 25.0e-6},
    {species::ETHANOLAMINE, 1.0e-9, 1.0, 0.5e-3, 0.5e-3, 0.0, false, 25.0e-6},
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

  cfg.hdf5.schedule.summary = 1;
  cfg.hdf5.schedule.agents = 5;
  cfg.hdf5.schedule.grid = 0;
  cfg.hdf5.schedule.lineage = 100;
  cfg.hdf5.schedule.genome = 100;

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

bool parse_bool_config(std::string_view val) {
  if (val == "1") return true;
  if (val == "0") return false;
  std::string lower;
  lower.reserve(val.size());
  for (char c : val) {
    lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  }
  return lower == "true" || lower == "yes" || lower == "on";
}

}  // namespace

void InputParser::finalize_config(SimulationConfig& cfg) {
  constexpr Real k_z_lambda = 25.0e-6;

  if (cfg.chem_env.oxygen.enabled && find_chemical_spec(cfg.chemicals, species::OXYGEN) < 0) {
    cfg.chemicals.emplace_back(
        species::OXYGEN, cfg.chem_env.oxygen.D_free, 1.0,
        cfg.chem_env.oxygen.epithelial_conc, cfg.chem_env.oxygen.epithelial_conc,
        0.0, true, k_z_lambda);
  }

  if (cfg.chem_env.acetate.enabled) {
    Int idx = find_chemical_spec(cfg.chemicals, species::ACETATE);
    if (idx < 0) {
      cfg.chemicals.emplace_back(
          species::ACETATE, cfg.chem_env.acetate.D_free, 1.0, 0.0, 0.0, 0.0, false, k_z_lambda);
    } else {
      auto& spec = cfg.chemicals[static_cast<size_t>(idx)];
      spec.diff_coeff = cfg.chem_env.acetate.D_free;
      spec.initial_conc = 0.0;
      spec.boundary_conc = 0.0;
    }
  }

  if (cfg.chem_env.mucin.enabled) {
    if (find_chemical_spec(cfg.chemicals, species::MUCIN) < 0) {
      cfg.chemicals.emplace_back(
          species::MUCIN, cfg.chem_env.mucin.D_free, cfg.chem_env.mucin.retardation,
          cfg.chem_env.mucin.initial_conc, cfg.chem_env.mucin.initial_conc,
          0.0, false, k_z_lambda);
    }
    cfg.vbf.use_dynamic_mucin = true;
  }

  if (cfg.chem_env.siderophore.enabled) {
    if (find_chemical_spec(cfg.chemicals, species::SIDEROPHORE) < 0) {
      cfg.chemicals.emplace_back(
          species::SIDEROPHORE, cfg.chem_env.siderophore.D_free, 1.0,
          0.0, 0.0, 0.0, false, k_z_lambda);
    }
  }
}

namespace {

using FlatKeyHandler = bool (*)(SimulationConfig&, const std::string&, const std::string&);

bool apply_time_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "total_time")           { cfg.time.total_time = parse_config_real(key, val); return true; }
  if (key == "bio_dt")               { cfg.time.bio_dt = parse_config_real(key, val); return true; }
  if (key == "output_interval")      { cfg.time.output_interval = parse_config_real(key, val); return true; }
  if (key == "seed")                 { cfg.seed = static_cast<uint64_t>(parse_config_int(key, val)); return true; }
  return false;
}

bool apply_domain_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "grid_dx")              { cfg.domain.grid_dx = parse_config_real(key, val); return true; }
  if (key == "domain_x")             { cfg.domain.hi[0] = parse_config_real(key, val); return true; }
  if (key == "domain_y")             { cfg.domain.hi[1] = parse_config_real(key, val); return true; }
  if (key == "domain_z")             { cfg.domain.hi[2] = parse_config_real(key, val); return true; }
  if (key == "hash_cell_size")       { cfg.domain.hash_cell_size = parse_config_real(key, val); return true; }
  if (key == "ghost_width")          { cfg.domain.ghost_width = parse_config_real(key, val); return true; }
  return false;
}

bool apply_advection_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "mucus_thickness")      { cfg.advection.mucus_thickness = parse_config_real(key, val); return true; }
  if (key == "radial_turnover")      { cfg.advection.radial_turnover = parse_config_real(key, val); return true; }
  if (key == "distal_transit")       { cfg.advection.distal_transit_time = parse_config_real(key, val); return true; }
  if (key == "distal_length")        { cfg.advection.distal_length = parse_config_real(key, val); return true; }
  if (key == "profile_alpha")        { cfg.advection.profile_alpha = parse_config_real(key, val); return true; }
  if (key == "taylor_aris_enabled")  { cfg.advection.taylor_aris_enabled = parse_bool_config(val); return true; }
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
  if (key == "colicin_release_rate") { cfg.qssa.colicin_release_rate = parse_config_real(key, val); return true; }
  if (key == "microcin_secretion")   { cfg.qssa.microcin_secretion = parse_config_real(key, val); return true; }
  if (key == "use_fmm")              { cfg.qssa.use_fmm = (val == "true" || val == "1"); return true; }
  if (key == "fmm_theta")            { cfg.qssa.fmm_theta = parse_config_real(key, val); return true; }
  if (key == "fmm_expansion_order")  { cfg.qssa.fmm_expansion_order = parse_config_int(key, val); return true; }
  return false;
}

bool apply_vbf_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "vbf_density")          { cfg.vbf.density = parse_config_real(key, val); return true; }
  if (key == "vbf_viscosity")        { cfg.vbf.viscosity = parse_config_real(key, val); return true; }
  if (key == "vbf_drag_coeff")       { cfg.vbf.drag_coeff = parse_config_real(key, val); return true; }
  if (key == "vbf_nutrient_sink")    { cfg.vbf.nutrient_sink = parse_config_real(key, val); return true; }
  if (key == "vbf_mucin_liberation") { cfg.vbf.mucin_liberation = parse_config_real(key, val); return true; }
  if (key == "vbf_carrying_cap")     { cfg.vbf.carrying_cap = parse_config_real(key, val); return true; }
  if (key == "vbf_mucin_z_gradient") { cfg.vbf.mucin_z_gradient_enabled = (val == "true" || val == "1"); return true; }
  if (key == "vbf_mucin_z_lambda")   { cfg.vbf.mucin_z_gradient_lambda = parse_config_real(key, val); return true; }
  return false;
}

bool apply_chemical_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "carbon_z_gradient") {
    for (auto& c : cfg.chemicals) {
      if (c.name == species::CARBON) { c.z_gradient_enabled = (val == "true" || val == "1"); return true; }
    }
    return true;
  }
  if (key == "carbon_z_lambda") {
    for (auto& c : cfg.chemicals) {
      if (c.name == species::CARBON) { c.z_gradient_lambda = parse_config_real(key, val); return true; }
    }
    return true;
  }
  if (key == "sos_lysis_prob")       { cfg.fixes.bacteriocin.sos_lysis_prob = parse_config_real(key, val); return true; }
  if (key == "sos_basal_rate")       { cfg.fixes.bacteriocin.sos_basal_rate = parse_config_real(key, val); return true; }
  if (key == "sos_cross_induction_rate") {
    cfg.fixes.bacteriocin.sos_cross_induction_rate = parse_config_real(key, val);
    return true;
  }
  if (key == "retardation_basic")    { cfg.fixes.bacteriocin.retardation_basic = parse_config_real(key, val); return true; }
  if (key == "retardation_acidic")   { cfg.fixes.bacteriocin.retardation_acidic = parse_config_real(key, val); return true; }
  if (key == "retardation_neutral")  { cfg.fixes.bacteriocin.retardation_neutral = parse_config_real(key, val); return true; }
  if (key == "D_free_colicin")       { cfg.fixes.bacteriocin.D_free_colicin = parse_config_real(key, val); return true; }
  if (key == "burst_molecules")      { cfg.fixes.bacteriocin.burst_molecules = parse_config_real(key, val); return true; }
  if (key == "microcin_mu_penalty")  { cfg.fixes.bacteriocin.microcin_mu_penalty = parse_config_real(key, val); return true; }
  return false;
}

bool apply_receptor_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "kd_b12_btuB")          { cfg.fixes.receptor.kd_b12_btuB = parse_config_real(key, val); return true; }
  if (key == "kd_colicinE_btuB")     { cfg.fixes.receptor.kd_colicinE_btuB = parse_config_real(key, val); return true; }
  if (key == "kd_enterobactin")       { cfg.fixes.receptor.kd_enterobactin = parse_config_real(key, val); return true; }
  if (key == "kd_colicinB_fepA")      { cfg.fixes.receptor.kd_colicinB_fepA = parse_config_real(key, val); return true; }
  if (key == "kd_lin_enterobactin")   { cfg.fixes.receptor.kd_lin_enterobactin = parse_config_real(key, val); return true; }
  if (key == "kd_colicinIa_cirA")    { cfg.fixes.receptor.kd_colicinIa_cirA = parse_config_real(key, val); return true; }
  if (key == "kill_rate_colicin")     { cfg.fixes.receptor.kill_rate_colicin = parse_config_real(key, val); return true; }
  if (key == "kill_rate_microcin")    { cfg.fixes.receptor.kill_rate_microcin = parse_config_real(key, val); return true; }
  if (key == "immunity_factor")       { cfg.fixes.receptor.immunity_factor = parse_config_real(key, val); return true; }
  return false;
}

bool apply_conjugation_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "pili_length")           { cfg.fixes.conjugation.pili_length = parse_config_real(key, val); return true; }
  if (key == "base_transfer_rate")    { cfg.fixes.conjugation.base_transfer_rate = parse_config_real(key, val); return true; }
  if (key == "shear_critical")          { cfg.fixes.conjugation.shear_critical = parse_config_real(key, val); return true; }
  if (key == "plasmid_copy_cost")     { cfg.fixes.conjugation.plasmid_copy_cost = parse_config_real(key, val); return true; }
  if (key == "pili_heterogeneity")    { cfg.fixes.conjugation.pili_heterogeneity = (val == "true" || val == "1"); return true; }
  if (key == "pili_length_min")       { cfg.fixes.conjugation.pili_length_min = parse_config_real(key, val); return true; }
  if (key == "pili_length_max")       { cfg.fixes.conjugation.pili_length_max = parse_config_real(key, val); return true; }
  return false;
}

bool apply_mutation_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "bi_duplication_rate")     { cfg.fixes.mutation.bi_duplication_rate = parse_config_real(key, val); return true; }
  if (key == "bi_recombination_rate")   { cfg.fixes.mutation.bi_recombination_rate = parse_config_real(key, val); return true; }
  if (key == "receptor_mutation_rate")  { cfg.fixes.mutation.receptor_mutation_rate = parse_config_real(key, val); return true; }
  if (key == "super_killer_rate")       { cfg.fixes.mutation.super_killer_rate = parse_config_real(key, val); return true; }
  if (key == "partial_resistance_rate") { cfg.fixes.mutation.partial_resistance_rate = parse_config_real(key, val); return true; }
  if (key == "receptor_reduction")      { cfg.fixes.mutation.receptor_reduction = parse_config_real(key, val); return true; }
  if (key == "max_bi_loci")             { cfg.fixes.mutation.max_bi_loci = parse_config_int(key, val); return true; }
  if (key == "immunity_escape_prob")    { cfg.fixes.mutation.immunity_escape_prob = parse_config_real(key, val); return true; }
  if (key == "escape_affinity_lo")      { cfg.fixes.mutation.escape_affinity_lo = parse_config_real(key, val); return true; }
  if (key == "escape_affinity_hi")      { cfg.fixes.mutation.escape_affinity_hi = parse_config_real(key, val); return true; }
  if (key == "compensatory_rate")       { cfg.fixes.mutation.compensatory_rate = parse_config_real(key, val); return true; }
  if (key == "compensatory_reduction")  { cfg.fixes.mutation.compensatory_reduction = parse_config_real(key, val); return true; }
  return false;
}

bool apply_io_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "hdf5_file" || key == "hdf5.file") {
    validate_path_syntax(val);
    cfg.hdf5.filename = val;
    return true;
  }
  if (key == "checkpoint_file") {
    validate_path_syntax(val);
    cfg.checkpoint.file = val;
    return true;
  }
  if (key == "checkpoint_step")        { cfg.checkpoint.step = val; return true; }
  return false;
}

bool apply_hdf5_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "hdf5.enabled" || key == "hdf5_enabled") {
    cfg.hdf5.enabled = parse_bool_config(val); return true;
  }
  if (key == "hdf5.compression") {
    cfg.hdf5.compression = val; return true;
  }
  if (key == "hdf5.compression_level" || key == "hdf5_compression_level") {
    cfg.hdf5.compression_level = parse_config_int(key, val); return true;
  }
  if (key == "hdf5.schedule.summary" || key == "hdf5_schedule_summary") {
    cfg.hdf5.schedule.summary = parse_config_int(key, val); return true;
  }
  if (key == "hdf5.schedule.agents" || key == "hdf5_schedule_agents") {
    cfg.hdf5.schedule.agents = parse_config_int(key, val); return true;
  }
  if (key == "hdf5.schedule.grid" || key == "hdf5_schedule_grid") {
    cfg.hdf5.schedule.grid = parse_config_int(key, val); return true;
  }
  if (key == "hdf5.schedule.lineage" || key == "hdf5_schedule_lineage") {
    cfg.hdf5.schedule.lineage = parse_config_int(key, val); return true;
  }
  if (key == "hdf5.schedule.genome" || key == "hdf5_schedule_genome") {
    cfg.hdf5.schedule.genome = parse_config_int(key, val); return true;
  }
  return false;
}

bool apply_siderophore_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "siderophore.enabled" || key == "siderophore_enabled") {
    cfg.chem_env.siderophore.enabled = parse_bool_config(val); return true;
  }
  if (key == "siderophore.secretion_rate" || key == "siderophore_secretion_rate") {
    cfg.chem_env.siderophore.secretion_rate = parse_config_real(key, val); return true;
  }
  if (key == "siderophore.D_free" || key == "siderophore_D_free") {
    cfg.chem_env.siderophore.D_free = parse_config_real(key, val); return true;
  }
  if (key == "siderophore.chelation_rate" || key == "siderophore_chelation_rate") {
    cfg.chem_env.siderophore.chelation_rate = parse_config_real(key, val); return true;
  }
  if (key == "siderophore.Km_reimport" || key == "siderophore_Km_reimport") {
    cfg.chem_env.siderophore.Km_reimport = parse_config_real(key, val); return true;
  }
  if (key == "siderophore.recapture_fraction" || key == "siderophore_recapture_fraction") {
    cfg.chem_env.siderophore.recapture_fraction = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_dt_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "adaptive_dt_enabled")  { cfg.adaptive_dt.enabled = (val == "true" || val == "1"); return true; }
  if (key == "dt_min")               { cfg.adaptive_dt.min = parse_config_real(key, val); return true; }
  if (key == "dt_max")               { cfg.adaptive_dt.max = parse_config_real(key, val); return true; }
  if (key == "dt_safety")            { cfg.adaptive_dt.safety = parse_config_real(key, val); return true; }
  if (key == "dt_growth_limit")      { cfg.adaptive_dt.growth_limit = parse_config_real(key, val); return true; }
  return false;
}

bool apply_metabolism_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "division_threshold")      { cfg.fixes.metabolism.division_threshold = parse_config_real(key, val); return true; }
  if (key == "maintenance_rate")        { cfg.fixes.metabolism.maintenance_rate = parse_config_real(key, val); return true; }
  if (key == "metE_penalty")            { cfg.fixes.metabolism.metE_penalty = parse_config_real(key, val); return true; }
  if (key == "metE_acetate_km")         { cfg.fixes.metabolism.metE_acetate_km = parse_config_real(key, val); return true; }
  if (key == "metE_acetate_max_factor") { cfg.fixes.metabolism.metE_acetate_max_factor = parse_config_real(key, val); return true; }
  if (key == "eut_km")                  { cfg.fixes.metabolism.eut_km = parse_config_real(key, val); return true; }
  if (key == "eut_max_penalty")         { cfg.fixes.metabolism.eut_max_penalty = parse_config_real(key, val); return true; }
  if (key == "km_iron_primary")         { cfg.fixes.metabolism.km_iron_primary = parse_config_real(key, val); return true; }
  if (key == "km_iron_iroN")            { cfg.fixes.metabolism.km_iron_iroN = parse_config_real(key, val); return true; }
  if (key == "km_iron_iutA")            { cfg.fixes.metabolism.km_iron_iutA = parse_config_real(key, val); return true; }
  if (key == "km_iron_fiu")             { cfg.fixes.metabolism.km_iron_fiu = parse_config_real(key, val); return true; }
  return false;
}

bool apply_mechanics_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "hertz_k")            { cfg.fixes.mechanics.hertz_k = parse_config_real(key, val); return true; }
  if (key == "hertzian_enabled")   { cfg.fixes.mechanics.hertzian_enabled = parse_bool_config(val); return true; }
  if (key == "adhesion_enabled")   { cfg.fixes.mechanics.adhesion_enabled = parse_bool_config(val); return true; }
  if (key == "adhesion_strength")  { cfg.fixes.mechanics.adhesion_strength = parse_config_real(key, val); return true; }
  if (key == "adhesion_range")     { cfg.fixes.mechanics.adhesion_range = parse_config_real(key, val); return true; }
  return false;
}

bool apply_misc_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "gpu_enabled")          { cfg.gpu.enabled = parse_bool_config(val); return true; }
  if (key == "gpu_device_id")        { cfg.gpu.device_id = parse_config_int(key, val); return true; }
  if (key == "profile_steps")        { cfg.profile_steps = (val == "true" || val == "1"); return true; }
  return false;
}

bool apply_oxygen_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "oxygen.enabled" || key == "oxygen_enabled") {
    cfg.chem_env.oxygen.enabled = parse_bool_config(val); return true;
  }
  if (key == "oxygen.epithelial_conc" || key == "oxygen_epithelial_conc") {
    cfg.chem_env.oxygen.epithelial_conc = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.D_free" || key == "oxygen_D_free") {
    cfg.chem_env.oxygen.D_free = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.Km" || key == "oxygen_Km") {
    cfg.chem_env.oxygen.Km = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.boost_max" || key == "oxygen_boost_max") {
    cfg.chem_env.oxygen.boost_max = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.q_consumption" || key == "oxygen_q_consumption") {
    cfg.chem_env.oxygen.q_consumption = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.vbf_sink" || key == "oxygen_vbf_sink") {
    cfg.chem_env.oxygen.vbf_sink = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.k_ROS" || key == "oxygen_k_ROS") {
    cfg.chem_env.oxygen.k_ROS = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_acetate_dyn_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "acetate.enabled" || key == "acetate_enabled") {
    cfg.chem_env.acetate.enabled = parse_bool_config(val); return true;
  }
  if (key == "acetate.D_free" || key == "acetate_D_free") {
    cfg.chem_env.acetate.D_free = parse_config_real(key, val); return true;
  }
  if (key == "acetate.vbf_production" || key == "acetate_vbf_production") {
    cfg.chem_env.acetate.vbf_production = parse_config_real(key, val); return true;
  }
  if (key == "acetate.vbf_consumption" || key == "acetate_vbf_consumption") {
    cfg.chem_env.acetate.vbf_consumption = parse_config_real(key, val); return true;
  }
  if (key == "acetate.overflow_threshold" || key == "acetate_overflow_threshold") {
    cfg.chem_env.acetate.overflow_threshold = parse_config_real(key, val); return true;
  }
  if (key == "acetate.overflow_rate" || key == "acetate_overflow_rate") {
    cfg.chem_env.acetate.overflow_rate = parse_config_real(key, val); return true;
  }
  if (key == "acetate.scavenge_Km" || key == "acetate_scavenge_Km") {
    cfg.chem_env.acetate.scavenge_Km = parse_config_real(key, val); return true;
  }
  if (key == "acetate.scavenge_rate" || key == "acetate_scavenge_rate") {
    cfg.chem_env.acetate.scavenge_rate = parse_config_real(key, val); return true;
  }
  if (key == "acetate.epithelial_uptake" || key == "acetate_epithelial_uptake") {
    cfg.chem_env.acetate.epithelial_uptake = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_mucin_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "mucin.enabled" || key == "mucin_enabled") {
    cfg.chem_env.mucin.enabled = parse_bool_config(val); return true;
  }
  if (key == "mucin.secretion_rate" || key == "mucin_secretion_rate") {
    cfg.chem_env.mucin.secretion_rate = parse_config_real(key, val); return true;
  }
  if (key == "mucin.Km_degradation" || key == "mucin_Km_degradation") {
    cfg.chem_env.mucin.Km_degradation = parse_config_real(key, val); return true;
  }
  if (key == "mucin.k_liberation" || key == "mucin_k_liberation") {
    cfg.chem_env.mucin.k_liberation = parse_config_real(key, val); return true;
  }
  if (key == "mucin.initial_conc" || key == "mucin_initial_conc") {
    cfg.chem_env.mucin.initial_conc = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_protease_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "protease.enabled" || key == "protease_enabled") {
    cfg.chem_env.protease.enabled = parse_bool_config(val); return true;
  }
  if (key == "protease.default_half_life" || key == "protease_default_half_life") {
    cfg.chem_env.protease.default_half_life = parse_config_real(key, val); return true;
  }
  if (key == "protease.dilution_rate" || key == "protease_dilution_rate") {
    cfg.chem_env.protease.dilution_rate = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_fur_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "fur.enabled" || key == "fur_enabled") {
    cfg.cell_bio.fur.enabled = parse_bool_config(val); return true;
  }
  if (key == "fur.Km" || key == "fur_Km") {
    cfg.cell_bio.fur.Km = parse_config_real(key, val); return true;
  }
  if (key == "fur.upregulation_max" || key == "fur_upregulation_max") {
    cfg.cell_bio.fur.upregulation_max = parse_config_real(key, val); return true;
  }
  if (key == "fur.receptor_max" || key == "fur_receptor_max") {
    cfg.cell_bio.fur.receptor_max = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_cdi_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "cdi.enabled" || key == "cdi_enabled") {
    cfg.cell_bio.cdi.enabled = parse_bool_config(val); return true;
  }
  if (key == "cdi.kill_rate" || key == "cdi_kill_rate") {
    cfg.cell_bio.cdi.kill_rate = parse_config_real(key, val); return true;
  }
  if (key == "cdi.contact_radius" || key == "cdi_contact_radius") {
    cfg.cell_bio.cdi.contact_radius = parse_config_real(key, val); return true;
  }
  if (key == "cdi.corpse_persistence" || key == "cdi_corpse_persistence") {
    cfg.cell_bio.cdi.corpse_persistence = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_motility_key(SimulationConfig& cfg, const std::string& key, const std::string& val) {
  if (key == "motility.enabled" || key == "motility_enabled") {
    cfg.cell_bio.motility.enabled = parse_bool_config(val); return true;
  }
  if (key == "motility.swim_speed" || key == "motility_swim_speed") {
    cfg.cell_bio.motility.swim_speed = parse_config_real(key, val); return true;
  }
  if (key == "motility.run_mean_duration" || key == "motility_run_mean_duration") {
    cfg.cell_bio.motility.run_mean_duration = parse_config_real(key, val); return true;
  }
  if (key == "motility.stop_probability" || key == "motility_stop_probability") {
    cfg.cell_bio.motility.stop_probability = parse_config_real(key, val); return true;
  }
  if (key == "motility.stop_duration" || key == "motility_stop_duration") {
    cfg.cell_bio.motility.stop_duration = parse_config_real(key, val); return true;
  }
  if (key == "motility.chemotaxis_enabled" || key == "motility_chemotaxis_enabled") {
    cfg.cell_bio.motility.chemotaxis_enabled = parse_bool_config(val); return true;
  }
  if (key == "motility.chi_carbon" || key == "motility_chi_carbon") {
    cfg.cell_bio.motility.chi_carbon = parse_config_real(key, val); return true;
  }
  if (key == "motility.chi_oxygen" || key == "motility_chi_oxygen") {
    cfg.cell_bio.motility.chi_oxygen = parse_config_real(key, val); return true;
  }
  if (key == "motility.cluster_suppress_radius" || key == "motility_cluster_suppress_radius") {
    cfg.cell_bio.motility.cluster_suppress_radius = parse_config_real(key, val); return true;
  }
  if (key == "motility.cluster_suppress_threshold" || key == "motility_cluster_suppress_threshold") {
    cfg.cell_bio.motility.cluster_suppress_threshold = parse_config_int(key, val); return true;
  }
  if (key == "motility.cluster_tumble_factor" || key == "motility_cluster_tumble_factor") {
    cfg.cell_bio.motility.cluster_tumble_factor = parse_config_real(key, val); return true;
  }
  return false;
}

constexpr std::array<FlatKeyHandler, 23> k_flat_key_handlers = {
  apply_time_key,
  apply_domain_key,
  apply_advection_key,
  apply_qssa_key,
  apply_vbf_key,
  apply_chemical_key,
  apply_metabolism_key,
  apply_receptor_key,
  apply_conjugation_key,
  apply_mutation_key,
  apply_mechanics_key,
  apply_io_key,
  apply_hdf5_key,
  apply_dt_key,
  apply_misc_key,
  apply_oxygen_key,
  apply_acetate_dyn_key,
  apply_mucin_key,
  apply_protease_key,
  apply_siderophore_key,
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

  std::string key = trim_config(std::string_view(trimmed).substr(0, colon));
  std::string val = trim_config(std::string_view(trimmed).substr(colon + 1));
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
    if (!InputParser::apply_flat_key(cfg, key, val)) {
      std::cerr << "Warning: unknown config key '" << key << "' ignored\n";
    }
  }
}

}  // namespace

bool InputParser::apply_flat_key(SimulationConfig& cfg,
                                 const std::string& key,
                                 const std::string& val) {
  for (FlatKeyHandler handler : k_flat_key_handlers) {
    if (handler(cfg, key, val)) return true;
  }
  return false;
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
