/* -----------------------------------------------------------------------
   GutIBM – Input parser implementation
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "config_json.h"
#include "plasmid.h"
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <sstream>
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
    throw std::runtime_error("invalid config value for key '" + key + "'");
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

  return cfg;
}

void InputParser::apply_flat_key(SimulationConfig& cfg,
                                 const std::string& key,
                                 const std::string& val) {
  if (key == "total_time")             cfg.total_time = parse_real(key, val);
  else if (key == "bio_dt")            cfg.bio_dt = parse_real(key, val);
  else if (key == "output_interval")   cfg.output_interval = parse_real(key, val);
  else if (key == "seed")              cfg.seed = static_cast<uint64_t>(parse_int(key, val));
  else if (key == "grid_dx")           cfg.domain.grid_dx = parse_real(key, val);
  else if (key == "domain_x")          cfg.domain.hi[0] = parse_real(key, val);
  else if (key == "domain_y")          cfg.domain.hi[1] = parse_real(key, val);
  else if (key == "domain_z")          cfg.domain.hi[2] = parse_real(key, val);
  else if (key == "mucus_thickness")   cfg.advection.mucus_thickness = parse_real(key, val);
  else if (key == "radial_turnover")   cfg.advection.radial_turnover = parse_real(key, val);
  else if (key == "distal_transit")    cfg.advection.distal_transit_time = parse_real(key, val);
  else if (key == "peristaltic_enabled")   cfg.advection.peristaltic_enabled = (val == "true" || val == "1");
  else if (key == "peristaltic_period")    cfg.advection.peristaltic_period = parse_real(key, val);
  else if (key == "peristaltic_amplitude") cfg.advection.peristaltic_amplitude = parse_real(key, val);
  else if (key == "peristaltic_wavelength") cfg.advection.peristaltic_wavelength = parse_real(key, val);
  else if (key == "toxin_cutoff")      cfg.qssa.toxin_cutoff = parse_real(key, val);
  else if (key == "nutrient_cutoff")   cfg.qssa.nutrient_cutoff = parse_real(key, val);
  else if (key == "use_fmm")           cfg.qssa.use_fmm = (val == "true" || val == "1");
  else if (key == "fmm_theta")         cfg.qssa.fmm_theta = parse_real(key, val);
  else if (key == "fmm_expansion_order") cfg.qssa.fmm_expansion_order = parse_int(key, val);
  else if (key == "vbf_density")       cfg.vbf.density = parse_real(key, val);
  else if (key == "vbf_viscosity")     cfg.vbf.viscosity = parse_real(key, val);
  else if (key == "vbf_mucin_z_gradient")  cfg.vbf.mucin_z_gradient_enabled = (val == "true" || val == "1");
  else if (key == "vbf_mucin_z_lambda")    cfg.vbf.mucin_z_gradient_lambda = parse_real(key, val);
  else if (key == "carbon_z_gradient")     {
    for (auto& c : cfg.chemicals) {
      if (c.name == "carbon") { c.z_gradient_enabled = (val == "true" || val == "1"); break; }
    }
  }
  else if (key == "carbon_z_lambda")       {
    for (auto& c : cfg.chemicals) {
      if (c.name == "carbon") { c.z_gradient_lambda = parse_real(key, val); break; }
    }
  }
  else if (key == "sos_lysis_prob")    cfg.bacteriocin.sos_lysis_prob = parse_real(key, val);
  // Receptor Fix tunables
  else if (key == "kd_b12_btuB")        cfg.receptor.kd_b12_btuB = parse_real(key, val);
  else if (key == "kd_colicinE_btuB")  cfg.receptor.kd_colicinE_btuB = parse_real(key, val);
  else if (key == "kd_enterobactin")    cfg.receptor.kd_enterobactin = parse_real(key, val);
  else if (key == "kd_colicinB_fepA")   cfg.receptor.kd_colicinB_fepA = parse_real(key, val);
  else if (key == "kd_lin_enterobactin") cfg.receptor.kd_lin_enterobactin = parse_real(key, val);
  else if (key == "kd_colicinIa_cirA") cfg.receptor.kd_colicinIa_cirA = parse_real(key, val);
  else if (key == "kill_rate_colicin")  cfg.receptor.kill_rate_colicin = parse_real(key, val);
  else if (key == "kill_rate_microcin") cfg.receptor.kill_rate_microcin = parse_real(key, val);
  else if (key == "immunity_factor")    cfg.receptor.immunity_factor = parse_real(key, val);
  // Conjugation Fix tunables
  else if (key == "pili_length")        cfg.conjugation.pili_length = parse_real(key, val);
  else if (key == "base_transfer_rate") cfg.conjugation.base_transfer_rate = parse_real(key, val);
  else if (key == "shear_critical")     cfg.conjugation.shear_critical = parse_real(key, val);
  else if (key == "plasmid_copy_cost")  cfg.conjugation.plasmid_copy_cost = parse_real(key, val);
  else if (key == "pili_heterogeneity") cfg.conjugation.pili_heterogeneity = (val == "true" || val == "1");
  else if (key == "pili_length_min")    cfg.conjugation.pili_length_min = parse_real(key, val);
  else if (key == "pili_length_max")    cfg.conjugation.pili_length_max = parse_real(key, val);
  // Mutation Fix tunables
  else if (key == "bi_duplication_rate")    cfg.mutation.bi_duplication_rate = parse_real(key, val);
  else if (key == "bi_recombination_rate")  cfg.mutation.bi_recombination_rate = parse_real(key, val);
  else if (key == "receptor_mutation_rate") cfg.mutation.receptor_mutation_rate = parse_real(key, val);
  else if (key == "super_killer_rate")      cfg.mutation.super_killer_rate = parse_real(key, val);
  else if (key == "partial_resistance_rate") cfg.mutation.partial_resistance_rate = parse_real(key, val);
  else if (key == "receptor_reduction")     cfg.mutation.receptor_reduction = parse_real(key, val);
  else if (key == "max_bi_loci")            cfg.mutation.max_bi_loci = parse_int(key, val);
  else if (key == "immunity_escape_prob")   cfg.mutation.immunity_escape_prob = parse_real(key, val);
  else if (key == "escape_affinity_lo")     cfg.mutation.escape_affinity_lo = parse_real(key, val);
  else if (key == "escape_affinity_hi")     cfg.mutation.escape_affinity_hi = parse_real(key, val);
  else if (key == "compensatory_rate")      cfg.mutation.compensatory_rate = parse_real(key, val);
  else if (key == "compensatory_reduction") cfg.mutation.compensatory_reduction = parse_real(key, val);
  else if (key == "crypts_enabled")     cfg.advection.crypts_enabled = (val == "true" || val == "1");
  else if (key == "crypt_depth")       cfg.advection.crypt_depth = parse_real(key, val);
  else if (key == "crypt_exit_rate")   cfg.advection.crypt_exit_rate = parse_real(key, val);
  else if (key == "crypt_entry_rate")  cfg.advection.crypt_entry_rate = parse_real(key, val);
  else if (key == "crypt_carrying_capacity") cfg.advection.crypt_carrying_capacity = parse_int(key, val);
  else if (key == "hdf5_file")         cfg.hdf5.filename = val;
  else if (key == "hdf5_every")        cfg.hdf5.dump_every = parse_int(key, val);
  else if (key == "checkpoint_file")   cfg.checkpoint.file = val;
  else if (key == "checkpoint_step")   cfg.checkpoint.step = val;
  else if (key == "adaptive_dt_enabled") cfg.adaptive_dt_enabled = (val == "true" || val == "1");
  else if (key == "dt_min")            cfg.dt_min = parse_real(key, val);
  else if (key == "dt_max")            cfg.dt_max = parse_real(key, val);
  else if (key == "dt_safety")         cfg.dt_safety = parse_real(key, val);
  else if (key == "dt_growth_limit")   cfg.dt_growth_limit = parse_real(key, val);
  else if (key == "gpu_enabled")       cfg.gpu.enabled = (val == "true" || val == "1");
  else if (key == "gpu_device_id")     cfg.gpu.device_id = parse_int(key, val);
  else if (key == "profile_steps")   cfg.profile_steps = (val == "true" || val == "1");
}

SimulationConfig InputParser::parse(const std::string& filename) {
  SimulationConfig cfg = default_config();

  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    std::cerr << "Warning: cannot open config file '" << filename
              << "', using defaults\n";
    return cfg;
  }

  std::ostringstream content_stream;
  content_stream << ifs.rdbuf();
  const std::string content = content_stream.str();

  if (ConfigJson::parse_document(cfg, content)) {
    return cfg;
  }

  // Legacy line-oriented fallback for non-JSON configs.
  std::istringstream lines(content);
  std::string line;
  while (std::getline(lines, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#' || line[0] == '/' || line[0] == '{' || line[0] == '}')
      continue;

    auto colon = line.find(':');
    if (colon == std::string::npos) continue;

    std::string key = trim(line.substr(0, colon));
    std::string val = trim(line.substr(colon + 1));

    if (!key.empty() && key.front() == '_') continue;

    if (key.size() >= 2 && key.front() == '"' && key.back() == '"') {
      key = key.substr(1, key.size() - 2);
      if (!key.empty() && key.front() == '_') continue;
    }

    if (!val.empty() && val.back() == ',') val.pop_back();
    val = trim(val);
    if (val.size() >= 2 && val.front() == '"' && val.back() == '"') {
      val = val.substr(1, val.size() - 2);
    }

    apply_flat_key(cfg, key, val);
  }

  InitialStrainsParseResult strains = ConfigJson::parse_initial_strains(content);
  if (strains.found) {
    cfg.initial_strains = std::move(strains.strains);
  }

  EnabledFixesParseResult fixes = ConfigJson::parse_enabled_fixes(content);
  if (fixes.found) {
    cfg.enabled_fixes = std::move(fixes.names);
  }

  return cfg;
}

std::string InputParser::trim(const std::string& s) {
  auto start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) return "";
  auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end - start + 1);
}

Real InputParser::parse_real(const std::string& key, const std::string& val) {
  const std::string trimmed = trim(val);
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
  } catch (...) {
    warn_parse_failure("numeric", key, val);
    return 0.0;
  }
}

Int InputParser::parse_int(const std::string& key, const std::string& val) {
  const std::string trimmed = trim(val);
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
  } catch (...) {
    warn_parse_failure("integer", key, val);
    return 0;
  }
}

}  // namespace gutibm
