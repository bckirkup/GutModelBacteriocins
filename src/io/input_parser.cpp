/* -----------------------------------------------------------------------
   GutIBM – Input parser implementation
   ----------------------------------------------------------------------- */

#include "input_parser.h"
#include "plasmid.h"
#include <iostream>
#include <algorithm>
#include <stdexcept>

namespace gutibm {

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

SimulationConfig InputParser::parse(const std::string& filename) {
  SimulationConfig cfg = default_config();

  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    std::cerr << "Warning: cannot open config file '" << filename
              << "', using defaults\n";
    return cfg;
  }

  std::string line;
  while (std::getline(ifs, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#' || line[0] == '/' || line[0] == '{' || line[0] == '}')
      continue;

    // Simple key: value parsing
    auto colon = line.find(':');
    if (colon == std::string::npos) continue;

    std::string key = trim(line.substr(0, colon));
    std::string val = trim(line.substr(colon + 1));

    // Remove trailing comma and quotes
    if (!val.empty() && val.back() == ',') val.pop_back();
    val = trim(val);
    if (val.size() >= 2 && val.front() == '"' && val.back() == '"') {
      val = val.substr(1, val.size() - 2);
    }

    // Map keys to config
    if (key == "total_time")             cfg.total_time = parse_real(val);
    else if (key == "bio_dt")            cfg.bio_dt = parse_real(val);
    else if (key == "output_interval")   cfg.output_interval = parse_real(val);
    else if (key == "seed")              cfg.seed = static_cast<uint64_t>(parse_int(val));
    else if (key == "grid_dx")           cfg.domain.grid_dx = parse_real(val);
    else if (key == "domain_x")          cfg.domain.hi[0] = parse_real(val);
    else if (key == "domain_y")          cfg.domain.hi[1] = parse_real(val);
    else if (key == "domain_z")          cfg.domain.hi[2] = parse_real(val);
    else if (key == "mucus_thickness")   cfg.advection.mucus_thickness = parse_real(val);
    else if (key == "radial_turnover")   cfg.advection.radial_turnover = parse_real(val);
    else if (key == "distal_transit")    cfg.advection.distal_transit_time = parse_real(val);
    else if (key == "vbf_density")       cfg.vbf.density = parse_real(val);
    else if (key == "vbf_viscosity")     cfg.vbf.viscosity = parse_real(val);
    else if (key == "vbf_mucin_z_gradient")  cfg.vbf.mucin_z_gradient_enabled = (val == "true" || val == "1");
    else if (key == "vbf_mucin_z_lambda")    cfg.vbf.mucin_z_gradient_lambda = parse_real(val);
    else if (key == "carbon_z_gradient")     {
      for (auto& c : cfg.chemicals) {
        if (c.name == "carbon") { c.z_gradient_enabled = (val == "true" || val == "1"); break; }
      }
    }
    else if (key == "carbon_z_lambda")       {
      for (auto& c : cfg.chemicals) {
        if (c.name == "carbon") { c.z_gradient_lambda = parse_real(val); break; }
      }
    }
    else if (key == "sos_lysis_prob")    cfg.bacteriocin.sos_lysis_prob = parse_real(val);
    else if (key == "crypts_enabled")     cfg.advection.crypts_enabled = (val == "true" || val == "1");
    else if (key == "crypt_depth")       cfg.advection.crypt_depth = parse_real(val);
    else if (key == "crypt_exit_rate")   cfg.advection.crypt_exit_rate = parse_real(val);
    else if (key == "crypt_entry_rate")  cfg.advection.crypt_entry_rate = parse_real(val);
    else if (key == "crypt_carrying_capacity") cfg.advection.crypt_carrying_capacity = parse_int(val);
    else if (key == "hdf5_file")         cfg.hdf5.filename = val;
    else if (key == "hdf5_every")        cfg.hdf5.dump_every = parse_int(val);
  }

  return cfg;
}

std::string InputParser::trim(const std::string& s) {
  auto start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) return "";
  auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end - start + 1);
}

Real InputParser::parse_real(const std::string& val) {
  try { return std::stod(val); }
  catch (...) { return 0.0; }
}

Int InputParser::parse_int(const std::string& val) {
  try { return std::stoi(val); }
  catch (...) { return 0; }
}

}  // namespace gutibm
