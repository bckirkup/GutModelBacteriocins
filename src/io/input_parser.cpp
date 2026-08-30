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
#include <ranges>
#include "error.h"
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <cstdlib>
#include <array>
#include <cctype>
#include <cmath>
#include <format>
#include <limits>

namespace gutibm {

namespace {

bool strict_config_enabled() {
  const char* env = std::getenv("GUTIBM_STRICT_CONFIG");
  if (!env || env[0] == '\0') return true;
  return env[0] != '0';
}

bool strict_config_requested() {
  const char* env = std::getenv("GUTIBM_STRICT_CONFIG");
  return env != nullptr && env[0] != '\0' && env[0] != '0';
}

void warn_parse_failure(const char* kind,
                        std::string_view key,
                        const std::string& val) {
  if (strict_config_enabled()) {
    throw ConfigError(std::string("invalid config value for key '")
                      + std::string(key) + "'");
  }
  std::cerr << "Warning: config key '" << key << "' has invalid " << kind
            << " value '" << val << "' — using 0\n";
}

void warn_unknown_config_key(std::string_view key) {
  InputParser::handle_unknown_config_key(key);
}

std::string trim_config(std::string_view s) {
  auto start = s.find_first_not_of(" \t\r\n");
  if (start == std::string_view::npos) return "";
  auto end = s.find_last_not_of(" \t\r\n");
  return std::string(s.substr(start, end - start + 1));
}

Real parse_config_real(std::string_view key, const std::string& val) {
  const std::string trimmed = trim_config(val);
  if (trimmed.empty()) {
    warn_parse_failure("numeric", key, val);
    return 0.0;
  }
  if (trimmed == "inf" || trimmed == "+inf"
      || trimmed == "infinity" || trimmed == "+infinity") {
    return std::numeric_limits<Real>::infinity();
  }
  if (trimmed == "-inf" || trimmed == "-infinity") {
    return -std::numeric_limits<Real>::infinity();
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

Real parse_positive_config_real(std::string_view key, const std::string& val) {
  const Real result = parse_config_real(key, val);
  if (!std::isfinite(result) || result <= 0.0) {
    throw ConfigError("config key '" + std::string(key)
                      + "' must be a finite positive real");
  }
  return result;
}

Int parse_config_int(std::string_view key, const std::string& val) {
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

uint64_t parse_config_uint64(std::string_view key, const std::string& val) {
  const std::string trimmed = trim_config(val);
  if (trimmed.empty() || trimmed.front() == '-') {
    warn_parse_failure("unsigned integer", key, val);
    return 0;
  }

  try {
    size_t consumed = 0;
    const uint64_t result = std::stoull(trimmed, &consumed);
    if (consumed != trimmed.size()) {
      warn_parse_failure("unsigned integer", key, val);
      return 0;
    }
    return result;
  } catch (const std::invalid_argument&) {
    warn_parse_failure("unsigned integer", key, val);
    return 0;
  } catch (const std::out_of_range&) {
    warn_parse_failure("unsigned integer", key, val);
    return 0;
  }
}

Int parse_positive_config_int(std::string_view key, const std::string& val) {
  const Int result = parse_config_int(key, val);
  if (result < 1) {
    throw ConfigError("config key '" + std::string(key)
                      + "' must be an integer at least 1");
  }
  return result;
}

EpithelialBoundaryMode parse_epithelial_boundary_mode(
    std::string_view key, const std::string& value) {
  using enum EpithelialBoundaryMode;
  if (value == "dirichlet") return Dirichlet;
  if (value == "robin") return Robin;
  if (value == "flux") return Flux;
  throw ConfigError(
      "invalid " + std::string(key)
      + ": expected 'dirichlet', 'robin', or 'flux', got '" + value + "'");
}

void configure_epithelial_boundary(
    ChemicalSpec& spec, std::string_view species_name,
    std::string_view boundary_key, const std::string& boundary,
    Real transfer_coeff, Real flux) {
  spec.epithelial_transfer_coeff = transfer_coeff;
  spec.epithelial_flux = flux;
  spec.epithelial_boundary_mode =
      parse_epithelial_boundary_mode(boundary_key, boundary);
  if (spec.epithelial_boundary_mode == EpithelialBoundaryMode::Robin
      && transfer_coeff <= 0.0) {
    throw ConfigError(
        std::string(species_name) + ".epithelial_transfer_coeff must be "
        "positive for Robin epithelial boundary mode");
  }
  if (spec.epithelial_boundary_mode != EpithelialBoundaryMode::Dirichlet
      && spec.z_gradient_enabled) {
    throw ConfigError(
        std::string(species_name)
        + " z-gradient cannot be combined with Robin or flux "
          "epithelial boundary modes");
  }
}

}  // namespace

SimulationConfig InputParser::default_config() {
  SimulationConfig cfg;

  // Default chemical species. Nutrients and small molecules use stable
  // implicit diffusion; bacteriocins remain on the QSSA Green's-function path.
  cfg.chemicals = {
    {species::CARBON,      5.0e-10, 1.0, 5.0e-3, cfg.carbon_boundary_conc, 0.0, true,  25.0e-6, true},
    {species::IRON,        7.0e-10, 1.0, 1.0e-4, 1.0e-4, 0.0, false, 25.0e-6, true},
    {species::B12,         5.0e-10, 1.0, cfg.b12_initial_conc, cfg.b12_initial_conc, 0.0, false, 25.0e-6, true},
    {species::BACTERIOCIN_BTUB, 4.0e-11, 10.0, 0.0, 0.0, 1.0e-4, false, 25.0e-6, false},
    {species::BACTERIOCIN_FEPA, 4.0e-11, 10.0, 0.0, 0.0, 1.0e-4, false, 25.0e-6, false},
    {species::BACTERIOCIN_CIRA, 4.0e-11, 10.0, 0.0, 0.0, 1.0e-4, false, 25.0e-6, false},
    {species::BACTERIOCIN_FHUA, 4.0e-11, 10.0, 0.0, 0.0, 1.0e-4, false, 25.0e-6, false},
    {species::ACETATE,     1.2e-9, 1.0, 80.0,   80.0,   0.0, false, 25.0e-6, true},
    {species::ETHANOLAMINE, 1.0e-9, 1.0, 0.5e-3, 0.5e-3, 0.0, false, 25.0e-6, true},
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
  cfg.hdf5.schedule.provenance = 0;

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

bool matching_toxin_transport(const ChemicalSpec& lhs,
                              const ChemicalSpec& rhs) {
  // Exact equality is intentional: lumping one field requires identical
  // transport and decay coefficients, not merely numerically close values.
  return lhs.diff_coeff == rhs.diff_coeff
      && lhs.retardation == rhs.retardation
      && lhs.decay_rate == rhs.decay_rate;
}

void configure_toxin_species(SimulationConfig& cfg) {
  constexpr std::array<const char*, 4> receptor_names = {
      species::BACTERIOCIN_BTUB, species::BACTERIOCIN_FEPA,
      species::BACTERIOCIN_CIRA, species::BACTERIOCIN_FHUA};
  if (cfg.qssa.toxin_lumping != "lumped") return;

  std::array<Int, receptor_names.size()> receptor_indices{};
  Int receptor_count = 0;
  for (size_t i = 0; i < receptor_names.size(); ++i) {
    receptor_indices[i] =
        find_chemical_spec(cfg.chemicals, receptor_names[i]);
    if (receptor_indices[i] >= 0) ++receptor_count;
  }
  if (const Int lumped_idx =
          find_chemical_spec(cfg.chemicals, species::BACTERIOCIN_LUMPED);
      lumped_idx >= 0) {
    if (receptor_count != 0) {
      throw ConfigError(
          "toxin_lumping=lumped cannot combine lumped and receptor fields");
    }
    return;
  }
  if (receptor_count != static_cast<Int>(receptor_names.size())) {
    throw ConfigError(
        "toxin_lumping=lumped requires the four bacteriocin specifications");
  }

  const ChemicalSpec& reference_spec =
      cfg.chemicals[static_cast<size_t>(receptor_indices[0])];
  for (size_t i = 1; i < receptor_indices.size(); ++i) {
    const auto& spec =
        cfg.chemicals[static_cast<size_t>(receptor_indices[i])];
    if (!matching_toxin_transport(spec, reference_spec)) {
      throw ConfigError(
          "toxin_lumping=lumped requires matching diff_coeff, retardation, "
          "and decay_rate for all bacteriocin specifications");
    }
  }

  ChemicalSpec lumped = reference_spec;
  lumped.name = species::BACTERIOCIN_LUMPED;
  std::vector<ChemicalSpec> configured;
  configured.reserve(cfg.chemicals.size() - receptor_names.size() + 1);
  bool inserted = false;
  for (const auto& spec : cfg.chemicals) {
    if (const bool is_receptor =
            std::ranges::find(receptor_names, spec.name)
            != std::end(receptor_names);
        is_receptor) {
      if (!inserted) {
        configured.push_back(lumped);
        inserted = true;
      }
      continue;
    }
    configured.push_back(spec);
  }
  cfg.chemicals = std::move(configured);
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

bool is_toxin_species(std::string_view name) {
  return name == species::BACTERIOCIN_BTUB
      || name == species::BACTERIOCIN_FEPA
      || name == species::BACTERIOCIN_CIRA
      || name == species::BACTERIOCIN_FHUA
      || name == species::BACTERIOCIN_LUMPED;
}

void disable_subset_mechanism(SimulationConfig& cfg,
                              std::string_view fix_name,
                              std::string_view audit_name) {
  if (!fix_name.empty()) {
    cfg.disabled_fixes.emplace_back(fix_name);
  }
  cfg.disabled_mechanisms.emplace_back(audit_name);
}

void apply_species_subset(SimulationConfig& cfg) {
  if (cfg.species_subset == "full") {
    cfg.disabled_fixes.clear();
    cfg.disabled_mechanisms.clear();
    return;
  }
  if (cfg.qssa.toxin_lumping == "lumped") {
    throw ConfigError(
        "chemistry.species_subset cannot drop toxin fields while "
        "toxin_lumping=lumped");
  }

  cfg.disabled_fixes.clear();
  cfg.disabled_mechanisms.clear();
  disable_subset_mechanism(cfg, "receptor", "receptor_killing");
  disable_subset_mechanism(cfg, "bacteriocin", "bacteriocin_path");
  disable_subset_mechanism(cfg, "", "qssa_toxin_path");
  if (cfg.species_subset == "carbon_only") {
    disable_subset_mechanism(cfg, "", "oxygen");
    disable_subset_mechanism(cfg, "", "acetate");
    disable_subset_mechanism(cfg, "", "ethanolamine");
    disable_subset_mechanism(cfg, "", "mucin");
    disable_subset_mechanism(cfg, "", "siderophore");
    disable_subset_mechanism(cfg, "", "ferrichrome");
    disable_subset_mechanism(cfg, "", "quorum_sensing");
    disable_subset_mechanism(cfg, "", "motility_aerotaxis");
    disable_subset_mechanism(cfg, "", "motility_mucin_drag");
    disable_subset_mechanism(cfg, "", "motility_ai2_chemotaxis");
    disable_subset_mechanism(cfg, "", "fur_regulation");
    disable_subset_mechanism(cfg, "", "iron_uptake");
    disable_subset_mechanism(cfg, "", "b12_uptake");
    disable_subset_mechanism(cfg, "", "eut_ethanolamine_uptake");
    disable_subset_mechanism(cfg, "", "vbf_iron_sink");
    disable_subset_mechanism(cfg, "", "vbf_dynamic_mucin");
    cfg.chem_env.oxygen.enabled = false;
    cfg.chem_env.acetate.enabled = false;
    cfg.chem_env.mucin.enabled = false;
    cfg.chem_env.siderophore.enabled = false;
    cfg.chem_env.ferrichrome.enabled = false;
    cfg.quorum_sensing.enabled = false;
    cfg.quorum_sensing.ai2_chemotaxis_enabled = false;
    cfg.cell_bio.fur.enabled = false;
    cfg.cell_bio.motility.aerotaxis_enabled = false;
    cfg.cell_bio.motility.mucin_drag_enabled = false;
    cfg.fixes.metabolism.iron_uptake_enabled = false;
    cfg.fixes.metabolism.b12_uptake_enabled = false;
    cfg.fixes.metabolism.eut_enabled = false;
    cfg.vbf.nutrient_sink = 0.0;
    cfg.vbf.use_dynamic_mucin = false;
  }

  std::vector<ChemicalSpec> retained;
  retained.reserve(cfg.chemicals.size());
  for (const auto& spec : cfg.chemicals) {
    const bool keep = cfg.species_subset == "nutrient_only"
        ? !is_toxin_species(spec.name)
        : spec.name == species::CARBON;
    if (keep) retained.push_back(spec);
  }
  cfg.chemicals = std::move(retained);
}

}  // namespace

void InputParser::finalize_config(SimulationConfig& cfg) {
  constexpr Real k_z_lambda = 25.0e-6;

  if (cfg.qssa.low_screening_policy != "warn"
      && cfg.qssa.low_screening_policy != "error"
      && cfg.qssa.low_screening_policy != "allow") {
    throw ConfigError(
        "invalid qssa.low_screening_policy: expected 'warn', 'error', "
        "or 'allow'");
  }

  if (cfg.closure.zero_realization_grace_steps < 0) {
    throw ConfigError(
        "closure.zero_realization_grace_steps must be nonnegative");
  }
  if (cfg.closure.reaction_clip_tolerance_fraction < 0.0) {
    throw ConfigError(
        "closure.reaction_clip_tolerance_fraction must be nonnegative");
  }

  if (cfg.initial_population.anatomic_exclusion_floor < 0.0) {
    throw ConfigError(
        "initial_population.anatomic_exclusion_floor must be nonnegative");
  }
  if (cfg.initial_population.anatomic_exponential_scale <= 0.0) {
    throw ConfigError(
        "initial_population.anatomic_exponential_scale must be positive");
  }
  if (cfg.initial_population.anatomic_outer_extent
      <= cfg.initial_population.anatomic_exclusion_floor) {
    throw ConfigError(
        "initial_population.anatomic_outer_extent must be greater than "
        "anatomic_exclusion_floor");
  }

  if (cfg.initial_population.placement == "anatomic") {
    if (cfg.initial_population.anatomic_exclusion_floor < cfg.domain.lo[2]
        || cfg.initial_population.anatomic_exclusion_floor >= cfg.domain.hi[2]) {
      throw ConfigError(
          "initial_population.anatomic_exclusion_floor must be inside the "
          "domain");
    }
    if (cfg.initial_population.anatomic_outer_extent <= cfg.domain.lo[2]
        || cfg.initial_population.anatomic_outer_extent > cfg.domain.hi[2]) {
      throw ConfigError(
          "initial_population.anatomic_outer_extent must be inside the "
          "domain");
    }
  }

  if (cfg.initial_population.placement == "z_slab") {
    if (cfg.initial_population.z_min < cfg.domain.lo[2]
        || cfg.initial_population.z_min >= cfg.domain.hi[2]) {
      throw ConfigError(
          "initial_population.z_min must be inside the domain");
    }
    if (cfg.initial_population.z_max <= cfg.domain.lo[2]
        || cfg.initial_population.z_max > cfg.domain.hi[2]) {
      throw ConfigError(
          "initial_population.z_max must be inside the domain");
    }
    if (cfg.initial_population.z_min >= cfg.initial_population.z_max) {
      throw ConfigError(
          "initial_population.z_min must be less than z_max");
    }
  }

  apply_species_subset(cfg);

  if (cfg.chem_env.oxygen.enabled) {
    const Int idx = find_chemical_spec(cfg.chemicals, species::OXYGEN);
    if (idx < 0) {
      cfg.chemicals.emplace_back(
          species::OXYGEN, cfg.chem_env.oxygen.D_free, 1.0,
          cfg.chem_env.oxygen.epithelial_conc, cfg.chem_env.oxygen.epithelial_conc,
          0.0, true, k_z_lambda, true);
      cfg.chemicals.back().z_gradient_enabled =
          cfg.oxygen_z_gradient_enabled;
    } else {
      auto& spec = cfg.chemicals[static_cast<size_t>(idx)];
      spec.diff_coeff = cfg.chem_env.oxygen.D_free;
      spec.diffusion_enabled = true;
      spec.z_gradient_enabled = cfg.oxygen_z_gradient_enabled;
    }
  }

  if (cfg.chem_env.acetate.enabled) {
    Int idx = find_chemical_spec(cfg.chemicals, species::ACETATE);
    if (idx < 0) {
      cfg.chemicals.emplace_back(
          species::ACETATE, cfg.chem_env.acetate.D_free, 1.0,
          0.0, 0.0, 0.0, false, k_z_lambda, true);
    } else {
      auto& spec = cfg.chemicals[static_cast<size_t>(idx)];
      spec.diff_coeff = cfg.chem_env.acetate.D_free;
      spec.initial_conc = 0.0;
      spec.boundary_conc = 0.0;
      spec.diffusion_enabled = true;
    }
  }

  if (cfg.chem_env.mucin.enabled) {
    if (find_chemical_spec(cfg.chemicals, species::MUCIN) < 0) {
      cfg.chemicals.emplace_back(
          species::MUCIN, cfg.chem_env.mucin.D_free, cfg.chem_env.mucin.retardation,
          cfg.chem_env.mucin.initial_conc, cfg.chem_env.mucin.initial_conc,
          0.0, false, k_z_lambda, false);
    }
    cfg.vbf.use_dynamic_mucin = true;
  }

  if (cfg.chem_env.siderophore.enabled) {
    if (find_chemical_spec(cfg.chemicals, species::SIDEROPHORE) < 0) {
      cfg.chemicals.emplace_back(
          species::SIDEROPHORE, cfg.chem_env.siderophore.D_free, 1.0,
          0.0, 0.0, 0.0, false, k_z_lambda, true);
    } else {
      auto& spec = cfg.chemicals[static_cast<size_t>(
          find_chemical_spec(cfg.chemicals, species::SIDEROPHORE))];
      spec.diff_coeff = cfg.chem_env.siderophore.D_free;
      spec.diffusion_enabled = true;
    }
  }

  if (cfg.chem_env.siderophore.enabled
      && find_chemical_spec(cfg.chemicals, species::FERRIC_ENTEROBACTIN) < 0) {
    cfg.chemicals.emplace_back(
        species::FERRIC_ENTEROBACTIN, cfg.chem_env.siderophore.D_free, 1.0,
        0.0, 0.0, 0.0, false, k_z_lambda, true);
  }

  if (cfg.chem_env.ferrichrome.enabled
      && find_chemical_spec(cfg.chemicals, species::FERRICHROME) < 0) {
    cfg.chemicals.emplace_back(
        species::FERRICHROME, 1.0e-10, 1.0,
        cfg.chem_env.ferrichrome.initial_conc,
        cfg.chem_env.ferrichrome.boundary_conc,
        0.0, false, k_z_lambda, true);
  }

  configure_toxin_species(cfg);

  const Real slab_height = cfg.domain.hi[2] - cfg.domain.lo[2];
  // This is the worst-case screening because both flow components vanish at
  // z_lo, so a source at the epithelial wall is unscreened by lumen flow.
  // Keep the parser guard on sqrt(decay_rate / D_eff) * H rather than using
  // distal or radial lumen velocities.
  for (const auto& spec : cfg.chemicals) {
    if (!is_toxin_species(spec.name) || spec.diffusion_enabled) continue;
    const Real d_eff = spec.diff_coeff / spec.retardation;
    const Real k_h = d_eff > 0.0 && spec.decay_rate > 0.0
        ? std::sqrt(spec.decay_rate / d_eff) * slab_height : 0.0;
    if (k_h >= 0.05) continue;
    std::cerr << "Warning: QSSA Green's-function species '" << spec.name
              << "' has negligible slab screening (kH=" << k_h
              << ", require kH >= 0.05); the bounded steady-state image "
                 "series may not converge\n";
    if (strict_config_requested()) {
      throw ConfigError(
          "QSSA Green's-function species '" + spec.name
          + "' has negligible slab screening (kH < 0.05)");
    }
  }

  // Spec 11 — AI-2 autoinducer (no z-gradient; agent-produced only)
  if (cfg.quorum_sensing.enabled) {
    const Int idx = find_chemical_spec(cfg.chemicals, species::AI2);
    if (idx < 0) {
      cfg.chemicals.emplace_back(
          species::AI2, cfg.quorum_sensing.ai2_D_free, 1.0,
          0.0, 0.0, cfg.quorum_sensing.ai2_decay_rate,
          false, k_z_lambda, true);
    } else {
      auto& spec = cfg.chemicals[static_cast<size_t>(idx)];
      spec.diff_coeff = cfg.quorum_sensing.ai2_D_free;
      spec.decay_rate = cfg.quorum_sensing.ai2_decay_rate;
      spec.diffusion_enabled = true;
      spec.z_gradient_enabled = false;
      spec.initial_conc = 0.0;
      spec.boundary_conc = 0.0;
    }
  }

  auto& metabolism = cfg.fixes.metabolism;
  if (metabolism.uptake_limit == "none") {
    metabolism.uptake_limit_mode = UptakeLimitMode::None;
  } else if (metabolism.uptake_limit == "sherwood") {
    metabolism.uptake_limit_mode = UptakeLimitMode::Sherwood;
  } else if (metabolism.uptake_limit == "voxel") {
    metabolism.uptake_limit_mode = UptakeLimitMode::Voxel;
  } else if (metabolism.uptake_limit == "delivery") {
    metabolism.uptake_limit_mode = UptakeLimitMode::Delivery;
  } else {
    throw ConfigError(
        "invalid uptake_limit: expected 'none', 'sherwood', 'voxel', or "
        "'delivery', got '"
        + metabolism.uptake_limit + "'");
  }
  if (metabolism.delivery_far_field_radius < 0.0) {
    throw ConfigError(
        "metabolism.delivery_far_field_radius must be non-negative");
  }
  if (cfg.chemistry_decomposition == "slab"
      && metabolism.delivery_far_field_radius > 0.0) {
    throw ConfigError(std::format(
        "metabolism.delivery_far_field_radius={:f}"
        " is unsupported for regularized delivery deposition in slab "
        "chemistry; slab configurations must set "
        "metabolism.delivery_far_field_radius = 0.0 to opt into the "
        "grid-dependent single-voxel delivery model",
        metabolism.delivery_far_field_radius));
  }
  if (cfg.chem_env.oxygen.respiration_driver != "ambient"
      && cfg.chem_env.oxygen.respiration_driver != "funded") {
    throw ConfigError(
        "invalid oxygen.respiration_driver: expected 'ambient' or "
        "'funded', got '"
        + cfg.chem_env.oxygen.respiration_driver + "'");
  }
  cfg.chem_env.oxygen.respiration_driver_mode =
      cfg.chem_env.oxygen.respiration_driver == "funded"
          ? RespirationDriver::Funded : RespirationDriver::Ambient;
  if (cfg.chem_env.oxygen.respiration_driver_mode == RespirationDriver::Funded
      && (!cfg.chem_env.oxygen.delivery_uptake_enabled
          || metabolism.uptake_limit_mode != UptakeLimitMode::Delivery)) {
    throw ConfigError(
        R"(oxygen.respiration_driver="funded" requires )"
        "oxygen.delivery_uptake_enabled=true and "
        "metabolism.uptake_limit=\"delivery\"");
  }
  if (cfg.chem_env.oxygen.ros_driver != "ambient"
      && cfg.chem_env.oxygen.ros_driver != "funded") {
    throw ConfigError(
        "invalid oxygen.ros_driver: expected 'ambient' or 'funded', got '"
        + cfg.chem_env.oxygen.ros_driver + "'");
  }
  cfg.chem_env.oxygen.ros_driver_mode =
      cfg.chem_env.oxygen.ros_driver == "funded"
          ? RosDriver::Funded : RosDriver::Ambient;
  if (cfg.chem_env.oxygen.ros_driver_mode == RosDriver::Funded
      && (!cfg.chem_env.oxygen.delivery_uptake_enabled
          || metabolism.uptake_limit_mode != UptakeLimitMode::Delivery)) {
    throw ConfigError(
        R"(oxygen.ros_driver="funded" requires )"
        "oxygen.delivery_uptake_enabled=true and "
        "metabolism.uptake_limit=\"delivery\"");
  }
  if (cfg.chem_env.oxygen.ros_driver_mode == RosDriver::Funded
      && (cfg.chem_env.oxygen.k_ROS_respiratory > 0.0)
          == (cfg.chem_env.oxygen.k_ROS_funded > 0.0)) {
    if (cfg.chem_env.oxygen.k_ROS_respiratory > 0.0) {
      throw ConfigError(
          "oxygen.ros_driver=\"funded\" requires exactly one of "
          "oxygen.k_ROS_respiratory or oxygen.k_ROS_funded to be positive; "
          "both are set and the normalization is ambiguous");
    }
    throw ConfigError(
        "oxygen.ros_driver=\"funded\" requires exactly one of "
        "oxygen.k_ROS_respiratory or oxygen.k_ROS_funded to be positive; "
        "neither key is positive");
  }
  if (cfg.chem_env.oxygen.delivery_uptake_enabled
      && metabolism.uptake_limit_mode != UptakeLimitMode::Delivery) {
    throw ConfigError(
        "oxygen.delivery_uptake_enabled=true requires "
        "metabolism.uptake_limit=\"delivery\"");
  }
  for (auto& spec : cfg.chemicals) {
    spec.delivery_enabled = false;
  }
  if (metabolism.uptake_limit_mode == UptakeLimitMode::Delivery) {
    if (const Int carbon =
            find_chemical_spec(cfg.chemicals, species::CARBON);
        carbon >= 0) {
      cfg.chemicals[static_cast<size_t>(carbon)].delivery_enabled = true;
    }
    if (cfg.chem_env.oxygen.delivery_uptake_enabled) {
      const Int oxygen = find_chemical_spec(cfg.chemicals, species::OXYGEN);
      if (oxygen >= 0) {
        cfg.chemicals[static_cast<size_t>(oxygen)].delivery_enabled = true;
      }
    }
  }
  if (const Int carbon_idx =
          find_chemical_spec(cfg.chemicals, species::CARBON);
      carbon_idx >= 0) {
    auto& carbon = cfg.chemicals[static_cast<size_t>(carbon_idx)];
    configure_epithelial_boundary(
        carbon, "carbon", "carbon.epithelial_boundary",
        cfg.carbon_epithelial_boundary,
        cfg.carbon_epithelial_transfer_coeff, cfg.carbon_epithelial_flux);
  }
  const Int oxygen_idx = find_chemical_spec(cfg.chemicals, species::OXYGEN);
  if (oxygen_idx >= 0) {
    auto& oxygen = cfg.chemicals[static_cast<size_t>(oxygen_idx)];
    const auto oxygen_boundary_mode = parse_epithelial_boundary_mode(
        "oxygen.epithelial_boundary", cfg.oxygen_epithelial_boundary);
    oxygen.initial_conc =
        oxygen_boundary_mode == EpithelialBoundaryMode::Dirichlet
        ? cfg.chem_env.oxygen.epithelial_conc : 0.0;
    oxygen.boundary_conc = cfg.chem_env.oxygen.epithelial_conc;
    configure_epithelial_boundary(
        oxygen, "oxygen", "oxygen.epithelial_boundary",
        cfg.oxygen_epithelial_boundary,
        cfg.oxygen_epithelial_transfer_coeff, cfg.oxygen_epithelial_flux);
  }
}

namespace {

using FlatKeyHandler = bool (*)(SimulationConfig&, std::string_view, const std::string&);

bool apply_time_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "total_time")           { cfg.time.total_time = parse_config_real(key, val); return true; }
  if (key == "bio_dt")               { cfg.time.bio_dt = parse_config_real(key, val); return true; }
  if (key == "output_interval")      { cfg.time.output_interval = parse_config_real(key, val); return true; }
  if (key == "seed")                 { cfg.seed = parse_config_uint64(key, val); return true; }
  return false;
}

bool apply_domain_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "grid_dx")              { cfg.domain.grid_dx = parse_config_real(key, val); return true; }
  if (key == "chemistry_stride_x" || key == "domain.chemistry_stride_x"
      || key == "domain.chemistry_stride.x") {
    cfg.domain.chemistry_stride[0] = parse_positive_config_int(key, val);
    return true;
  }
  if (key == "chemistry_stride_y" || key == "domain.chemistry_stride_y"
      || key == "domain.chemistry_stride.y") {
    cfg.domain.chemistry_stride[1] = parse_positive_config_int(key, val);
    return true;
  }
  if (key == "chemistry_stride_z" || key == "domain.chemistry_stride_z"
      || key == "domain.chemistry_stride.z") {
    cfg.domain.chemistry_stride[2] = parse_positive_config_int(key, val);
    return true;
  }
  if (key == "grid_halo_width" || key == "domain.grid_halo_width") {
    cfg.domain.grid_halo_width = parse_positive_config_int(key, val);
    return true;
  }
  if (key == "domain_x")             { cfg.domain.hi[0] = parse_config_real(key, val); return true; }
  if (key == "domain_y")             { cfg.domain.hi[1] = parse_config_real(key, val); return true; }
  if (key == "domain_z")             { cfg.domain.hi[2] = parse_config_real(key, val); return true; }
  if (key == "hash_cell_size")       { cfg.domain.hash_cell_size = parse_config_real(key, val); return true; }
  if (key == "ghost_width")          { cfg.domain.ghost_width = parse_config_real(key, val); return true; }
  return false;
}

bool apply_chemistry_key(SimulationConfig& cfg, std::string_view key,
                         const std::string& val) {
  if (key == "chemistry.species_subset" || key == "species_subset") {
    if (val == "full" || val == "nutrient_only" || val == "carbon_only") {
      cfg.species_subset = val;
      return true;
    }
    throw ConfigError(
        "invalid species_subset: expected 'full', 'nutrient_only', or "
        "'carbon_only', got '" + val + "'");
  }
  if (key == "chemistry.toxin_evaluation" || key == "toxin_evaluation") {
    if (val == "grid" || val == "agents") {
      cfg.qssa.toxin_evaluation = val;
      return true;
    }
    throw ConfigError(
        "invalid toxin_evaluation: expected 'grid' or 'agents', got '" + val
        + "'");
  }
  if (key == "chemistry_decomposition" || key == "chemistry.decomposition") {
    if (val == "replicated" || val == "slab") {
      cfg.chemistry_decomposition = val;
      return true;
    }
    if (val == "interface") {
      throw ConfigError(
          "invalid chemistry_decomposition: interface is not yet implemented");
    }
    throw ConfigError(
        "invalid chemistry_decomposition: unknown mode '" + val + "'");
  }
  if (key == "chemistry.toxin_lumping" || key == "toxin_lumping") {
    if (val == "per_receptor" || val == "lumped") {
      cfg.qssa.toxin_lumping = val;
      return true;
    }
    throw ConfigError(
        "invalid toxin_lumping: expected 'per_receptor' or 'lumped', got '"
        + val + "'");
  }
  if (key == "image_series_relative_tolerance") {
    cfg.qssa.image_series_relative_tolerance =
        parse_config_real(key, val);
    return true;
  }
  if (key == "image_series_max_shells") {
    cfg.qssa.image_series_max_shells = parse_config_int(key, val);
    cfg.qssa.image_series_max_shells_explicit = true;
    return true;
  }
  if (key == "image_series_mode") {
    if (val == "corrected"
        || val == "pre_fix_duplicated_reflection") {
      cfg.qssa.image_series_mode = val;
      return true;
    }
    throw ConfigError(
        "invalid image_series_mode: expected 'corrected' or "
        "'pre_fix_duplicated_reflection'");
  }
  return false;
}

bool apply_advection_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "mucus_thickness")      { cfg.advection.mucus_thickness = parse_config_real(key, val); return true; }
  if (key == "radial_turnover")      { cfg.advection.radial_turnover = parse_config_real(key, val); return true; }
  if (key == "washout.trap") {
    if (val == "emergent") {
      cfg.advection.washout_trap = WashoutTrapMode::EMERGENT;
      return true;
    }
    if (val == "imposed") {
      cfg.advection.washout_trap = WashoutTrapMode::IMPOSED;
      return true;
    }
    throw ConfigError(
        "invalid washout.trap: expected 'emergent' or 'imposed', got '"
        + val + "'");
  }
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

bool apply_qssa_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "toxin_cutoff") {
    cfg.qssa.toxin_cutoff = parse_positive_config_real(key, val);
    return true;
  }
  if (key == "toxin.lumen_transfer_length"
      || key == "lumen_transfer_length") {
    const Real transfer_length = parse_config_real(key, val);
    if (!(transfer_length > 0.0)) {
      throw ConfigError(
          "config key '" + std::string(key)
          + "' must be positive and finite, or +inf to disable Robin transfer");
    }
    cfg.qssa.lumen_transfer_length = transfer_length;
    return true;
  }
  if (key == "toxin.lumen_transfer_basis") {
    if (val != "effective" && val != "free") {
      throw ConfigError(
          "invalid toxin.lumen_transfer_basis: expected 'effective' or 'free'");
    }
    cfg.qssa.lumen_transfer_basis = val;
    return true;
  }
  if (key == "low_screening_policy" || key == "qssa.low_screening_policy") {
    if (val != "warn" && val != "error" && val != "allow") {
      throw ConfigError(
          "invalid qssa.low_screening_policy: expected 'warn', 'error', "
          "or 'allow'");
    }
    cfg.qssa.low_screening_policy = val;
    return true;
  }
  if (key == "nutrient_cutoff")      { cfg.qssa.nutrient_cutoff = parse_config_real(key, val); return true; }
  if (key == "colicin_release_rate") { cfg.qssa.colicin_release_rate = parse_config_real(key, val); return true; }
  if (key == "microcin_secretion")   { cfg.qssa.microcin_secretion = parse_config_real(key, val); return true; }
  if (key == "use_fmm")              { cfg.qssa.use_fmm = (val == "true" || val == "1"); return true; }
  if (key == "fmm_theta")            { cfg.qssa.fmm_theta = parse_config_real(key, val); return true; }
  if (key == "fmm_expansion_order")  { cfg.qssa.fmm_expansion_order = parse_config_int(key, val); return true; }
  return false;
}

bool apply_vbf_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "vbf_density")          { cfg.vbf.density = parse_config_real(key, val); return true; }
  if (key == "vbf_viscosity")        { cfg.vbf.viscosity = parse_config_real(key, val); return true; }
  if (key == "vbf_drag_coeff")       { cfg.vbf.drag_coeff = parse_config_real(key, val); return true; }
  if (key == "vbf_nutrient_sink")    { cfg.vbf.nutrient_sink = parse_config_real(key, val); return true; }
  if (key == "vbf_mucin_liberation") { cfg.vbf.mucin_liberation = parse_config_real(key, val); return true; }
  if (key == "vbf_carrying_cap")     { cfg.vbf.carrying_cap = parse_config_real(key, val); return true; }
  if (key == "vbf_mucin_z_gradient") { cfg.vbf.mucin_z_gradient_enabled = (val == "true" || val == "1"); return true; }
  if (key == "vbf_mucin_z_lambda")   { cfg.vbf.mucin_z_gradient_lambda = parse_config_real(key, val); return true; }
  if (key == "vbf_carbon_sink_vmax") { cfg.vbf.carbon_sink_vmax = parse_config_real(key, val); return true; }
  if (key == "vbf_carbon_sink_km")   { cfg.vbf.carbon_sink_km = parse_config_real(key, val); return true; }
  if (key == "vbf_agent_carbon_coupling"
      || key == "vbf.agent_carbon_coupling"
      || key == "agent_carbon_coupling") {
    cfg.vbf.agent_carbon_coupling = parse_config_real(key, val);
    return true;
  }
  return false;
}

bool apply_chemical_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "carbon.boundary_conc" || key == "carbon_boundary_conc") {
    cfg.carbon_boundary_conc = parse_config_real(key, val);
    for (auto& c : cfg.chemicals) {
      if (c.name == species::CARBON) {
        c.boundary_conc = cfg.carbon_boundary_conc;
        if (cfg.carbon_z_amplitude <= 0.0) {
          c.initial_conc = cfg.carbon_boundary_conc;
        }
        return true;
      }
    }
    return true;
  }
  if (key == "carbon.z_amplitude" || key == "carbon_z_amplitude") {
    cfg.carbon_z_amplitude = parse_config_real(key, val);
    for (auto& c : cfg.chemicals) {
      if (c.name == species::CARBON) {
        c.initial_conc = cfg.carbon_z_amplitude;
        return true;
      }
    }
    return true;
  }
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
  if (key == "b12.initial_conc" || key == "b12_initial_conc"
      || key == "corrinoid.initial_conc"
      || key == "corrinoid_initial_conc") {
    cfg.b12_initial_conc = parse_config_real(key, val);
    for (auto& c : cfg.chemicals) {
      if (c.name == species::B12) {
        c.initial_conc = cfg.b12_initial_conc;
        c.boundary_conc = cfg.b12_initial_conc;
        return true;
      }
    }
    return true;
  }
  if (key == "carbon.epithelial_boundary"
      || key == "carbon_epithelial_boundary") {
    (void)parse_epithelial_boundary_mode(key, val);
    cfg.carbon_epithelial_boundary = val;
    return true;
  }
  if (key == "carbon.epithelial_transfer_coeff"
      || key == "carbon_epithelial_transfer_coeff") {
    cfg.carbon_epithelial_transfer_coeff = parse_config_real(key, val);
    return true;
  }
  if (key == "carbon.epithelial_flux"
      || key == "carbon_epithelial_flux") {
    cfg.carbon_epithelial_flux = parse_config_real(key, val);
    return true;
  }
  if (key == "oxygen.z_gradient" || key == "oxygen_z_gradient") {
    cfg.oxygen_z_gradient_enabled = parse_bool_config(val);
    return true;
  }
  if (key == "sos_lysis_prob")       { cfg.fixes.bacteriocin.sos_lysis_prob = parse_config_real(key, val); return true; }
  if (key == "sos_basal_rate")       { cfg.fixes.bacteriocin.sos_basal_rate = parse_config_real(key, val); return true; }
  if (key == "sos_cross_induction_rate") {
    cfg.fixes.bacteriocin.sos_cross_induction_rate = parse_config_real(key, val);
    return true;
  }
  if (key == "D_free_colicin")       { cfg.fixes.bacteriocin.D_free_colicin = parse_config_real(key, val); return true; }
  if (key == "burst_release_tau")    { cfg.fixes.bacteriocin.burst_release_tau = parse_positive_config_real(key, val); return true; }
  if (key == "microcin_mu_penalty")  { cfg.fixes.bacteriocin.microcin_mu_penalty = parse_config_real(key, val); return true; }
  return false;
}

bool apply_bacteriocin_key(SimulationConfig& cfg, std::string_view key,
                           const std::string& val) {
  if (key == "bacteriocin.mucin_charge.r_min") {
    cfg.fixes.bacteriocin.mucin_charge.r_min =
        parse_config_real(key, val);
    return true;
  }
  if (key == "bacteriocin.mucin_charge.amplitude") {
    cfg.fixes.bacteriocin.mucin_charge.amplitude =
        parse_config_real(key, val);
    return true;
  }
  if (key == "bacteriocin.mucin_charge.dz_half") {
    cfg.fixes.bacteriocin.mucin_charge.dz_half =
        parse_config_real(key, val);
    return true;
  }
  if (key == "bacteriocin.mucin_charge.width") {
    cfg.fixes.bacteriocin.mucin_charge.width =
        parse_positive_config_real(key, val);
    return true;
  }
  if (key == "bacteriocin.mucin_charge.ph") {
    cfg.fixes.bacteriocin.mucin_charge.ph =
        parse_config_real(key, val);
    return true;
  }
  return false;
}

bool apply_plasmid_override_key(SimulationConfig& cfg,
                                std::string_view key,
                                const std::string& val) {
  constexpr std::string_view prefix = "plasmid_overrides.";
  if (!key.starts_with(prefix)) return false;

  const std::string_view remainder = key.substr(prefix.size());
  const size_t separator = remainder.find('.');
  if (separator == std::string_view::npos) {
    throw ConfigError("plasmid override key must name a field");
  }
  const std::string plasmid_name(remainder.substr(0, separator));
  const std::string field(remainder.substr(separator + 1));
  const PlasmidEntry* entry = PlasmidLibrary::find(plasmid_name);
  if (entry == nullptr) {
    throw ConfigError("unknown plasmid name in override: " + plasmid_name);
  }

  const Real value = parse_config_real(key, val);
  if (!std::isfinite(value) || value < 0.0
      || (field == "retardation" && value == 0.0)
      || (field == "diff_coeff" && value == 0.0)) {
    throw ConfigError("plasmid override '" + std::string(key)
                      + "' must be finite and nonnegative"
                      + (field == "retardation" || field == "diff_coeff"
                             ? " and positive"
                             : ""));
  }

  auto& override_values = cfg.plasmid_overrides[entry->name];
  if (field == "retardation") {
    override_values.retardation = value;
    return true;
  }
  if (field == "diff_coeff") {
    override_values.diff_coeff = value;
    return true;
  }
  if (field == "burst_size") {
    override_values.burst_size = value;
    return true;
  }
  throw ConfigError("unknown plasmid override field: " + field);
}

bool apply_receptor_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "kd_b12_btuB")          { cfg.fixes.receptor.kd_b12_btuB = parse_config_real(key, val); return true; }
  // Alias — Spec 6 / Receptor Ligand Parameterization: the BtuB ligand is the
  // dominant corrinoid analog, not true cobalamin. Same underlying Kd field.
  if (key == "kd_corrinoid_btuB")    { cfg.fixes.receptor.kd_b12_btuB = parse_config_real(key, val); return true; }
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

bool apply_conjugation_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "pili_length")           { cfg.fixes.conjugation.pili_length = parse_config_real(key, val); return true; }
  if (key == "base_transfer_rate")    { cfg.fixes.conjugation.base_transfer_rate = parse_config_real(key, val); return true; }
  if (key == "shear_critical")          { cfg.fixes.conjugation.shear_critical = parse_config_real(key, val); return true; }
  if (key == "plasmid_copy_cost")     { cfg.fixes.conjugation.plasmid_copy_cost = parse_config_real(key, val); return true; }
  if (key == "pili_heterogeneity")    { cfg.fixes.conjugation.pili_heterogeneity = (val == "true" || val == "1"); return true; }
  if (key == "pili_length_min")       { cfg.fixes.conjugation.pili_length_min = parse_config_real(key, val); return true; }
  if (key == "pili_length_max")       { cfg.fixes.conjugation.pili_length_max = parse_config_real(key, val); return true; }
  return false;
}

bool apply_mutation_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
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

bool apply_io_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
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
  if (key == "restart.enabled" || key == "restart_enabled") {
    cfg.restart.enabled = parse_bool_config(val); return true;
  }
  if (key == "restart.directory" || key == "restart_directory") {
    validate_path_syntax(val);
    cfg.restart.directory = val;
    return true;
  }
  if (key == "restart.interval_steps" || key == "restart_interval_steps") {
    cfg.restart.interval_steps = parse_config_int(key, val); return true;
  }
  return false;
}

bool apply_hdf5_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
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
  if (key == "hdf5.schedule.provenance" || key == "hdf5_schedule_provenance") {
    cfg.hdf5.schedule.provenance = parse_config_int(key, val); return true;
  }
  return false;
}

bool apply_siderophore_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
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
  if (key == "siderophore.Vmax_reimport" || key == "siderophore_Vmax_reimport") {
    cfg.chem_env.siderophore.Vmax_reimport = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_ferrichrome_key(SimulationConfig& cfg, std::string_view key,
                           const std::string& val) {
  if (key == "ferrichrome.enabled" || key == "ferrichrome_enabled") {
    cfg.chem_env.ferrichrome.enabled = parse_bool_config(val);
    return true;
  }
  if (key == "ferrichrome.initial_conc" || key == "ferrichrome_initial_conc") {
    cfg.chem_env.ferrichrome.initial_conc = parse_config_real(key, val);
    return true;
  }
  if (key == "ferrichrome.boundary_conc"
      || key == "ferrichrome_boundary_conc") {
    cfg.chem_env.ferrichrome.boundary_conc = parse_config_real(key, val);
    return true;
  }
  return false;
}

bool apply_dt_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "adaptive_dt_enabled")  { cfg.adaptive_dt.enabled = (val == "true" || val == "1"); return true; }
  if (key == "dt_min")               { cfg.adaptive_dt.min = parse_config_real(key, val); return true; }
  if (key == "dt_max")               { cfg.adaptive_dt.max = parse_config_real(key, val); return true; }
  if (key == "dt_safety")            { cfg.adaptive_dt.safety = parse_config_real(key, val); return true; }
  if (key == "dt_growth_limit")      { cfg.adaptive_dt.growth_limit = parse_config_real(key, val); return true; }
  return false;
}

bool apply_immigration_key(SimulationConfig& cfg, std::string_view key,
                           const std::string& val) {
  if (key == "immigration.enabled") {
    cfg.immigration.enabled = parse_bool_config(val);
    return true;
  }
  if (key == "immigration.count") {
    cfg.immigration.count = parse_config_int(key, val);
    return true;
  }
  if (key == "immigration.strain_index") {
    cfg.immigration.strain_index = parse_config_int(key, val);
    return true;
  }
  if (key == "immigration.placement") {
    if (val != "uniform" && val != "at_distance" && val != "z_slab") {
      throw ConfigError("invalid immigration.placement: " + val);
    }
    cfg.immigration.placement = val;
    return true;
  }
  if (key == "immigration.distance") {
    cfg.immigration.distance = parse_config_real(key, val);
    return true;
  }
  if (key == "immigration.distance_tolerance") {
    cfg.immigration.distance_tolerance = parse_config_real(key, val);
    return true;
  }
  if (key == "immigration.distance_reference") {
    if (val != "nearest_agent" && val != "centroid") {
      throw ConfigError("invalid immigration.distance_reference: " + val);
    }
    cfg.immigration.distance_reference = val;
    return true;
  }
  if (key == "immigration.z_min") {
    cfg.immigration.z_min = parse_config_real(key, val);
    return true;
  }
  if (key == "immigration.z_max") {
    cfg.immigration.z_max = parse_config_real(key, val);
    return true;
  }
  if (key == "immigration.schedule") {
    if (val != "pulse" && val != "continuous") {
      throw ConfigError("invalid immigration.schedule: " + val);
    }
    cfg.immigration.schedule = val;
    return true;
  }
  if (key == "immigration.step") {
    cfg.immigration.step = parse_config_int(key, val);
    return true;
  }
  if (key == "immigration.rate") {
    cfg.immigration.rate = parse_config_real(key, val);
    return true;
  }
  return false;
}

bool apply_initial_population_key(SimulationConfig& cfg,
                                  std::string_view key,
                                  const std::string& val) {
  if (key == "initial_population.placement") {
    if (val != "legacy" && val != "z_slab" && val != "anatomic") {
      throw ConfigError("invalid initial_population.placement: " + val);
    }
    cfg.initial_population.placement = val;
    return true;
  }
  if (key == "initial_population.z_min") {
    cfg.initial_population.z_min = parse_config_real(key, val);
    return true;
  }
  if (key == "initial_population.z_max") {
    cfg.initial_population.z_max = parse_config_real(key, val);
    return true;
  }
  if (key == "initial_population.anatomic_exclusion_floor"
      || key == "initial_population_anatomic_exclusion_floor") {
    cfg.initial_population.anatomic_exclusion_floor =
        parse_config_real(key, val);
    return true;
  }
  if (key == "initial_population.anatomic_exponential_scale"
      || key == "initial_population_anatomic_exponential_scale") {
    cfg.initial_population.anatomic_exponential_scale =
        parse_config_real(key, val);
    return true;
  }
  if (key == "initial_population.anatomic_outer_extent"
      || key == "initial_population_anatomic_outer_extent") {
    cfg.initial_population.anatomic_outer_extent =
        parse_config_real(key, val);
    return true;
  }
  return false;
}

bool apply_closure_key(SimulationConfig& cfg, std::string_view key,
                       const std::string& val) {
  if (key == "closure.enforce_delivery_realization"
      || key == "closure_enforce_delivery_realization") {
    cfg.closure.enforce_delivery_realization = parse_bool_config(val);
    return true;
  }
  if (key == "closure.zero_realization_grace_steps"
      || key == "closure_zero_realization_grace_steps") {
    cfg.closure.zero_realization_grace_steps = parse_config_int(key, val);
    return true;
  }
  if (key == "closure.enforce_reaction_clip"
      || key == "closure_enforce_reaction_clip") {
    cfg.closure.enforce_reaction_clip = parse_bool_config(val);
    return true;
  }
  if (key == "closure.reaction_clip_tolerance_fraction"
      || key == "closure_reaction_clip_tolerance_fraction") {
    cfg.closure.reaction_clip_tolerance_fraction =
        parse_config_real(key, val);
    return true;
  }
  return false;
}

bool apply_metabolism_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "division_threshold")      { cfg.fixes.metabolism.division_threshold = parse_config_real(key, val); return true; }
  if (key == "bacteriostasis_threshold") { cfg.fixes.metabolism.bacteriostasis_threshold = parse_config_real(key, val); return true; }
  if (key == "maintenance_rate")        { cfg.fixes.metabolism.maintenance_rate = parse_config_real(key, val); return true; }
  if (key == "carbon_maintenance_rate"
      || key == "metabolism.carbon_maintenance_rate") {
    cfg.fixes.metabolism.carbon_maintenance_rate =
        parse_config_real(key, val);
    return true;
  }
  if (key == "delivery_far_field_radius"
      || key == "metabolism.delivery_far_field_radius"
      || key == "metabolism_delivery_far_field_radius") {
    cfg.fixes.metabolism.delivery_far_field_radius =
        parse_config_real(key, val);
    return true;
  }
  if (key == "acid_inhibition_enabled"
      || key == "metabolism.acid_inhibition_enabled"
      || key == "metabolism_acid_inhibition_enabled") {
    cfg.fixes.metabolism.acid_inhibition_enabled = parse_bool_config(val); return true;
  }
  if (key == "acid_inhibition_max"
      || key == "metabolism.acid_inhibition_max"
      || key == "metabolism_acid_inhibition_max") {
    cfg.fixes.metabolism.acid_inhibition_max = parse_config_real(key, val); return true;
  }
  if (key == "acid_inhibition_Ki"
      || key == "metabolism.acid_inhibition_Ki"
      || key == "metabolism_acid_inhibition_Ki") {
    cfg.fixes.metabolism.acid_inhibition_Ki = parse_config_real(key, val); return true;
  }
  if (key == "acetate_pKa"
      || key == "metabolism.acetate_pKa"
      || key == "metabolism_acetate_pKa") {
    cfg.fixes.metabolism.acetate_pKa = parse_config_real(key, val); return true;
  }
  if (key == "metE_penalty")            { cfg.fixes.metabolism.metE_penalty = parse_config_real(key, val); return true; }
  if (key == "metE_acetate_km")         { cfg.fixes.metabolism.metE_acetate_km = parse_config_real(key, val); return true; }
  if (key == "metE_acetate_max_factor") { cfg.fixes.metabolism.metE_acetate_max_factor = parse_config_real(key, val); return true; }
  if (key == "eut_km")                  { cfg.fixes.metabolism.eut_km = parse_config_real(key, val); return true; }
  if (key == "eut_max_penalty")         { cfg.fixes.metabolism.eut_max_penalty = parse_config_real(key, val); return true; }
  if (key == "km_iron_primary")         { cfg.fixes.metabolism.km_iron_primary = parse_config_real(key, val); return true; }
  if (key == "km_iron_iroN")            { cfg.fixes.metabolism.km_iron_iroN = parse_config_real(key, val); return true; }
  if (key == "km_iron_iutA")            { cfg.fixes.metabolism.km_iron_iutA = parse_config_real(key, val); return true; }
  if (key == "km_iron_fiu")             { cfg.fixes.metabolism.km_iron_fiu = parse_config_real(key, val); return true; }
  if (key == "uptake_limit" || key == "metabolism.uptake_limit") {
    if (val == "none" || val == "sherwood" || val == "voxel"
        || val == "delivery") {
      cfg.fixes.metabolism.uptake_limit = val;
      return true;
    }
    throw ConfigError(
        "invalid uptake_limit: expected 'none', 'sherwood', 'voxel', or "
        "'delivery', got '"
        + val + "'");
  }
  return false;
}

bool apply_mechanics_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "hertz_k")            { cfg.fixes.mechanics.hertz_k = parse_config_real(key, val); return true; }
  if (key == "hertzian_enabled")   { cfg.fixes.mechanics.hertzian_enabled = parse_bool_config(val); return true; }
  if (key == "adhesion_enabled")   { cfg.fixes.mechanics.adhesion_enabled = parse_bool_config(val); return true; }
  if (key == "adhesion_strength")  { cfg.fixes.mechanics.adhesion_strength = parse_config_real(key, val); return true; }
  if (key == "adhesion_range")     { cfg.fixes.mechanics.adhesion_range = parse_config_real(key, val); return true; }
  return false;
}

bool apply_misc_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "gpu_enabled")          { cfg.gpu.enabled = parse_bool_config(val); return true; }
  if (key == "gpu_device_id")        { cfg.gpu.device_id = parse_config_int(key, val); return true; }
  if (key == "profile_steps")        { cfg.profile_steps = (val == "true" || val == "1"); return true; }
  if (key == "dysbiosis_threshold")  { cfg.dysbiosis_threshold = parse_config_real(key, val); return true; }
  if (key == "dysbiosis_sampling_interval") {
    cfg.dysbiosis_sampling_interval = parse_config_real(key, val);
    return true;
  }
  if (key == "dysbiosis_sample_count") {
    cfg.dysbiosis_sample_count = parse_config_int(key, val);
    return true;
  }
  return false;
}

bool apply_oxygen_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
  if (key == "oxygen.enabled" || key == "oxygen_enabled") {
    cfg.chem_env.oxygen.enabled = parse_bool_config(val); return true;
  }
  if (key == "oxygen.epithelial_conc" || key == "oxygen_epithelial_conc") {
    cfg.chem_env.oxygen.epithelial_conc = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.epithelial_boundary"
      || key == "oxygen_epithelial_boundary") {
    (void)parse_epithelial_boundary_mode(key, val);
    cfg.oxygen_epithelial_boundary = val;
    return true;
  }
  if (key == "oxygen.epithelial_transfer_coeff"
      || key == "oxygen_epithelial_transfer_coeff") {
    cfg.oxygen_epithelial_transfer_coeff = parse_config_real(key, val);
    return true;
  }
  if (key == "oxygen.epithelial_flux"
      || key == "oxygen_epithelial_flux") {
    cfg.oxygen_epithelial_flux = parse_config_real(key, val);
    return true;
  }
  if (key == "oxygen.delivery_uptake_enabled"
      || key == "oxygen_delivery_uptake_enabled") {
    cfg.chem_env.oxygen.delivery_uptake_enabled = parse_bool_config(val);
    return true;
  }
  if (key == "oxygen.respiration_driver"
      || key == "oxygen_respiration_driver") {
    if (val != "ambient" && val != "funded") {
      throw ConfigError(
          "invalid oxygen.respiration_driver: expected 'ambient' or "
          "'funded', got '" + val + "'");
    }
    cfg.chem_env.oxygen.respiration_driver = val;
    return true;
  }
  if (key == "oxygen.ros_driver" || key == "oxygen_ros_driver") {
    if (val != "ambient" && val != "funded") {
      throw ConfigError(
          "invalid oxygen.ros_driver: expected 'ambient' or 'funded', got '"
          + val + "'");
    }
    cfg.chem_env.oxygen.ros_driver = val;
    return true;
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
  if (key == "oxygen.metabolic_switch_enabled"
      || key == "oxygen_metabolic_switch_enabled") {
    cfg.chem_env.oxygen.metabolic_switch_enabled = parse_bool_config(val); return true;
  }
  if (key == "oxygen.mu_crit" || key == "oxygen_mu_crit") {
    cfg.chem_env.oxygen.mu_crit = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.aerobic_mu_factor" || key == "oxygen_aerobic_mu_factor") {
    cfg.chem_env.oxygen.aerobic_mu_factor = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.anaerobic_mu_factor" || key == "oxygen_anaerobic_mu_factor") {
    cfg.chem_env.oxygen.anaerobic_mu_factor = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.aerobic_carbon_cost_factor"
      || key == "oxygen_aerobic_carbon_cost_factor") {
    cfg.chem_env.oxygen.aerobic_carbon_cost_factor = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.anaerobic_carbon_cost_factor"
      || key == "oxygen_anaerobic_carbon_cost_factor") {
    cfg.chem_env.oxygen.anaerobic_carbon_cost_factor = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.tau_metabolic_switch"
      || key == "oxygen_tau_metabolic_switch") {
    cfg.chem_env.oxygen.tau_metabolic_switch = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.ferm_acid_yield" || key == "oxygen_ferm_acid_yield") {
    cfg.chem_env.oxygen.ferm_acid_yield = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.anaerobic_maintenance_factor"
      || key == "oxygen_anaerobic_maintenance_factor") {
    cfg.chem_env.oxygen.anaerobic_maintenance_factor = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.q_consumption" || key == "oxygen_q_consumption") {
    cfg.chem_env.oxygen.q_consumption = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.q_maintenance" || key == "oxygen_q_maintenance") {
    cfg.chem_env.oxygen.q_maintenance = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.vbf_sink" || key == "oxygen_vbf_sink") {
    cfg.chem_env.oxygen.vbf_sink = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.k_ROS" || key == "oxygen_k_ROS") {
    cfg.chem_env.oxygen.k_ROS = parse_config_real(key, val); return true;
  }
  if (key == "oxygen.k_ROS_respiratory"
      || key == "oxygen_k_ROS_respiratory") {
    cfg.chem_env.oxygen.k_ROS_respiratory = parse_config_real(key, val);
    return true;
  }
  if (key == "oxygen.k_ROS_funded"
      || key == "oxygen_k_ROS_funded") {
    cfg.chem_env.oxygen.k_ROS_funded = parse_config_real(key, val);
    return true;
  }
  return false;
}

bool apply_acetate_dyn_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
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

bool apply_mucin_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
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

bool apply_protease_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
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

bool apply_fur_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
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

bool apply_cdi_key(SimulationConfig& cfg, std::string_view key, const std::string& val) {
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

bool apply_motility_run_key(SimulationConfig& cfg, std::string_view key,
                            const std::string& val) {
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
  return false;
}

bool apply_motility_chemotaxis_key(SimulationConfig& cfg, std::string_view key,
                                   const std::string& val) {
  if (key == "motility.chemotaxis_enabled" || key == "motility_chemotaxis_enabled") {
    cfg.cell_bio.motility.chemotaxis_enabled = parse_bool_config(val); return true;
  }
  if (key == "motility.chi_carbon" || key == "motility_chi_carbon") {
    cfg.cell_bio.motility.chi_carbon = parse_config_real(key, val); return true;
  }
  if (key == "motility.chemotaxis_threshold" || key == "motility_chemotaxis_threshold") {
    cfg.cell_bio.motility.chemotaxis_threshold = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_motility_aerotaxis_key(SimulationConfig& cfg, std::string_view key,
                                  const std::string& val) {
  if (key == "motility.aerotaxis_enabled" || key == "motility_aerotaxis_enabled") {
    cfg.cell_bio.motility.aerotaxis_enabled = parse_bool_config(val); return true;
  }
  if (key == "motility.aerotaxis_sensitivity" || key == "motility_aerotaxis_sensitivity") {
    cfg.cell_bio.motility.aerotaxis_sensitivity = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_motility_energy_key(SimulationConfig& cfg, std::string_view key,
                               const std::string& val) {
  if (key == "motility.energy_taxis_enabled" || key == "motility_energy_taxis_enabled") {
    cfg.cell_bio.motility.energy_taxis_enabled = parse_bool_config(val); return true;
  }
  if (key == "motility.energy_taxis_floor" || key == "motility_energy_taxis_floor") {
    cfg.cell_bio.motility.energy_taxis_floor = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_motility_surface_key(SimulationConfig& cfg, std::string_view key,
                                const std::string& val) {
  if (key == "motility.surface_sensing_enabled" || key == "motility_surface_sensing_enabled") {
    cfg.cell_bio.motility.surface_sensing_enabled = parse_bool_config(val); return true;
  }
  if (key == "motility.surface_sensing_depth" || key == "motility_surface_sensing_depth") {
    cfg.cell_bio.motility.surface_sensing_depth = parse_config_real(key, val); return true;
  }
  if (key == "motility.surface_sensing_floor" || key == "motility_surface_sensing_floor") {
    cfg.cell_bio.motility.surface_sensing_floor = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_motility_mucin_key(SimulationConfig& cfg, std::string_view key,
                              const std::string& val) {
  if (key == "motility.mucin_drag_enabled" || key == "motility_mucin_drag_enabled") {
    cfg.cell_bio.motility.mucin_drag_enabled = parse_bool_config(val); return true;
  }
  if (key == "motility.mucin_drag_reference" || key == "motility_mucin_drag_reference") {
    cfg.cell_bio.motility.mucin_drag_reference = parse_config_real(key, val); return true;
  }
  return false;
}

bool apply_motility_cluster_key(SimulationConfig& cfg, std::string_view key,
                                const std::string& val) {
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

bool apply_motility_key(SimulationConfig& cfg, std::string_view key,
                        const std::string& val) {
  return apply_motility_run_key(cfg, key, val)
      || apply_motility_chemotaxis_key(cfg, key, val)
      || apply_motility_aerotaxis_key(cfg, key, val)
      || apply_motility_energy_key(cfg, key, val)
      || apply_motility_surface_key(cfg, key, val)
      || apply_motility_mucin_key(cfg, key, val)
      || apply_motility_cluster_key(cfg, key, val);
}

bool apply_quorum_sensing_key(SimulationConfig& cfg, std::string_view key,
                              const std::string& val) {
  if (key == "quorum_sensing.enabled" || key == "quorum_sensing_enabled") {
    cfg.quorum_sensing.enabled = parse_bool_config(val); return true;
  }
  if (key == "quorum_sensing.ai2_basal_rate" || key == "quorum_sensing_ai2_basal_rate") {
    cfg.quorum_sensing.ai2_basal_rate = parse_config_real(key, val); return true;
  }
  if (key == "quorum_sensing.ai2_growth_coupled"
      || key == "quorum_sensing_ai2_growth_coupled") {
    cfg.quorum_sensing.ai2_growth_coupled = parse_config_real(key, val); return true;
  }
  if (key == "quorum_sensing.lsr_vmax" || key == "quorum_sensing_lsr_vmax") {
    cfg.quorum_sensing.lsr_vmax = parse_config_real(key, val); return true;
  }
  if (key == "quorum_sensing.lsr_km" || key == "quorum_sensing_lsr_km") {
    cfg.quorum_sensing.lsr_km = parse_config_real(key, val); return true;
  }
  if (key == "quorum_sensing.ai2_D_free" || key == "quorum_sensing_ai2_D_free") {
    cfg.quorum_sensing.ai2_D_free = parse_config_real(key, val); return true;
  }
  if (key == "quorum_sensing.ai2_decay_rate" || key == "quorum_sensing_ai2_decay_rate") {
    cfg.quorum_sensing.ai2_decay_rate = parse_config_real(key, val); return true;
  }
  // Spec key name is ai2_chemotaxis (maps to ai2_chemotaxis_enabled)
  if (key == "quorum_sensing.ai2_chemotaxis" || key == "quorum_sensing_ai2_chemotaxis") {
    cfg.quorum_sensing.ai2_chemotaxis_enabled = parse_bool_config(val); return true;
  }
  if (key == "quorum_sensing.chi_ai2" || key == "quorum_sensing_chi_ai2") {
    cfg.quorum_sensing.chi_ai2 = parse_config_real(key, val); return true;
  }
  return false;
}

constexpr std::array<FlatKeyHandler, 31> k_flat_key_handlers = {
  apply_time_key,
  apply_domain_key,
  apply_chemistry_key,
  apply_advection_key,
  apply_qssa_key,
  apply_vbf_key,
  apply_chemical_key,
  apply_bacteriocin_key,
  apply_plasmid_override_key,
  apply_metabolism_key,
  apply_receptor_key,
  apply_conjugation_key,
  apply_mutation_key,
  apply_mechanics_key,
  apply_io_key,
  apply_hdf5_key,
  apply_dt_key,
  apply_immigration_key,
  apply_initial_population_key,
  apply_closure_key,
  apply_misc_key,
  apply_oxygen_key,
  apply_acetate_dyn_key,
  apply_mucin_key,
  apply_protease_key,
  apply_siderophore_key,
  apply_ferrichrome_key,
  apply_fur_key,
  apply_cdi_key,
  apply_motility_key,
  apply_quorum_sensing_key,
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
      warn_unknown_config_key(key);
    }
  }
}

}  // namespace

void InputParser::handle_unknown_config_key(std::string_view key) {
  if (strict_config_enabled()) {
    throw ConfigError("unknown config key '" + std::string(key) + "'");
  }
  std::cerr << "Warning: unknown config key '" << key << "' ignored\n";
}

bool InputParser::apply_flat_key(SimulationConfig& cfg,
                                 std::string_view key,
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
    throw IOError("cannot open config file '" + filename + "': " + ex.what());
  }

  std::ifstream ifs(config_path);
  if (!ifs.is_open()) {
    throw IOError("cannot open config file '" + config_path + "'");
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
