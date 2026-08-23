/* -----------------------------------------------------------------------
   GutIBM – Minimal JSON helpers for config fragments
   ----------------------------------------------------------------------- */

#include "config_json.h"

#include <cctype>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string_view>
#include "error.h"
#include <utility>

namespace gutibm {

namespace {

class JsonCursor;

void apply_json_scalar(SimulationConfig& cfg, const std::string& key, JsonCursor& cursor);

int hex_value(char c) {
  if (c >= '0' && c <= '9') return c - '0';
  if (c >= 'a' && c <= 'f') return c - 'a' + 10;
  if (c >= 'A' && c <= 'F') return c - 'A' + 10;
  throw ConfigError("invalid hex digit in JSON unicode escape");
}

void append_utf8(std::string& out, uint32_t codepoint) {
  if (codepoint <= 0x7f) {
    out.push_back(static_cast<char>(codepoint));
  } else if (codepoint <= 0x7ff) {
    out.push_back(static_cast<char>(0xc0 | (codepoint >> 6)));
    out.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
  } else if (codepoint <= 0xffff) {
    out.push_back(static_cast<char>(0xe0 | (codepoint >> 12)));
    out.push_back(static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f)));
    out.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
  } else {
    out.push_back(static_cast<char>(0xf0 | (codepoint >> 18)));
    out.push_back(static_cast<char>(0x80 | ((codepoint >> 12) & 0x3f)));
    out.push_back(static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f)));
    out.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
  }
}

class JsonCursor {
 public:
  explicit JsonCursor(std::string text) : text_(std::move(text)) {}

  void skip_ws() {
    while (pos_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_]))) {
      ++pos_;
    }
  }

  bool match(char expected) {
    skip_ws();
    if (pos_ < text_.size() && text_[pos_] == expected) {
      ++pos_;
      return true;
    }
    return false;
  }

  std::string parse_string() {
    skip_ws();
    if (pos_ >= text_.size() || text_[pos_] != '"') {
      throw ConfigError("expected JSON string");
    }
    ++pos_;
    std::string out;
    while (pos_ < text_.size()) {
      char c = text_[pos_++];
      if (c == '"') return out;
      if (c == '\\') {
        if (pos_ >= text_.size()) throw ConfigError("truncated JSON escape");
        char esc = text_[pos_++];
        if (esc == '"' || esc == '\\' || esc == '/') out.push_back(esc);
        else if (esc == 'b') out.push_back('\b');
        else if (esc == 'f') out.push_back('\f');
        else if (esc == 'n') out.push_back('\n');
        else if (esc == 'r') out.push_back('\r');
        else if (esc == 't') out.push_back('\t');
        else if (esc == 'u') append_unicode_escape(out);
        else throw ConfigError("unsupported JSON escape");
      } else {
        if (static_cast<unsigned char>(c) < 0x20) {
          throw ConfigError("unescaped control character in JSON string");
        }
        out.push_back(c);
      }
    }
    throw ConfigError("unterminated JSON string");
  }

  uint16_t parse_unicode_code_unit() {
    if (text_.size() - pos_ < 4) {
      throw ConfigError("truncated JSON unicode escape");
    }
    uint16_t code_unit = 0;
    for (int i = 0; i < 4; ++i) {
      code_unit = static_cast<uint16_t>(
          (code_unit << 4) | hex_value(text_[pos_++]));
    }
    return code_unit;
  }

  void append_unicode_escape(std::string& out) {
    const uint16_t first = parse_unicode_code_unit();
    if (first >= 0xd800 && first <= 0xdbff) {
      if (text_.size() - pos_ < 6 || text_[pos_] != '\\' ||
          text_[pos_ + 1] != 'u') {
        throw ConfigError("missing low surrogate in JSON unicode escape");
      }
      pos_ += 2;
      const uint16_t second = parse_unicode_code_unit();
      if (second < 0xdc00 || second > 0xdfff) {
        throw ConfigError("invalid low surrogate in JSON unicode escape");
      }
      const uint32_t codepoint =
          0x10000u + ((static_cast<uint32_t>(first) - 0xd800u) << 10) +
          (static_cast<uint32_t>(second) - 0xdc00u);
      append_utf8(out, codepoint);
      return;
    }
    if (first >= 0xdc00 && first <= 0xdfff) {
      throw ConfigError("unexpected low surrogate in JSON unicode escape");
    }
    append_utf8(out, first);
  }

  Real parse_number() {
    skip_ws();
    size_t start = pos_;
    if (pos_ < text_.size() && (text_[pos_] == '-' || text_[pos_] == '+')) ++pos_;
    while (pos_ < text_.size() &&
           (std::isdigit(static_cast<unsigned char>(text_[pos_])) || text_[pos_] == '.' ||
            text_[pos_] == 'e' || text_[pos_] == 'E' || text_[pos_] == '-' ||
            text_[pos_] == '+')) {
      ++pos_;
    }
    if (start == pos_) throw ConfigError("expected JSON number");
    return std::stod(text_.substr(start, pos_ - start));
  }

  bool parse_bool() {
    skip_ws();
    if (text_.compare(pos_, 4, "true") == 0) {
      pos_ += 4;
      return true;
    }
    if (text_.compare(pos_, 5, "false") == 0) {
      pos_ += 5;
      return false;
    }
    throw ConfigError("expected JSON boolean");
  }

  std::vector<std::string> parse_string_array() {
    std::vector<std::string> out;
    if (!match('[')) throw ConfigError("expected JSON array");
    skip_ws();
    if (match(']')) return out;
    while (true) {
      out.push_back(parse_string());
      skip_ws();
      if (match(']')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON array");
    }
    return out;
  }

  std::vector<SimulationConfig::InitialStrain> parse_strain_array() {
    std::vector<SimulationConfig::InitialStrain> out;
    if (!match('[')) throw ConfigError("expected JSON array");
    skip_ws();
    if (match(']')) return out;
    while (true) {
      out.push_back(parse_strain_object());
      skip_ws();
      if (match(']')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON array");
    }
    return out;
  }

  void parse_receptor_expression_object(
      std::map<std::string, Real>& expressions) {
    if (!match('{')) throw ConfigError("expected receptor expression object");
    skip_ws();
    if (match('}')) return;
    while (true) {
      const std::string receptor_name = parse_string();
      if (!match(':')) throw ConfigError("expected ':' in receptor expression");
      const auto receptor = receptor_type_from_name(receptor_name);
      if (!receptor.has_value()) {
        throw ConfigError("unknown receptor name: " + receptor_name);
      }
      const Real expression = parse_number();
      if (!std::isfinite(expression) || expression < 0.0 || expression > 1.0) {
        throw ConfigError("receptor expression must be finite and in [0, 1]");
      }
      expressions[receptor_name] = expression;
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in receptor expression");
    }
  }

  void parse_plasmid_overrides_object(
      SimulationConfig& cfg) {
    if (!match('{')) throw ConfigError("expected plasmid overrides object");
    skip_ws();
    if (match('}')) return;
    while (true) {
      const std::string plasmid_name = parse_string();
      if (!match(':')) throw ConfigError("expected ':' in plasmid override");
      if (!match('{')) throw ConfigError("expected plasmid override object");
      skip_ws();
      if (!match('}')) {
        while (true) {
          const std::string field = parse_string();
          if (!match(':')) throw ConfigError("expected ':' in plasmid override");
          const Real value = parse_number();
          std::ostringstream value_text;
          value_text << std::setprecision(17) << value;
          if (!InputParser::apply_flat_key(
                  cfg, "plasmid_overrides." + plasmid_name + "." + field,
                  value_text.str())) {
            throw ConfigError("unknown plasmid override field: " + field);
          }
          skip_ws();
          if (match('}')) break;
          if (!match(',')) throw ConfigError("expected ',' in plasmid override");
        }
      }
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in plasmid overrides");
    }
  }

  char peek() {
    skip_ws();
    return pos_ < text_.size() ? text_[pos_] : '\0';
  }

  void skip_value() {
    skip_ws();
    if (pos_ >= text_.size()) throw ConfigError("unexpected end of JSON");

    char c = text_[pos_];
    if (c == '"') {
      (void)parse_string();
      return;
    }
    if (c == '{') {
      skip_object();
      return;
    }
    if (c == '[') {
      skip_array();
      return;
    }
    if (std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+') {
      (void)parse_number();
      return;
    }
    if (text_.compare(pos_, 4, "true") == 0) {
      pos_ += 4;
      return;
    }
    if (text_.compare(pos_, 5, "false") == 0) {
      pos_ += 5;
      return;
    }
    if (text_.compare(pos_, 4, "null") == 0) {
      pos_ += 4;
      return;
    }
    throw ConfigError("unsupported JSON value");
  }

  SimulationConfig::InitialStrain parse_strain_object() {
    if (!match('{')) throw ConfigError("expected JSON object");

    SimulationConfig::InitialStrain strain{};

    skip_ws();
    if (match('}')) return strain;

    while (true) {
      std::string key = parse_string();
      if (!match(':')) throw ConfigError("expected ':' in JSON object");

      if (!key.empty() && key.front() == '_') {
        skip_value();
      } else if (key == "type") {
        strain.type = static_cast<Int>(parse_number());
      } else if (key == "count") {
        strain.count = static_cast<Int>(parse_number());
      } else if (key == "mu_max") {
        strain.mu_max = parse_number();
      } else if (key == "plasmids") {
        strain.plasmids = parse_string_array();
      } else if (key == "conjugative") {
        strain.conjugative = parse_bool();
      } else if (key == "cdi_type") {
        strain.cdi_type = static_cast<uint16_t>(parse_number());
      } else if (key == "cdi_immunity") {
        strain.cdi_immunity = static_cast<uint16_t>(parse_number());
      } else if (key == "receptor_expression"
                 || key == "receptor_genotype" || key == "receptors") {
        parse_receptor_expression_object(strain.receptor_expression);
      } else {
        skip_value();
      }

      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON object");
    }

    return strain;
  }

  void skip_object() {
    if (!match('{')) throw ConfigError("expected JSON object");
    skip_ws();
    if (match('}')) return;
    while (true) {
      (void)parse_string();
      if (!match(':')) throw ConfigError("expected ':' in JSON object");
      skip_value();
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON object");
    }
  }

  void skip_array() {
    if (!match('[')) throw ConfigError("expected JSON array");
    skip_ws();
    if (match(']')) return;
    while (true) {
      skip_value();
      skip_ws();
      if (match(']')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON array");
    }
  }

  void parse_hdf5_schedule(SimulationConfig& cfg) {
    if (!match('{')) throw ConfigError("expected JSON object for hdf5.schedule");
    skip_ws();
    if (match('}')) return;
    while (true) {
      std::string key = parse_string();
      if (!match(':')) throw ConfigError("expected ':' in JSON object");
      if (key == "grid_species") {
        cfg.hdf5.schedule.grid_species = parse_string_array();
      } else {
        const std::string flat_key = "hdf5.schedule." + key;
        apply_json_scalar(cfg, flat_key, *this);
      }
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON object");
    }
  }

  void parse_hdf5_object(SimulationConfig& cfg) {
    if (!match('{')) throw ConfigError("expected JSON object for hdf5");
    skip_ws();
    if (match('}')) return;
    while (true) {
      std::string key = parse_string();
      if (!match(':')) throw ConfigError("expected ':' in JSON object");
      if (key == "schedule") {
        parse_hdf5_schedule(cfg);
      } else {
        const std::string flat_key = "hdf5." + key;
        apply_json_scalar(cfg, flat_key, *this);
      }
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON object");
    }
  }

  void parse_restart_object(SimulationConfig& cfg) {
    if (!match('{')) throw ConfigError("expected JSON object for restart");
    skip_ws();
    if (match('}')) return;
    while (true) {
      std::string key = parse_string();
      if (!match(':')) throw ConfigError("expected ':' in JSON object");
      const std::string flat_key = "restart." + key;
      apply_json_scalar(cfg, flat_key, *this);
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON object");
    }
  }

  void parse_prefixed_object(SimulationConfig& cfg,
                             std::string_view prefix,
                             std::string_view name) {
    if (!match('{')) {
      throw ConfigError("expected JSON object for " + std::string(name));
    }
    skip_ws();
    if (match('}')) return;
    while (true) {
      std::string key = parse_string();
      if (!match(':')) throw ConfigError("expected ':' in JSON object");
      const std::string flat_key = std::string(prefix) + "." + key;
      apply_json_scalar(cfg, flat_key, *this);
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON object");
    }
  }

  void parse_immigration_object(SimulationConfig& cfg) {
    parse_prefixed_object(cfg, "immigration", "immigration");
  }

  void parse_initial_population_object(SimulationConfig& cfg) {
    parse_prefixed_object(cfg, "initial_population", "initial_population");
  }

  void parse_advection_object(SimulationConfig& cfg) {
    if (!match('{')) throw ConfigError("expected JSON object for advection");
    skip_ws();
    if (match('}')) return;
    while (true) {
      std::string key = parse_string();
      if (!match(':')) throw ConfigError("expected ':' in JSON object");
      if (key == "washout") {
        parse_prefixed_object(cfg, "washout", "advection.washout");
      } else {
        apply_json_scalar(cfg, "advection." + key, *this);
      }
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON object");
    }
  }

  void parse_chemistry_object(SimulationConfig& cfg) {
    if (!match('{')) throw ConfigError("expected JSON object for chemistry");
    skip_ws();
    if (match('}')) return;
    while (true) {
      std::string key = parse_string();
      if (!match(':')) throw ConfigError("expected ':' in JSON object");
      const std::string flat_key = "chemistry." + key;
      apply_json_scalar(cfg, flat_key, *this);
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON object");
    }
  }

  void parse_bacteriocin_object(SimulationConfig& cfg) {
    if (!match('{')) {
      throw ConfigError("expected JSON object for bacteriocin");
    }
    skip_ws();
    if (match('}')) return;
    while (true) {
      const std::string key = parse_string();
      if (!match(':')) throw ConfigError("expected ':' in JSON object");
      if (key == "mucin_charge") {
        parse_prefixed_object(cfg, "bacteriocin.mucin_charge",
                              "bacteriocin.mucin_charge");
      } else {
        apply_json_scalar(cfg, "bacteriocin." + key, *this);
      }
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON object");
    }
  }

  void parse_domain_object(SimulationConfig& cfg) {
    if (!match('{')) throw ConfigError("expected JSON object for domain");
    skip_ws();
    if (match('}')) return;
    while (true) {
      std::string key = parse_string();
      if (!match(':')) throw ConfigError("expected ':' in JSON object");
      if (key == "chemistry_stride") {
        if (!match('{')) {
          throw ConfigError("expected JSON object for domain.chemistry_stride");
        }
        skip_ws();
        if (match('}')) {
          throw ConfigError("domain.chemistry_stride requires x, y, and z");
        }
        while (true) {
          const std::string axis = parse_string();
          if (!match(':')) throw ConfigError("expected ':' in JSON object");
          apply_json_scalar(
              cfg, "domain.chemistry_stride." + axis, *this);
          skip_ws();
          if (match('}')) break;
          if (!match(',')) throw ConfigError("expected ',' in JSON object");
        }
      } else {
        apply_json_scalar(cfg, "domain." + key, *this);
      }
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw ConfigError("expected ',' in JSON object");
    }
  }

  std::string text_;
  size_t pos_ = 0;
};

void apply_json_scalar(SimulationConfig& cfg, const std::string& key, JsonCursor& cursor) {
  const char c = cursor.peek();
  bool handled = false;
  if (c == '"') {
    handled = InputParser::apply_flat_key(cfg, key, cursor.parse_string());
  } else if (c == 't' || c == 'f') {
    handled = InputParser::apply_flat_key(cfg, key, cursor.parse_bool() ? "true" : "false");
  } else if (std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+') {
    const Real value = cursor.parse_number();
    std::ostringstream oss;
    oss << value;
    handled = InputParser::apply_flat_key(cfg, key, oss.str());
  } else {
    cursor.skip_value();
    return;
  }
  // Keys beginning with '_' are documented as comments and handled by the
  // caller; every other scalar key that no handler recognized is a typo or a
  // version mismatch, so surface it instead of silently ignoring it.
  if (!handled && !key.empty() && key.front() != '_') {
    std::cerr << "Warning: unknown config key '" << key << "' ignored\n";
    const char* strict = std::getenv("GUTIBM_STRICT_CONFIG");
    if (strict != nullptr && strict[0] != '\0' && strict[0] != '0') {
      throw ConfigError("unknown config key '" + key + "'");
    }
  }
}

size_t find_json_object(const std::string& content) {
  const size_t pos = content.find_first_not_of(" \t\r\n");
  if (pos == std::string::npos || content[pos] != '{') return std::string::npos;
  return pos;
}

size_t find_initial_strains_array(std::string_view content) {
  const std::string key = "\"initial_strains\":";
  size_t key_pos = content.find(key);
  if (key_pos == std::string::npos) return std::string::npos;

  size_t bracket = content.find('[', key_pos + key.size());
  if (bracket == std::string::npos) return std::string::npos;
  return bracket;
}

size_t find_fixes_array(std::string_view content) {
  const std::string key = "\"fixes\":";
  size_t key_pos = content.find(key);
  if (key_pos == std::string::npos) return std::string::npos;

  size_t bracket = content.find('[', key_pos + key.size());
  if (bracket == std::string::npos) return std::string::npos;
  return bracket;
}

}  // namespace

InitialStrainsParseResult ConfigJson::parse_initial_strains(const std::string& content) {
  InitialStrainsParseResult result;
  size_t array_pos = find_initial_strains_array(content);
  if (array_pos == std::string::npos) {
    return result;
  }

  result.found = true;
  JsonCursor cursor(content.substr(array_pos));

  if (!cursor.match('[')) {
    std::cerr << "Warning: malformed initial_strains array — using default strains\n";
    result.found = false;
    return result;
  }

  cursor.skip_ws();
  if (cursor.match(']')) {
    return result;
  }

  try {
    while (true) {
      result.strains.push_back(cursor.parse_strain_object());
      cursor.skip_ws();
      if (cursor.match(']')) break;
      if (!cursor.match(',')) {
        throw ConfigError("expected ',' between strain objects");
      }
    }
  } catch (const ConfigError& ex) {
    const std::string message = ex.what();
    if (message.find("unknown receptor name") != std::string::npos) {
      throw;
    }
    std::cerr << "Warning: failed to parse initial_strains: " << ex.what()
              << " — using default strains\n";
    result.found = false;
    result.strains.clear();
  }

  return result;
}

EnabledFixesParseResult ConfigJson::parse_enabled_fixes(const std::string& content) {
  EnabledFixesParseResult result;
  size_t array_pos = find_fixes_array(content);
  if (array_pos == std::string::npos) {
    return result;
  }

  result.found = true;
  JsonCursor cursor(content.substr(array_pos));

  try {
    result.names = cursor.parse_string_array();
  } catch (const ConfigError& ex) {
    std::cerr << "Warning: failed to parse fixes: " << ex.what()
              << " — using default fix list\n";
    result.found = false;
    result.names.clear();
  }

  return result;
}

bool ConfigJson::parse_document(SimulationConfig& cfg, const std::string& content) {
  const size_t object_pos = find_json_object(content);
  if (object_pos == std::string::npos) {
    return false;
  }

  try {
    JsonCursor cursor(content.substr(object_pos));
    if (!cursor.match('{')) return false;

    cursor.skip_ws();
    if (cursor.match('}')) return true;

    while (true) {
      const std::string key = cursor.parse_string();
      if (!cursor.match(':')) throw ConfigError("expected ':' after key");

      if (!key.empty() && key.front() == '_') {
        cursor.skip_value();
      } else if (key == "initial_strains") {
        cfg.initial_strains = cursor.parse_strain_array();
      } else if (key == "plasmid_overrides") {
        cursor.parse_plasmid_overrides_object(cfg);
      } else if (key == "fixes") {
        cfg.enabled_fixes = cursor.parse_string_array();
      } else if (key == "hdf5") {
        cursor.parse_hdf5_object(cfg);
      } else if (key == "restart") {
        cursor.parse_restart_object(cfg);
      } else if (key == "immigration") {
        cursor.parse_immigration_object(cfg);
      } else if (key == "initial_population") {
        cursor.parse_initial_population_object(cfg);
      } else if (key == "closure") {
        cursor.parse_prefixed_object(cfg, "closure", "closure");
      } else if (key == "domain") {
        cursor.parse_domain_object(cfg);
      } else if (key == "advection") {
        cursor.parse_advection_object(cfg);
      } else if (key == "washout") {
        cursor.parse_prefixed_object(cfg, "washout", "washout");
      } else if (key == "chemistry") {
        cursor.parse_chemistry_object(cfg);
      } else if (key == "bacteriocin") {
        cursor.parse_bacteriocin_object(cfg);
      } else {
        apply_json_scalar(cfg, key, cursor);
      }

      cursor.skip_ws();
      if (cursor.match('}')) break;
      if (!cursor.match(',')) throw ConfigError("expected ',' between object fields");
    }

    return true;
  } catch (const ConfigError& ex) {
    const std::string message = ex.what();
    if (message.find("invalid immigration.") == 0
        || message.find("invalid initial_population.") == 0
        || message.find("invalid chemistry_decomposition") == 0
        || message.find("invalid species_subset") == 0
        || message.find("invalid toxin_evaluation") == 0
        || message.find("invalid toxin_lumping") == 0
        || message.find("invalid washout.trap") == 0
        || message.find("chemistry stride") != std::string::npos
        || message.find("chemistry_stride") != std::string::npos
        || message.find("grid_halo_width") != std::string::npos
        || message.find("unknown receptor name") != std::string::npos
        || message.find("unknown plasmid name") != std::string::npos
        || message.find("plasmid override") != std::string::npos
        || message.find("unknown config key") != std::string::npos) {
      throw;
    }
    std::cerr << "Warning: JSON config parse failed: " << ex.what()
              << " — falling back to legacy parser\n";
    return false;
  }
}

std::string ConfigJson::serialize_document(const SimulationConfig& cfg) {
  std::ostringstream out;
  out << std::setprecision(17);
  bool first = true;
  const auto key = [&out, &first](std::string_view name) {
    if (!first) out << ',';
    first = false;
    out << '"' << name << "\":";
  };
  const auto string_value = [&out](const std::string& value) {
    out << '"';
    for (const char c : value) {
      if (c == '\\' || c == '"') out << '\\';
      if (c == '\n') out << "\\n";
      else if (c == '\r') out << "\\r";
      else if (c == '\t') out << "\\t";
      else out << c;
    }
    out << '"';
  };
  const auto string_key = [&key, &string_value](std::string_view name,
                                                 const std::string& value) {
    key(name);
    string_value(value);
  };
  const auto real_key = [&key, &out](std::string_view name, Real value) {
    key(name);
    out << value;
  };
  const auto int_key = [&key, &out](std::string_view name, auto value) {
    key(name);
    out << value;
  };
  const auto bool_key = [&key, &out](std::string_view name, bool value) {
    key(name);
    out << (value ? "true" : "false");
  };

  out << '{';
  real_key("total_time", cfg.time.total_time);
  real_key("bio_dt", cfg.time.bio_dt);
  real_key("output_interval", cfg.time.output_interval);
  int_key("seed", cfg.seed);
  real_key("grid_dx", cfg.domain.grid_dx);
  real_key("domain_x", cfg.domain.hi[0]);
  real_key("domain_y", cfg.domain.hi[1]);
  real_key("domain_z", cfg.domain.hi[2]);
  real_key("hash_cell_size", cfg.domain.hash_cell_size);
  real_key("ghost_width", cfg.domain.ghost_width);
  int_key("domain.chemistry_stride_x", cfg.domain.chemistry_stride[0]);
  int_key("domain.chemistry_stride_y", cfg.domain.chemistry_stride[1]);
  int_key("domain.chemistry_stride_z", cfg.domain.chemistry_stride[2]);
  int_key("domain.grid_halo_width", cfg.domain.grid_halo_width);
  string_key("hdf5.file", cfg.hdf5.filename);
  if (!cfg.checkpoint.file.empty()) {
    string_key("checkpoint_file", cfg.checkpoint.file);
  }
  if (!cfg.checkpoint.step.empty()) {
    string_key("checkpoint_step", cfg.checkpoint.step);
  }
  string_key("initial_population.placement", cfg.initial_population.placement);
  real_key("initial_population.z_min", cfg.initial_population.z_min);
  real_key("initial_population.z_max", cfg.initial_population.z_max);
  bool_key("immigration.enabled", cfg.immigration.enabled);
  int_key("immigration.count", cfg.immigration.count);
  int_key("immigration.strain_index", cfg.immigration.strain_index);
  string_key("immigration.placement", cfg.immigration.placement);
  real_key("immigration.distance", cfg.immigration.distance);
  real_key("immigration.distance_tolerance", cfg.immigration.distance_tolerance);
  string_key("immigration.distance_reference", cfg.immigration.distance_reference);
  real_key("immigration.z_min", cfg.immigration.z_min);
  real_key("immigration.z_max", cfg.immigration.z_max);
  string_key("immigration.schedule", cfg.immigration.schedule);
  int_key("immigration.step", cfg.immigration.step);
  real_key("immigration.rate", cfg.immigration.rate);
  string_key("chemistry_decomposition", cfg.chemistry_decomposition);
  string_key("species_subset", cfg.species_subset);
  string_key("uptake_limit", cfg.fixes.metabolism.uptake_limit);
  string_key("chemistry.toxin_evaluation", cfg.qssa.toxin_evaluation);
  string_key("chemistry.toxin_lumping", cfg.qssa.toxin_lumping);
  real_key("mucus_thickness", cfg.advection.mucus_thickness);
  real_key("radial_turnover", cfg.advection.radial_turnover);
  string_key("washout.trap", cfg.advection.washout_trap
      == WashoutTrapMode::IMPOSED ? "imposed" : "emergent");
  real_key("distal_transit", cfg.advection.distal_transit_time);
  real_key("distal_length", cfg.advection.distal_length);
  real_key("profile_alpha", cfg.advection.profile_alpha);
  bool_key("taylor_aris_enabled", cfg.advection.taylor_aris_enabled);
  bool_key("peristaltic_enabled", cfg.advection.peristaltic_enabled);
  real_key("peristaltic_period", cfg.advection.peristaltic_period);
  real_key("peristaltic_amplitude", cfg.advection.peristaltic_amplitude);
  real_key("peristaltic_wavelength", cfg.advection.peristaltic_wavelength);
  bool_key("crypts_enabled", cfg.advection.crypts_enabled);
  real_key("crypt_depth", cfg.advection.crypt_depth);
  real_key("crypt_exit_rate", cfg.advection.crypt_exit_rate);
  real_key("crypt_entry_rate", cfg.advection.crypt_entry_rate);
  int_key("crypt_carrying_capacity", cfg.advection.crypt_carrying_capacity);
  real_key("vbf_density", cfg.vbf.density);
  real_key("vbf_viscosity", cfg.vbf.viscosity);
  real_key("vbf_drag_coeff", cfg.vbf.drag_coeff);
  real_key("vbf_nutrient_sink", cfg.vbf.nutrient_sink);
  real_key("vbf_mucin_liberation", cfg.vbf.mucin_liberation);
  real_key("vbf_carrying_cap", cfg.vbf.carrying_cap);
  bool_key("vbf_mucin_z_gradient", cfg.vbf.mucin_z_gradient_enabled);
  real_key("vbf_mucin_z_lambda", cfg.vbf.mucin_z_gradient_lambda);
  real_key("vbf_carbon_sink_vmax", cfg.vbf.carbon_sink_vmax);
  real_key("vbf_carbon_sink_km", cfg.vbf.carbon_sink_km);
  real_key("vbf.agent_carbon_coupling", cfg.vbf.agent_carbon_coupling);
  real_key("carbon_boundary_conc", cfg.carbon_boundary_conc);
  bool_key("oxygen.enabled", cfg.chem_env.oxygen.enabled);
  real_key("oxygen.epithelial_conc", cfg.chem_env.oxygen.epithelial_conc);
  bool_key("oxygen.delivery_uptake_enabled",
           cfg.chem_env.oxygen.delivery_uptake_enabled);
  string_key("oxygen.respiration_driver",
             cfg.chem_env.oxygen.respiration_driver);
  real_key("oxygen.D_free", cfg.chem_env.oxygen.D_free);
  real_key("oxygen.Km", cfg.chem_env.oxygen.Km);
  bool_key("oxygen.metabolic_switch_enabled",
           cfg.chem_env.oxygen.metabolic_switch_enabled);
  real_key("oxygen.mu_crit", cfg.chem_env.oxygen.mu_crit);
  real_key("oxygen.aerobic_mu_factor", cfg.chem_env.oxygen.aerobic_mu_factor);
  real_key("oxygen.anaerobic_mu_factor",
           cfg.chem_env.oxygen.anaerobic_mu_factor);
  real_key("oxygen.aerobic_carbon_cost_factor",
           cfg.chem_env.oxygen.aerobic_carbon_cost_factor);
  real_key("oxygen.anaerobic_carbon_cost_factor",
           cfg.chem_env.oxygen.anaerobic_carbon_cost_factor);
  real_key("oxygen.tau_metabolic_switch",
           cfg.chem_env.oxygen.tau_metabolic_switch);
  real_key("oxygen.ferm_acid_yield", cfg.chem_env.oxygen.ferm_acid_yield);
  real_key("oxygen.anaerobic_maintenance_factor",
           cfg.chem_env.oxygen.anaerobic_maintenance_factor);
  real_key("oxygen.boost_max", cfg.chem_env.oxygen.boost_max);
  real_key("oxygen.q_consumption", cfg.chem_env.oxygen.q_consumption);
  real_key("oxygen.q_maintenance", cfg.chem_env.oxygen.q_maintenance);
  real_key("oxygen.vbf_sink", cfg.chem_env.oxygen.vbf_sink);
  real_key("oxygen.k_ROS", cfg.chem_env.oxygen.k_ROS);
  bool_key("acetate.enabled", cfg.chem_env.acetate.enabled);
  real_key("acetate.D_free", cfg.chem_env.acetate.D_free);
  real_key("acetate.vbf_production", cfg.chem_env.acetate.vbf_production);
  real_key("acetate.vbf_consumption", cfg.chem_env.acetate.vbf_consumption);
  real_key("acetate.overflow_threshold", cfg.chem_env.acetate.overflow_threshold);
  real_key("acetate.overflow_rate", cfg.chem_env.acetate.overflow_rate);
  real_key("acetate.scavenge_rate", cfg.chem_env.acetate.scavenge_rate);
  real_key("acetate.scavenge_Km", cfg.chem_env.acetate.scavenge_Km);
  real_key("acetate.epithelial_uptake", cfg.chem_env.acetate.epithelial_uptake);
  bool_key("mucin.enabled", cfg.chem_env.mucin.enabled);
  real_key("mucin.initial_conc", cfg.chem_env.mucin.initial_conc);
  real_key("mucin.secretion_rate", cfg.chem_env.mucin.secretion_rate);
  real_key("mucin.Km_degradation", cfg.chem_env.mucin.Km_degradation);
  real_key("mucin.k_liberation", cfg.chem_env.mucin.k_liberation);
  bool_key("protease.enabled", cfg.chem_env.protease.enabled);
  real_key("protease.default_half_life", cfg.chem_env.protease.default_half_life);
  real_key("protease.dilution_rate", cfg.chem_env.protease.dilution_rate);
  bool_key("siderophore.enabled", cfg.chem_env.siderophore.enabled);
  real_key("siderophore.secretion_rate", cfg.chem_env.siderophore.secretion_rate);
  real_key("siderophore.D_free", cfg.chem_env.siderophore.D_free);
  real_key("siderophore.chelation_rate", cfg.chem_env.siderophore.chelation_rate);
  real_key("siderophore.Km_reimport", cfg.chem_env.siderophore.Km_reimport);
  real_key("siderophore.Vmax_reimport", cfg.chem_env.siderophore.Vmax_reimport);
  bool_key("ferrichrome.enabled", cfg.chem_env.ferrichrome.enabled);
  real_key("ferrichrome.initial_conc", cfg.chem_env.ferrichrome.initial_conc);
  real_key("ferrichrome.boundary_conc", cfg.chem_env.ferrichrome.boundary_conc);
  bool_key("quorum_sensing.enabled", cfg.quorum_sensing.enabled);
  real_key("quorum_sensing.ai2_basal_rate", cfg.quorum_sensing.ai2_basal_rate);
  real_key("quorum_sensing.ai2_growth_coupled", cfg.quorum_sensing.ai2_growth_coupled);
  real_key("quorum_sensing.lsr_vmax", cfg.quorum_sensing.lsr_vmax);
  real_key("quorum_sensing.lsr_km", cfg.quorum_sensing.lsr_km);
  real_key("quorum_sensing.ai2_D_free", cfg.quorum_sensing.ai2_D_free);
  real_key("quorum_sensing.ai2_decay_rate", cfg.quorum_sensing.ai2_decay_rate);
  bool_key("quorum_sensing.ai2_chemotaxis", cfg.quorum_sensing.ai2_chemotaxis_enabled);
  real_key("quorum_sensing.chi_ai2", cfg.quorum_sensing.chi_ai2);
  real_key("toxin_cutoff", cfg.qssa.toxin_cutoff);
  real_key("nutrient_cutoff", cfg.qssa.nutrient_cutoff);
  real_key("colicin_release_rate", cfg.qssa.colicin_release_rate);
  real_key("microcin_secretion", cfg.qssa.microcin_secretion);
  bool_key("use_fmm", cfg.qssa.use_fmm);
  real_key("fmm_theta", cfg.qssa.fmm_theta);
  int_key("fmm_expansion_order", cfg.qssa.fmm_expansion_order);
  real_key("division_threshold", cfg.fixes.metabolism.division_threshold);
  real_key("bacteriostasis_threshold", cfg.fixes.metabolism.bacteriostasis_threshold);
  real_key("maintenance_rate", cfg.fixes.metabolism.maintenance_rate);
  real_key("carbon_maintenance_rate",
           cfg.fixes.metabolism.carbon_maintenance_rate);
  bool_key("acid_inhibition_enabled",
           cfg.fixes.metabolism.acid_inhibition_enabled);
  real_key("acid_inhibition_max", cfg.fixes.metabolism.acid_inhibition_max);
  real_key("acid_inhibition_Ki", cfg.fixes.metabolism.acid_inhibition_Ki);
  real_key("acetate_pKa", cfg.fixes.metabolism.acetate_pKa);
  real_key("metE_penalty", cfg.fixes.metabolism.metE_penalty);
  real_key("metE_acetate_km", cfg.fixes.metabolism.metE_acetate_km);
  real_key("metE_acetate_max_factor", cfg.fixes.metabolism.metE_acetate_max_factor);
  real_key("eut_km", cfg.fixes.metabolism.eut_km);
  real_key("eut_max_penalty", cfg.fixes.metabolism.eut_max_penalty);
  real_key("km_iron_primary", cfg.fixes.metabolism.km_iron_primary);
  real_key("km_iron_iroN", cfg.fixes.metabolism.km_iron_iroN);
  real_key("km_iron_iutA", cfg.fixes.metabolism.km_iron_iutA);
  real_key("km_iron_fiu", cfg.fixes.metabolism.km_iron_fiu);
  real_key("kd_b12_btuB", cfg.fixes.receptor.kd_b12_btuB);
  real_key("kd_colicinE_btuB", cfg.fixes.receptor.kd_colicinE_btuB);
  real_key("kd_enterobactin", cfg.fixes.receptor.kd_enterobactin);
  real_key("kd_colicinB_fepA", cfg.fixes.receptor.kd_colicinB_fepA);
  real_key("kd_lin_enterobactin", cfg.fixes.receptor.kd_lin_enterobactin);
  real_key("kd_colicinIa_cirA", cfg.fixes.receptor.kd_colicinIa_cirA);
  real_key("kill_rate_colicin", cfg.fixes.receptor.kill_rate_colicin);
  real_key("kill_rate_microcin", cfg.fixes.receptor.kill_rate_microcin);
  real_key("immunity_factor", cfg.fixes.receptor.immunity_factor);
  real_key("sos_lysis_prob", cfg.fixes.bacteriocin.sos_lysis_prob);
  real_key("sos_basal_rate", cfg.fixes.bacteriocin.sos_basal_rate);
  real_key("D_free_colicin", cfg.fixes.bacteriocin.D_free_colicin);
  real_key("burst_release_tau", cfg.fixes.bacteriocin.burst_release_tau);
  real_key("microcin_mu_penalty", cfg.fixes.bacteriocin.microcin_mu_penalty);
  real_key("sos_cross_induction_rate", cfg.fixes.bacteriocin.sos_cross_induction_rate);
  real_key("bacteriocin.mucin_charge.r_min",
           cfg.fixes.bacteriocin.mucin_charge.r_min);
  real_key("bacteriocin.mucin_charge.amplitude",
           cfg.fixes.bacteriocin.mucin_charge.amplitude);
  real_key("bacteriocin.mucin_charge.dz_half",
           cfg.fixes.bacteriocin.mucin_charge.dz_half);
  real_key("bacteriocin.mucin_charge.width",
           cfg.fixes.bacteriocin.mucin_charge.width);
  real_key("bacteriocin.mucin_charge.ph",
           cfg.fixes.bacteriocin.mucin_charge.ph);
  real_key("pili_length", cfg.fixes.conjugation.pili_length);
  real_key("base_transfer_rate", cfg.fixes.conjugation.base_transfer_rate);
  real_key("shear_critical", cfg.fixes.conjugation.shear_critical);
  real_key("plasmid_copy_cost", cfg.fixes.conjugation.plasmid_copy_cost);
  bool_key("pili_heterogeneity", cfg.fixes.conjugation.pili_heterogeneity);
  real_key("pili_length_min", cfg.fixes.conjugation.pili_length_min);
  real_key("pili_length_max", cfg.fixes.conjugation.pili_length_max);
  real_key("bi_duplication_rate", cfg.fixes.mutation.bi_duplication_rate);
  real_key("bi_recombination_rate", cfg.fixes.mutation.bi_recombination_rate);
  real_key("receptor_mutation_rate", cfg.fixes.mutation.receptor_mutation_rate);
  real_key("super_killer_rate", cfg.fixes.mutation.super_killer_rate);
  real_key("partial_resistance_rate", cfg.fixes.mutation.partial_resistance_rate);
  real_key("receptor_reduction", cfg.fixes.mutation.receptor_reduction);
  int_key("max_bi_loci", cfg.fixes.mutation.max_bi_loci);
  real_key("immunity_escape_prob", cfg.fixes.mutation.immunity_escape_prob);
  real_key("escape_affinity_lo", cfg.fixes.mutation.escape_affinity_lo);
  real_key("escape_affinity_hi", cfg.fixes.mutation.escape_affinity_hi);
  real_key("compensatory_rate", cfg.fixes.mutation.compensatory_rate);
  real_key("compensatory_reduction", cfg.fixes.mutation.compensatory_reduction);
  real_key("hertz_k", cfg.fixes.mechanics.hertz_k);
  bool_key("hertzian_enabled", cfg.fixes.mechanics.hertzian_enabled);
  bool_key("adhesion_enabled", cfg.fixes.mechanics.adhesion_enabled);
  real_key("adhesion_strength", cfg.fixes.mechanics.adhesion_strength);
  real_key("adhesion_range", cfg.fixes.mechanics.adhesion_range);
  bool_key("gpu_enabled", cfg.gpu.enabled);
  int_key("gpu_device_id", cfg.gpu.device_id);
  bool_key("adaptive_dt_enabled", cfg.adaptive_dt.enabled);
  real_key("dt_min", cfg.adaptive_dt.min);
  real_key("dt_max", cfg.adaptive_dt.max);
  real_key("dt_safety", cfg.adaptive_dt.safety);
  real_key("dt_growth_limit", cfg.adaptive_dt.growth_limit);
  bool_key("profile_steps", cfg.profile_steps);
  real_key("dysbiosis_threshold", cfg.dysbiosis_threshold);
  real_key("dysbiosis_sampling_interval", cfg.dysbiosis_sampling_interval);
  int_key("dysbiosis_sample_count", cfg.dysbiosis_sample_count);
  bool_key("closure.enforce_delivery_realization",
           cfg.closure.enforce_delivery_realization);
  int_key("closure.zero_realization_grace_steps",
          cfg.closure.zero_realization_grace_steps);
  bool_key("closure.enforce_reaction_clip",
           cfg.closure.enforce_reaction_clip);
  real_key("closure.reaction_clip_tolerance_fraction",
           cfg.closure.reaction_clip_tolerance_fraction);
  bool_key("restart.enabled", cfg.restart.enabled);
  string_key("restart.directory", cfg.restart.directory);
  int_key("restart.interval_steps", cfg.restart.interval_steps);
  bool_key("hdf5.enabled", cfg.hdf5.enabled);
  string_key("hdf5.compression", cfg.hdf5.compression);
  int_key("hdf5.compression_level", cfg.hdf5.compression_level);
  int_key("hdf5.schedule.summary", cfg.hdf5.schedule.summary);
  int_key("hdf5.schedule.agents", cfg.hdf5.schedule.agents);
  int_key("hdf5.schedule.grid", cfg.hdf5.schedule.grid);
  int_key("hdf5.schedule.lineage", cfg.hdf5.schedule.lineage);
  int_key("hdf5.schedule.genome", cfg.hdf5.schedule.genome);
  int_key("hdf5.schedule.provenance", cfg.hdf5.schedule.provenance);
  real_key("b12.initial_conc", cfg.b12_initial_conc);
  out << ",\"hdf5.schedule.grid_species\":[";
  for (size_t i = 0; i < cfg.hdf5.schedule.grid_species.size(); ++i) {
    if (i != 0) out << ',';
    string_value(cfg.hdf5.schedule.grid_species[i]);
  }
  out << ']';
  out << ",\"initial_strains\":[";
  for (size_t i = 0; i < cfg.initial_strains.size(); ++i) {
    if (i != 0) out << ',';
    const auto& strain = cfg.initial_strains[i];
    out << "{\"type\":" << strain.type
        << ",\"count\":" << strain.count
        << ",\"mu_max\":" << strain.mu_max
        << ",\"conjugative\":" << (strain.conjugative ? "true" : "false")
        << ",\"cdi_type\":" << strain.cdi_type
        << ",\"cdi_immunity\":" << strain.cdi_immunity
        << ",\"plasmids\":[";
    for (size_t j = 0; j < strain.plasmids.size(); ++j) {
      if (j != 0) out << ',';
      string_value(strain.plasmids[j]);
    }
    out << "],\"receptor_expression\":{";
    size_t receptor_index = 0;
    for (const auto& [name, expression] : strain.receptor_expression) {
      if (receptor_index++ != 0) out << ',';
      string_value(name);
      out << ':' << expression;
    }
    out << "}}";
  }
  out << "],\"plasmid_overrides\":{";
  size_t plasmid_index = 0;
  for (const auto& [name, values] : cfg.plasmid_overrides) {
    if (plasmid_index++ != 0) out << ',';
    string_value(name);
    out << ":{";
    bool first_value = true;
    const auto emit_override = [&out, &first_value, &string_value](
                                   std::string_view field,
                                   const std::optional<Real>& value) {
      if (!value.has_value()) return;
      if (!first_value) out << ',';
      first_value = false;
      string_value(std::string(field));
      out << ':' << *value;
    };
    emit_override("retardation", values.retardation);
    emit_override("diff_coeff", values.diff_coeff);
    emit_override("burst_size", values.burst_size);
    out << '}';
  }
  out << "},\"fixes\":[";
  for (size_t i = 0; i < cfg.enabled_fixes.size(); ++i) {
    if (i != 0) out << ',';
    string_value(cfg.enabled_fixes[i]);
  }
  out << "],\"disabled_fixes\":[";
  for (size_t i = 0; i < cfg.disabled_fixes.size(); ++i) {
    if (i != 0) out << ',';
    string_value(cfg.disabled_fixes[i]);
  }
  out << "],\"disabled_mechanisms\":[";
  for (size_t i = 0; i < cfg.disabled_mechanisms.size(); ++i) {
    if (i != 0) out << ',';
    string_value(cfg.disabled_mechanisms[i]);
  }
  out << "]}";
  return out.str();
}

}  // namespace gutibm
