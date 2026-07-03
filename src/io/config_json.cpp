/* -----------------------------------------------------------------------
   GutIBM – Minimal JSON helpers for config fragments
   ----------------------------------------------------------------------- */

#include "config_json.h"

#include <cctype>
#include <iostream>
#include <sstream>
#include <string_view>
#include "error.h"
#include <utility>

namespace gutibm {

namespace {

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
        else if (esc == 'n') out.push_back('\n');
        else if (esc == 't') out.push_back('\t');
        else throw ConfigError("unsupported JSON escape");
      } else {
        out.push_back(c);
      }
    }
    throw ConfigError("unterminated JSON string");
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
    strain.mu_max = 5.0e-4;
    strain.conjugative = false;

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

  std::string text_;
  size_t pos_ = 0;
};

void apply_json_scalar(SimulationConfig& cfg, const std::string& key, JsonCursor& cursor) {
  const char c = cursor.peek();
  if (c == '"') {
    InputParser::apply_flat_key(cfg, key, cursor.parse_string());
  } else if (c == 't' || c == 'f') {
    InputParser::apply_flat_key(cfg, key, cursor.parse_bool() ? "true" : "false");
  } else if (std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+') {
    const Real value = cursor.parse_number();
    std::ostringstream oss;
    oss << value;
    InputParser::apply_flat_key(cfg, key, oss.str());
  } else {
    cursor.skip_value();
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
      } else if (key == "fixes") {
        cfg.enabled_fixes = cursor.parse_string_array();
      } else {
        apply_json_scalar(cfg, key, cursor);
      }

      cursor.skip_ws();
      if (cursor.match('}')) break;
      if (!cursor.match(',')) throw ConfigError("expected ',' between object fields");
    }

    return true;
  } catch (const ConfigError& ex) {
    std::cerr << "Warning: JSON config parse failed: " << ex.what()
              << " — falling back to legacy parser\n";
    return false;
  }
}

}  // namespace gutibm
