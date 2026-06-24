/* -----------------------------------------------------------------------
   GutIBM – Minimal JSON helpers for config fragments
   ----------------------------------------------------------------------- */

#include "config_json.h"

#include <cctype>
#include <iostream>
#include <stdexcept>

namespace gutibm {

namespace {

class JsonCursor {
 public:
  explicit JsonCursor(const std::string& text) : text_(text) {}

  void skip_ws() {
    while (pos_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_]))) {
      ++pos_;
    }
  }

  bool eof() const { return pos_ >= text_.size(); }

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
      throw std::runtime_error("expected JSON string");
    }
    ++pos_;
    std::string out;
    while (pos_ < text_.size()) {
      char c = text_[pos_++];
      if (c == '"') return out;
      if (c == '\\') {
        if (pos_ >= text_.size()) throw std::runtime_error("truncated JSON escape");
        char esc = text_[pos_++];
        if (esc == '"' || esc == '\\' || esc == '/') out.push_back(esc);
        else if (esc == 'n') out.push_back('\n');
        else if (esc == 't') out.push_back('\t');
        else throw std::runtime_error("unsupported JSON escape");
      } else {
        out.push_back(c);
      }
    }
    throw std::runtime_error("unterminated JSON string");
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
    if (start == pos_) throw std::runtime_error("expected JSON number");
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
    throw std::runtime_error("expected JSON boolean");
  }

  std::vector<std::string> parse_string_array() {
    std::vector<std::string> out;
    if (!match('[')) throw std::runtime_error("expected JSON array");
    skip_ws();
    if (match(']')) return out;
    while (true) {
      out.push_back(parse_string());
      skip_ws();
      if (match(']')) break;
      if (!match(',')) throw std::runtime_error("expected ',' in JSON array");
    }
    return out;
  }

  SimulationConfig::InitialStrain parse_strain_object() {
    if (!match('{')) throw std::runtime_error("expected JSON object");

    SimulationConfig::InitialStrain strain{};
    strain.mu_max = 5.0e-4;
    strain.conjugative = false;

    skip_ws();
    if (match('}')) return strain;

    while (true) {
      std::string key = parse_string();
      if (!match(':')) throw std::runtime_error("expected ':' in JSON object");

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
      } else {
        skip_value();
      }

      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw std::runtime_error("expected ',' in JSON object");
    }

    return strain;
  }

  void skip_object() {
    if (!match('{')) throw std::runtime_error("expected JSON object");
    skip_ws();
    if (match('}')) return;
    while (true) {
      (void)parse_string();
      if (!match(':')) throw std::runtime_error("expected ':' in JSON object");
      skip_value();
      skip_ws();
      if (match('}')) break;
      if (!match(',')) throw std::runtime_error("expected ',' in JSON object");
    }
  }

  void skip_array() {
    if (!match('[')) throw std::runtime_error("expected JSON array");
    skip_ws();
    if (match(']')) return;
    while (true) {
      skip_value();
      skip_ws();
      if (match(']')) break;
      if (!match(',')) throw std::runtime_error("expected ',' in JSON array");
    }
  }

 private:
  void skip_value() {
    skip_ws();
    if (pos_ >= text_.size()) throw std::runtime_error("unexpected end of JSON");

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
    throw std::runtime_error("unsupported JSON value");
  }

  const std::string& text_;
  size_t pos_ = 0;
};

size_t find_initial_strains_array(const std::string& content) {
  const std::string key = "\"initial_strains\":";
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
  std::string array_text = content.substr(array_pos);
  JsonCursor cursor(array_text);

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
        throw std::runtime_error("expected ',' between strain objects");
      }
    }
  } catch (const std::exception& ex) {
    std::cerr << "Warning: failed to parse initial_strains: " << ex.what()
              << " — using default strains\n";
    result.found = false;
    result.strains.clear();
  }

  return result;
}

}  // namespace gutibm
