/* -----------------------------------------------------------------------
   GutIBM – Minimal JSON helpers for config fragments
   Parses the initial_strains array from GutIBM input files.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_CONFIG_JSON_H
#define GUTIBM_CONFIG_JSON_H

#include "input_parser.h"
#include <string>
#include <vector>

namespace gutibm {

struct InitialStrainsParseResult {
  bool found = false;
  std::vector<SimulationConfig::InitialStrain> strains;
};

struct EnabledFixesParseResult {
  bool found = false;
  std::vector<std::string> names;
};

class ConfigJson {
 public:
  // Parse initial_strains from full file content. Returns found=false when
  // the key is absent; found=true with strains (possibly empty) on success.
  static InitialStrainsParseResult parse_initial_strains(const std::string& content);

  // Parse fixes from full file content. Returns found=false when absent.
  static EnabledFixesParseResult parse_enabled_fixes(const std::string& content);

  // Parse a strict JSON config document into cfg. Returns true on success.
  static bool parse_document(SimulationConfig& cfg, const std::string& content);
};

}  // namespace gutibm

#endif  // GUTIBM_CONFIG_JSON_H
