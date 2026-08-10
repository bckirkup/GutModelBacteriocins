#ifndef GUTIBM_IMMIGRATION_CONFIG_H
#define GUTIBM_IMMIGRATION_CONFIG_H

#include "types.h"

#include <string>

namespace gutibm {

struct ImmigrationConfig {
  bool enabled = false;
  Int count = 1;
  Int strain_index = 0;
  std::string placement = "uniform";
  Real distance = 0.0;
  Real distance_tolerance = 0.0;
  std::string distance_reference = "nearest_agent";
  Real z_min = 0.0;
  Real z_max = 0.0;
  std::string schedule = "pulse";
  Int step = 0;
  Real rate = 0.0;
};

}  // namespace gutibm

#endif  // GUTIBM_IMMIGRATION_CONFIG_H
