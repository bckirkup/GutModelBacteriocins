#ifndef GUTIBM_GPU_DIAGNOSTIC_FORMAT_H
#define GUTIBM_GPU_DIAGNOSTIC_FORMAT_H

#include "types.h"
#include "chemical_field.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <charconv>
#include <iostream>
#include <limits>
#include <string>

namespace gutibm::gpu_diagnostic {

inline std::string format_real(Real value) {
  std::array<char, 64> buffer{};
  const auto result = std::to_chars(
      buffer.data(), buffer.data() + buffer.size(), value,
      std::chars_format::general, std::numeric_limits<Real>::max_digits10);
  assert(result.ec == std::errc{});
  return std::string(buffer.data(), result.ptr);
}

struct ConcentrationSummary {
  Real minimum = std::numeric_limits<Real>::infinity();
  Real mean = 0.0;
  Real maximum = -std::numeric_limits<Real>::infinity();
};

inline ConcentrationSummary summarize_species(const ChemicalField& chem,
                                              Int species) {
  ConcentrationSummary summary;
  const Int cells = chem.global_ncells();
  for (Int cell = 0; cell < cells; ++cell) {
    const Real value = chem.conc_global(species, cell);
    summary.minimum = std::min(summary.minimum, value);
    summary.maximum = std::max(summary.maximum, value);
    summary.mean += value;
  }
  summary.mean /= static_cast<Real>(cells);
  return summary;
}

inline void print_concentration_diagnostics(
    const char* test_name, const char* label, const ChemicalField& cpu,
    const ChemicalField& gpu) {
  std::cerr << "[gpu_diag][" << test_name << "]";
  if (label != nullptr && label[0] != '\0') {
    std::cerr << "[" << label << "]";
  }
  std::cerr << "[species]\n";
  for (Int species = 0; species < cpu.num_species(); ++species) {
    const ConcentrationSummary cpu_summary = summarize_species(cpu, species);
    const ConcentrationSummary gpu_summary = summarize_species(gpu, species);
    std::cerr << "  name=" << cpu.spec(species).name
              << " cpu_min=" << format_real(cpu_summary.minimum)
              << " cpu_mean=" << format_real(cpu_summary.mean)
              << " cpu_max=" << format_real(cpu_summary.maximum)
              << " gpu_min=" << format_real(gpu_summary.minimum)
              << " gpu_mean=" << format_real(gpu_summary.mean)
              << " gpu_max=" << format_real(gpu_summary.maximum) << "\n";
  }
}

}  // namespace gutibm::gpu_diagnostic

#endif  // GUTIBM_GPU_DIAGNOSTIC_FORMAT_H
