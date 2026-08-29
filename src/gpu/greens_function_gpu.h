#ifndef GUTIBM_GREENS_FUNCTION_GPU_H
#define GUTIBM_GREENS_FUNCTION_GPU_H

#include "types.h"
#include "greens_function.h"
#include "robin_correction_table.h"
#include <cstddef>
#include <memory>
#include <vector>

namespace gutibm {

class Domain;
class AdvectionField;

constexpr size_t kMaximumRobinDeviceTables = 256;

std::vector<int> make_robin_launch_table_indices(
    const std::vector<std::shared_ptr<const robin::Table>>& source_tables,
    std::vector<std::shared_ptr<const robin::Table>>& launch_tables);

std::vector<size_t> robin_host_fallback_sources(
    const Domain& domain,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params);

// GPU-accelerated Green's function superposition. Returns false if GPU unavailable.
bool gpu_superpose_to_grid(
    const Domain& domain,
    const AdvectionField& adv,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    std::vector<Real>& grid_conc,
    Real cutoff_radius, uint64_t* cap_hits = nullptr,
    uint64_t* kernel_evaluations = nullptr,
    uint64_t* low_screening_evaluations = nullptr,
    uint64_t* negative_field_count = nullptr,
    Real* most_negative_field = nullptr);

// Deposit near-field superposition directly into a device concentration buffer.
bool gpu_superpose_to_device(
    const Domain& domain,
    const AdvectionField& adv,
    const std::vector<Vec3>& sources,
    const std::vector<GreensFunctionParams>& params,
    const std::vector<Real>& strength_factors,
    double* d_grid_conc,
    Real cutoff_radius, uint64_t* cap_hits = nullptr,
    uint64_t* kernel_evaluations = nullptr,
    uint64_t* low_screening_evaluations = nullptr,
    uint64_t* negative_field_count = nullptr,
    Real* most_negative_field = nullptr);

}  // namespace gutibm

#endif  // GUTIBM_GREENS_FUNCTION_GPU_H
