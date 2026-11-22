/* -----------------------------------------------------------------------
   GutIBM – ChemicalField multi-rank host memory estimator, pre-allocation
   OOM guard, auto-slab policy, and PinnedHostBuffer unit tests (#480).
   ----------------------------------------------------------------------- */

#include "chemical_field.h"
#include "config_json.h"
#include "device_memory.h"
#include "error.h"
#include "gpu_profile.h"
#include "input_parser.h"
#include "simulation.h"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>

using namespace gutibm;

namespace {

constexpr double kBytesPerGiB = 1024.0 * 1024.0 * 1024.0;
constexpr Int kShippedFullSpecies = 11;

void test_memory_footprint_exact_arithmetic() {
  constexpr uint64_t kCells1000 = 1000ULL * 1000ULL * 50ULL;
  constexpr uint64_t kCells500 = 500ULL * 500ULL * 50ULL;
  constexpr uint64_t kExpectedPerRank1000 =
      2ULL * static_cast<uint64_t>(kShippedFullSpecies) * kCells1000
      * sizeof(Real);
  constexpr uint64_t kExpectedPerRank500 =
      2ULL * static_cast<uint64_t>(kShippedFullSpecies) * kCells500
      * sizeof(Real);

  // Verify 1000x1000x50 analytical footprint: ~8.20 GiB / rank, ~16.39 GiB / 2,
  // ~32.78 GiB / 4, ~65.57 GiB / 8 in replicated mode without delivery.
  const std::array<Int, 4> rank_counts = {1, 2, 4, 8};
  for (const Int nprocs : rank_counts) {
    const auto rep1000 = ChemicalField::estimate_memory_footprint(
        1000, 1000, 50, kShippedFullSpecies, nprocs, nprocs, false, false);
    assert(rep1000.global_cells == kCells1000);
    assert(rep1000.storage_cells_per_rank == kCells1000);
    assert(rep1000.conc_bytes_per_rank
           == static_cast<uint64_t>(kShippedFullSpecies) * kCells1000
                  * sizeof(Real));
    assert(rep1000.reac_bytes_per_rank == rep1000.conc_bytes_per_rank);
    assert(rep1000.delivery_bytes_per_rank == 0ULL);
    assert(rep1000.scratch_bytes_per_rank == 0ULL);
    assert(rep1000.total_bytes_per_rank == kExpectedPerRank1000);
    assert(rep1000.total_bytes_per_node
           == kExpectedPerRank1000 * static_cast<uint64_t>(nprocs));

    const double gib_per_rank =
        static_cast<double>(rep1000.total_bytes_per_rank) / kBytesPerGiB;
    const double gib_per_node =
        static_cast<double>(rep1000.total_bytes_per_node) / kBytesPerGiB;
    assert(gib_per_rank > 8.0 && gib_per_rank < 8.4);
    if (nprocs == 2) {
      assert(gib_per_node > 16.0 && gib_per_node < 16.8);
    }
    if (nprocs == 8) {
      assert(gib_per_node > 65.0 && gib_per_node < 66.0);
    }

    // 500x500x50 is strictly 1/4 of 1000x1000x50.
    const auto rep500 = ChemicalField::estimate_memory_footprint(
        500, 500, 50, kShippedFullSpecies, nprocs, nprocs, false, false);
    assert(rep500.total_bytes_per_rank == kExpectedPerRank500);
    assert(rep500.total_bytes_per_rank * 4ULL == rep1000.total_bytes_per_rank);
    assert(rep500.total_bytes_per_node * 4ULL == rep1000.total_bytes_per_node);

    // Slab mode with zero halo divides per-rank footprint cleanly by nprocs
    // (1000 and 500 are divisible by 1, 2, 4, 8 except 500/8=62 rem 4 -> max 63).
    const auto slab1000 = ChemicalField::estimate_memory_footprint(
        1000, 1000, 50, kShippedFullSpecies, nprocs, nprocs, true, false, 0);
    assert(slab1000.total_bytes_per_rank
           == kExpectedPerRank1000 / static_cast<uint64_t>(nprocs));
    assert(slab1000.total_bytes_per_node == kExpectedPerRank1000);

    // Slab with halo=2 still yields strong savings for nprocs >= 2.
    const auto slab1000_halo = ChemicalField::estimate_memory_footprint(
        1000, 1000, 50, kShippedFullSpecies, nprocs, nprocs, true, false, 2);
    if (nprocs > 1) {
      assert(slab1000_halo.total_bytes_per_rank < rep1000.total_bytes_per_rank);
      assert(slab1000_halo.total_bytes_per_node < rep1000.total_bytes_per_node);
      assert(rep1000.slab_savings_factor < 1.0);
      assert(rep1000.slab_savings_factor > 0.0);
    }
  }

  // Delivery mode accounts for VBF (2 arrays) + delivery (4 arrays) + scratch.
  const auto deliv1000 = ChemicalField::estimate_memory_footprint(
      1000, 1000, 50, kShippedFullSpecies, 2, 2, false, true);
  constexpr uint64_t kExpectedVbf =
      2ULL * static_cast<uint64_t>(kShippedFullSpecies) * kCells1000
      * sizeof(Real);
  constexpr uint64_t kExpectedDelivery =
      4ULL * static_cast<uint64_t>(kShippedFullSpecies) * kCells1000
      * sizeof(Real);
  constexpr uint64_t kExpectedScratch =
      kCells1000 * (sizeof(char) + sizeof(Int));
  assert(deliv1000.vbf_bytes_per_rank == kExpectedVbf);
  assert(deliv1000.delivery_bytes_per_rank == kExpectedDelivery);
  assert(deliv1000.scratch_bytes_per_rank == kExpectedScratch);
  assert(deliv1000.total_bytes_per_rank
         == kExpectedPerRank1000 + kExpectedVbf + kExpectedDelivery
                + kExpectedScratch);
  assert(deliv1000.total_bytes_per_node == 2ULL * deliv1000.total_bytes_per_rank);

  std::cout << "  test_memory_footprint_exact_arithmetic: PASSED\n";
}

void test_proc_meminfo_parser() {
  const std::string sample =
      "MemTotal:       65800000 kB\n"
      "MemFree:        12000000 kB\n"
      "MemAvailable:   32768000 kB\n"
      "Buffers:          500000 kB\n";
  const auto parsed = ChemicalField::parse_proc_meminfo_available_bytes(sample);
  assert(parsed.has_value());
  assert(*parsed == 32768000ULL * 1024ULL);

  assert(!ChemicalField::parse_proc_meminfo_available_bytes(
              "MemTotal: 1000 kB\nMemFree: 500 kB\n").has_value());
  assert(!ChemicalField::parse_proc_meminfo_available_bytes(
              "MemAvailable: notanumber kB\n").has_value());
  assert(!ChemicalField::parse_proc_meminfo_available_bytes("").has_value());

  std::cout << "  test_proc_meminfo_parser: PASSED\n";
}

void test_guard_fail_fast_and_auto_slab_pure() {
  ChemicalFieldMemoryGuardInput base;
  base.nx = 1000;
  base.ny = 1000;
  base.nz = 50;
  base.num_species = kShippedFullSpecies;
  base.nprocs = 2;
  base.ranks_per_node = 2;
  base.slab_mode = false;
  base.delivery_any = false;
  base.halo_width = 0;
  base.required_slab_halo_width = 5;
  base.slab_feature_eligible = true;
  base.max_host_chemical_gb = 12.0;  // Below ~16.39 GiB two-rank replicated.

  // 1. Default "error" policy fails fast before allocation and recommends slab.
  {
    auto input = base;
    input.replicated_oom_policy = "error";
    const auto decision = ChemicalField::evaluate_memory_guard(input);
    assert(decision.has_budget);
    assert(decision.exceeds_budget);
    assert(decision.emit_replicated_diagnostic);
    assert(decision.should_throw);
    assert(!decision.switch_to_slab);
    assert(decision.diagnostic_message.find("\"slab\"") != std::string::npos);
  }

  // 2. "auto_slab" switches cleanly when slab fits inside budget (~8.36 GiB <= 12 GiB).
  {
    auto input = base;
    input.replicated_oom_policy = "auto_slab";
    const auto decision = ChemicalField::evaluate_memory_guard(input);
    assert(decision.exceeds_budget);
    assert(decision.switch_to_slab);
    assert(!decision.should_throw);
    assert(decision.diagnostic_message.find("auto_slab") != std::string::npos);
  }

  // 3. "auto_slab" refuses to switch and fails fast when a feature blocks slab.
  {
    auto input = base;
    input.replicated_oom_policy = "auto_slab";
    input.slab_feature_eligible = false;
    const auto decision = ChemicalField::evaluate_memory_guard(input);
    assert(decision.exceeds_budget);
    assert(!decision.switch_to_slab);
    assert(decision.should_throw);
  }

  // 4. "warn" policy emits diagnostic without throwing or switching.
  {
    auto input = base;
    input.replicated_oom_policy = "warn";
    const auto decision = ChemicalField::evaluate_memory_guard(input);
    assert(decision.exceeds_budget);
    assert(!decision.switch_to_slab);
    assert(!decision.should_throw);
    assert(!decision.diagnostic_message.empty());
  }

  // 5. Auto-budget from /proc/meminfo (90% of MemAvailable): 16 GiB RAM triggers
  //    at 14.4 GiB budget; 64 GiB RAM passes while emitting the replicated notice.
  {
    auto input = base;
    input.max_host_chemical_gb = 0.0;
    input.replicated_oom_policy = "error";
    input.available_host_memory_bytes =
        static_cast<uint64_t>(16.0 * kBytesPerGiB);
    const auto tight = ChemicalField::evaluate_memory_guard(input);
    assert(tight.exceeds_budget);
    assert(tight.should_throw);

    input.available_host_memory_bytes =
        static_cast<uint64_t>(64.0 * kBytesPerGiB);
    const auto roomy = ChemicalField::evaluate_memory_guard(input);
    assert(!roomy.exceeds_budget);
    assert(!roomy.should_throw);
    assert(roomy.emit_replicated_diagnostic);
  }

  std::cout << "  test_guard_fail_fast_and_auto_slab_pure: PASSED\n";
}

void test_simulation_pre_allocation_fail_fast() {
  SimulationConfig cfg = InputParser::default_config();
  cfg.hdf5.enabled = false;
  cfg.gpu.enabled = false;
  cfg.domain.lo = {0.0, 0.0, 0.0};
  cfg.domain.hi = {100.0e-6, 100.0e-6, 50.0e-6};
  cfg.domain.grid_dx = 5.0e-6;  // 20x20x10 = 4000 cells (~704 KiB for 11 species)
  cfg.domain.hash_cell_size = 10.0e-6;
  cfg.memory.max_host_chemical_gb = 1.0e-5;  // ~10 KiB -> must fail before alloc
  cfg.memory.replicated_oom_policy = "error";

  Simulation sim;
  bool caught_memory_limit = false;
  try {
    sim.init(cfg);
  } catch (const MemoryLimitError& ex) {
    caught_memory_limit = true;
    const std::string msg = ex.what();
    assert(msg.find("budget exceeded") != std::string::npos);
  }
  assert(caught_memory_limit);
  // Pre-allocation guarantee: ChemicalField never allocated any cells.
  assert(sim.chemical_field().ncells() == 0);

  // Confirm MemoryLimitError is also catchable as ConfigError.
  bool caught_as_config_error = false;
  try {
    Simulation sim2;
    sim2.init(cfg);
  } catch (const ConfigError&) {
    caught_as_config_error = true;
  }
  assert(caught_as_config_error);

  // Policy "warn" allows initialization to proceed on the small grid.
  cfg.memory.replicated_oom_policy = "warn";
  Simulation sim_warn;
  sim_warn.init(cfg);
  assert(sim_warn.chemical_field().ncells() == 4000);
  // Non-delivery run must NOT allocate delivery arrays on host.
  assert(!sim_warn.chemical_field().delivery_storage_allocated());

  std::cout << "  test_simulation_pre_allocation_fail_fast: PASSED\n";
}

void test_pinned_host_buffer_and_profile_fields() {
  PinnedHostBuffer<double> buf;
  assert(buf.empty());
  assert(buf.size() == 0);

  buf.ensure_capacity(64);
  assert(buf.size() >= 64);
  assert(buf.data() != nullptr);
  for (size_t i = 0; i < buf.size(); ++i) {
    buf.data()[i] = static_cast<double>(i) * 1.5;
  }

  // Requesting a smaller capacity preserves existing storage without realloc.
  const double* ptr_before = buf.data();
  const size_t size_before = buf.size();
  buf.ensure_capacity(16);
  assert(buf.size() == size_before);
  assert(buf.data() == ptr_before);
  assert(std::abs(buf.data()[10] - 15.0) <= 1.0e-12);

  // Move construction transfers ownership cleanly.
  PinnedHostBuffer<double> moved(std::move(buf));
  assert(moved.size() == size_before);
  assert(std::abs(moved.data()[10] - 15.0) <= 1.0e-12);
  assert(buf.empty());

  // Move assignment transfers ownership cleanly.
  PinnedHostBuffer<double> assigned;
  assigned = std::move(moved);
  assert(assigned.size() == size_before);
  assert(std::abs(assigned.data()[10] - 15.0) <= 1.0e-12);

  // GpuProfile exposes h2d_bytes, d2h_bytes, and pinned_staging_active.
  gpu_transfer_profile_reset();
  gpu_transfer_profile_set_enabled(true);
  gpu_transfer_set_pinned_staging_active(true);
  const GpuProfile prof = gpu_transfer_profile_snapshot();
  assert(prof.h2d_bytes == 0ULL);
  assert(prof.d2h_bytes == 0ULL);
  assert(prof.pinned_staging_active);
  assert(gpu_transfer_pinned_staging_active());
  gpu_transfer_set_pinned_staging_active(false);
  gpu_transfer_profile_set_enabled(false);

  std::cout << "  test_pinned_host_buffer_and_profile_fields: PASSED\n";
}

void test_memory_config_json_parsing() {
  SimulationConfig cfg = InputParser::default_config();
  const std::string doc =
      R"({"memory": {"max_host_chemical_gb": 14.5, "replicated_oom_policy": "auto_slab"}})";
  assert(ConfigJson::parse_document(cfg, doc));
  InputParser::finalize_config(cfg);
  assert(std::abs(cfg.memory.max_host_chemical_gb - 14.5) <= 1.0e-12);
  assert(cfg.memory.replicated_oom_policy == "auto_slab");

  const std::string serialized = ConfigJson::serialize_document(cfg);
  SimulationConfig roundtrip = InputParser::default_config();
  assert(ConfigJson::parse_document(roundtrip, serialized));
  assert(std::abs(roundtrip.memory.max_host_chemical_gb - 14.5) <= 1.0e-12);
  assert(roundtrip.memory.replicated_oom_policy == "auto_slab");

  bool rejected_neg = false;
  try {
    SimulationConfig bad = InputParser::default_config();
    ConfigJson::parse_document(bad, R"({"memory.max_host_chemical_gb": -1.0})");
  } catch (const ConfigError&) {
    rejected_neg = true;
  }
  assert(rejected_neg);

  bool rejected_policy = false;
  try {
    SimulationConfig bad = InputParser::default_config();
    ConfigJson::parse_document(
        bad, R"({"memory.replicated_oom_policy": "silent_oom"})");
  } catch (const ConfigError&) {
    rejected_policy = true;
  }
  assert(rejected_policy);

  std::cout << "  test_memory_config_json_parsing: PASSED\n";
}

}  // namespace

int main() {
  std::cout << "=== ChemicalField Memory Guard & Pinned Staging Tests ===\n";
  test_memory_footprint_exact_arithmetic();
  test_proc_meminfo_parser();
  test_guard_fail_fast_and_auto_slab_pure();
  test_simulation_pre_allocation_fail_fast();
  test_pinned_host_buffer_and_profile_fields();
  test_memory_config_json_parsing();
  std::cout << "All ChemicalField memory guard tests passed.\n";
  return 0;
}
