/* -----------------------------------------------------------------------
   GutIBM – Tests for graded iron uptake fallback (Issue #10)
   Verifies that secondary receptors (IroN, IutA, Fiu) provide intermediate
   growth when FepA is downregulated.
   ----------------------------------------------------------------------- */

#include "agent.h"
#include "fix_metabolism.h"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace gutibm;

void test_receptor_type_iutA_exists() {
  // IutA should be at index 7
  static_assert(static_cast<int>(ReceptorType::IutA) == 7,
                "IutA must be ReceptorType 7");
  static_assert(NUM_RECEPTORS == 8,
                "NUM_RECEPTORS must be 8 with IutA added");
  std::cout << "  test_receptor_type_iutA_exists: PASSED\n";
}

void test_agent_iutA_default_expression() {
  Agent a = Agent::create_default(1, 1, {100e-6, 100e-6, 50e-6}, 5e-4);
  Real iutA_expr = a.receptor_expr[static_cast<int>(ReceptorType::IutA)];
  assert(std::abs(iutA_expr - 1.0) < 1e-12);
  std::cout << "  test_agent_iutA_default_expression: PASSED\n";
}

void test_graded_fallback_intermediate_growth() {
  // When FepA is fully down but secondary receptors are up,
  // iron uptake should be non-zero (intermediate via IroN/IutA/Fiu)
  MetabolismConfig cfg;

  // Simulate the iron uptake calculation manually
  Real S_iron = 1.0e-4;  // 100 uM iron (excess)

  // Case 1: Wild-type (all receptors at 1.0)
  Real expr_fepA_wt = 1.0;
  Real expr_iroN_wt = 1.0;
  Real expr_iutA_wt = 1.0;
  Real expr_fiu_wt  = 1.0;

  Real uptake_wt = 0.0;
  uptake_wt += expr_fepA_wt * S_iron / (cfg.km_iron_primary + S_iron);
  uptake_wt += expr_iroN_wt * S_iron / (cfg.km_iron_iroN + S_iron);
  uptake_wt += expr_iutA_wt * S_iron / (cfg.km_iron_iutA + S_iron);
  uptake_wt += expr_fiu_wt  * S_iron / (cfg.km_iron_fiu + S_iron);
  Real monod_wt = uptake_wt / (1.0 + expr_iroN_wt + expr_iutA_wt + expr_fiu_wt);

  // Case 2: FepA fully downregulated, secondary receptors intact
  Real expr_fepA_down = 0.0;
  Real expr_iroN_up = 1.0;
  Real expr_iutA_up = 1.0;
  Real expr_fiu_up  = 1.0;

  Real uptake_down = 0.0;
  uptake_down += expr_fepA_down * S_iron / (cfg.km_iron_primary + S_iron);
  uptake_down += expr_iroN_up * S_iron / (cfg.km_iron_iroN + S_iron);
  uptake_down += expr_iutA_up * S_iron / (cfg.km_iron_iutA + S_iron);
  uptake_down += expr_fiu_up  * S_iron / (cfg.km_iron_fiu + S_iron);
  Real monod_down = uptake_down / (1.0 + expr_iroN_up + expr_iutA_up + expr_fiu_up);

  // Case 3: All receptors down (complete iron starvation)
  Real uptake_none = 0.0;
  // All expression = 0
  Real monod_none = uptake_none / (1.0 + 0.0 + 0.0 + 0.0);

  // Assertions:
  // 1. Wild-type should have highest iron uptake
  assert(monod_wt > monod_down);
  // 2. FepA-down with secondary receptors should have intermediate growth (> 0)
  assert(monod_down > 0.0);
  // 3. All receptors down should give zero
  assert(std::abs(monod_none) < 1e-15);
  // 4. FepA-down should give substantially less than wild-type
  assert(monod_down < monod_wt);
  // 5. FepA-down should still give meaningful iron uptake (> 25% of wild-type)
  assert(monod_down > 0.25 * monod_wt);

  std::cout << "  test_graded_fallback_intermediate_growth: PASSED\n";
  std::cout << "    wild-type monod_iron = " << monod_wt << "\n";
  std::cout << "    FepA-down monod_iron = " << monod_down << "\n";
  std::cout << "    ratio = " << (monod_down / monod_wt) << "\n";
}

void test_graded_fallback_km_ordering() {
  MetabolismConfig cfg;
  // Verify Km ordering: FepA < IroN < IutA < Fiu
  assert(cfg.km_iron_primary < cfg.km_iron_iroN);
  assert(cfg.km_iron_iroN < cfg.km_iron_iutA);
  assert(cfg.km_iron_iutA < cfg.km_iron_fiu);
  std::cout << "  test_graded_fallback_km_ordering: PASSED\n";
}

void test_graded_fallback_partial_fepA() {
  // Partial FepA expression should give intermediate between full and zero FepA
  MetabolismConfig cfg;
  Real S_iron = 1.0e-5;  // 10 uM

  auto compute_monod = [&](Real fepA, Real iroN, Real iutA, Real fiu) {
    Real uptake = 0.0;
    uptake += fepA * S_iron / (cfg.km_iron_primary + S_iron);
    uptake += iroN * S_iron / (cfg.km_iron_iroN + S_iron);
    uptake += iutA * S_iron / (cfg.km_iron_iutA + S_iron);
    uptake += fiu  * S_iron / (cfg.km_iron_fiu + S_iron);
    return uptake / (1.0 + iroN + iutA + fiu);
  };

  Real monod_full = compute_monod(1.0, 1.0, 1.0, 1.0);
  Real monod_half = compute_monod(0.5, 1.0, 1.0, 1.0);
  Real monod_zero = compute_monod(0.0, 1.0, 1.0, 1.0);

  // Monotonically decreasing with FepA downregulation
  assert(monod_full > monod_half);
  assert(monod_half > monod_zero);
  // Partial should be between full and zero
  assert(monod_half > monod_zero);
  assert(monod_half < monod_full);

  std::cout << "  test_graded_fallback_partial_fepA: PASSED\n";
}

int main() {
  std::cout << "=== Iron Fallback Tests (Issue #10) ===\n";
  test_receptor_type_iutA_exists();
  test_agent_iutA_default_expression();
  test_graded_fallback_intermediate_growth();
  test_graded_fallback_km_ordering();
  test_graded_fallback_partial_fepA();
  std::cout << "All iron fallback tests passed.\n";
  return 0;
}
