/* -----------------------------------------------------------------------
   GutIBM – Spatial hash tests
   ----------------------------------------------------------------------- */

#include "spatial_hash.h"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace gutibm;

void test_insert_and_query() {
  SpatialHash hash;
  hash.init({0, 0, 0}, {1e-3, 1e-3, 1e-4}, 10e-6);

  hash.insert(0, {5e-6, 5e-6, 5e-6});
  hash.insert(1, {8e-6, 5e-6, 5e-6});
  hash.insert(2, {500e-6, 500e-6, 50e-6});

  // Query near agent 0 — should find 0 and 1
  auto result = hash.query_radius({5e-6, 5e-6, 5e-6}, 15e-6);
  assert(result.size() >= 2);

  // Query near agent 2 — should find only 2
  auto result2 = hash.query_neighbors({500e-6, 500e-6, 50e-6});
  bool found2 = false;
  for (auto idx : result2) {
    if (idx == 2) found2 = true;
  }
  assert(found2);

  std::cout << "  test_insert_and_query: PASSED\n";
}

void test_empty_query() {
  SpatialHash hash;
  hash.init({0, 0, 0}, {1e-3, 1e-3, 1e-4}, 10e-6);

  auto result = hash.query_radius({500e-6, 500e-6, 50e-6}, 5e-6);
  assert(result.empty());

  std::cout << "  test_empty_query: PASSED\n";
}

void test_clear() {
  SpatialHash hash;
  hash.init({0, 0, 0}, {1e-3, 1e-3, 1e-4}, 10e-6);

  hash.insert(0, {5e-6, 5e-6, 5e-6});
  hash.clear();

  auto result = hash.query_neighbors({5e-6, 5e-6, 5e-6});
  assert(result.empty());

  std::cout << "  test_clear: PASSED\n";
}

int main() {
  std::cout << "=== Spatial Hash Tests ===\n";
  test_insert_and_query();
  test_empty_query();
  test_clear();
  std::cout << "All spatial hash tests passed.\n";
  return 0;
}
