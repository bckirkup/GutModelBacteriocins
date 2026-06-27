/* -----------------------------------------------------------------------
   GutIBM – Path validation unit tests
   ----------------------------------------------------------------------- */

#include "path_utils.h"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>

using namespace gutibm;

namespace {

void expect_throw(void (*fn)(), const char* label) {
  bool threw = false;
  try {
    fn();
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
  std::cout << "  " << label << ": PASSED\n";
}

void test_rejects_parent_traversal() {
  expect_throw([] { validate_path_syntax("../etc/passwd"); },
               "test_rejects_parent_traversal");
}

void test_rejects_null_byte() {
  expect_throw([] { validate_path_syntax(std::string("safe\0path", 9)); },
               "test_rejects_null_byte");
}

void test_secure_temp_file_is_private() {
  const std::string path = secure_temp_file("gutibm_test_");
  struct stat st {};
  assert(stat(path.c_str(), &st) == 0);
  assert(S_ISREG(st.st_mode));
  assert((st.st_mode & S_IRWXG) == 0);
  assert((st.st_mode & S_IRWXO) == 0);
  std::remove(path.c_str());
  std::cout << "  test_secure_temp_file_is_private: PASSED\n";
}

void test_validate_input_file() {
  const std::string path = secure_temp_file("gutibm_input_");
  {
    std::ofstream out(path);
    out << "ok";
  }

  const std::string canonical = validate_input_file_path(path);
  assert(!canonical.empty());
  std::remove(path.c_str());
  std::cout << "  test_validate_input_file: PASSED\n";
}

void test_resolve_test_h5_path_env_override() {
  const std::string path = secure_temp_file("gutibm_env_");
  setenv("GUTIBM_TEST_H5", path.c_str(), 1);
  const std::string resolved = resolve_test_h5_path("GUTIBM_TEST_H5", "unit");
  assert(resolved == path);
  unsetenv("GUTIBM_TEST_H5");
  std::remove(path.c_str());
  std::cout << "  test_resolve_test_h5_path_env_override: PASSED\n";
}

}  // namespace

int main() {
  std::cout << "=== Path Validation Tests ===\n";
  test_rejects_parent_traversal();
  test_rejects_null_byte();
  test_secure_temp_file_is_private();
  test_validate_input_file();
  test_resolve_test_h5_path_env_override();
  std::cout << "All path validation tests passed.\n";
  return 0;
}
