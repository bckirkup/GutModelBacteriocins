/* -----------------------------------------------------------------------
   Shared named-shard dispatch for multi-sub-test binaries registered via
   gutibm_add_test(... SHARDS ...). Each binary defines a kTests table;
   run with shard-name arguments it executes just those sub-tests, and
   with no arguments it runs all of them.
   ----------------------------------------------------------------------- */

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <string>

namespace testshards {

struct NamedTest {
  const char* name;
  void (*run)();
};

template <std::size_t N>
int shard_main(int argc, char** argv,
               const std::array<NamedTest, N>& tests,
               const char* banner, const char* pass_message) {
  std::cout << banner;
  if (argc > 1) {
    for (int i = 1; i < argc; ++i) {
      const std::string shard = argv[i];
      const auto it = std::find_if(
          tests.begin(), tests.end(),
          [&shard](const NamedTest& test) { return shard == test.name; });
      if (it == tests.end()) {
        std::cerr << "unknown test shard: " << shard << "\n";
        return 2;
      }
      it->run();
    }
    return 0;
  }
  for (const NamedTest& test : tests) {
    test.run();
  }
  std::cout << pass_message;
  return 0;
}

}  // namespace testshards
