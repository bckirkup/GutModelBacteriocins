---
name: sonarqube-cpp
description: SonarQube C++ coding patterns for GutIBM. Use when writing or fixing C++ in src/ or tests/ to avoid recurring Sonar findings (exceptions, paths, RNG, modernization, complexity).
---

# SonarQube C++ Patterns — GutIBM

Rule-keyed reference for `src/` and `tests/`. Project: `bckirkup_GutModelBacteriocins`.

See `.agents/skills/sonarqube-gutibm/SKILL.md` for batch PR workflow and priority order.

## Exceptions (S112, S1181)

Use domain exceptions from `src/core/error.h`:

```cpp
#include "error.h"

throw gutibm::ConfigError("missing key: " + key);
throw gutibm::HDF5Error("dataset not found");
throw gutibm::PathError("path traversal detected");
throw gutibm::SimulationError("invalid state");
```

Hierarchy: `Error` ← `ConfigError`, `IOError` ← `HDF5Error`/`PathError`, `CheckpointError`, `SimulationError`.

Catch specific types — never bare `catch (...)`:

```cpp
catch (const ConfigError& e) { ... }
catch (const IOError& e) { ... }
catch (const std::exception& e) { ... }  // last resort in parse helpers only
```

## Enum conversion (S7035)

Use `to_underlying()` from `src/core/types.h`:

```cpp
int idx = to_underlying(ReceptorType::BtuB);           // Good
int idx = static_cast<int>(ReceptorType::BtuB);        // Bad — S7035
```

Exception: compile-time constants like `NUM_RECEPTORS` may keep `static_cast<int>()`.

## Variable declarations (S1659)

One variable per statement:

```cpp
std::vector<double> x(n);
std::vector<double> y(n);
int lo = 0;
int hi = 0;
```

## Loops and ranges (S5566, S6197)

Prefer range-for over index loops when iterating containers:

```cpp
for (const auto& agent : agents) { ... }              // Good
for (const auto& part : fs::path(path)) { ... }       // Good

if (std::ranges::contains(path.parts(), "..")) { ... } // Good — parent traversal check
```

For erase-remove, prefer C++20 `std::erase` / `std::ranges::remove` where available.

## Path and temp file safety (S5443)

All user-supplied paths go through `src/io/path_utils.{h,cpp}`:

| Function | Use |
|----------|-----|
| `validate_path_syntax()` | Reject `..`, null bytes |
| `validate_input_file_path()` | Read existing files |
| `validate_output_file_path()` | Write to existing parent dir |
| `secure_temp_file()` | Private temp via `mkstemp` |
| `resolve_test_h5_path()` | Test HDF5 output with env override |

`temp_directory()` validates `TMPDIR` / system temp:

- Must exist and be a directory
- Rejects symlinked temp dirs
- World-writable dirs must have the sticky bit (`/tmp` is OK)

Never call `mkdir` / `open` on a user path before validation.

## Pseudorandom numbers (S2245)

| Context | Pattern |
|---------|---------|
| Simulation / Fixes | `gutibm::RNG` with explicit seed (`src/core/random.h`) |
| Unit tests | `RNG rng(42);` + `rng.uniform(lo, hi)` — not `std::mt19937` |
| Crypto / secrets | **Not applicable** — do not use `mt19937` for security |

`random.h` documents simulation-only PRNG:

```cpp
std::mt19937_64 gen_;  // NOSONAR cpp:S2245 — reproducible IbM stochasticity via explicit seed
```

## C++17/20 modernization

### Init-statement in if (S6004)

```cpp
if (auto it = map.find(key); it != map.end()) {   // Good
  use(it->second);
}
```

Skip when the variable must be `static` (init-statement cannot replace function-scope static).

### string_view for read-only strings (S6009)

```cpp
void parse_key(std::string_view key);   // Good — no copy
void parse_key(const std::string& key); // Bad for read-only params
```

### Emplacement (S6003, S6011)

```cpp
vec.emplace_back(x, y, z);                    // Good
vec.push_back(Vec3{x, y, z});                 // S6011 if temporary
```

### Template argument deduction (S6012)

```cpp
std::vector<int> v;                             // Good
std::vector<int> v(10);                         // Bad if deducible — S6012
```

**Caution:** S6012 false-positives when the template arg cannot be deduced (e.g. `MyClass<int>(10)` when default ctor arg is required). Compile after each change.

## Complexity and nesting (S134, S3776, S107)

### Deep nesting (S134) — CRITICAL severity smell

Extract helpers, use early return, invert guards:

```cpp
// Bad — 4 levels
for (...) {
  if (a) {
    for (...) {
      if (b) { do_work(); }
    }
  }
}

// Good
for (...) {
  if (!a) continue;
  process_outer(...);
}
```

Hotspots: `src/diffusion/fmm_kernel.cpp`, `tests/test_config_diversity.cpp`.

### Cognitive complexity (S3776)

Split functions above 15 complexity into named helpers. Prefer one responsibility per function.

### Too many parameters (S107)

Group related parameters into a struct (e.g. `FMMBuildContext`).

## General hygiene

| Rule | Pattern |
|------|---------|
| S4962 | `nullptr`, not `NULL` or `0` for pointers |
| S5350 | Pass non-trivial objects by `const&` when read-only |
| S5272 | Don't read values after `std::move` |
| S108 | Empty blocks need a comment explaining why |
| S1481 | Remove unused vars; `(void)param;` for intentional unused |
| S5827 | Mark noexcept move operations where appropriate |
| S6007 | `[[nodiscard]]` on functions whose return must not be ignored |

## Do not batch blindly

| Rule | Risk |
|------|------|
| S6012 | May break compilation when deduction fails |
| S6004 | Invalid on function-scope `static` locals |
| S134 / S3776 | Behavior change risk — refactor one file, run full tests |
| S6197 | Requires C++20 ranges support in toolchain |

## Reference implementations

- Paths: `src/io/path_utils.cpp`
- RNG: `src/core/random.h`, `tests/test_fmm.cpp` (post-#102 pattern)
- Exceptions: `src/core/error.h`
- Enum helper: `src/core/types.h` (`to_underlying`)
