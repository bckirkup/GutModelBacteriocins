---
name: testing-gutibm
description: Test the GutIBM C++/Python simulation end-to-end. Covers build verification, CTest execution, Python pytest, and known preexisting failures. Use when verifying code changes or running post-merge validation.
---

# Testing GutIBM

Shell-based testing only (no GUI/browser). Do not record — collect command outputs as text evidence.

## Quick Verification Sequence

```bash
# 1. Build with warnings (zero new warnings required)
cd build && cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON \
  -DCMAKE_CXX_FLAGS="-Wall -Wextra" && make -j$(nproc) 2>&1

# 2. Run C++ tests
cd build && ctest --output-on-failure

# 3. Python lint + tests
cd python && ruff check .
cd python && pytest tests/ -v -m "not integration"
```

## Preexisting Test Failures

The Devin VM has a broken `libhwloc` that crashes MPI initialization. These tests fail on ALL branches — do not attribute them to code changes:

| Test | Failure Mode | Root Cause |
|------|-------------|------------|
| `scaling_benchmark` | MPI init crash (SIGFPE) | `libhwloc` integer divide-by-zero in `hwloc_topology_load` |
| `hdf5_roundtrip` | MPI init crash | Same libhwloc issue |
| `hdf5_checkpoint` | MPI init crash | Same libhwloc issue |
| `hdf5_roundtrip_parallel` | Numerical exception | MPI init failure cascade |
| `mpi_multi_rank` | Numerical exception | MPI init failure cascade |
| `config_diversity` | `fp_subset != fp_full` assertion | Preexisting test logic issue (not MPI-related) |

**Expected baseline:** 26/32 CTest pass, 6 fail. If you see fewer than 26 passing, investigate — a code change may have introduced a regression.

To confirm a failure is preexisting, checkout the base branch, rebuild, and run the specific test:

```bash
git stash && git checkout <base_commit>
cd build && cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON && make -j$(nproc)
ctest -R <test_name> --output-on-failure
git checkout - && git stash pop
```

## Checking for New Warnings

After building with `-Wall -Wextra`, filter warnings to only your modified files:

```bash
make -j$(nproc) 2>&1 | grep 'warning:' | grep -E '(file1|file2|file3)' | wc -l
```

Result must be 0. There are ~150 preexisting warnings in base class virtual functions (`vbf.h`, `fix.h`) and test files — these are expected.

## Python Test Environment

```bash
cd python && pip install -e ".[dev]"  # installs ruff, pytest, etc.
```

Python tests: 46 non-integration tests + 1 deselected (integration marker). All should pass.

## SonarQube-sensitive test patterns

When writing or reviewing Python tests, follow `.agents/skills/sonarqube-python/SKILL.md`:

- Float comparisons: `pytest.approx()` or `np.testing.assert_allclose` — never `== 0.0`
- Paths in tests: use `tmp_path` fixtures; production path helpers live in `gut_ibm_tools.path_utils`

When writing C++ tests with random fixtures, use `gutibm::RNG` (see `.agents/skills/sonarqube-cpp/SKILL.md`), not `std::mt19937`.

## Golden File Considerations

When modifying RNG usage (e.g., switching from `np.random.seed()` to `np.random.default_rng()`), the deterministic stream changes. This requires updating golden files:

- `python/tests/fixtures/sample_hdf5_golden.json` — local test golden values
- `python/tests/fixtures/eari_vadi_ci_golden.json` — CI regression golden values (value comes from CI run, not local)

The EARI/VADI golden must match what CI computes. Push the code change first, let CI fail, read the actual value from the CI log, then update the golden file.

## Adversarial Tests for Exception/Type Changes

When modifying exception types or enum conversions, verify at runtime, not just at compile time:

- **Exception hierarchy:** Run `test_input_parser`, `test_smoke`, `test_agent` — these exercise throw/catch paths
- **`to_underlying()` polyfill:** Run `test_receptor`, `test_iron_fallback` — these use enum-to-int in biological computations where wrong values produce wrong outcomes
- **Python RNG:** Verify determinism with same seed and divergence with different seeds

## Devin Secrets Needed

None — all testing is local. No external services or credentials required.
