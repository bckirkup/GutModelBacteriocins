---
name: sonarqube-gutibm
description: Plan and execute SonarQube remediation for GutIBM (bckirkup_GutModelBacteriocins). Use when fixing Sonar issues, batching cleanup PRs, or preventing new findings in C++/Python/shell code.
---

# SonarQube Remediation — GutIBM

Orchestration playbook for the SonarCloud project **`bckirkup_GutModelBacteriocins`**
([dashboard](https://sonarcloud.io/project/overview?id=bckirkup_GutModelBacteriocins)).

## When to use this skill

- Triaging or fixing open SonarQube issues
- Planning batch remediation PRs
- Reviewing PRs that touch I/O, tests, or randomness
- Writing new code that should not reintroduce cleared rule families

## Related skills

| Skill | Scope |
|-------|-------|
| `.agents/skills/sonarqube-cpp/SKILL.md` | C++ rule patterns and examples |
| `.agents/skills/sonarqube-python/SKILL.md` | Python rule patterns and examples |
| `.agents/skills/gut-ibm/SKILL.md` | Build, test, Fix architecture |
| `.agents/skills/testing-gutibm/SKILL.md` | Validation commands after fixes |

## Query open issues

SonarCloud API (no auth needed for public project metrics):

```bash
# Counts by type and severity
curl -s "https://sonarcloud.io/api/issues/search?componentKeys=bckirkup_GutModelBacteriocins&resolved=false&facets=types,severities,rules&ps=1" | python3 -m json.tool

# Issues for a specific rule
curl -s "https://sonarcloud.io/api/issues/search?componentKeys=bckirkup_GutModelBacteriocins&resolved=false&rules=cpp:S5566&ps=50"
```

With `sonarqube-cli` installed and authenticated:

```bash
sonar list issues -p bckirkup_GutModelBacteriocins --format table
sonar list issues -p bckirkup_GutModelBacteriocins --types CODE_SMELL --severities CRITICAL
```

## Fix priority order

Work in this order — same sequence used for PRs #101–#102:

1. **BUG** — reliability defects (e.g. `python:S1244` float equality)
2. **VULNERABILITY** — security (path traversal, unsafe temp dirs, weak PRNG in sensitive paths)
3. **CODE_SMELL — CRITICAL severity** — high-impact maintainability (`S134` nesting, `S3776` complexity)
4. **Mechanical smells** — safe batch transforms (`S5566`, `S1659`, `S6012`, `S1192`)
5. **Modernization** — C++17/20 (`S6009`, `S6004`, `S6197`)
6. **Judgment / suppression** — documented `NOSONAR` only where fixing is wrong

> **Severity ≠ type.** A "critical code smell" is still a maintainability finding, not a security bug.

## Batch PR workflow

Each remediation PR should:

1. Target **one rule family** or **one subsystem** (not random files).
2. Branch: `cursor/sonar-<phase>-<topic>-75ac`
3. Run validation before push:

```bash
cd build && cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON \
  -DCMAKE_CXX_FLAGS="-Wall -Wextra" && make -j$(nproc)
cd build && ctest -L unit -LE slow --output-on-failure
cd python && ruff check . && pytest tests/ -v -m "not integration"
```

4. For diffusion/simulation refactors, also run integration targets:

```bash
cd build && ctest -R "smoke|config_diversity" --output-on-failure
```

5. Commit message cites rule keys and approximate issue count cleared.
6. Re-scan on SonarCloud after merge before starting the next batch.

## Remaining work map (post #101 bugs, #102 vulnerabilities)

~280 code smells remain. Top rule families (batch in this order):

| Phase | Rules | ~Count | Risk |
|-------|-------|--------|------|
| 1 Mechanical | `cpp:S5566`, `cpp:S1659`, `cpp:S6012`, `cpp:S6011`, `python:S1192`, `shelldre:S7679` | ~120 | Low |
| 2 Modernize | `cpp:S6009`, `cpp:S6004`, `cpp:S6197`, `cpp:S6003` | ~50 | Low–medium |
| 3 Structural | `cpp:S134`, `cpp:S3776`, `cpp:S107` | ~45 | High — one subsystem per PR |
| 4 Stragglers | single-count rules, justified `NOSONAR` | ~20 | Case-by-case |

## NOSONAR policy

Use `// NOSONAR cpp:RULE — <reason>` sparingly. Approved cases:

| Rule | When to suppress |
|------|------------------|
| `cpp:S2245` | `gutibm::RNG` / seed-pinned test randomness — simulation only, not crypto |
| `cpp:S7035` | `to_underlying()` polyfill in `types.h` (intentional C++17 shim) |
| Other | **Never** for path I/O, exception handling, or float equality in tests |

Do **not** modify tests to make Sonar pass — fix the implementation or use the correct assertion pattern.

## New-code checklist

Before committing C++/Python changes in areas Sonar scans:

```
□ C++ exceptions from src/core/error.h (not std::runtime_error)
□ to_underlying() for enum→int in runtime code
□ gutibm::RNG for simulation randomness (not raw std::mt19937 in src/)
□ path_utils / prepare_output_file for any user-supplied filesystem path
□ pytest.approx / np.isclose for float assertions in Python tests
□ One variable per declaration (S1659)
□ Range-for or std::ranges when iterating containers (S5566)
□ np.random.default_rng() for Python randomness (S6709)
```

## Issue budget for PRs

- **Block merge:** new BUG or VULNERABILITY on the PR
- **Prefer fix:** new CODE_SMELL on changed lines
- **Target:** quality gate GREEN after Phase 3 structural work
