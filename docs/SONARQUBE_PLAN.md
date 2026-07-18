# SonarQube Remediation Plan — GutIBM

**Project:** [bckirkup_GutModelBacteriocins](https://sonarcloud.io/project/overview?id=bckirkup_GutModelBacteriocins)

**Status (Jul 2026):** Quality gate **OK** (0 BUG, 0 VULNERABILITY). Pragmatic
clear-to-zero:

| Batch | Action | Status |
|-------|--------|--------|
| **A** | Mechanical code fixes (~25 smells) | Done (`cursor/sonar-mechanical-cleanup-aead` / PR) |
| **B** | Won’t Fix accepted complexity/architecture debt (~54) | Script ready — run after merge re-scan |
| **C** | This doc + skill remaining-work map | Done |

## Policy

| Category | Action |
|----------|--------|
| **BUG** | Fix immediately — blocks merge |
| **VULNERABILITY** | Fix immediately — blocks merge |
| **New smells on changed lines** | Fix opportunistically when touching the file |
| **Accepted debt (Batch B rules)** | Multicriteria for scanner runs + SonarCloud Won’t Fix for auto-analysis dashboard |

The quality gate should pass on **reliability + security** (and current
maintainability-on-new-code conditions). Maintainability smells in diffusion
kernels, GPU headers, and legacy Fix modules are not worth batch-refactoring for
a research prototype.

## Why multicriteria alone is not enough

[`sonar-project.properties`](../sonar-project.properties) ignores accepted smell
families via `sonar.issue.ignore.multicriteria`. **SonarCloud automatic analysis
does not apply those exclusions**, so ignored rules can still appear as open
issues on the dashboard. Clearing the dashboard requires either a code fix or a
SonarCloud **Won’t Fix** (or False Positive) resolution.

## Batch A — Mechanical fixes (done)

Cleared in code (do not re-suppress these rules project-wide):

| Rule | Fix |
|------|-----|
| `python:S1192` | Path traversal message constant |
| `cpp:S6009` | `std::string_view` for flat keys / compression helpers |
| `cpp:S5566` | `std::ranges::any_of` in QSSA GPU parity test |
| `cpp:S3358` | Nested ternary split in `simulation.cpp` |
| `cpp:S1854` | Dead stores to `stopped_for_population` |
| `cpp:S1905` | Redundant cast removed |
| `cpp:S125` | Comment reworded (not “commented-out code”) |
| `cpp:S5827` | Redundant type → `auto` |
| `cpp:S5421` | Mutable global replaced with function-local counter |
| `cpp:S1188` | Long test lambdas → named helpers |
| `cpp:S5812` | Flattened `namespace gutibm::test` |
| `cpp:S6177` | `using enum` |
| `cpp:S5945` | C arrays → `std::array` for HDF5 dims |
| `cpp:S6022` | `std::byte` in agent-transfer append |

## Batch B — Accepted debt (Won’t Fix)

| Rule family | Approx. count | Reason |
|-------------|---------------|--------|
| `cpp:S134` nesting | 11 | Hot kernels / receptor / GPU |
| `cpp:S107` param count | 10 | Diffusion APIs need context-struct redesign |
| `cpp:S6004` init-if | 10 | Low-value modernization |
| `cpp:S3776` / `python:S3776` | 11 | Parser, HDF5, GPU, batch CLI |
| `cpp:S995` const ptr | 4 | GPU buffer mutability |
| `cpp:S5008` `void*` | 2 | HDF5 C API buffers |
| `cpp:S1820` / `cpp:S1448` | 3 | `Simulation` / GPU types size |
| `cpp:S3656` protected | 1 | NUFEB-style `Fix` base |
| `cpp:S924` nested break | 1 | Coupled to `Simulation::run` |
| `cpp:S7034` `contains` | 1 | `string_view::contains` is C++23; project is C++20 |

### Clear dashboard after Batch A merges

1. Wait for SonarCloud automatic analysis on `main` (or the PR).
2. With a token that can **Administer Issues**:

```bash
export SONAR_TOKEN=...   # SonarCloud user token
python3 scripts/sonar_wont_fix_debt.py --dry-run   # preview
python3 scripts/sonar_wont_fix_debt.py             # resolve
```

3. Or in the UI: Issues → filter by the rules above → Bulk Change → Won’t Fix,
   comment pointing at this doc.
4. Keep multicriteria entries in `sonar-project.properties` for CI/scanner runs.
   Do **not** add new exclusions for rules fixed in Batch A.

## Security (fixed — do not regress)

### `pythonsecurity:S8707` — path traversal on CLI output

**Problem:** Sonar's taint engine does not recognize custom `path_utils` validators, so
`prepare_output_file` / `write_text_file` kept flagging `mkdir` and `open` sinks even after
validation.

**Fix (two parts):**

1. **Code** (`python/gut_ibm_tools/path_utils.py`):
   - Run full parent/output validation **before** any `mkdir`
   - `_trusted_output_path()` rebuilds cwd-relative paths from regex-sanitized
     segments (breaks S8707 taint from CLI args on SonarCloud automatic analysis)
   - Use `open()` on the trusted path instead of passing user input through
   - `_mkdir_validated_parents()` creates only pre-checked directory segments

2. **SAST config** (`sonar/pythonsecurity-s8707.json` + `sonar-project.properties`):
   - Registers `validate_path_syntax`, `validate_output_path`, `prepare_output_file`, and
     `_ensure_output_within_cwd` as validators for rule S8707
   - Upload the same JSON in SonarCloud → Project Settings → SAST Engine → Python custom
     configuration if the scanner property is not picked up automatically

**All user-supplied paths** must go through `path_utils` — see
`.agents/skills/sonarqube-python/SKILL.md`.

### C++ path I/O (`cpp:S5443`, `cpp:S2083`)

Already handled in `src/io/path_utils.cpp` with `validate_temp_directory()` + `mkstemp`.
Approved `NOSONAR` only where env temp dirs are validated before use.

## Suppressed code smells (scanner multicriteria)

Configured in `sonar-project.properties` via `sonar.issue.ignore.multicriteria`
(accepted debt families, including Batch B rules). Re-scan with the scanner
honors these; automatic analysis still needs Won’t Fix (script above).

## If a suppressed rule fires on new code

1. Fix it if the change is small (same PR).
2. If the rule is wrong for the pattern, add a documented `NOSONAR` (C++) — never for path I/O.
3. Do **not** add new multicriteria exclusions without updating this doc.

## Monitoring

```bash
# Open issue count
curl -s "https://sonarcloud.io/api/issues/search?\
componentKeys=bckirkup_GutModelBacteriocins&resolved=false&ps=1" \
  | python3 -c "import json,sys; print('open:', json.load(sys.stdin)['total'])"

# By type
curl -s "https://sonarcloud.io/api/issues/search?\
componentKeys=bckirkup_GutModelBacteriocins&resolved=false&facets=types&ps=1" \
  | python3 -m json.tool
```

**Target:** 0 open BUG/VULNERABILITY at all times; 0 total open issues after
Batch A merge + Batch B Won’t Fix.

## Related docs

- `.agents/skills/sonarqube-gutibm/SKILL.md` — agent workflow
- `.agents/skills/sonarqube-cpp/SKILL.md` — C++ patterns
- `.agents/skills/sonarqube-python/SKILL.md` — Python patterns
- `scripts/sonar_wont_fix_debt.py` — Batch B dashboard clear
