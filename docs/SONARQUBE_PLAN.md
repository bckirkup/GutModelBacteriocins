# SonarQube Remediation Plan — GutIBM

**Project:** [bckirkup_GutModelBacteriocins](https://sonarcloud.io/project/overview?id=bckirkup_GutModelBacteriocins)

**Status (Jun 2026):** Phases 1–4 cleared ~200 findings (bugs, vulns, mechanical, structural).
**Remaining after this PR:** target **0 open issues** — 2 security fixes + 83 smells suppressed.

## Policy

| Category | Action |
|----------|--------|
| **BUG** | Fix immediately — blocks merge |
| **VULNERABILITY** | Fix immediately — blocks merge |
| **CODE_SMELL** | Accepted technical debt — suppressed via `sonar-project.properties` |
| **New smells on changed lines** | Fix opportunistically when touching the file |

The quality gate should pass on **reliability + security** only. Maintainability smells in
diffusion kernels, GPU headers, and legacy Fix modules are not worth batch-refactoring for
a research prototype.

## Security (fixed — do not regress)

### `pythonsecurity:S8707` — path traversal on CLI output

**Problem:** Sonar's taint engine does not recognize custom `path_utils` validators, so
`prepare_output_file` / `write_text_file` kept flagging `mkdir` and `open` sinks even after
validation.

**Fix (two parts):**

1. **Code** (`python/gut_ibm_tools/path_utils.py`):
   - Run full parent/output validation **before** any `mkdir`
   - Use `open()` with a re-validated path instead of `Path.write_text`
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

## Suppressed code smells (83 issues, 38 rule families)

Configured in `sonar-project.properties` via `sonar.issue.ignore.multicriteria`.
Re-scan after merge should close these without code changes.

| Rule family | Count | Why suppressed |
|-------------|-------|----------------|
| `cpp:S3608` | 11 | FMM kernel loop style — refactor risk in hot path |
| `cpp:S107` | 7 | Long parameter lists in diffusion APIs — needs struct grouping per subsystem |
| `cpp:S134` | 5 | Deep nesting in FMM/QSSA — partial Phase 3 work remains |
| `cpp:S5817` / `S6177` / `S1820` | 12 | Const-correctness / rule-of-five in core headers |
| `cpp:S5213` / `S5812` | 6 | GPU header include guards / export macros |
| `cpp:S6004` / `S995` / `S6022` | 7 | C++17 init-if / cast / emplace — low value |
| `cpp:S6045` / `S6231` / `S6495` | 7 | Template deduction / structured bindings — compile-risky batch |
| `cpp:S1121` | 2 | Assignment in sub-expression in FMM |
| `cpp:S1144` | 2 | Unused private methods — intentional stubs |
| `cpp:S6185` | 2 | `std::make_unique` vs `new` in lineage tracker |
| `cpp:S7034` | 1 | `string_view::contains` — C++23, toolchain is C++17 |
| Single-count stragglers | 17 | One-off style in tests/GPU — fix when file is touched |
| `python:S117` | 1 | PEP8 naming in `validation.py` CLI helper |

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

**Target:** 0 open BUG/VULNERABILITY at all times; 0 total open issues after re-scan.

## Related docs

- `.agents/skills/sonarqube-gutibm/SKILL.md` — agent workflow
- `.agents/skills/sonarqube-cpp/SKILL.md` — C++ patterns
- `.agents/skills/sonarqube-python/SKILL.md` — Python patterns
