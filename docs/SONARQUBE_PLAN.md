# SonarQube Remediation Plan — GutIBM

**Project:** [bckirkup_GutModelBacteriocins](https://sonarcloud.io/project/overview?id=bckirkup_GutModelBacteriocins)

**Status (Aug 2026):** 541 open issues (0 BUG, 0 VULNERABILITY after PR #370),
under a clear-to-zero sweep. The standard is that an open finding must be
either fixed in code or resolved with a recorded, defensible reason — nothing
is left sitting on the dashboard as "probably fine".

| Batch | Action | Status |
|-------|--------|--------|
| **A** | Mechanical code fixes (~25 smells) | Done (Jul 2026) |
| **B** | Won’t Fix accepted complexity/architecture debt | Reclassified — see below |
| **C** | This doc + skill remaining-work map | Done |
| **D** | `pythonsecurity:S6549` manifest path taint (8) | Done in code (PR #370) |
| **E** | Concurrency triage of `cpp:S8379` (13) | 1 real race fixed (PR #371); 2 under audit; rest accepted |
| **F** | Type/template modernization sweep (~137) | In progress |
| **G** | Algorithms, control flow, and the `cpp:S1669` BLOCKER | Queued |
| **H** | `cpp:S6004` init-if (56) | Queued — reclassified from accepted debt to code fix |
| **I** | Python / Docker / shell findings (9) | Queued |

## Policy

| Category | Action |
|----------|--------|
| **BUG** | Fix immediately — blocks merge |
| **VULNERABILITY** | Fix immediately — blocks merge. Never Won’t Fix, never suppress, even when the finding is unexploitable — a reader who has to reason about whether it is safe has already paid the cost |
| **New smells on changed lines** | Fix opportunistically when touching the file |
| **Mechanical / modernization smells** | Fix in code, provided the fix is available in C++20 on the CI toolchain and changes no behaviour |
| **Accepted debt** | Only the four categories below, each with its reason recorded per resolution by `scripts/sonar_wont_fix_debt.py` |

Not permitted as a way to clear a finding: `// NOSONAR` / `# NOSONAR`, new
`sonar.issue.ignore.multicriteria` entries, and Won’t Fix without a reason that
survives being read out loud.

### The four accepted-debt categories

1. **Toolchain** — the fix needs a C++23 library feature. The project is
   `CMAKE_CXX_STANDARD 20` and CI builds with GCC 11, whose libstdc++ has no
   `<format>` at all: `cpp:S7034` (`std::string::contains`), `cpp:S7035`
   (`std::to_underlying`), `cpp:S6185` and `cpp:S6484` (`std::format`).
   Revisit as one batch when the toolchain moves, not case by case.
2. **Numerical reproducibility** — `cpp:S6179` wants `std::lerp`/`std::midpoint`
   in the Robin correction-table interpolation and the metabolic-mode blend.
   Neither is bit-identical to the current arithmetic, and both sites are
   validated against Python oracles at ~1e-9 and regression-guarded. Trading
   reproducibility of the scientific output for a style rule is not a trade.
3. **Synchronization the rule cannot see** — `cpp:S8379` wants a mutex on
   `mutable` members that are in fact protected by OpenMP `atomic update`,
   per-thread slots, or serial-only mutation. A mutex in the QSSA/Green's
   function hot loop would serialize it for no correctness gain. This category
   is earned per finding, not per rule: triaging it turned up one genuine race
   (below), so the rule stays open on the dashboard until the last two findings
   are audited.
4. **Architecture of a research prototype** — parameter counts, nesting,
   cognitive complexity, and type size in the diffusion kernels, the NUFEB-style
   `Fix` base, and the config parser (`cpp:S107`, `S134`, `S3776`,
   `python:S3776`, `S1820`, `S1448`, `S995`, `S5008`, `S3656`, `S924`).
   Addressing these is a redesign that would put working scientific code at
   risk to move a maintainability rating.

`cpp:S6004` (init-if, 56 findings) was previously in this list as "low-value
modernization" and has been moved out: it is a mechanical, behaviour-preserving
C++17 fix, so it gets fixed.

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

Counts as of the Aug 2026 inventory. `scripts/sonar_wont_fix_debt.py` holds the
authoritative rule list and records the reason on each individual resolution,
so a reader of the dashboard sees why without finding this file.

| Rule family | Count | Category | Reason |
|-------------|-------|----------|--------|
| `cpp:S7034` `contains` | 61 | Toolchain | `std::string::contains` is C++23 |
| `cpp:S6185` / `cpp:S6484` `std::format` | 11 | Toolchain | libstdc++ 11 has no `<format>` |
| `cpp:S7035` `to_underlying` | 5 | Toolchain | C++23 |
| `cpp:S6179` `std::lerp` | 11 | Numerical | not bit-identical; FP reproducibility |
| `cpp:S107` param count | 45 | Architecture | diffusion/GPU APIs need a context-struct redesign |
| `cpp:S134` nesting | 41 | Architecture | hot kernels / receptor / GPU |
| `cpp:S3776` / `python:S3776` | 39 | Architecture | parser, HDF5, GPU, batch CLI |
| `cpp:S1820` / `cpp:S1448` | 17 | Architecture | `Simulation` / GPU type size |
| `cpp:S995` const ptr | 5 | Architecture | GPU buffer mutability |
| `cpp:S5008` `void*` | 4 | Architecture | HDF5 C API buffers |
| `cpp:S3656` protected | 1 | Architecture | NUFEB-style `Fix` base contract |
| `cpp:S924` nested break | 1 | Architecture | coupled to `Simulation::run` |

### `cpp:S8379` — held open on purpose

The rule flags 13 `mutable` members as needing a mutex. Twelve are protected by
something the rule does not model: `#pragma omp atomic update` and a per-thread
`kernel_evaluations_by_thread_` vector in `GreensFunction`, `omp critical` in
`FixMetabolism`, and serial-only mutation of `Simulation::fixes_` and the HDF5
`run_provenance_written_` flag.

The thirteenth was real. `QSSASolver::sampled_toxin_conc()` lazily called
`field.samples.resize()` on a shared vector from inside `fix_receptor`'s
`omp parallel for`, reachable whenever an agent count grew since the sampling
pass — routine, because `FixMetabolism::compute()` divides earlier in the same
biology phase. Fixed in PR #371 by making the out-of-range case non-mutating,
which also let `sampled_fields_`/`sampled_nuclease_fields_` drop `mutable` so
the pattern cannot come back silently.

Because of that, the family is **not** in the Won't Fix rule list: the two
`FixMetabolism` findings (`delivery_support_cache_`, `delivery_support_stencil_`)
are still under audit — the fast-path cache `find()` is unsynchronized while a
miss inserts under `omp critical` — and a wholesale resolution would bury them.
Add `cpp:S8379` to the script only once that audit lands.

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

### `pythonsecurity:S8707` / `pythonsecurity:S6549` — path traversal and tainted input paths

**Problem:** Sonar's taint engine does not recognize custom `path_utils` validators, so
`prepare_output_file` / `write_text_file` kept flagging `mkdir` and `open` sinks even after
validation. `S6549` additionally flags paths taken from parsed JSON/HDF5 content
when they reach filesystem calls.

**Fix (two parts):**

1. **Code** (`python/gut_ibm_tools/path_utils.py`):
   - Run full parent/output validation **before** any `mkdir`
   - `_trusted_output_path()` rebuilds cwd-relative paths from regex-sanitized
     segments (breaks S8707 taint from CLI args on SonarCloud automatic analysis)
   - Use `open()` on the trusted path instead of passing user input through
   - `_mkdir_validated_parents()` creates only pre-checked directory segments
   - `prepare_output_directory()` for CLI work dirs (`aws_batch_qa`) — same trusted rebuild
   - `validate_input_path_within()` rebuilds paths declared in parsed file content
     from an explicit trusted root and regex-allowlisted segments before reading them;
     this enforces manifest containment rather than trusting a serialized path

2. **SAST config** (`sonar/pythonsecurity-s8707.json` + `sonar-project.properties`):
   - Registers `validate_path_syntax`, `validate_output_path`, `prepare_output_file`,
     `prepare_output_directory`, `_ensure_output_within_cwd`, and
     `validate_input_path_within` as validators for S8707/S6549
   - Upload the same JSON in SonarCloud → Project Settings → SAST Engine → Python custom
     configuration if the scanner property is not picked up automatically. This
     file is scanner-only; SonarCloud automatic analysis clears these findings
     through the code-structure rebuild in `validate_input_path_within()`, not
     through this configuration.

### GitHub Actions pip (`githubactions:S8541` / `S8544`)

Install third-party deps from the hash-locked
[`.github/requirements-ci.txt`](../.github/requirements-ci.txt) with
`pip install --only-binary=:all: --require-hashes -r ...`. Load the local package via
`PYTHONPATH` (no editable install in CI). Regenerate from
`.github/requirements-ci.in` with `pip-compile --generate-hashes`.

### Docker COPY ownership (`docker:S6504`)

Runtime `COPY` of `gut_ibm` / `entry.sh` uses `--chown=root:root --chmod=755` so the
non-root `gutibm` user can execute but not rewrite the binaries.

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
