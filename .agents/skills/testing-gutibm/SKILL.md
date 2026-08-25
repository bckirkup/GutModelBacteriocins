---
name: testing-gutibm
description: Test the GutIBM C++/Python simulation end-to-end. Covers build verification, CTest execution, Python pytest, and known preexisting failures. Use when verifying code changes or running post-merge validation.
---

# Testing GutIBM

Shell-based testing only (no GUI/browser). Do not record — collect command outputs as text evidence.

## Quick Verification Sequence

```bash
# 1. Build with warnings (zero new warnings required)
#    g++-13 is REQUIRED: src/core/simulation.cpp includes <format>, and the
#    default /usr/bin/c++ on this VM is too old -> "fatal error: format: No
#    such file or directory". Pass -DCMAKE_CXX_COMPILER=g++-13 for EVERY
#    build dir, including MPI ones (see the MPICH recipe below).
cd build && cmake .. -DCMAKE_CXX_COMPILER=g++-13 \
  -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON \
  -DCMAKE_CXX_FLAGS="-Wall -Wextra" && make -j$(nproc) 2>&1

# 2. Run C++ tests
cd build && ctest --output-on-failure

# 3. Python lint + tests
cd python && ruff check .
cd python && pytest tests/ -v -m "not integration"
```

### Never test a shared, dirty checkout

Multiple agents may share one clone, and a tree that moves mid-build produces
fake errors that look like real code defects (e.g. a header edited between the
compile and your `grep` yields `no declaration matches ...`, and stale build
dirs yield phantom "struct Agent has no member" errors). Before building, run
`git status --short`; if there are unexpected modifications, clone fresh and
check out the exact commit under test:

```bash
git clone https://github.com/bckirkup/GutModelBacteriocins /home/ubuntu/repos/gutibm-test
cd /home/ubuntu/repos/gutibm-test && git checkout <commit>
```

For a merge commit `M`, `M^1` is the pre-merge `main` baseline — use it to
rebuild and rerun any fingerprint-style target you need to classify. To decide
whether a `-Wall -Wextra` warning is new, `git blame` the line and test
`git merge-base --is-ancestor <blame_sha> <M^1>`: false means the PR introduced
it.

Always create a separate build dir per configuration; reconfiguring an existing
CMake cache with a different `CMAKE_CXX_COMPILER` or `HDF5_ROOT` silently keeps
stale objects.

## Preexisting Test Failures

The Devin VM's OpenMPI cannot initialize: `libhwloc` hits an integer
divide-by-zero (SIGFPE) inside `hwloc_topology_load` before `MPI_Init` returns.
Every target that touches MPI therefore fails in an OpenMPI-linked build on ALL
branches — do not attribute these to code changes:

| Test | Failure Mode | Root Cause |
|------|-------------|------------|
| `scaling_benchmark` | MPI init crash (SIGFPE) | `libhwloc` integer divide-by-zero in `hwloc_topology_load` |
| `hdf5_roundtrip` | MPI init crash | Same libhwloc issue |
| `hdf5_checkpoint` | MPI init crash | Same libhwloc issue |
| `hdf5_restart` | MPI init crash | Same libhwloc issue |
| `hdf5_roundtrip_parallel` | Numerical exception | MPI init failure cascade |
| `hdf5_checkpoint_parallel` | Numerical exception | MPI init failure cascade |
| `mpi_multi_rank` | Numerical exception | MPI init failure cascade |
| `mpi_four_rank` | Numerical exception | MPI init failure cascade |

As of the merge of PR #308 the observed OpenMPI-build baseline is **86/94 CTest
passed, 8 failed, 18 `gpu`-labelled Skipped**, with exactly the eight targets
above failing. `config_diversity` and `openmp_parity` now PASS — treat any
movement in those two (or any other fingerprint-style target) as a candidate
regression and confirm against `M^1` before judging. The test count grows over
time; compare the *set of failing target names*, not the raw pass count.

The OpenMPI breakage also means the OpenMPI-linked `gut_ibm` binary cannot run
**even on a single rank** (`orte_ess_init` / "Unable to start a daemon on the
local node"). Any end-to-end simulation run must use the non-MPI build below.

## Running MPI targets: use the MPICH build

MPICH works on this VM when pointed at the minimal hwloc XML. Build a separate
MPICH build dir (note `MPICH_CXX=g++-13`, and the MPICH HDF5 root — mixing the
OpenMPI HDF5 into an MPICH binary is a silent hazard):

```bash
mkdir build-mpich && cd build-mpich
MPICH_CXX=g++-13 cmake .. \
  -DCMAKE_CXX_COMPILER=mpicxx.mpich \
  -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON \
  -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec.hydra \
  -DHDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/mpich -DHDF5_PREFER_PARALLEL=ON
MPICH_CXX=g++-13 make -j$(nproc)

export HWLOC_XMLFILE=$HOME/.config/gutibm/minimal_hwloc.xml
ctest --output-on-failure -R 'mpi_multi_rank|hdf5_roundtrip_parallel'
mpiexec.hydra -n 4 ./tests/test_mpi_four_rank   # see caveat below
```

**Caveat:** the `mpi_four_rank` (and `cuda_aware_mpi_reaction`) CTest commands
hardcode OpenMPI's `--oversubscribe` flag in `tests/CMakeLists.txt`. MPICH's
hydra rejects it with `parse_args ... error parsing input array`, so the target
fails under `ctest` in an MPICH build for tooling reasons only. Invoke the
executable directly with `mpiexec.hydra -n 4` and check for the
`All MPI four-rank tests passed.` line before calling it a pass or a failure.

**Some CTest MPI targets have no same-named executable.** `hdf5_roundtrip_parallel`
is not a build target — it is the *serial* `test_hdf5_roundtrip` binary launched on
2 ranks (see `tests/CMakeLists.txt`). `make test_hdf5_roundtrip_parallel` fails with
`No rule to make target`. Read the `add_test(...)` block to find the real
`$<TARGET_FILE:...>` before concluding a target is missing or broken:

```bash
grep -n -A6 'NAME hdf5_roundtrip_parallel' tests/CMakeLists.txt
MPICH_CXX=g++-13 make -j$(nproc) test_hdf5_roundtrip
mpiexec.hydra -n 2 ./tests/test_hdf5_roundtrip   # expect: All HDF5 round-trip tests passed.
```

Under MPICH, `test_mpi_multi_rank` (2), `test_mpi_four_rank` (4) and
`test_hdf5_roundtrip` (1 and 2 ranks) all pass by hand, which is how you show the
OpenMPI-build CTest failures are environmental rather than code defects.

## Serial end-to-end simulation runs (non-MPI build)

```bash
mkdir build-nompi && cd build-nompi
cmake .. -DCMAKE_CXX_COMPILER=g++-13 -DGUTIBM_USE_MPI=OFF -DGUTIBM_USE_HDF5=ON \
  -DHDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/serial
make -j$(nproc)
./gut_ibm /path/to/input.json
```

`-DHDF5_ROOT=.../hdf5/serial` is mandatory with `GUTIBM_USE_MPI=OFF`: CMake
otherwise picks the OpenMPI HDF5, whose `H5public.h` includes `<mpi.h>` and
fails to compile.

**Runtime budget:** `examples/single_colony/input.json` uses a 1 mm x 1 mm x
0.1 mm domain at `grid_dx=2e-6` (12.5 M cells) and its default
`total_time=86400` is 1440 steps — on the order of hours serially, and it slows
down as the population grows. For verification, shorten `total_time` (1800 s /
30 steps is enough for flux-ledger and inertness properties; 7200 s / 120 steps
if you need a longer trajectory) and assert full-horizon completion against the
*shortened* horizon.

**Memory:** one single-colony run holds roughly 2.7 GB RSS on a ~7.9 GB VM. Run
E2E arms strictly **one at a time** — two concurrent arms invite the kernel OOM
killer (`dmesg | grep -i oom-kill`), which silently truncates a run.

**Never trust stdout for horizon completion.** A run that was OOM-killed,
stop-requested or otherwise cut short still prints `Simulation complete.` with a
`Steps taken:` line. Only the HDF5 provenance distinguishes them:

```python
/run_provenance/halt_reason_code          # 0 = no dysbiosis guard halt
/run_provenance/termination_reason_code   # 0 = reached horizon, 1 = dysbiosis, 2 = cut short
/run_provenance/termination_step          # must equal your configured step count
/run_provenance/termination_time          # must equal your configured total_time
/run_provenance/completed_total_time      # 1 = horizon reached
```

Assert `termination_reason_code == 0` **and** `termination_time == total_time`;
report both. A process that exits 0 has not necessarily reached its horizon.

**`/run_provenance/git_sha` is stamped at CMake *configure* time, not build
time.** If you `git checkout` a new commit and only re-run `make`, the binary
contains the new code but the HDF5 records the *old* SHA, which makes output
look like it came from the wrong commit. Re-run `cmake ..` after every checkout
in a reused build dir, and cross-check the recorded SHA against
`git rev-parse HEAD` before trusting a provenance-based claim.

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

Result must be 0. A clean `g++-13 -Wall -Wextra` build of the whole tree emits
only a handful of warnings (unused function/variable in `src/gpu/fmm_gpu.cpp`,
`tests/test_greens_function_gpu.cpp`, `tests/test_mechanics.cpp`), so any
warning in a file your change touches stands out. A common new-code finding is
`-Wsign-compare` from comparing `sim.agents().size()` (returns `Int`) against a
`std::vector::size()` (`size_t`) — cast one side.

## Python Test Environment

```bash
cd python && pip install -e ".[dev]"  # installs ruff, pytest, etc.
```

If a venv already exists in another clone (e.g. `python/.venv311`), you can
reuse its interpreter against a fresh clone without reinstalling — `PYTHONPATH`
wins over the editable install's path hook:

```bash
VENV=/path/to/other/clone/python/.venv311/bin/python
$VENV -m ruff check python/
cd python && PYTHONPATH=$PWD $VENV -m pytest tests/ -m "not integration" -q
# verify you imported the right tree:
PYTHONPATH=$PWD $VENV -c "import gut_ibm_tools; print(gut_ibm_tools.__file__)"
```

Python tests: 131 non-integration tests pass, 2 deselected (integration marker)
as of PR #308/#309.

## Language standard: C++20, not C++23

`CMakeLists.txt` sets `CMAKE_CXX_STANDARD 20`. C++23-only library APIs compile in
no build dir here — the recurring offender is `std::string::contains` in tests
(`error: 'const std::string' has no member named 'contains'`). Use
`s.find(x) != std::string::npos`. A single such line makes the whole target
unbuildable, so a "passing" CTest run that simply never built the target can hide
it: after building, confirm the executables you care about actually exist.

## Proving a refactor is numerically inert

For a refactor claimed to be behaviour-preserving, build the **same**
configuration at the merge base and at the head, run the identical input under
both, and diff the whole per-step summary series rather than the final line:

```bash
git worktree add /home/ubuntu/repos/gutibm-base <merge_base>
# identical non-MPI configure in each, then run the same config JSON into two .h5
```

Walk every dataset under each `/summary/step_XXXXXX/` group (`n_total`,
`num_agents`, `chem/mean_*`, all `nutrient_flux/*_cumulative`) and report the
first diverging step plus field. Decide **before** running which fields an
intended behaviour change is allowed to move, and treat everything else as a
finding. Re-running the same binary twice is a cheap determinism control and also
re-stamps provenance (see `git_sha` caveat above).

## SonarQube-sensitive test patterns

When writing or reviewing Python tests, follow `.agents/skills/sonarqube-python/SKILL.md`:

- Float comparisons: `pytest.approx()` or `np.testing.assert_allclose` — never `== 0.0`
- Paths in tests: use `tmp_path` fixtures; production path helpers live in `gut_ibm_tools.path_utils`

When writing C++ tests with random fixtures, use `gutibm::RNG` (see `.agents/skills/sonarqube-cpp/SKILL.md`), not `std::mt19937`.

## Validation output considerations

The EARI/VADI CI scenario uses finite-value, bounds, liveness, and FISH
relationship checks rather than stored trajectory baselines. Keep documented
full-run biological targets in `validation_regression.py`; do not turn a short
run's observed values into new thresholds.

## GPU Runtime and Memory Validation

A passing hosted `cuda-compile` job proves CUDA compilation, not physical GPU
execution. Hosted runners may have the CUDA toolkit without a visible device,
and GPU tests can return success after their no-device skip.

Before claiming GPU runtime behavior passed:

```bash
command -v nvcc
nvidia-smi
ctest --test-dir build-cuda -L gpu --output-on-failure -V
```

- Require `nvcc` and a device listed by `nvidia-smi`.
- Use verbose CTest output or run the executable directly so `SKIPPED` messages
  are visible.
- Report CUDA compilation, no-device skips, and physical GPU execution as three
  separate assertions.
- For memory-sensitive grid tests, compute
  `nx × ny × nz × species × 2 × sizeof(double)` for concentration and reaction
  arrays, then account for simultaneous CPU and GPU field owners.
- Preserve GPU-specific coverage when shrinking fixtures. For directional
  diffusion, retain any intentionally long axis while reducing transverse
  dimensions.

## Verifying nutrient-flux / uptake accounting changes

Flux ledgers are written per summary step under
`/summary/step_XXXXXX/nutrient_flux/`, with a `species_names` dataset to index
by species (`carbon`, ...). Useful datasets:
`reaction_clip_cumulative`, `agent_uptake_cumulative`,
`uptake_demand_cumulative`, `maintenance_cumulative`,
`uptake_shortfall_cumulative`.

To prove an uptake-accounting change actually altered the physics rather than
renaming a counter, run **two arms of the same config** differing only in
`metabolism.uptake_limit` and contrast the ledgers — e.g. `delivery` should show
`reaction_clip_cumulative[carbon]` exactly `0.0` while the default `none` arm
still shows a nonzero clip. Also assert
`agent_uptake_cumulative <= uptake_demand_cumulative` and that the population is
neither instantly extinct nor unbounded. A single arm cannot distinguish a real
fix from a relabelled counter.

Config-rejection checks are worth running too: an invalid or unsupported
combination should abort with a message naming the offending keys. Note the
binary reports these as an uncaught `gutibm::ConfigError` via `terminate`
(exit 134 / SIGABRT), not a graceful exit 1 — assert on the message text plus a
nonzero exit rather than on a specific code.

## Adversarial Tests for Exception/Type Changes

When modifying exception types or enum conversions, verify at runtime, not just at compile time:

- **Exception hierarchy:** Run `test_input_parser`, `test_smoke`, `test_agent` — these exercise throw/catch paths
- **`to_underlying()` polyfill:** Run `test_receptor`, `test_iron_fallback` — these use enum-to-int in biological computations where wrong values produce wrong outcomes
- **Python RNG:** Verify determinism with same seed and divergence with different seeds

## Devin Secrets Needed

None — all testing is local. No external services or credentials required.

## Testing delivery-mode uptake from a JSON config (end-to-end)

`metabolism.uptake_limit="delivery"` plus `metabolism.delivery_far_field_radius`
(physical metres; `0.0` = historical single-voxel deposition) needs several
non-obvious config prerequisites before any assertion is meaningful.

### Config keys that are easy to get wrong

- Configs use **flat dotted keys** (`"metabolism.uptake_limit"`, `"grid_dx"`,
  `"domain_x"`), but a few settings are only reachable through **nested objects**
  handled by `ConfigJson` — notably `"domain": {"grid_halo_width": N}` and
  `"hdf5": {"schedule": {...}}`. A flat `"grid_halo_width"` is silently ignored by
  the JSON path even though `input_parser.cpp` accepts it on the legacy path.
- **`hdf5.schedule.grid_species` defaults to empty, so NO grid arrays are
  written.** Set `"grid_species": ["all"]` (with `"grid": 1`) or you cannot check
  for negative concentrations at all. Summary only exposes `chem/mean_*` — there
  is no min. Grid arrays live at `grid/step_XXXXXX/<species>`.
- Prove the *default* of a parameter actually shipped by omitting the key and
  reading `/run_provenance/resolved_config` from the HDF5 file, which serializes
  every resolved value. This is the only trustworthy evidence that a default
  took effect rather than a config value.
- `oxygen.delivery_uptake_enabled=true` **requires** `uptake_limit="delivery"`
  (hard `ConfigError`), so it cannot be combined with a legacy control arm.
- Slab chemistry requires `grid_halo_width >= ceil(ghost_width / dx_x)`
  (`ghost_width` default 10 µm). Beware floating point: a 1.8e-4 domain at
  "2 µm" gives dx slightly under 2e-6, so the requirement is **6**, not 5. When
  testing that a config is refused for reason X, make it otherwise valid first,
  or you will attribute an unrelated refusal to X.

### Choosing a regime, or you will measure nothing

The delivery sink is **first-order and implicitly integrated** in the
backward-Euler solves, so it cannot by itself drive a cell negative. Consequences:

- `delivery_retry_events_*` is usually **0**, and under stress it **saturates**
  at ~3 events/step (`kMaxDeliveryRetries = 2` per solve × ~3 directional
  solves) — identical at any radius. Retry counts are therefore a poor
  discriminator; always normalize per step and state the ceiling.
- The useful observables are the **funded fraction**
  `(agent_uptake_cumulative + maintenance_cumulative) / uptake_demand_cumulative`,
  `delivery_reduction_cumulative`, and the **local minimum concentration** from
  the grid arrays.
- Regimes are bimodal and most configs are degenerate. Measured with default
  carbon (`carbon_boundary_conc=5.0e-3`), 180–200 µm domain, dx 2 µm, 10 steps:
  at ≤500 agents with `mu_max=5.5e-4` the funded fraction is exactly 1.0000 at
  *every* radius (contrast inert); only around 1000 agents with `mu_max≈5.5e-3`
  does it fall below 1 enough to compare. Sweep agent count and `mu_max` to find
  a non-degenerate regime **before** asserting a contrast.
- Do **not** create scarcity by lowering `carbon_boundary_conc` far below the
  hard-coded 5.0e-3 initial concentration (`input_parser.cpp:163`, there is no
  carbon initial-conc key). That imposes a huge artificial boundary/initial
  gradient, fills the field with negatives and can invert the result you are
  testing. Raise demand instead.
- Compare radii against the **voxel vs support-sphere mass scale**: a 2 µm voxel
  holds ~4e-20 mol carbon at 5e-3, while a 10 µm sphere holds ~2.1e-17 mol
  (~520 voxels). Per-agent per-step demand must sit between those to be
  delivery-limited at `0.0` but not at 10 µm.

### Caveat when comparing two radii (or any two physical settings)

Changing the radius changes uptake, hence growth, hence the population, so the
two arms are **different trajectories with different demand** — not a controlled
measurement at fixed demand. Small funded-fraction gaps are partly dynamical
divergence. More robust comparisons: the minimum concentration field, and
**grid-independence** (run each setting at two `grid_dx` values over the same
physical domain and compare how much each moves). A physical-radius support
should move less across resolutions than voxel-local deposition.

### Known pre-existing negative-concentration behaviour (do not blame a PR)

- With `oxygen.delivery_uptake_enabled=true`, oxygen goes substantially negative
  (order −1e-5 against a 5.5e-5 epithelial scale, over half the cells) at
  **both** radius 10 µm and radius 0.0, and not at all when the flag is off.
- Under high carbon stress, carbon also goes negative at both radii.
- Retries fire and `delivery_reduction_*` is nonzero yet negatives persist, so
  the retry/reduction machinery does not guarantee positivity.
Always run the radius-0.0 (or feature-off) control before reporting negativity
as a regression.

### Config rejection ergonomics

Invalid delivery configurations (`gpu_enabled` + delivery, bogus `uptake_limit`,
positive radius + slab chemistry) surface as an **uncaught `gutibm::ConfigError`
→ `terminate` → SIGABRT, exit 134**, not a graceful exit with a clean message on
stderr only. The message text is correct and names the offending keys; the exit
path is ugly. Report it as ergonomics, and note that when the same message is
raised from two enforcement sites (parser `finalize_config` and
`Simulation::init`), stdout alone cannot tell you which one fired.
