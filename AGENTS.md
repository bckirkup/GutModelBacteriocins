# AGENTS.md — AI Agent Guidelines for GutIBM

## Repository Purpose

Massively parallel 3D Individual-based Model (IbM) for simulating
*Enterobacteriaceae* population dynamics and genomic diversity in the colonic
mucus layer. C++17 core with Python analysis tools. Built on NUFEB-2
framework philosophy with MPI domain decomposition and HDF5 I/O.

**Version:** 0.1.0 (early research prototype — not production-HPC-ready)

## First Steps

1. Read `.agents/skills/gut-ibm/SKILL.md` for build commands, test inventory, config keys, and Fix registration workflow.
2. Skim `docs/MECHANISMS.md` before editing any Fix module.
3. Before editing I/O, tests, or randomness, read `.agents/skills/sonarqube-gutibm/SKILL.md` (and the C++/Python sub-skills as needed).
4. Check **Known Bugs & Landmines** below before debugging unexpected behavior.

## Setup

```bash
mkdir -p build && cd build
cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON
make -j$(nproc)
```

Python analysis tools:

```bash
cd python && pip install -e ".[viz]"
pip install ruff pytest   # dev tooling
```

## Validation Commands

Run these before committing:

```bash
# C++ build (with warnings)
cd build && cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON \
  -DCMAKE_CXX_FLAGS="-Wall -Wextra" && make -j$(nproc)

# C++ tests
cd build && ctest --output-on-failure

# Python tools (if modified)
cd python && ruff check .
cd python && pytest tests/ -v -m "not integration"
```

## Architecture Rules

- **NUFEB-2 Fix architecture** — biological rules are modular Fix plugins (`FixRegistry` in `src/fixes/fix_registry.cpp`)
- **QSSA for chemistry** — no explicit PDE solvers; analytical Green's function kernels
- **VBF for anaerobes** — 99% of background flora is a continuum field, not discrete agents
- **Spatial hashing** — O(N) neighbor lookups, not O(N²)
- **MPI domain decomposition** — 1D slab along x-axis; ghost exchange + agent migration
- **Never modify tests to make them pass** — fix the implementation

### Timestep modules (do not reorder casually)

```
pre_step → biology (all Fixes compute) → chemistry (QSSA) → physics (advection + mechanics) → post_step → MPI migrate → washout → cleanup
```

Chemistry is instantaneous. Bio timestep (`bio_dt` = 60 s default) is decoupled from chemistry.

## Key Files

| Path | Purpose |
|------|---------|
| `src/core/simulation.cpp` | Main loop, MPI migration, washout |
| `src/core/domain.cpp` | Slab decomposition, `owner_rank()`, ghost bounds |
| `src/core/spatial_hash.*` | O(N) neighbor queries |
| `src/fields/` | ChemicalField, AdvectionField, VBF |
| `src/diffusion/` | Green's function, QSSA solver, Barnes-Hut octree |
| `src/fixes/fix_registry.cpp` | Fix plugin registration and factory |
| `src/genome/plasmid.cpp` | Plasmid library (`ColE1`, `ColB`, …) |
| `src/io/input_parser.cpp` | JSON + legacy flat-key config parser |
| `src/io/hdf5_writer.cpp` | Parallel HDF5 output + genome checkpoint groups |
| `src/io/hdf5_reader.cpp` | Checkpoint restart snapshots |
| `python/gut_ibm_tools/` | HDF5 reader, analysis, validation, visualization |
| `examples/` | `single_colony/`, `diversity_paradox/`, `eari_vadi_validation/` |
| `tests/` | 30 CTest targets (see test map below) |
| `.agents/skills/gut-ibm/SKILL.md` | Hands-on development reference |
| `.agents/skills/sonarqube-gutibm/SKILL.md` | SonarQube remediation workflow |
| `.agents/skills/sonarqube-cpp/SKILL.md` | C++ SonarQube patterns |
| `.agents/skills/sonarqube-python/SKILL.md` | Python SonarQube patterns |

## Key Concepts

- **Advective Double-Bind (EARI)**: Growth/resistance trade-off leading to washout
- **Combinatorial Washout Trap (VADI)**: Physical expulsion when μ_realized < γ_flow
- **Lethal Core/Halo**: Spatial toxin zones based on isoelectric point (pI)
- **Comet-tails**: Downstream-elongated inhibitory zones from mucus flow
- **BI Locus**: Bacteriocin-Immunity genetic region subject to recombination
- **Method of Images**: Boundary condition technique for diffusion kernels
- **MPS (Mating-Pair Stabilization)**: Shear-dependent HGT probability

## Known Bugs & Landmines

| Issue | Status | Notes |
|-------|--------|-------|
| **#40 Metabolic washout** | Fixed | `check_washout()` uses `mu_realized < gamma_flow` |
| **#41 MPI state loss** | Fixed | `agent_transfer.cpp` serializes crypt, affinities, immunity escape |
| **#42 Plasmid names** | Fixed | `PlasmidLibrary::find()` + aliases; warn on unknown names |
| **#43 Multi-rank tests** | Fixed | `mpi_multi_rank` + `hdf5_roundtrip_parallel` CTest targets |
| **#78 parse_real() silent zero** | Fixed | Invalid numerics log warnings; `GUTIBM_STRICT_CONFIG=1` aborts |
| GPU portability | Open | CUDA CI compiles/tests; no multi-GPU or production HPC path |
| Large-scale MPI scaling | Open | 2-rank tests only; manual `mpirun -np 4+` for HPC validation |

When writing tests that involve plasmids, use **`ColE1`/`ColB`** (legacy `colicin_E1` aliases still resolve) and assert `agent.genome.bi_loci.size() > 0`.

## Test Coverage Map

### C++ (CTest — 30 targets)

**Fast unit (`ctest -L unit -LE slow`):** spatial hash, Green's functions, agent/plasmid, iron fallback, octree, FMM, conjugation, z-gradient, domain decomp, acetate/MetE, peristaltic advection, ethanolamine, adaptive dt, agent transfer pack/unpack, fix registry, input parser, bacteriocin, receptor, mutation.

**Slow unit:** mechanics, immunity escape.

**Integration (`ctest -L integration`):** smoke (end-to-end biology), config diversity (fixture/example fingerprints must differ), HDF5 round-trip, HDF5 checkpoint restart.

**MPI:** `mpi_multi_rank` (`mpirun -np 2`), `hdf5_roundtrip_parallel`.

**OpenMP:** `openmp_parity` + `scripts/compare_openmp_parity.sh` (serial vs OpenMP fingerprint).

**GPU (CUDA job):** `greens_function_gpu`, `gpu_smoke` + `scripts/compare_gpu_parity.sh`.

**Benchmark:** `scaling_benchmark` (issue #55 smoke counts).

**Config diversity guardrail:** `test_config_diversity` runs short simulations from parser fixtures and example JSON files and asserts distinct deterministic fingerprints — catches configs silently reverting to defaults.

### Python (pytest in CI)

`python-lint` job runs `ruff check`, import smoke, and `pytest tests/ -m "not integration"`. Coverage includes HDF5 reader, analysis helpers, validation regression helpers, and FISH observation models (#25).

### CI jobs (`.github/workflows/ci.yml`)

| Job | What it exercises |
|-----|-------------------|
| `unit-tests` | Fast CTest unit shard |
| `integration-tests` | Smoke, config diversity, HDF5, slow/benchmark tests |
| `openmp-parity` | Serial vs OpenMP build fingerprints |
| `cuda-compile` | GPU build + parity script |
| `eari-vadi-validation` | Short EARI/VADI + FISH golden regression (#56, #25) |
| `python-lint` | JSON syntax, ruff, pytest (fast) |

### Gaps (add tests when touching these areas)

- Multi-rank MPI beyond 2 processes (`mpirun -np 4+`)
- Python integration pytest in CI (marked `integration` today)
- Metabolic washout trap as a dedicated long-horizon regression
- OpenMP equivalence on stochastic (toxin-kill) scenarios

## Adding Features — Agent Checklist

### New Fix module
1. `src/fixes/fix_<name>.{h,cpp}` inheriting `Fix`
2. Config struct in `input_parser.h` + defaults in `default_config()`
3. Register in `FixRegistry::register_fix()` inside `fix_registry.cpp` (or call from a static initializer)
4. `tests/test_<name>.cpp` + entry in `tests/CMakeLists.txt`
5. Update `docs/MECHANISMS.md` if biological behavior changes

### New diffusion / QSSA kernel
1. `src/diffusion/` — QSSA-compatible, Method of Images
2. Wire through `QSSASolver`
3. Analytical verification test

### MPI-sensitive changes
- Guard all MPI calls with rank checks
- Ensure collectives are called by all ranks
- Update `pack_agent`/`unpack_agent` if adding agent or genome fields
- Extend `test_mpi_multi_rank.cpp` when adding migration-sensitive state

### Python changes
- Import from submodules (`gut_ibm_tools.hdf5_reader`, not top-level)
- Add pytest tests; run `ruff check python/`

### New config keys
- Add to `InputParser::apply_flat_key()` and/or `config_json.cpp`
- Add parser fixture under `tests/fixtures/` + assertion in `test_input_parser.cpp`
- Extend `test_config_diversity.cpp` if the key should change simulation outcomes

## Configuration Quick Reference

| What | How |
|------|-----|
| Run with file | `./gut_ibm examples/single_colony/input.json` |
| Default strains | Resident (`ColE1`+`ColB`, 500) + immigrant (100), no plasmids |
| Strain setup in tests | `cfg.initial_strains` on `SimulationConfig` or `initial_strains` JSON array |
| Fix selection | `fixes` JSON array or `cfg.enabled_fixes` (empty = all registered) |
| Fix tunables | Flat keys (`kd_colicinE_btuB`, `bi_duplication_rate`, …) — see `docs/PARAMETERS.md` |
| Plasmid names | `ColE1`, `ColE2`, `ColB`, `ColIa`, `ColM`, `MccV` |
| HDF5 interval | `hdf5_every` (input file) or `cfg.hdf5.dump_every` (code) |
| Checkpoint restart | `checkpoint_file` + optional `checkpoint_step` in input JSON |
| Disable HDF5 in tests | `cfg.hdf5.enabled = false` |
| Strict config | `GUTIBM_STRICT_CONFIG=1` aborts on invalid numerics |
| MPI decomp axis | `cfg.domain.mpi_decomp_axis` (default 0 = x) |
| Barnes-Hut FMM | `use_fmm`, `fmm_theta`, `fmm_expansion_order` in input JSON |
| Peristaltic mixing | `peristaltic_*` keys in input JSON |
| Chemical environment (Spec 1) | `oxygen.enabled`, `acetate.enabled`, `mucin.enabled`, `protease.enabled` + nested keys in `docs/PARAMETERS.md` |
| Cell biology (Spec 3) | `fur.enabled`, `cdi.enabled`, `motility.enabled` + nested keys; per-strain `cdi_type`, `cdi_immunity` |
| GPU | `gpu_enabled` in input JSON (CUDA build required) |

Full parameter docs: `docs/PARAMETERS.md`.

## Spec Documents

- `EARI.md` — Eco-Advective Receptor Interference framework
- `VADI.md` — Viscous Advective-Diffusion Interference framework
- `docs/MECHANISMS.md` — per-Fix biological mechanisms
- `docs/API.md` — class reference
- `docs/PARAMETERS.md` — configurable parameters
- `docs/SCALING.md` — agent-count benchmarks and profiling

## Code Conventions

- C++17 with modern idioms (smart pointers, RAII, move semantics)
- CMake build system; sources GLOB'd — no manual source list for new `.cpp` in `src/`
- MPI for parallelism (guard all MPI calls with rank checks)
- HDF5 for I/O (parallel HDF5 when MPI enabled)
- Python: NumPy, SciPy, h5py; matplotlib optional via `[viz]` extra
- Fix hooks: `init()`, `pre_step()`, `compute()`, `post_step()` — not `pre_force`/`post_force`

## PR Requirements

- Clean build with `-Wall -Wextra` (no new warnings)
- All CTest tests pass
- New features include tests with **biological outcome assertions**
- MPI-safe (no deadlocks, proper collective operations)
- Update `examples/` if adding user-facing config
- Update this file and `SKILL.md` if changing architecture, config keys, or known bugs

## Open Issue Tracker (Jun 2026)

Recently closed critical path: #40–#43, #56, #75–#81, #25 (FISH models), #29 (higher-order FMM), #55 (scaling benchmarks).

Remaining long-horizon: #33 (GPU production path), larger-scale MPI/HPC validation.

**Project board:** [docs/PROJECT_BOARD.md](docs/PROJECT_BOARD.md) — kanban layout, PR bundles, merge order. Run `./scripts/setup_project_board.sh` to create GitHub labels, milestones, and a Projects v2 board.
