---
name: gut-ibm-development
description: Build, test, and develop the GutIBM C++/MPI simulation. Covers CMake configuration, running tests, adding Fix modules, and using the Python analysis tools.
---

# GutIBM Development Skill

Hands-on reference for building, testing, and extending GutIBM. Read `AGENTS.md` for
architecture rules, known bugs, and landmines.

## Build

```bash
mkdir -p build && cd build
cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON
make -j$(nproc)
```

Debug build with warnings (run before every commit):

```bash
cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-Wall -Wextra"
make -j$(nproc)
```

OpenMP (optional, off by default):

```bash
cmake .. -DGUTIBM_USE_OPENMP=ON ...
```

CUDA GPU acceleration (optional, off by default):

```bash
cmake .. -DGUTIBM_USE_CUDA=ON ...
# Runtime: gpu_enabled true in input file, or cfg.gpu.enabled = true in code
```

Prerequisites (Ubuntu):

```bash
sudo apt-get install cmake libopenmpi-dev openmpi-bin libhdf5-mpi-dev
# Optional GPU build:
sudo apt-get install nvidia-cuda-toolkit
```

Sources are collected via `file(GLOB_RECURSE ...)` in `CMakeLists.txt` — new `.cpp`
files under `src/{core,fields,diffusion,fixes,genome,io}/` are picked up automatically.
No manual source-list edit needed for new Fix or diffusion files.

## Run Tests

```bash
cd build && ctest --output-on-failure
```

CI shards tests by label (see `.github/workflows/ci.yml`):

```bash
ctest -L unit -LE slow --output-on-failure      # fast unit gate
ctest -L 'integration|slow|benchmark' -LE gpu   # integration job
ctest -L gpu                                     # CUDA job only
```

30 CTest targets (2-rank MPI tests included; no `mpirun -np 4+` gate yet):

| Test | File | Focus |
|------|------|-------|
| `spatial_hash` | `test_spatial_hash.cpp` | Insert, query, clear |
| `greens_function` | `test_greens_function.cpp` | Radial symmetry, comet-tail, pI classes |
| `agent` | `test_agent.cpp` | Agent pool, plasmid library |
| `iron_fallback` | `test_iron_fallback.cpp` | Secondary iron receptors |
| `octree` | `test_octree.cpp` | Barnes-Hut FMM vs exact GF |
| `fmm` | `test_fmm.cpp` | Higher-order FMM accuracy |
| `conjugation` | `test_conjugation.cpp` | Pili length heterogeneity |
| `smoke` | `test_smoke.cpp` | End-to-end mini simulation |
| `config_diversity` | `test_config_diversity.cpp` | Distinct fingerprints across configs |
| `z_gradient` | `test_z_gradient.cpp` | Z-dependent nutrient gradients |
| `domain_decomp` | `test_domain_decomp.cpp` | Slab logic (single rank) |
| `acetate_mete` | `test_acetate_mete.cpp` | Acetate inhibition of MetE |
| `advection_peristaltic` | `test_advection_peristaltic.cpp` | Peristaltic velocity modulation |
| `immunity_escape` | `test_immunity_escape.cpp` | Affinity-neutralization matrix |
| `mechanics` | `test_mechanics.cpp` | Hertzian contact, adhesion |
| `ethanolamine` | `test_ethanolamine.cpp` | Nutrient-conditional eut penalty |
| `openmp_parity` | `test_openmp_parity.cpp` | Determinism + cross-build fingerprint |
| `adaptive_dt` | `test_adaptive_dt.cpp` | Adaptive timestep selection |
| `agent_transfer` | `test_agent_transfer.cpp` | MPI pack/unpack round-trip |
| `fix_registry` | `test_fix_registry.cpp` | Default Fix plugin registration |
| `input_parser` | `test_input_parser.cpp` | Example JSON config files |
| `bacteriocin` | `test_bacteriocin.cpp` | SOS lysis / secretion unit tests |
| `receptor` | `test_receptor.cpp` | TBDT binding / killing unit tests |
| `mutation` | `test_mutation.cpp` | BI locus evolution unit tests |
| `hdf5_roundtrip` | `test_hdf5_roundtrip.cpp` | Writer/reader schema parity |
| `hdf5_checkpoint` | `test_hdf5_checkpoint.cpp` | Checkpoint restart |
| `mpi_multi_rank` | `test_mpi_multi_rank.cpp` | 2-rank agent migration |
| `scaling_benchmark` | `test_scaling_benchmark.cpp` | Agent-count timing smoke |
| `greens_function_gpu` | `test_greens_function_gpu.cpp` | GPU vs CPU GF parity (CUDA build) |
| `gpu_smoke` | `test_gpu_smoke.cpp` | Short CPU vs GPU simulation fingerprint |

## Run Simulation

```bash
cd build
mpirun -np 4 ./gut_ibm ../examples/single_colony/input.json
```

Single-process (debugging):

```bash
./gut_ibm ../examples/single_colony/input.json
```

No config argument → `InputParser::default_config()` (resident + immigrant strains).

Examples: `examples/single_colony/input.json`, `examples/diversity_paradox/input.json`.

## Python Analysis Tools

```bash
cd python
pip install -e ".[viz]"
pip install ruff pytest   # dev tools (not in setup.py extras yet)
```

Import from package root (re-exported in `__init__.py`):

```python
from gut_ibm_tools import GutIBMData, analysis, validation, visualization
```

CI runs `ruff check python/`, fast pytest (`-m "not integration"`), JSON config validation (`scripts/validate_config_json.sh`), and EARI/VADI golden regression (`scripts/validate_eari_vadi.sh`).

## Simulation Timestep (read before touching `simulation.cpp`)

Each biological step in `Simulation::step()`:

1. Clear ghosts, zero reaction accumulators
2. `exchange_ghost_agents()` — MPI boundary ghosts for neighbor queries
3. `rebuild_spatial_hash()`, `update_grid_coupling()`
4. Fix `pre_step(dt)` hooks
5. **Biology** — all Fixes `compute(dt)` (metabolism, bacteriocin, receptor, conjugation, mutation, mechanics)
6. Clear ghosts
7. **Chemistry** — QSSA bacteriocin field, nutrient depletion, VBF coupling, grid update
8. **Physics** — mucus advection, VBF drag, mechanics repulsion
9. Fix `post_step(dt)` hooks
10. `migrate_agents()` — MPI cross-rank transfer
11. `check_washout()`, `remove_dead_agents()`, `allreduce_global_stats()`

Chemistry is instantaneous (QSSA); bio timestep (`bio_dt`, default 60 s) is decoupled.

## Adding a New Fix Module

1. Create `src/fixes/fix_<name>.h` and `src/fixes/fix_<name>.cpp`
2. Inherit from `Fix` (`src/fixes/fix.h`)
3. Implement `compute(Real dt)` (required); override `init()`, `pre_step()`, `post_step()` as needed
4. **Register in `FixRegistry`** — add `register_fix("my_fix", ...)` in `src/fixes/fix_registry.cpp` (or call `FixRegistry::register_fix` from module init)
5. Add config struct to `src/io/input_parser.h` (`SimulationConfig`) with sensible defaults in `default_config()`
6. Add test in `tests/test_<name>.cpp` and register in `tests/CMakeLists.txt` (copy an existing test block)
7. Update `examples/` if the feature is user-facing

Current Fix modules (hardcoded order in `simulation.cpp`):

| Fix | Config struct | Role |
|-----|---------------|------|
| `FixMetabolism` | `MetabolismConfig` | Triple Monod growth, iron fallback, division/death |
| `FixBacteriocin` | `BacteriocinConfig` | SOS lysis, microcin secretion |
| `FixReceptor` | `ReceptorConfig` | TBDT competitive binding, toxin killing |
| `FixConjugation` | `ConjugationConfig` | F-pili HGT, shear-dependent MPS |
| `FixMutation` | `MutationConfig` | BI locus evolution, receptor downregulation |
| `FixMechanics` | `MechanicsConfig` | Hertzian repulsion, EPS adhesion |

## Adding a New Diffusion Kernel

1. Add kernel logic in `src/diffusion/` (see `greens_function.cpp`, `octree.cpp`)
2. Use Method of Images for boundary conditions; must be QSSA-compatible (steady-state, no explicit dt)
3. Wire into `QSSASolver` (`src/diffusion/qssa_solver.cpp`) or field update in `src/fields/`
4. Add unit test in `tests/` comparing against analytical solution
5. For large N, consider Barnes-Hut FMM (`QSSAConfig::use_fmm`, `fmm_theta`)

## Configuration

### Two ways to configure

1. **C++ tests / programmatic** — build `SimulationConfig` directly or start from `InputParser::default_config()` and override fields.
2. **Input file** — line-oriented pseudo-JSON parsed by `InputParser::parse()`. Supports flat keys plus JSON arrays for `initial_strains` and `fixes`.

### Input file keys (parsed today)

| Key | Maps to |
|-----|---------|
| `total_time`, `bio_dt`, `output_interval`, `seed` | Time control |
| `domain_x`, `domain_y`, `domain_z`, `grid_dx` | Domain size / resolution |
| `mucus_thickness`, `radial_turnover`, `distal_transit` | Advection |
| `vbf_density`, `vbf_viscosity`, `vbf_mucin_z_gradient`, `vbf_mucin_z_lambda` | VBF |
| `carbon_z_gradient`, `carbon_z_lambda` | Carbon z-gradient |
| `sos_lysis_prob` | Bacteriocin Fix |
| `crypts_enabled`, `crypt_depth`, `crypt_exit_rate`, `crypt_entry_rate`, `crypt_carrying_capacity` | Crypt refugia |
| `hdf5_file`, `hdf5_every` | HDF5 output |
| `profile_steps` | Per-step profiling (`docs/SCALING.md`) |
| `checkpoint_file`, `checkpoint_step` | Checkpoint restart |
| `adaptive_dt_enabled`, `dt_min`, `dt_max`, `dt_safety`, `dt_growth_limit` | Adaptive timestep |
| `initial_strains` | JSON array of strain objects |
| `fixes` | JSON array of Fix plugin names (execution order) |
| `kd_colicinE_btuB`, `kill_rate_colicin`, … | Receptor Fix tunables |
| `base_transfer_rate`, `pili_heterogeneity`, … | Conjugation Fix tunables |
| `bi_duplication_rate`, `max_bi_loci`, … | Mutation Fix tunables |
| `use_fmm`, `fmm_theta`, `fmm_expansion_order` | Barnes-Hut FMM |
| `peristaltic_enabled`, `peristaltic_period`, … | Peristaltic advection |
| `gpu_enabled` | GPU acceleration (CUDA build) |

Invalid numerics log warnings; `GUTIBM_STRICT_CONFIG=1` aborts. Unknown keys are ignored.

### Initial strains (JSON array or code)

```cpp
SimulationConfig::InitialStrain s;
s.type = 1;
s.count = 500;
s.mu_max = 5.5e-4;
s.plasmids = {"ColE1", "ColB"};  // exact match required — see below
s.conjugative = true;
cfg.initial_strains.push_back(s);
```

### Fix selection (JSON array or code)

```json
"fixes": ["metabolism", "bacteriocin", "receptor", "conjugation", "mutation", "mechanics"]
```

```cpp
cfg.enabled_fixes = {"metabolism", "mechanics"};  // empty = all defaults
```

Register new Fix factories in `src/fixes/fix_registry.cpp` — no `simulation.cpp` edits needed.

### Plasmid library names (canonical)

Use **exact** names from `PlasmidLibrary::entries()` in `src/genome/plasmid.cpp`:

| Name | Receptor | Notes |
|------|----------|-------|
| `ColE1` | BtuB | Lethal Core (pI ~9) |
| `ColE2` | BtuB | Lethal Halo |
| `ColB` | FepA | Lethal Halo, conjugative |
| `ColIa` | CirA | Neutral, conjugative |
| `ColM` | FhuA | Lethal Core |
| `MccV` | CirA | Microcin V, conjugative |

**Do not use** unknown plasmid names — `PlasmidLibrary::find()` resolves canonical names (`ColE1`, …) and legacy aliases (`colicin_E1`, …). Unknown names log a warning and spawn agents without BI loci.

### MPI decomposition

1D slab decomposition along `DomainConfig::mpi_decomp_axis` (default `0` = x, distal flow). Set in code, not input file. Ghost agents exchanged at slab boundaries; agents migrate via `migrate_agents()`.

**Caveat:** MPI serialization lives in `src/core/agent_transfer.cpp` — update pack/unpack when adding agent or genome fields.

### HDF5 output

- Writer: `src/io/hdf5_writer.cpp` — parallel when MPI enabled
- **No C++ checkpoint restart** — write-only; Python `GutIBMData` reads output
- Layout per step: `atoms/`, `grid/`, `metadata/`, `lineage/`
- Disable in tests: `cfg.hdf5.enabled = false`

## Writing Good Tests

```cpp
SimulationConfig cfg = InputParser::default_config();
cfg.domain.hi = {100e-6, 100e-6, 50e-6};
cfg.hdf5.enabled = false;
cfg.initial_strains.clear();
// ... configure strains with ColE1 etc. ...
Simulation sim;
sim.init(cfg);
sim.run();
// Assert biological outcomes, not just crash-free execution
```

Always assert mechanism outcomes when testing biology (e.g. `bi_loci.size() > 0`, kill counts, washout events).

For config keys and parser fixtures, extend `test_config_diversity.cpp` so distinct settings produce distinct simulation fingerprints — this catches silent overrides to `default_config()`.

## Spec Documents

- `EARI.md` — Eco-Advective Receptor Interference (advective double-bind)
- `VADI.md` — Viscous Advective-Diffusion Interference (combinatorial washout trap)
- `docs/MECHANISMS.md` — per-Fix biological detail
- `docs/PARAMETERS.md` — full parameter reference (some keys not yet in parser)
- `docs/SCALING.md` — 10⁶–10⁷ agent benchmarks, memory, FMM tuning
- `docs/API.md` — class/function reference

## GPU

CUDA acceleration lives in `src/gpu/`. Enable with `-DGUTIBM_USE_CUDA=ON` and `gpu_enabled true` at runtime. Falls back to OpenMP/serial CPU when CUDA is unavailable. See `src/gpu/README.md` and issue #33.
