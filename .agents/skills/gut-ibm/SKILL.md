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

Serial Release + HDF5 build instructions, the mandatory `g++-13`, and the
serial `-DHDF5_ROOT=.../hdf5/serial` setting are maintained in
`.agents/skills/testing-gutibm/SKILL.md`; use that guide rather than duplicating
the environment-specific recipe here.

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

CTest targets (including custom script and MPI targets; inventory from
`tests/CMakeLists.txt`):

| Test | Labels | Focus |
|------|--------|-------|
| `metabolic_mode` | unit, metabolism | Metabolism mode selection |
| `assertions_enabled` | unit | Test assertions remain enabled |
| `progress_report` | unit | Progress output |
| `dysbiosis_guardrail` | unit | Dysbiosis guard thresholds |
| `spatial_hash` | unit | Insert, query, clear |
| `greens_function` | unit | Radial symmetry and boundary kernels |
| `agent` | unit | Agent pool and plasmid library |
| `iron_fallback` | unit | Secondary iron receptors |
| `octree` | unit | Barnes-Hut FMM versus exact Green's function |
| `fmm` | unit | Higher-order FMM accuracy |
| `qssa_stoichiometry` | unit | QSSA reaction stoichiometry |
| `conjugation` | unit | Pili length heterogeneity and transfer |
| `z_gradient` | unit | Z-dependent nutrient gradients |
| `nutrient_diffusion` | unit | Implicit diffusion residuals and invariants |
| `domain_decomp` | unit | Slab decomposition logic |
| `acetate_mete` | unit | Acetate inhibition of MetE |
| `protease_decay` | unit | Protease half-life decay |
| `oxygen_gradient` | unit | Oxygen registration and z-gradient |
| `O2_growth_boost` | unit | Aerobic metabolism boost |
| `mucin_liberation` | unit | Mucin-to-carbon liberation |
| `mechanism_wiring` | integration | Cross-mechanism mass balance and coupling |
| `advection_peristaltic` | unit | Peristaltic velocity modulation |
| `ethanolamine` | unit | Nutrient-conditional EUT penalty |
| `adaptive_dt` | unit | Adaptive timestep selection |
| `agent_transfer` | unit | Agent pack/unpack round-trip |
| `gpu_receptor_layout` | unit | GPU receptor layout |
| `fix_registry` | unit | Fix plugin registration |
| `path_utils` | unit | Safe path handling |
| `input_parser` | unit | Example and parser configurations |
| `config_ingestion` | unit | Every parser key reaches `SimulationConfig` |
| `kill_provenance` | unit | Kill provenance metadata |
| `counter_resume` | unit | Counter persistence and resume |
| `vbf_accounting` | unit | VBF source/sink accounting |
| `data_defaults` | unit | Analysis data defaults |
| `bacteriocin` | unit | SOS lysis and secretion |
| `per_receptor_toxin` | unit | Per-receptor toxin behavior |
| `toxin_lumping` | unit | Toxin lumping behavior |
| `species_subset` | unit | Species-subset configuration |
| `agent_toxin_sampling` | unit | Agent toxin sampling |
| `siderophore_depletion` | unit | Siderophore depletion |
| `receptor` | unit | TBDT binding and killing |
| `mutation` | unit | BI-locus evolution |
| `fur` | unit | Fur regulation |
| `uptake_limit` | unit | Agent uptake limitation |
| `uptake_limit_population_resolution` | integration, slow | Population-scale delivery resolution |
| `maintenance_carbon` | unit, integration | Maintenance carbon funding |
| `cdi` | unit | Contact-dependent inhibition |
| `motility` | unit | Motility and taxis |
| `quorum_sensing` | unit | AI-2 quorum sensing |
| `immigration` | unit | Immigration schedules |
| `scaling_benchmark` | benchmark, slow | Agent-count timing smoke |
| `mechanics` | unit, slow | Hertzian contact and adhesion |
| `immunity_escape` | unit, slow | Affinity-neutralization matrix |
| `smoke` | integration | End-to-end mini simulation |
| `toxin_sentinel` | integration, science | Toxin sentinel invariants |
| `washout_trap` | integration, slow, science | Long-horizon metabolic washout |
| `feature_combinations` | integration | Feature-combination scenarios |
| `config_diversity` | integration | Distinct config fingerprints |
| `openmp_parity` | openmp | Serial/OpenMP fingerprints |
| `hdf5_roundtrip` | integration, hdf5 | HDF5 writer/reader parity |
| `hdf5_checkpoint` | integration, hdf5 | HDF5 checkpoint writing |
| `hdf5_restart` | integration, hdf5 | C++ checkpoint restart |
| `hdf5_schedule` | integration, hdf5 | HDF5 schedule handling |
| `hdf5_compression` | integration, hdf5 | HDF5 compression |
| `hdf5_live_agents` | integration, hdf5 | Live-agent HDF5 output |
| `summary_events` | integration, hdf5 | Summary event output |
| `hdf5_graceful_shutdown` | integration, hdf5 | Graceful HDF5 shutdown |
| `hdf5_roundtrip_parallel` | integration, mpi, hdf5 | Parallel HDF5 round-trip |
| `hdf5_checkpoint_parallel` | integration, mpi, hdf5 | Parallel HDF5 checkpoint |
| `mpi_multi_rank` | integration, mpi | Two-rank migration |
| `mpi_four_rank` | integration, mpi | Four-rank slab and migration |
| `mpi_gpu_multi_rank` | integration, mpi, gpu | MPI plus GPU chemistry |
| `cuda_aware_mpi_reaction` | integration, mpi, gpu | CUDA-aware MPI reaction reduction |
| `greens_function_gpu` | gpu | GPU/CPU Green's-function parity |
| `gpu_diffusion` | gpu | GPU/CPU nutrient diffusion parity |
| `gpu_chemical_field` | gpu | GPU chemical-field facade parity |
| `gpu_feature_combinations` | gpu, integration | GPU feature combinations |
| `gpu_production_path` | gpu, integration | Production GPU chemistry path |
| `gpu_smoke` | gpu | CPU/GPU smoke fingerprint |
| `gpu_reproducibility` | gpu, integration | GPU reproducibility |
| `gpu_scaling_benchmark` | gpu, benchmark | GPU scaling smoke |
| `qssa_gpu_parity` | gpu | GPU QSSA parity |
| `gpu_toxin_burst_parity` | gpu, integration | GPU toxin burst parity |
| `gpu_nutrient_feedback` | gpu, integration | GPU nutrient feedback |
| `gpu_metabolism_fur` | gpu, integration | GPU metabolism/Fur |
| `gpu_uptake_limit` | gpu, integration | GPU uptake limitation |
| `spatial_hash_gpu_csr` | gpu | GPU spatial hash CSR |
| `mechanics_gpu_parity` | gpu | GPU mechanics parity |
| `gpu_kernel_units` | gpu | GPU kernel units |
| `gpu_target_manifest` | unit | GPU target manifest |
| `rebuild_and_run_script` | unit | Rebuild-and-run script |
| `aws_capacity_script` | unit | AWS capacity script |
| `checkpoint_retention` | unit | Checkpoint retention |
| `aws_entrypoint_runtime_files` | unit | AWS entrypoint runtime files |
| `process_group_signal` | unit | Process-group signal handling |

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

CI runs `ruff check python/`, fast pytest (`-m "not integration"`), JSON config validation (`scripts/validate_config_json.sh`), and EARI/VADI invariant/bounds plus FISH checks (`scripts/validate_eari_vadi.sh`).

## Simulation Timestep (read before touching `simulation.cpp`)

Each biological step in `Simulation::step()`:

1. Clear ghosts, zero reaction accumulators
2. `exchange_ghost_agents()` — MPI boundary ghosts for neighbor queries
3. `rebuild_spatial_hash()`, `update_grid_coupling()`
4. Fix `pre_step(dt)` hooks
5. **Biology** — all Fixes `compute(dt)` (metabolism, bacteriocin, receptor, conjugation, mutation, mechanics)
6. Clear ghosts
7. **Chemistry** — QSSA bacteriocin field, rank-summed nutrient reactions, VBF coupling, grid update, implicit nutrient diffusion, boundaries
8. **Physics** — mucus advection, VBF drag, mechanics repulsion
9. Fix `post_step(dt)` hooks
10. `migrate_agents()` — MPI cross-rank transfer
11. `check_washout()`, `remove_dead_agents()`, `allreduce_global_stats()`

Toxin transport is instantaneous QSSA. Nutrient and small-molecule transport uses L-stable backward-Euler directional diffusion once per `bio_dt` (default 60 s), with no explicit CFL substeps.

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
| `FixQuorumSensing` | `QuorumSensingConfig` | AI-2 production / Lsr import (Spec 11) |
| `FixBacteriocin` | `BacteriocinConfig` | SOS lysis, microcin secretion |
| `FixReceptor` | `ReceptorConfig` | TBDT competitive binding, toxin killing |
| `FixMotility` | `MotilityConfig` | Run-and-reverse + Weber–Fechner taxis (Spec 10v2/11) |
| `FixConjugation` | `ConjugationConfig` | F-pili HGT, shear-dependent MPS |
| `FixCdi` | `CdiConfig` | Contact-dependent inhibition |
| `FixMutation` | `MutationConfig` | BI locus evolution, receptor downregulation |
| `FixMechanics` | `MechanicsConfig` | Hertzian repulsion, EPS adhesion |

## Adding a New Diffusion Kernel

1. Bacteriocins: add QSSA Green's-function logic in `src/diffusion/` and wire it through `QSSASolver`; use Method of Images and consider FMM for large source counts.
2. Nutrients/small molecules: add stable implicit field logic in `src/fields/chemical_field.cpp`; do not use an explicit stencil at `bio_dt`.
3. Preserve periodic x/y and luminal zero-flux z conditions. Epithelial z=0
   defaults to Dirichlet, with runtime-selectable Robin/flux delivery per
   configured species; keep each mode's accounting and solver form consistent.
4. Sum rank-local reaction grids before global VBF coupling and diffusion.
5. Add analytic residual/invariant checks plus enable/coefficient sensitivity, positivity, MPI-equality, and CPU/GPU parity coverage.

## Configuration

### Two ways to configure

1. **C++ tests / programmatic** — build `SimulationConfig` directly or start from `InputParser::default_config()` and override fields.
2. **Input file** — line-oriented pseudo-JSON parsed by `InputParser::parse()`. Supports flat keys plus JSON arrays for `initial_strains` and `fixes`.

### Input file keys (parsed today)

| Key | Maps to |
|-----|---------|
| `total_time`, `bio_dt`, `output_interval` (console/lineage only), `seed` | Time control |
| `domain_x`, `domain_y`, `domain_z`, `grid_dx` | Domain size / resolution |
| `mucus_thickness`, `radial_turnover`, `distal_transit` | Advection |
| `vbf_density`, `vbf_viscosity`, `vbf_mucin_z_gradient`, `vbf_mucin_z_lambda` | VBF |
| `carbon_z_gradient`, `carbon_z_lambda` | Carbon z-gradient |
| `sos_lysis_prob` | Bacteriocin Fix |
| `crypts_enabled`, `crypt_depth`, `crypt_exit_rate`, `crypt_entry_rate`, `crypt_carrying_capacity` | Crypt refugia |
| `hdf5_file`, `hdf5.schedule.*`, `hdf5.compression` | HDF5 output (Spec 4 layered schema) |
| `siderophore.enabled`, `siderophore.secretion_rate`, … | Siderophore dynamics (Spec 4) |
| `profile_steps` | Per-step profiling (`docs/SCALING.md`) |
| `checkpoint_file`, `checkpoint_step` | Checkpoint restart |
| `restart.enabled`, `restart.directory`, `restart.interval_steps` | Closed midstream restarts (Tier 2) |
| `adaptive_dt_enabled`, `dt_min`, `dt_max`, `dt_safety`, `dt_growth_limit` | Adaptive timestep |
| `initial_strains` | JSON array of strain objects |
| `fixes` | JSON array of Fix plugin names (execution order) |
| `kd_colicinE_btuB`, `kill_rate_colicin`, … | Receptor Fix tunables |
| `base_transfer_rate`, `pili_heterogeneity`, … | Conjugation Fix tunables |
| `bi_duplication_rate`, `max_bi_loci`, … | Mutation Fix tunables |
| `use_fmm`, `fmm_theta`, `fmm_expansion_order` | Barnes-Hut FMM |
| `peristaltic_enabled`, `peristaltic_period`, … | Peristaltic advection |
| `oxygen.enabled`, `oxygen.epithelial_conc`, … | O₂ field and aerobic boost (Spec 1) |
| `oxygen.delivery_uptake_enabled` | Oxygen delivery uptake; default `false`, requires delivery uptake mode |
| `oxygen.ros_driver` | ROS driver; default `ambient`, supported values `ambient` and `funded` |
| `oxygen.k_ROS` | Ambient ROS coefficient; default `0.0` (set positive to opt in) |
| `oxygen.k_ROS_respiratory` | Funded specific-flux ROS coefficient; default `0.0` |
| `oxygen.k_ROS_funded` | Funded absolute-flux ROS coefficient in mol⁻¹; default `0.0` |
| `carbon.epithelial_boundary` | Carbon z=0 boundary; default `dirichlet`, with `robin`/`flux` supported |
| `carbon.epithelial_transfer_coeff` | Carbon Robin transfer coefficient; default `0.0`, positive when Robin is selected |
| `carbon.epithelial_flux` | Carbon epithelial flux; default `0.0` |
| `oxygen.epithelial_boundary` | Oxygen z=0 boundary; default `dirichlet`, with `robin`/`flux` supported |
| `oxygen.epithelial_transfer_coeff` | Oxygen Robin transfer coefficient; default `0.0`, positive when Robin is selected |
| `oxygen.epithelial_flux` | Oxygen epithelial flux; default `0.0` |
| `metabolism.uptake_limit` | Agent uptake cap; default `none`, supported `none`, `voxel`, `sherwood`, `delivery` |
| `metabolism.delivery_far_field_radius` | Physical delivery support; default `10 µm`; positive values are refused for slab chemistry, which requires explicit `0.0` |
| `initial_population.placement` | Founder placement; default `legacy`, with `anatomic` and `z_slab` also supported |
| `acetate.enabled`, `acetate.vbf_production`, … | Dynamic acetate (Spec 1) |
| `mucin.enabled`, `mucin.initial_conc`, … | Mucin polymer field (Spec 1) |
| `protease.enabled`, `protease.default_half_life`, … | Protease toxin decay (Spec 1) |
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

Cell-aligned 1D slab decomposition along `DomainConfig::mpi_decomp_axis` (default `0` = x, distal flow). Set in code, not input file. Ghost agents exchanged at slab boundaries; agents migrate via `migrate_agents()`.

**Caveat:** MPI serialization lives in `src/core/agent_transfer.cpp` — update pack/unpack when adding agent or genome fields.

ChemicalField `replicated` storage is the default. In `slab` mode, concentration
and reaction arrays are rank-local owned x-slabs with configured concentration
halos; global agent cell indices map explicitly into storage. Periodic-x
diffusion uses an exact full-line exchange, while y/z operators remain local.
HDF5 grid output and checkpoint/restart support slab mode: owned x-blocks are
gathered for global grid datasets, and checkpoint grids restore only owned
cells before halo refresh. The GPU mirror still rejects slab mode pending a
follow-up implementation.

### HDF5 output

- Writer: `src/io/hdf5_writer.cpp` — parallel when MPI enabled
- Reader: `src/io/hdf5_reader.cpp`; `Simulation::init_from_checkpoint` restores
  C++ checkpoints selected through `checkpoint.*` or closed `restart.*`
  configuration. The `hdf5_checkpoint` and `hdf5_restart` CTest targets cover
  these paths.
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

## Batch runs

Resumable parameter scans: see [docs/BATCH_RUNNER.md](../../docs/BATCH_RUNNER.md).
For **AWS Batch Spot + CUDA** campaign runs (Stage 3 scale), see
[docs/AWS_BATCH.md](../../docs/AWS_BATCH.md) and draft bits under `deploy/aws/`.
Progress / Spot / cost helpers (local clone with AWS profile):

```bash
gut-ibm-aws-status <jobId> --checkpoint-prefix s3://…/ckpt --array-index 0
gut-ibm-aws-estimate --instance-type g5.2xlarge --wall-hours 24 --array-size 12
gut-ibm-aws-qa --output-prefix s3://…/out --input-prefix s3://…/jobs
```

Quick commands:

```bash
cd python && pip install -e ".[dev]"
python -m gut_ibm_tools.batch_runner examples/batch_scan/batch.json --dry-run
python -m gut_ibm_tools.batch_runner examples/batch_scan/batch.json
python -m gut_ibm_tools.batch_runner examples/batch_scan/batch.json --resume
python -m gut_ibm_tools.batch_runner examples/batch_scan/batch.json --status
```

**Cartesian sweep** (`sweep` keys are multiplied):

```json
{
  "output_dir": "batch_results/eari_kd_seed",
  "base_config": "examples/eari_vadi_validation/input.json",
  "binary": "build/gut_ibm",
  "sweep": {
    "seed": [4092, 4093, 4094],
    "kd_colicinE_btuB": [2e-10, 8e-10]
  }
}
```

**Explicit run list**:

```json
{
  "output_dir": "batch_results/custom",
  "base_config": "examples/diversity_paradox/input.json",
  "binary": "build/gut_ibm",
  "runs": [
    { "id": "baseline", "overrides": {} },
    { "id": "short", "overrides": { "total_time": 3600 } }
  ]
}
```

Overrides support dot paths for nested fields (e.g. `"initial_strains.0.count": 200`). Optional `validate` block runs `validation_regression` after each successful simulation. See `examples/batch_scan/batch.json` for a starter manifest; `examples/batch_scan/batch_ci.json` is the single-job smoke used in CI (`scripts/smoke_batch_runner.sh`).

## Spec Documents

- `EARI.md` — Eco-Advective Receptor Interference (advective double-bind)
- `VADI.md` — Viscous Advective-Diffusion Interference (combinatorial washout trap)
- `docs/MECHANISMS.md` — per-Fix biological detail
- `docs/PARAMETERS.md` — full parameter reference (some keys not yet in parser)
- `docs/SCALING.md` — 10⁶–10⁷ agent benchmarks, memory, FMM tuning
- `docs/API.md` — class/function reference

Immigration keys: `immigration.enabled`, `immigration.count`,
`immigration.strain_index`, `immigration.placement`, `immigration.distance`,
`immigration.distance_tolerance`, `immigration.distance_reference`,
`immigration.z_min`, `immigration.z_max`, `immigration.schedule`,
`immigration.step`, and `immigration.rate`.

## GPU

CUDA acceleration lives in `src/gpu/`. Enable with `-DGUTIBM_USE_CUDA=ON` and `gpu_enabled true` at runtime. Falls back to OpenMP/serial CPU when CUDA is unavailable. See `src/gpu/README.md` and issue #33.

## SonarQube

The project is scanned on SonarCloud as **`bckirkup_GutModelBacteriocins`**.

Before fixing or writing code that touches I/O, tests, randomness, or exceptions, read:

1. `.agents/skills/sonarqube-gutibm/SKILL.md` — remediation workflow and PR batching
2. `.agents/skills/sonarqube-cpp/SKILL.md` — C++ rule patterns
3. `.agents/skills/sonarqube-python/SKILL.md` — Python rule patterns

Quick rules: domain exceptions from `error.h`, `to_underlying()` for enums, `gutibm::RNG` for simulation PRNG, `path_utils` / `prepare_output_file` for user paths, `pytest.approx` for float test assertions.
