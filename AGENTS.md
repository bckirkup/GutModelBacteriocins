# AGENTS.md — AI Agent Guidelines for GutIBM

## Repository Purpose

Massively parallel 3D Individual-based Model (IbM) for simulating
*Enterobacteriaceae* population dynamics and genomic diversity in the colonic
mucus layer. C++20 core with Python analysis tools. Built on NUFEB-2
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
- **Hybrid chemical transport** — analytical QSSA Green's functions for bacteriocins; stable implicit directional diffusion for nutrients and small molecules
- **VBF for anaerobes** — 99% of background flora is a continuum field, not discrete agents
- **Spatial hashing** — O(N) neighbor lookups, not O(N²)
- **MPI domain decomposition** — cell-aligned 1D slab along x-axis; ghost exchange + agent migration
- **Chemical decomposition** — `replicated` is the default; `slab` stores owned x-cells
  with concentration halos, local y/z operators, and an exact periodic-x exchange.
  Slab HDF5 grid output gathers owned cells into global datasets, and slab
  checkpoint/restart restores owned cells before refreshing halos. The GPU
  mirror uses local x storage and owned-cell accounting in slab mode; periodic
  x diffusion takes an exact host transpose round trip.
  Integer per-axis chemistry strides are a scientific model variant, not a
  performance switch: carbon diffuses approximately 77 µm in 60 s versus the
  default 2 µm lateral cell, so x/y coarsening tests lateral microgradient
  effects while z remains fine for epithelial oxygen/mucin stratification.
- **Never modify tests to make them pass** — fix the implementation

### Timestep modules (do not reorder casually)

```
pre_step → biology (all Fixes compute) → chemistry (QSSA toxins + implicit nutrient diffusion) → physics (advection + mechanics) → post_step → MPI migrate → washout → cleanup
```

Chemical transport is applied once per biological step. Toxins use instantaneous QSSA; nutrient reactions are followed by L-stable backward-Euler directional diffusion at `bio_dt` (60 s default).

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
| `docs/DELIVERY_ROS_CAMPAIGN.md` | Consolidated delivery/ROS campaign findings and retired claims |
| `python/gut_ibm_tools/` | HDF5 reader, analysis, validation, visualization |
| `python/gut_ibm_tools/colony.py` and `spatial_stats.py` | Colony catalogs and 3-D spatial observables |
| `examples/` | `single_colony/`, `diversity_paradox/`, `eari_vadi_validation/`, `cell_biology/`, `batch_scan/`, `scaling_benchmark/` |
| `tests/` | CTest targets (see test map below) |
| `.agents/skills/gut-ibm/SKILL.md` | Hands-on development reference |
| `.agents/skills/gutibm-campaign-ops/SKILL.md` | AWS campaign operations and checkpoint safety |
| `.agents/skills/sonarqube-gutibm/SKILL.md` | SonarQube remediation workflow |
| `.agents/skills/sonarqube-cpp/SKILL.md` | C++ SonarQube patterns |
| `.agents/skills/sonarqube-python/SKILL.md` | Python SonarQube patterns |
| `docs/AWS_BATCH.md` | AWS Batch Spot + CUDA deployment plan + Phase Observability |
| `deploy/aws/` | Draft Dockerfile / Batch entrypoint / status heartbeat for that plan |

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
| **#40 Metabolic washout** | Fixed | `washout.trap=emergent` is the default and actual transport determines `outflow_boundary`; `imposed` retains the explicit comparison variant as `outflow_washout`. |
| **#41 MPI state loss** | Fixed | `agent_transfer.cpp` serializes crypt, affinities, immunity escape |
| **#42 Plasmid names** | Fixed | `PlasmidLibrary::find()` + aliases; warn on unknown names |
| **#43 Multi-rank tests** | Fixed | `mpi_multi_rank` + `hdf5_roundtrip_parallel` CTest targets |
| **#78 parse_real() silent zero** | Fixed | Invalid numerics log warnings; `GUTIBM_STRICT_CONFIG=1` aborts |
| **Dysbiosis guard noise sensitivity** | Fixed | The former per-sample strict-increase and non-deceleration criterion was effectively unfirable on stochastic trajectories: one density dip reset the seven-sample window. The guard now uses a positive net rise and aggregate half-window increment means while retaining the above-threshold requirement. |
| **Dysbiosis halt artifact visibility** | Fixed | Run-level termination metadata is written under `/run_provenance/`, so guard-halted runs are distinguishable from full-horizon completions even when the final summary was written before the guard evaluation. |
| **Carbon VBF overdraw** | Fixed | The Monod carbon sink is implicitly integrated on both host and CUDA paths from one shared closed-form helper, and realized removal is reported; reaction clips are accounted. Agent-side uptake overdraw is closed by the `uptake_limit` model: `delivery` integrates the per-agent first-order sink with diffusion and funds growth only from realized removal; `sherwood` caps per-agent per-step uptake at the quasi-steady diffusive delivery `4*pi*D_eff*r*C_local`, funds the biomass increment and `mu_realized` at the same fraction, and reports demanded uptake plus bound-agent counts; `voxel` is a diagnostic-only grid-artifact probe and `none` (default) keeps the historical unfunded behaviour. The legacy zero-order agent sink could not brake density because it was clipped before diffusion in a single voxel, so nearly all requested carbon was never removed from the field. Under gradient-preserving transport, the diffusion solve subtracts the reference profile to solve the perturbation `P`, includes the `-sG` sink RHS, and restores `G` afterward; biology-time delivery reads the stored full field with `conc_global`, not `P+G` a second time. |
| **Ambient ROS mortality** | Fixed | ROS-driven SOS induction was priced as `k_ROS × [O₂]_ambient × mu_realized` — oxygen no cell acquired — and with the uncited `k_ROS = 1e2` it accounted for 99.5% of cumulative SOS hazard and for every population loss in the oxygen series (#332). The shipped `k_ROS` default is now `0.0`, retiring ambient ROS mortality by default. A funded driver prices induction from each agent's funded respiratory oxygen flux, in either a specific (kg/mol) or absolute (mol⁻¹, per-generation calibratable) normalization (#339), and the four SOS hazard components are emitted per step. The ambient path remains reachable by explicitly setting a positive `oxygen.k_ROS`; see `docs/DELIVERY_ROS_CAMPAIGN.md`. |
| **Delivery rationing accepted negative fields** | Fixed | Prescribed delivery mass was removed unconditionally and, after two retries, a negative solve was accepted, so the field borrowed against undelivered supply (min oxygen `-9.9e-6` at 500 agents) while remaining mass-exact — a positivity failure, not an accounting one (#338). Rationing now iterates to feasibility, local reductions first and a global bisection guarantee second, never accepts a negative solve, and reports a per-species rationing factor. |
| **Shipped ROS/uptake defaults after the delivery/ROS campaign** | Partially decided | `oxygen.k_ROS` now ships as `0.0`, so the ambient ROS mortality path is retired by default while remaining available through an explicit positive coefficient. `oxygen.ros_driver=ambient` and `metabolism.uptake_limit=none` remain shipped because the parser refuses `uptake_limit=delivery` with `gpu_enabled=true` until CUDA delivery parity exists. The remaining open decision is whether and how to enable funded uptake/ROS defaults after that parity work. |
| **Truncated-run indistinguishability** | Fixed | `/run_provenance/termination_cause_code` is pessimistically written as `incomplete_unknown` and overwritten only at an explicit run termination; stdout names the cause and distinguishes early exits. |
| GPU portability | Open | Production chemistry + mechanics on GPU; FMM M2L tree-walk on CPU; multi-GPU NCCL not wired |
| Large-scale MPI scaling | Partial | `mpi_four_rank` CTest (`mpirun -np 4`, includes periodic-x ring); manual `mpirun -np 8+` on HPC |
| **MPI hang np>2 (periodic x)** | Fixed | Ghost/migrate used sequential `Sendrecv(lo)` then `Sendrecv(hi)`, deadlocking the default periodic-x ring for `np>2`. Now non-blocking `Isend`/`Irecv`+`Waitall`. `np=2` used the collapsed path and was unaffected. |
| **Checkpoint resume mu_max wipe** | Fixed | Closed restarts wrote `mu_realized` but restore set `mu_max = mu_realized`. After Spot resume from a stressed snapshot (`mu≈5e-6`), growth could not recover and combinatorial washout killed the population in one step. Now persist/restore `/mu_max` (+ `/in_crypt`); legacy files fall back to strain `mu_max`. `init_from_checkpoint` also prints `GPU: ON` so `REQUIRE_GPU=1` no longer false-fails clean resumes. |
| **MPI ghost reaction/uptake double count** | Fixed | Ghost agents are marked explicitly and excluded from owner-committed reaction writes, uptake accounting, oxygen depletion, and siderophore biomass aggregation before rank reductions. |
| **Density-coupled VBF sink saw only rank-local agents** | Fixed | Spec 12 Change 1 builds a per-voxel occupancy histogram from owned agents. In `replicated` chemical decomposition (the default) every rank applies VBF over the whole global grid, so with `nprocs > 1` and `vbf.agent_carbon_coupling != 0` each rank felt only its own agents' carbon competition and the ranks' fields diverged from each other and from serial — silently, since agent reactions are already reduced before the VBF term. The histogram is now globally summed in replicated mode, and the device VBF path defers to the reduced host path in that configuration. Inert at the default coupling of `0.0`; it would have corrupted any multi-rank carbon-competition sweep. |
| **MPI ghost division duplication** | Fixed | Ghost agents are excluded from division, so only the owning rank creates daughters before migration. |
| **Uncounted SOS/phage lysis deaths** | Fixed | `mortality_lysis` records actual delayed SOS/phage deaths separately from induction counters and participates in global interval/cumulative event accounting and population closure. |
| **Population-ledger semantics** | Fixed | Starvation does not kill: `bacteriostasis_threshold` classifies viable, non-reproducing bacteriostatic cells and kills nothing. `bacteriostatic_live_agents` plus `washout_trapped_live_agents` are instantaneous stocks, never closure counters. `washout.trap=emergent` is the default. Persisted `ProvenanceCause` values remain fixed; value 4, formerly starvation, is retired. |
| **Founder placement scales with domain height** | Fixed | `initial_population.placement=z_slab` with explicit `z_min`/`z_max` keeps founders in a scenario-defined band; the default `legacy` policy preserves the historical band for compatibility. |
| **Exploding inertial mechanics** | Fixed | The previous update used `Δx = F·dt/(2m)` as a displacement, so a 0.1 µm overlap could move a cell approximately 1.7 m in one 60 s step. Contact could eject cells from the domain, making every spatial observable and boundary/washout attribution produced before the overdamped fix suspect. Mechanics now uses Stokes mobility with `vbf_viscosity` and caps displacement at `0.1·r`. |
| **Cross-boundary HGT to ghost recipient** | Open | Plasmid transfer from a local donor to a ghost recipient is lost when ghosts are cleared, under-counting cross-boundary conjugation. |
| **Dynamic mucin carbon liberation is dimensionally wrong** | Open | `dynamic_mucin_liberation()` returns `k_liberation [1/s] * vbf.density [cells/m^3] * Monod(M)`, i.e. `1e-4 * 1e11 = 1e7` mol/m^3/s at defaults — eleven orders above the static `vbf.mucin_liberation` term it replaces. Measured with only `mucin.enabled=true`: carbon reaches `3.4e7` mol/m^3 at 26 µm after one 60 s step against a `5e-3` physiological scale, then washes out over roughly ten steps. Mucin is also emptied outside the epithelial layer in one step while carbon receives the full unlimited amount, so the two species are not stoichiometrically coupled, and `mucin.secretion_rate` is inert because liberation is Monod-saturated (`M=1e-2` vs `Km_degradation=1e-3`). `mucin.enabled` is `false` by default, so no shipped result is affected; a mucin-driven carbon arm must not be run until `k_liberation` is redefined as a per-organism rate (mol/(cell·s)) with mucin-limited release. See `docs/CARBON_LADDER_CAMPAIGN.md`. |
| **Carbon reference amplitude has no config key** | Fixed | `carbon.boundary_conc` now sets both the epithelial boundary and the prescribed z-gradient amplitude by default, closing the inversion measured at `carbon.boundary_conc=1e-3`: C(2 µm)=1.00e-3 below C(26 µm)=1.75e-3. Use `carbon.z_amplitude` or the legacy `carbon_z_amplitude` spelling to explicitly override the amplitude while retaining the boundary value; see `docs/CARBON_LADDER_CAMPAIGN.md`. |
| **Shipped per-generation lysis is ~2.8%, not the documented 1%** | Open | `sos_lysis_prob = 0.01` is documented as "1% per division", but `fix_metabolism.cpp` marks `just_divided` on *both* mother and daughter, so each division draws the hazard twice and a cell is marked twice per lifetime (born, then divides) — ≈2% per generation. The separate `sos_basal_rate = 1e-6/s` adds `1 - exp(-1e-6·T_gen)` = 0.82% at the measured uncrowded `T_gen ≈ 2.3 h` and 1.2% at the 3.3 h measured at 80 founders. Total ≈2.8% per generation for a producer with `oxygen.k_ROS = 0.0` and no funded uptake, i.e. in the shipped configuration; nuclease cross-induction adds a density-dependent term on top. Consequence for calibration: the `k_ROS_funded = 6.2e11 / 1.24e12 / 3.1e12 mol⁻¹` coefficients are *increments* on that base, so an arm labelled "1% per generation" actually lyses at ≈3.8%. Report realized lysis as `mortality_lysis / divisions`, never as configured. See `docs/SPEC13_LYSIS_SELECTION.md`. |
| **Colicin retardation and burst size are hardcoded per plasmid** | Open | `plasmid.cpp` sets ColE1/E2 `retardation = 50`, copied into the burst by `fix_bacteriocin.cpp` and consumed as `D_eff = diff_coeff/retardation = 8e-13 m²/s` — 72× below the literature analogue `7e-11/1.2`, i.e. ~8.5× kill-zone radius and ~72× area. The `retardation_basic/neutral/acidic` config keys have no live consumer (`retardation_for_pI` is dead code), so sweeping them produces identical runs, and `burst_size` has no key at all. These are larger levers on colicin efficacy than the lysis prior, so any lysis coefficient calibrated against a displacement time series at the current values absorbs the transport error. See `experiments/rps_campaign/rps_spec_audit.md` §B and `docs/SPEC13_LYSIS_SELECTION.md`. |
| **No configurable resistant strain** | Open | `SimulationConfig::InitialStrain` has no receptor-genotype field and `parse_strain_object` ignores unknown keys, so `receptor_mutations` parses silently and does nothing; every founder gets `receptor_expr_base.fill(1.0)` and resistance arises only from `FixMutation` at 1e-7/step. The model can express colicin producer vs sensitive, but not the three-strain non-transitive (RPS) system, so it cannot yet be compared to Kirkup & Riley 2004 as an RPS result. When the genotype field lands, note that a knockout already costs growth through `Km_b12 = km_b12/(expr·affinity)` (9% at 1 µM corrinoid, 25% at 0.3 µM), so an additional `mu_max` cost double-charges resistance. See `experiments/rps_campaign/rps_spec_audit.md` §C. |

When writing tests that involve plasmids, use **`ColE1`/`ColB`** (legacy `colicin_E1` aliases still resolve) and assert `agent.genome.bi_loci.size() > 0`.

## Test Coverage Map

### C++ (CTest)

**Fast unit (`ctest -L unit -LE slow`):** spatial hash, Green's functions, agent/plasmid, iron fallback, octree, FMM, conjugation, z-gradient, nutrient diffusion, domain decomp, acetate/MetE, protease decay, oxygen gradient, O₂ growth boost, mucin liberation, peristaltic advection, ethanolamine, adaptive dt, agent transfer pack/unpack, fix registry, input parser, config ingestion (every parser key is tracked into `SimulationConfig`), qssa stoichiometry, bacteriocin, receptor, mutation, uptake limitation.

**Slow unit:** mechanics, immunity escape.

**Integration (`ctest -L integration`):** smoke (end-to-end biology), config diversity (fixture/example fingerprints must differ), mechanism wiring (cross-spec mass-balance + directional coupling; see `docs/WIRING_AUDIT.md`), population-scale delivery resolution (`uptake_limit_population_resolution`), HDF5 round-trip, HDF5 checkpoint restart, washout trap long-horizon regression (`washout_trap`, issue #160).

**MPI:** `mpi_multi_rank` (`mpirun -np 2`), `mpi_four_rank` (`mpirun -np 4`), `mpi_gpu_multi_rank` (2-rank + `gpu_enabled`), `cuda_aware_mpi_reaction` (CUDA-aware device Allreduce when available), `hdf5_roundtrip_parallel`.

**OpenMP:** `openmp_parity` + `scripts/compare_openmp_parity.sh` (serial vs OpenMP deterministic + stochastic toxin-kill fingerprints)

**GPU (CUDA job):** `greens_function_gpu`, `gpu_diffusion`,
`gpu_chemical_field`, `gpu_feature_combinations`, `gpu_production_path`,
`gpu_smoke`, `gpu_reproducibility`, `gpu_scaling_benchmark`,
`qssa_gpu_parity`, `gpu_toxin_burst_parity`, `gpu_nutrient_feedback`,
`gpu_metabolism_fur`,
`spatial_hash_gpu_csr`, `mechanics_gpu_parity`, `gpu_kernel_units`,
`mpi_gpu_multi_rank`, and `cuda_aware_mpi_reaction`, plus
`scripts/compare_gpu_parity.sh`. GPU-labelled CTest targets return CTest's skip
code `77` when CUDA is not compiled or no physical device is available, so
CTest reports them as `Skipped`, not `Passed`. Setting `REQUIRE_GPU=1` converts
either condition into a named test failure; this variable is shared with
`deploy/aws/entry.sh`, where it means the run must actually report GPU use.
CUDA compilation alone is not equivalent to physical GPU execution. New
`gpu`-labelled CTest targets must also be added to the explicit target list in
`.github/workflows/ci.yml`; otherwise the CUDA job reports them as `Not Run`.

**Benchmark:** `scaling_benchmark` (issue #55 smoke counts). The larger
scaling sweeps are on-demand only via
`scripts/run_scaling_benchmark.sh` and
`scripts/run_gpu_scaling_benchmark.sh`; the nightly workflow was removed after
OOM failures at the `1000x1000x50` grid in the two-rank case. Chemical-field
storage is confirmed by code and a local OOM measurement to use the global
grid per rank, requiring approximately 16 GiB for two ranks before other
allocations. This memory behavior is not fixed.

**Config diversity guardrail:** `test_config_diversity` runs short simulations from parser fixtures and example JSON files and asserts distinct deterministic fingerprints — catches configs silently reverting to defaults.

### Python (pytest in CI)

`python-lint` job runs `ruff check`, import smoke, and `pytest tests/ -m "not integration"`. Coverage includes HDF5 reader, analysis helpers, validation regression helpers, and FISH observation models (#25).

### CI jobs (`.github/workflows/ci.yml`)

| Job | What it exercises |
|-----|-------------------|
| `unit-tests` | Fast CTest unit shard |
| `integration-tests` | Smoke, config diversity, HDF5, batch runner smoke, Python integration pytest (#155) |
| `gpu-parity` | CPU vs GPU fingerprint check via `compare_gpu_parity.sh` (#158) |
| `openmp-parity` | Serial vs OpenMP deterministic + stochastic toxin-kill fingerprints |
| `cuda-compile` | CUDA compile + GPU test targets (single arch, no duplicate parity rebuild) |
| `gpu-device-tests` | Physical T4 GPU verification-of-record via AWS Batch (`REQUIRE_GPU=1 ctest -L gpu`) |
| `eari-vadi-validation` | Short EARI/VADI + FISH invariant/bounds checks (#56, #25) |
| `python-lint` | JSON syntax, ruff, pytest (fast), batch runner dry-run |

### Gaps (add tests when touching these areas)

- Multi-rank MPI beyond 4 processes (`mpirun -np 8+` on HPC)
- Multi-GPU NCCL

## Adding Features — Agent Checklist

### New Fix module
1. `src/fixes/fix_<name>.{h,cpp}` inheriting `Fix`
2. Config struct in `input_parser.h` + defaults in `default_config()`
3. Register in `FixRegistry::register_fix()` inside `fix_registry.cpp` (or call from a static initializer)
4. `tests/test_<name>.cpp` + entry in `tests/CMakeLists.txt`
5. Update `docs/MECHANISMS.md` if biological behavior changes

### New diffusion / QSSA kernel
1. Bacteriocins: add QSSA-compatible Green's-function logic in `src/diffusion/` and wire through `QSSASolver`
2. Nutrients/small molecules: add stable implicit field logic in `src/fields/chemical_field.cpp`; never add an explicit biological-timestep stencil
3. Preserve periodic x/y and luminal zero-flux z conditions. Epithelial z=0
   defaults to Dirichlet, with runtime-selectable Robin/flux delivery for
   configured species; keep each mode's accounting and solver form consistent.
4. Add analytic/invariant, coefficient/enable sensitivity, positivity, MPI, and GPU-parity tests

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
- **Add an ingestion probe in `tests/test_config_ingestion.cpp`** (`build_probes()` for flat keys, `array_and_strain_keys()` for `config_json.cpp` array/object keys). This test tracks *every* parser key: its completeness guard scans the parser sources for `key == "..."` literals and fails CI if any parsed key lacks a probe (or a probe references a key no longer parsed).
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
| HDF5 schedule | Nested `hdf5.schedule.*` in input JSON or `cfg.hdf5.schedule` (summary/agents/grid/lineage/genome/provenance intervals) |
| Checkpoint restart | `checkpoint_file` + optional `checkpoint_step` in input JSON |
| Closed midstream restarts | `restart.enabled` + `restart.directory` + `restart.interval_steps` (Tier 2: agents+grid; AWS uploads immutable `step_*.h5` + `latest.json`) |
| Disable HDF5 in tests | `cfg.hdf5.enabled = false` |
| Strict config | `GUTIBM_STRICT_CONFIG=1` aborts on invalid numerics |
| MPI decomp axis | `cfg.domain.mpi_decomp_axis` (default 0 = x) |
| Barnes-Hut FMM | `use_fmm`, `fmm_theta`, `fmm_expansion_order` in input JSON |
| Peristaltic mixing | `peristaltic_*` keys in input JSON |
| Chemical environment (Spec 1) | `oxygen.enabled`, `acetate.enabled`, `mucin.enabled`, `protease.enabled` + nested keys in `docs/PARAMETERS.md` |
| Cell biology (Spec 3 / 10v2) | `fur.enabled`, `cdi.enabled`, `motility.enabled` + aerotaxis/energy/surface/mucin keys; per-strain `cdi_type`, `cdi_immunity` |
| Quorum sensing (Spec 11) | `quorum_sensing.enabled` + AI-2 production/import/chemotaxis keys in `docs/PARAMETERS.md` |
| GPU | `gpu_enabled` in input JSON (CUDA build required) |
| Immigration | `immigration.enabled` plus `immigration.count`, `strain_index`, placement and schedule keys; pulse `step` is relative to run start |
| Initial founder placement | `initial_population.placement=z_slab` plus `initial_population.z_min`/`z_max`; scenario-wide founder band, independent of `domain_z` when configured |

Full parameter docs: `docs/PARAMETERS.md`.

## Spec Documents

- `EARI.md` — Eco-Advective Receptor Interference framework
- `VADI.md` — Viscous Advective-Diffusion Interference framework
- `docs/MECHANISMS.md` — per-Fix biological mechanisms
- `docs/WIRING_AUDIT.md` — cross-spec mechanism wiring map (species mass balance, coupling points, open questions) + wiring test strategy
- `docs/API.md` — class reference
- `docs/CONFIG_FORMAT.md` — strict JSON input format
- `docs/BATCH_RUNNER.md` — resumable parameter-scan CLI
- `docs/PARAMETERS.md` — configurable parameters
- `docs/SCALING.md` — agent-count benchmarks and profiling
- `docs/AWS_BATCH.md` — AWS Batch Spot + CUDA deployment + Phase Observability
- `docs/MULTI_SCALE_EXPERIMENTATION.md` — nested short/small vs host-scale experiment ladder
- `docs/BRANCHING_FROM_CHECKPOINTS.md` — measured calibration checkpoint fork runbook
- `docs/PRE_SUBMISSION_CHECKLIST.md` — operator checks before campaign submission

## Code Conventions

- C++20 with modern idioms (smart pointers, RAII, move semantics)
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

## Open Issue Tracker (Jul 2026)

Jun 2026 queue closed (#40–#81, #25, #29, #33, #55). Post–GPU ROI backlog **#154–#161** complete via #162–#167.

| Issue | Topic | Status |
|-------|-------|--------|
| #154 | MPI four-rank validation | Done (#163) |
| #155 | Python integration pytest in CI | Done (#162) |
| #156 | CUDA-aware MPI reaction reduce | Done (#163) |
| #157 | GPU mechanics force kernel | Done (#164) |
| #158 | GPU parity CI | Done (#162) |
| #159 | Sub-quadratic FMM M2L | Done (#164) |
| #160 | Metabolic washout trap regression | Done (`washout_trap` CTest) |
| #161 | OpenMP stochastic parity | Done (`FINGERPRINT_STOCHASTIC` in openmp parity CI) |

Remaining long-horizon: `mpirun -np 8+` on HPC, GPU FMM octree traversal, Fur on device.

**Project board:** [docs/PROJECT_BOARD.md](docs/PROJECT_BOARD.md) — kanban layout, PR bundles, merge order. Run `./scripts/setup_project_board.sh` to create GitHub labels, milestones, and a Projects v2 board.

## Cursor Cloud specific instructions

System deps (`cmake`, `build-essential`, MPI, parallel HDF5) and Python deps are
provisioned by the startup update script; the notes below are non-obvious caveats
for building/running/testing in this environment.

- **Build with GCC, not the default `c++`.** The default `cc`/`c++` on this VM is
  Clang 18, which selects the gcc-14 toolchain but there is no `libstdc++-14-dev`,
  so any Clang link fails with `cannot find -lstdc++`. Always configure with
  `CC=gcc CXX=g++` (gcc-13, which has full dev libs and matches CI):
  `CC=gcc CXX=g++ cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON`. Standard
  build/test/run commands are in `.agents/skills/gut-ibm/SKILL.md`.
- **All 41 CTest targets pass here.** The `libhwloc`/MPI "preexisting failures"
  table in `.agents/skills/testing-gutibm/SKILL.md` describes a *different* (Devin)
  VM; on Cursor Cloud MPI/HDF5 init works and the full suite (incl. `mpi_multi_rank`,
  `hdf5_*`, `scaling_benchmark`) passes.
- **`ruff`/`pytest` install to `~/.local/bin`**, which is not on `PATH` by default —
  prefix with `PATH="$HOME/.local/bin:$PATH"` or add it, or call `python3 -m pytest`.
- **Running the simulation:** binary is `build/gut_ibm <config.json>`. For
  `mpirun -np N` with N greater than the available cores, add `--oversubscribe`.
  The full `examples/single_colony/input.json` (86400 s / 1440 steps) is compute-
  heavy and writes a very large HDF5; use a short/small config for quick checks.
