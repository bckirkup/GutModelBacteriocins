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
cd python && pip install -e ".[viz,dev]"   # viz (matplotlib) + dev (ruff, pytest)
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
| `src/layer2/` | Spec 13 Layer-2 `MucusSegment`, patch table, live-patch re-entry, perturbation pulses |
| `src/layer3/` | Spec 13 Layer-3 `ColonicChain` and `LumenCompartment` multi-segment colon dynamics |
| `src/io/input_parser.cpp` | JSON + legacy flat-key config parser |
| `src/io/hdf5_writer.cpp` | Parallel HDF5 output + genome/RNG checkpoint groups |
| `src/io/hdf5_reader.cpp` | Checkpoint restart snapshots + RNG replay |
| `docs/DELIVERY_ROS_CAMPAIGN.md` | Consolidated delivery/ROS campaign findings and retired claims |
| `docs/EXTERNAL_AUDIT_2026-08.md` | External audit claims, per-claim verdicts, and the decisions taken |
| `docs/GPU_COST_BENEFIT_2026-08.md` | Measured T4 cost/benefit per uptake mode and scale, with raw records under `bench_results/` |
| `docs/GPU_PRECISION_CAMPAIGN_2026-09.md` | Host/device outcome interchangeability (p1 precision + 1440-step agent-dense science arm) and speedup vs agent load, records under `bench_results/gpu_precision_2026-09/` |
| `docs/SPEC13_MULTISCALE.md` | Layer-2/Layer-3 multiscale architecture, live physics, and occupancy readouts |
| `docs/SPEC14_PROPHAGE_MODULE.md` | Spec 14 temperate prophage carriage, Hill induction, virion PDE, and presets |
| `docs/QSSA_TRANSIENT_BOUNDS.md` | QSSA vs transient PDE error bounds (`qssa_transient_bounds` CTest) |
| `python/gut_ibm_tools/` | HDF5 reader, analysis, validation, visualization, pellet/FISH observation, checkpoint analyze CLI |
| `python/gut_ibm_tools/colony.py` and `spatial_stats.py` | Colony catalogs and 3-D spatial observables |
| `examples/` | `single_colony/`, `diversity_paradox/`, `eari_vadi_validation/`, `cell_biology/`, `oxygen_boundary_presets/`, `prophage_presets/`, `immigration_fork/`, `fork_profile/`, `batch_scan/`, `gpu_cost_benefit_benchmark/`, `scaling_benchmark/` |
| `tests/` | CTest targets (see test map below) |
| `.agents/skills/gut-ibm/SKILL.md` | Hands-on development reference |
| `.agents/skills/gutibm-campaign-ops/SKILL.md` | AWS campaign operations and checkpoint safety |
| `.agents/skills/sonarqube-gutibm/SKILL.md` | SonarQube remediation workflow |
| `.agents/skills/sonarqube-cpp/SKILL.md` | C++ SonarQube patterns |
| `.agents/skills/sonarqube-python/SKILL.md` | Python SonarQube patterns |
| `docs/AWS_BATCH.md` | AWS Batch Spot + CUDA deployment plan + Phase Observability |
| `deploy/aws/` | Dockerfile / Batch entrypoint / status heartbeat and AWS ops scripts |
| `docs/LANDMINES_ARCHIVE.md` | Resolved landmines (same table format, `Fixed in` PR column) |

## Key Concepts

- **Advective Double-Bind (EARI)**: Growth/resistance trade-off leading to washout
- **Combinatorial Washout Trap (VADI)**: Physical expulsion when μ_realized < γ_flow
- **Lethal Core/Halo**: Spatial toxin zones based on isoelectric point (pI)
- **Comet-tails**: Downstream-elongated inhibitory zones from mucus flow
- **BI Locus**: Bacteriocin-Immunity genetic region subject to recombination
- **Method of Images**: Boundary condition technique for diffusion kernels
- **MPS (Mating-Pair Stabilization)**: Shear-dependent HGT probability

## Known Bugs & Landmines

Resolved landmines live in `docs/LANDMINES_ARCHIVE.md` (same table format plus a `Fixed in` column); only open and partial rows stay here.

| Issue | Status | Notes |
|-------|--------|-------|
| **Shipped ROS/uptake defaults after the delivery/ROS campaign** | Partially decided | `oxygen.k_ROS` now ships as `0.0`, so the ambient ROS mortality path is retired by default while remaining available through an explicit positive coefficient. `oxygen.ros_driver=ambient` and `metabolism.uptake_limit=none` remain shipped; explicit GPU delivery is permitted via device delivery (`device_delivery`) or host-forced delivery (`host_forced_delivery`) recorded in `/run_provenance/chemistry_placement`, and agent-side delivery prepare/commit + funded respiration now run on CUDA (PR #490, `docs/CUDA_DELIVERY_PARITY.md`). The remaining open decision is whether and how to enable funded uptake/ROS defaults after benchmark evidence. |
| GPU portability | Open, cost/benefit measured | Production chemistry + mechanics + delivery metabolism + prophage hazard/adsorption + device-resident FMM far-field (both local-expansion and direct GPU octree multipole traversal across replicated and slab layouts, #159) on GPU; CDI corpse participation is parity-preserving through a shared host/device predicate; multi-GPU NCCL not wired. Device coverage is evidenced only by the T4 CI job. Post-fix T4 at 500x500x50 with OpenMP: device chemistry is 1.91x faster than host for `uptake_limit=none` and 2.86x end-to-end at 1e5 agents — the earlier 1.29x 'regression' was host-device marshalling, since removed. Device-resident stepping (#408/#409/#480, PR #493) keeps `d_conc_` and `d_reac_` on device across consecutive GPU steps (eliminating full-grid `step_conc_upload`, `chem_diffusion_result`, `qssa_toxin_host_sync`, and `chem_reactions` transfers between output dumps) and routes unavoidable staging through RAII `PinnedHostBuffer<T>` (`cudaMallocHost` with pageable fallback), reporting `h2d_bytes`, `d2h_bytes`, and `pinned_staging_active` in `GpuProfile` and `/run_provenance/`. Full-matrix T4 wall-clock remeasurement awaits the next `g4dn.xlarge` Batch pass. See `docs/GPU_COST_BENEFIT_2026-08.md`. |
| Large-scale MPI scaling | Partial | `mpi_four_rank` CTest (`mpirun -np 4`, includes periodic-x ring); manual `mpirun -np 8+` on HPC |
| Image series with wall-normal flow | Mitigated (`auto` + hybrid CPU/GPU kernel, PR #495); legacy default retained | Translations + z-reversed reflections satisfy zero diffusive flux at both walls rather than physical zero total flux `-D dC/dz + U_z C = 0` when `U_z != 0`; exact only for wall-parallel flow. The exact three-branch hybrid Sturm-Liouville / Hankel kernel (`robin::evaluate_physical_sealed_concentration`) now runs natively on both CPU and CUDA GPU (`rel L2 <= 2.5e-12` vs 500-mode SL and `<= 1.3e-3` vs axisymmetric FVM across ColE1 `Pe_z = 0.822–2.325` and ColE2; zero sealed GPU host fallback), selectable globally via `qssa.drift_correction=true` or automatically per high-`Pe_z` source via `qssa.drift_envelope_policy=auto` above `qssa.peclet_z_threshold` (default `0.05`). Shipped defaults retain `qssa.drift_envelope_policy=warn` and `qssa.drift_correction=false` for historical bit-reproducibility (so the warning still fires at shipped ColE1 `Pe_z=0.822/2.325`), and killing/induction outcome movement between legacy and `auto` modes remains to be measured in a campaign. Device parity is evidenced by the T4 CI job (`greens_function_gpu`, `qssa_gpu_parity`). See `docs/EXTERNAL_AUDIT_2026-08.md` and `docs/NEUMANN_WALL_NORMAL_DRIFT.md`. |
| **Colicin retardation and burst size are hardcoded per plasmid** | Open, partly corrected | `plasmid.cpp` derives retardation from pI (`retardation_from_pI`): 50.2 for ColE1 (pI 9.0, basic, mucin-bound) but **2.04** for ColE2 (complex pI 6.5, near-neutral) — the nuclease colicin is nearly unretarded (`D_eff = 1.7e-11 m²/s`), so transport is *not* what limits it; ColE1 remains 72x below the literature analogue `7e-11/1.2`. The dead global `retardation_basic/neutral/acidic` keys and `retardation_for_pI` were removed in `0fdf590` (PR #300; recorded in `docs/LANDMINES_ARCHIVE.md`), and per-plasmid `retardation`, `diff_coeff`, and `burst_size` are configurable via `plasmid_overrides.<plasmid>.<field>` — what stays open is whether the *shipped* values are right, since a lysis coefficient calibrated at current retardation absorbs the transport error. See `experiments/rps_campaign/rps_spec_audit.md` §B and `docs/SPEC13_LYSIS_SELECTION.md`. |
| **~~No configurable resistant strain~~ — withdrawn; the cost is the real gap** | Open (cost unmeasured) | The withdrawn claim was wrong: `initial_strains[].receptor_expression` (aliases `receptor_genotype`, `receptors`) is parsed, range-validated, applied to `receptor_expr`, inherited, MPI-shipped, and checkpointed — a `{"BtuB": 0.0}` founder is a complete *btuB* R-strain phenocopy, so C/S/R is constructible today. The open gap is the **cost**: `Km_b12 = km_b12/(expr·affinity)` clamps at a 100x Km inflation — a ~9% growth penalty at shipped `b12_initial_conc=1e-3` — still unmeasured at population scale (the probe that was to demonstrate it is the next row). Classify strains by `receptor_expr[BtuB]` or agent `type`, never `state`. See `docs/SPEC13_LYSIS_SELECTION.md` §4a. |
| **Strain composition is not observable at patch scale** | Open | The 12-arm `experiments/rps_probe/` run at `3c166c1` found a 200x200x100 µm patch cannot report *which* strain wins, only how many cells there are: founders collapse to ~30 live against hundreds of boundary exports, then whichever few lineages survived the crash sweep it — identical-treatment seeds gave R-only, S-only, and C+S finals; one strain was extinct in 8 of 12 arms. RPS/displacement contrasts need Layer 2 reseeding or Layer 3, >=5 seeds with paired within-seed statistics, never between-arm mean ratios; the producer treatment was nearly inert in the same run (11 receptor-mediated kills vs 7629 divisions). See `docs/SPEC13_LYSIS_SELECTION.md` §4a-result and `docs/CARBON_LADDER_CAMPAIGN.md`. |
| **Colicin E potency is set by the corrinoid pool, and it ships where colicins do not matter** | Open (decoupling toggle added in #488) | `FixReceptor::toxin_occupancy` is competitive (`apparent_kd = kd_colicinE_btuB * (1 + [corrinoid]/kd_b12_btuB)`); shipped `b12_initial_conc = 1e-3` leaves BtuB corrinoid-occupied 1000-fold, so colicin E faces an apparent Kd of `5e-4` rather than `5e-7`. A five-level corrinoid ladder (`experiments/spec14_sec8/`, 50 arms, paired nulls, five seeds) swings kills per producer lysis 61x (0.121 to 7.38); PR #488 (#464) added `metabolism.b12_growth_essential` (default `true`) so experiments can decouple BtuB receptor shielding from obligate B12 Monod growth arrest, plus carriage-matched ColE2/ColIa/ColB controls. Efficacy results at the shipped pool still describe a near-inert toxin; whether ~1 µM free mucus corrinoid is right remains a parameter decision. See `docs/SPEC14_SECTION8_VALIDATION.md`. |
| **A ratio read off the final output step measures the division wave, not selection** | Open (analysis practice; canonical estimator added in #484) | Dense cultures divide in alternating synchronized waves, so a last-sample `log10(n1/n2)` oscillates ±0.2 decades with the phase — at `b12=1e-5` the final sample reported +0.16 decades against a +0.065 wave-averaged value (2.5x overstatement). Report a tail median and a fitted slope over multiple waves using the canonical estimators in `gut_ibm_tools.analysis` (PR #484, #481) and `experiments/spec14_sec8/analyze_corrinoid.py`; this is the second final-sample statistic after `rps_probe` to produce a wrong campaign number. |
| **Spec 14's motivating in-vitro exclusion does not occur in GutIBM** | Open (campaign re-evaluation; `FixProphage` landed default-off in #494) | The §8 no-code validation (`docs/SPEC14_SECTION8_VALIDATION.md`, run at `5fb4250`) found **V1 FAIL at every corrinoid level**: the Spec 14 target is 100x exclusion in 6 h (0.33 decades/h) and the ladder's best was 0.090 decades/h — parametric (corrinoid competition) and structural (a lysis forfeits a whole lineage while a killed sensitive is largely replaced; break-even sits between 1.17 and 7.38 kills per lysis). PR #488 (#464) landed `b12_growth_essential` decoupling and carriage-matched controls, and PR #494 (#475) landed the default-off temperate prophage module (`prophage.enabled=false`, `docs/SPEC14_PROPHAGE_MODULE.md`); empirical §8 V1 re-evaluation and *in vivo* lysogen-displacement testing remain open campaign work. |
| **A sequence census is an upper bound on a functionally-defined parameter, not its central value** | Open (parameterization practice) | When a functionally-defined parameter (inducible, viable, active) is sourced from a sequence census, the census value is an upper bound — three instances: mucosal density (16S qPCR 100-1000x CFU; `docs/gutibm_density_reconciliation.md`), lysogen prevalence (`lysogen_prevalence_init` counts *inducible* prophage but surveys count integrated sequences), induction rate. Record the census value as a bound and sweep the functional fraction; prevalence and per-lysogen rate are not independent axes since the observable is their product. See `docs/SPEC14_PRIOR_REVIEW.md`. |

`uptake_limit=voxel` is retained as the grid-artifact probe behind the Spec 12 findings (exercised by the `tests/test_uptake_limit.cpp` resolution probes; rationale in `docs/SPEC12_DENSITY_LIMITATION.md`), not as a production uptake mode.

When writing tests that involve plasmids, use **`ColE1`/`ColB`** (legacy `colicin_E1` aliases still resolve) and assert `agent.genome.bi_loci.size() > 0`.

## Test Coverage Map

### C++ (CTest)

**Fast unit (`ctest -L unit -LE slow`):** 118 targets covering spatial hash, Green's functions, Neumann image series & wall-normal drift (13 shards + `wall_normal_peclet`), Robin lumen BC (18 shards), agent/plasmid, iron fallback, octree, FMM (11 fast shards + `benchmark_100k_sources`), qssa stoichiometry, QSSA transient bounds (`qssa_transient_bounds`), metabolism, chemical field memory guard (`chemical_field_memory_guard`), conjugation, z-gradient, nutrient diffusion, domain decomp, acetate/MetE, protease decay, oxygen gradient, O₂ growth boost, mucin liberation, units & stoichiometry (`vbf_mucin_stoichiometry`), peristaltic advection, ethanolamine, adaptive dt, agent transfer pack/unpack, fix registry, input parser, config ingestion (every parser key is tracked into `SimulationConfig`), bacteriocin, receptor, temperate prophage (`prophage`), mutation, fur, cdi, motility, quorum sensing, immigration, uptake limitation, delivery support cache concurrency, mechanics participation, HDF5 fail-closed, Spec 13 Layer-2 (`layer2_segment`, `layer2_reentry` 3 shards) and Layer-3 (`layer3_chain`), and AWS/ops helper scripts.

**Slow unit:** mechanics, immunity escape.

**Integration (`ctest -L integration`):** 34 targets covering smoke (end-to-end biology), config diversity (fixture/example fingerprints must differ), feature combinations, mechanism wiring (cross-spec mass-balance + directional coupling; see `docs/WIRING_AUDIT.md`), toxin sentinel, temperate prophage (`prophage`), population-scale delivery resolution (`uptake_limit_population_resolution`), maintenance carbon funding, HDF5 round-trip/checkpoint/restart/RNG-replay (`hdf5_rng_replay`)/schedule/compression/live-agents/graceful-shutdown/summary-events, MPI, GPU integration, and washout trap long-horizon regression (`washout_trap`, issue #160).

**MPI:** `mpi_multi_rank` (`mpirun -np 2`), `vbf_mucin_stoichiometry_mpi` (`mpirun -np 2`), `mpi_four_rank` (`mpirun -np 4`), `mpi_gpu_multi_rank` (2-rank + `gpu_enabled`), `cuda_aware_mpi_reaction` (CUDA-aware device Allreduce when available), `hdf5_roundtrip_parallel`, `hdf5_checkpoint_parallel`, `vbf_density_replicated_mpi`, `vbf_density_replicated_serial`.

**OpenMP:** `openmp_parity` + `delivery_support_cache_race` + `scripts/compare_openmp_parity.sh` (serial vs OpenMP deterministic + stochastic toxin-kill fingerprints)

**GPU (CUDA job — 25 `gpu`-labelled targets):** `greens_function_gpu`, `gpu_diffusion`,
`gpu_chemical_field`, `gpu_residency_audit`, `gpu_reaction_residency`,
`gpu_feature_combinations`, `gpu_production_path`,
`gpu_smoke`, `gpu_reproducibility`, `gpu_scaling_benchmark`,
`qssa_gpu_parity`, `gpu_fmm_parity`, `neumann_image_series_gpu`, `robin_lumen_gpu`,
`gpu_toxin_burst_parity`, `gpu_toxin_source_clearing`, `gpu_nutrient_feedback`,
`gpu_metabolism_fur`, `gpu_uptake_limit`, `gpu_prophage_parity`,
`spatial_hash_gpu_csr`, `mechanics_gpu_parity`, `gpu_kernel_units`,
`mpi_gpu_multi_rank`, and `cuda_aware_mpi_reaction`, plus
`scripts/compare_gpu_parity.sh`. GPU-labelled CTest targets return CTest's skip
code `77` when CUDA is not compiled or no physical device is available, so
CTest reports them as `Skipped`, not `Passed`. Setting `REQUIRE_GPU=1` converts
either condition into a named test failure; this variable is shared with
`deploy/aws/entry.sh`, where it means the run must actually report GPU use.
CUDA compilation alone is not equivalent to physical GPU execution. New
`gpu`-labelled CTest targets must also be added to the explicit target list in
`.github/workflows/ci.yml`; otherwise the CUDA job reports them as `Not Run`
(enforced by `gpu_target_manifest` CTest).

**Benchmark:** `scaling_benchmark` (issue #55 smoke counts). The larger
scaling sweeps are on-demand only via
`scripts/run_scaling_benchmark.sh` and
`scripts/run_gpu_scaling_benchmark.sh`; the nightly workflow was removed after
OOM failures at the `1000x1000x50` grid in the two-rank case. In `replicated`
mode, `ChemicalField` allocates the global grid per rank (`~8.20 GiB / rank`,
`~16.39 GiB / 2 ranks`, `~65.57 GiB / 8 ranks` at `1000x1000x50` across 11
species after #480 lazy delivery-vector allocation). `Simulation::init()` now
runs `ChemicalField::evaluate_memory_guard()` before allocating, enforcing
`memory.max_host_chemical_gb` (or 90% of `/proc/meminfo` `MemAvailable`) across
`ranks_per_node` with `memory.replicated_oom_policy` (`"error"`, `"auto_slab"`,
or `"warn"`) and emitting a rank-0 diagnostic recommending
`"chemistry_decomposition": "slab"` (PR #493).

**Config diversity guardrail:** `test_config_diversity` runs short simulations from parser fixtures and example JSON files and asserts distinct deterministic fingerprints — catches configs silently reverting to defaults.

### Python (pytest in CI)

`python-lint` job runs `ruff check`, import smoke, and `pytest tests/ -m "not integration"`. Coverage (280 test functions across `python/tests/`) includes HDF5 reader, analysis & division-wave estimators, validation regression helpers, FISH and fecal-pellet observation models (#25, #474), colony/spatial stats, transport metrics, batch runner, and AWS Batch helpers.

### CI jobs (`.github/workflows/ci.yml` + `gpu-device-tests.yml`)

| Job | What it exercises |
|-----|-------------------|
| `unit-tests` | Fast CTest unit shard (`ctest -L unit -LE slow`) |
| `serial-build` | MPI/HDF5-off build with `-Wall -Wextra` |
| `sanitizers` | ASan+UBSan unit tests (serial HDF5, excludes `species_subset` OOM) |
| `tsan` | ThreadSanitizer unit tests (serial HDF5, excludes `species_subset` OOM) |
| `integration-tests` | Smoke, config diversity, HDF5, MPI, slow unit, batch runner smoke, Python integration pytest (#155) |
| `gpu-parity` | CPU vs GPU fingerprint check via `compare_gpu_parity.sh` (#158) |
| `openmp-parity` | Serial vs OpenMP deterministic + stochastic toxin-kill fingerprints |
| `cuda-compile` | CUDA compile + GPU test targets (single arch, no duplicate parity rebuild) |
| `gpu-device-tests` | Physical T4 GPU verification-of-record via AWS Batch (`REQUIRE_GPU=1 ctest -L gpu`, in `.github/workflows/gpu-device-tests.yml`) |
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
| Plasmid names | `ColE1`, `ColE2`, `ColB`, `ColIa`, `ColM`, `MccV`; per-plasmid overrides via `plasmid_overrides.<plasmid>.<field>` |
| HDF5 schedule | Nested `hdf5.schedule.*` in input JSON or `cfg.hdf5.schedule` (summary/agents/grid/lineage/genome/provenance intervals) |
| Checkpoint restart | `checkpoint_file` + optional `checkpoint_step` in input JSON (bitwise RNG replay on same-seed/same-rank continuation; reseeded on seed forks) |
| Closed midstream restarts | `restart.enabled` + `restart.directory` + `restart.interval_steps` (Tier 2: agents+grid+RNG; AWS uploads immutable `step_*.h5` + `latest.json`) |
| Disable HDF5 in tests | `cfg.hdf5.enabled = false` |
| Strict config | Unset/`GUTIBM_STRICT_CONFIG=1` aborts on invalid numerics and unknown keys; `GUTIBM_STRICT_CONFIG=0` is the explicit lenient escape hatch |
| MPI decomp axis | `cfg.domain.mpi_decomp_axis` (default 0 = x) |
| Barnes-Hut FMM | `use_fmm`, `fmm_theta`, `fmm_expansion_order` in input JSON |
| Robin lumen transfer | `toxin.lumen_transfer_length` in input JSON |
| Peristaltic mixing | `peristaltic_*` keys in input JSON |
| Chemical environment (Spec 1) | `oxygen.enabled`, `oxygen.epithelial_boundary` (`dirichlet`/`robin`/`flux`, presets in `examples/oxygen_boundary_presets/`), `acetate.enabled`, `mucin.enabled`, `protease.enabled` + nested keys in `docs/PARAMETERS.md` |
| Cell biology (Spec 3 / 10v2) | `fur.enabled`, `cdi.enabled`, `motility.enabled` + aerotaxis/energy/surface/mucin keys; per-strain `cdi_type`, `cdi_immunity` |
| Quorum sensing (Spec 11) | `quorum_sensing.enabled` + AI-2 production/import/chemotaxis keys in `docs/PARAMETERS.md` |
| BtuB/B12 decoupling | `metabolism.b12_growth_essential` (default `true`; set `false` to decouple B12 growth limitation from BtuB receptor competition) |
| Wall-normal Peclet QSSA | `qssa.drift_envelope_policy` (`warn`/`error`/`allow`/`auto`), `qssa.peclet_z_threshold` (`0.05`), `qssa.drift_correction` — see `docs/NEUMANN_WALL_NORMAL_DRIFT.md` |
| Chemical memory guard | `memory.max_host_chemical_gb` (`0.0` = 90% available RAM) and `memory.replicated_oom_policy` (`error`/`auto_slab`/`warn`) |
| Multiscale Layer 2 / Layer 3 (Spec 13) | `layer2.*` (including `live_physics_mode`, `live_inherit_oxygen`, `live_align_patch_table_delivery`, `audit.patch_indices`, perturbation pulses) and `layer3.*` — see `docs/SPEC13_MULTISCALE.md` |
| Temperate prophage (Spec 14) | `prophage.enabled` + carriage/Hill induction/virion keys in `docs/SPEC14_PROPHAGE_MODULE.md` and `docs/PARAMETERS.md` (presets in `examples/prophage_presets/`) |
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
- `docs/CHECKPOINT_ANALYSIS.md` — stream Spec-4 `step_*.h5` checkpoints → CSV/PNG
- `docs/PARAMETERS.md` — configurable parameters
- `docs/SCALING.md` — agent-count benchmarks and profiling
- `docs/AWS_BATCH.md` — AWS Batch Spot + CUDA deployment + Phase Observability
- `docs/MULTI_SCALE_EXPERIMENTATION.md` — nested short/small vs host-scale experiment ladder
- `docs/SPEC13_MULTISCALE.md` — Layer-2/Layer-3 multiscale colon simulation reference
- `docs/SPEC14_PROPHAGE_MODULE.md` — Spec 14 temperate prophage module and dual-regime presets
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

## Open Issue Tracker (Nov 2026)

Jun–Aug 2026 queues (#40–#81, #154–#161 except #159, #256–#306) and the entire Nov 2026 Cycle 2 queue (#464–#483 via PRs #484–#495) are **closed**. Open GitHub issues as of 2026-11-22:

| Issue | Topic | Status |
|-------|-------|--------|
| #446 | Campaign storage index (automated S3 collector) | Open |
| #159 | FMM: sub-quadratic CPU M2L landed (#164); full GPU far-field octree traversal remains open | Open |

Remaining long-horizon: `mpirun -np 8+` on HPC, multi-GPU NCCL, GPU FMM octree traversal (#159).

**Project board:** [docs/PROJECT_BOARD.md](docs/PROJECT_BOARD.md) — kanban layout, completed Cycle 2 wave, and open queue. Run `./scripts/setup_project_board.sh` to create GitHub labels, milestones, and a Projects v2 board.

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
- **Non-GPU CTest targets pass here (171 total with MPI+HDF5 enabled; 25 `gpu`-labelled targets skip without CUDA/GPU).** The `libhwloc`/MPI "preexisting failures"
  table in `.agents/skills/testing-gutibm/SKILL.md` describes a *different* (Devin)
  VM; on Cursor Cloud MPI/HDF5 init works and the full non-GPU suite (incl. `mpi_multi_rank`,
  `hdf5_*`, `scaling_benchmark`) passes.
- **`ruff`/`pytest` install to `~/.local/bin`**, which is not on `PATH` by default —
  prefix with `PATH="$HOME/.local/bin:$PATH"` or add it, or call `python3 -m pytest`.
- **Running the simulation:** binary is `build/gut_ibm <config.json>`. For
  `mpirun -np N` with N greater than the available cores, add `--oversubscribe`.
  The full `examples/single_colony/input.json` (86400 s / 1440 steps) is compute-
  heavy and writes a very large HDF5; use a short/small config for quick checks.
