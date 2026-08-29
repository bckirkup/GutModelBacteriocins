# Independent Audit of GutModelBacteriocins

**For:** Benjamin Kirkup and Devin  
**Audit date:** 2026-08-29 UTC  
**Audit type:** Independent scientific software audit

## 1. Executive verdict

GutModelBacteriocins, also called GutIBM, is a serious and unusually transparent research prototype. Its strongest features are substantive rather than cosmetic: modular mechanism fixes with visible execution order; current-step toxin refresh before receptor killing; species-specific receptor and ligand competition; evolvable partial resistance; timestep-aware stochastic hazards; detailed mass, event, and provenance ledgers; assertions retained in Release tests; broad serial, OpenMP, Message Passing Interface (MPI), HDF5, Python, and CUDA test coverage; true CUDA kernels rather than a nominal wrapper; a real Amazon Web Services (AWS) T4 device gate; and documentation that records important known limitations rather than hiding them. These features make the model inspectable and make decisive repair and ablation work feasible. Evidence is consolidated from the scientific-validity, CUDA, software-engineering, and claim-evidence audits. (task:cdd046b0-61fe-459b-9725-e2b2364d6654; task:180ea40c-de07-4423-ada9-8675f07c6e84; task:3969f802-a51c-4b43-8a0d-d2ec122df88a; task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae)

The dissent is also clear. The current repository is suitable for **hypothesis generation, mechanism exploration, and relative comparisons within a declared operating envelope**. It is **not yet defensible as a validated causal explanation of Enterobacteriaceae diversity**. Two critical defects affect the physical and resource foundations of the principal claim: the reflecting-wall toxin kernel does not implement the complete two-wall Neumann image series, and the default uptake mode allows biomass growth that is not constrained by realized carbon removal. Fail-open configuration handling adds a separate scientific control-plane hazard. Several high-severity backend, restart, event-semantics, and validation defects further prevent unrestricted pooling or causal interpretation. (task:cdd046b0-61fe-459b-9725-e2b2364d6654; task:3969f802-a51c-4b43-8a0d-d2ec122df88a; task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae)

This is not a judgment that the EARI/VADI mechanism is false. Eco-Advective Receptor Interference (EARI) and Viscous Advective-Diffusion Interference (VADI) are represented in code through toxin release and transport, receptor competition, metabolic penalties, mutation, flow, and washout. The judgment is narrower: current outputs cannot identify that causal explanation over resource-only, toxin-only, neutral priority-effect, or imposed-threshold alternatives, and several advertised observables are presently invalid or untested. A methods paper can claim a **testable implementation of a coupled mechanistic hypothesis**. It should not yet claim a validated explanation of the diversity paradox. (task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae)

> ### Version and execution scope
> - **Authoritative repository:** `https://github.com/bckirkup/GutModelBacteriocins`
> - **Audited code snapshot:** `main` at `7bfa7c644de47c2875e58cc6971691d29b2fe856`, recorded in `.AUDIT_SOURCE.json:2-6`.
> - **Main at report time:** `89c40bf896e947ef7bc8ea1cee59e6fa5fe62390`. Comparison showed changes only to `AGENTS.md` and Spec 14/density/reference documentation, not executable source. The code findings therefore still apply to latest `main` at report time.
> - **Latest-main continuous integration:** eight reported jobs succeeded, covering unit, integration, serial, OpenMP parity, Python, EARI/VADI regression, CUDA compile/tests, and GPU parity.
> - **Local execution:** Python package installation passed. `pytest -s` reported **133 passed, 2 deselected**. Native C++ compilation stopped on missing `<format>` under GCC 12. This is an environment limit plus a README toolchain-documentation gap, not a failure of the audited scientific code. Local hardware was CPU-only and lacked usable MPI/GCC 13, so GPU, MPI, and C++ binaries were not executed locally and are not penalized.
> - **GPU evidence:** ordinary GPU parity may pass by an explicit no-device skip. The separate AWS Batch T4 workflow is real, requires `nvidia-smi` and non-skipped GPU CTests, and has multiple recent successful pull-request runs. It did not run at the exact audited merge SHA, so exact-SHA physical-device execution remains unverified.
> - **Source integrity:** repository source was not modified.

## 2. Strengths to preserve

These should be retained while defects are repaired.

| Strength | Why it matters | Evidence |
|---|---|---|
| Modular `Fix` architecture and explicit order | Mechanisms can be ablated and their ordering can be audited. The sequence is documented and visible in the simulation loop. | `docs/MECHANISMS.md:89-103`; `src/core/simulation.cpp:1645-1682,1713-1768` (task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae) |
| Current-step toxin refresh | Toxin sources are refreshed immediately before receptor killing, avoiding a simple one-step exposure lag. | `src/core/simulation.cpp:1658-1682,1739-1747` (task:cdd046b0-61fe-459b-9725-e2b2364d6654) |
| Species-specific receptor competition | BtuB/B12, FepA/ferric-enterobactin, CirA/linearized ferric-enterobactin, and FhuA/ferrichrome are distinct. Competitive occupancy and partial-affinity resistance are explicit. | `src/fixes/fix_receptor.h:27-41`; `src/fixes/fix_receptor.cpp:114-159,194-210` (task:cdd046b0-61fe-459b-9725-e2b2364d6654) |
| Evolvable partial resistance | Receptor downregulation, separate toxin/ligand affinity changes, and compensatory mutations create testable escape routes rather than a binary declared phenotype. | `src/fixes/fix_mutation.cpp:95-131,155-169` (task:cdd046b0-61fe-459b-9725-e2b2364d6654) |
| Timestep-aware hazards | SOS and receptor death use exponential event conversion rather than a fixed per-step probability. Receptor random draws remain serial on the host in the GPU path. | `src/fixes/fix_bacteriocin.cpp:99-107`; `src/fixes/fix_receptor.cpp:49-82,141-159` (task:cdd046b0-61fe-459b-9725-e2b2364d6654; task:180ea40c-de07-4423-ada9-8675f07c6e84) |
| Mass, event, and provenance ledgers | Boundary fluxes, VBF source/sinks, uptake, maintenance, clipping, lysis, colicin death, washout, HGT, mutation, and termination state are observable. These outputs can support causal tests after semantics are corrected. | `docs/CARBON_BUDGET.md:114-143`; `src/io/hdf5_writer.cpp:834-966,1006-1051` (task:cdd046b0-61fe-459b-9725-e2b2364d6654; task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae) |
| Release assertions and broad CI | Tests retain assertions in Release. The suite checks conservation, directional behavior, parity, restart content, and failure paths, not only process exit. | `tests/CMakeLists.txt:34-46`; `.github/workflows/ci.yml:13-219` (task:3969f802-a51c-4b43-8a0d-d2ec122df88a) |
| Real CUDA implementation and device gate | Metabolism, receptor probability, chemistry, diffusion, Green deposition, VBF coupling, and mechanics have device kernels. AWS T4 testing requires a physical device and non-skipped GPU tests. | `src/gpu/*.cu`; `.github/workflows/gpu-device-tests.yml:23-53`; `scripts/run_gpu_device_tests.sh:218-241,341-357` (task:180ea40c-de07-4423-ada9-8675f07c6e84) |
| Candid limitations | Unfunded uptake, approximation limits, GPU boundaries, and long-horizon uncertainty are documented. This lowers the cost of correction and should not be replaced with stronger claims. | `AGENTS.md:127-131`; `docs/SCALING.md:164-186`; `src/gpu/README.md:100,133-143` (task:cdd046b0-61fe-459b-9725-e2b2364d6654; task:180ea40c-de07-4423-ada9-8675f07c6e84) |

## 3. Ranked verified defects

“Verified defect” means source logic, shipped configuration, or analysis behavior contradicts its stated semantics, governing construction, or requested scientific record. Risks that require empirical discrimination are separated in Section 4.

| Rank | Severity | Finding | Scientific consequence | Evidence | Cheapest repair and acceptance test |
|---:|---|---|---|---|---|
| 1 | **P0 / Critical** | **Incomplete and duplicated two-wall Neumann image series.** For a slab bounded at 0 and `H`, the required source families are `z_s + 2nH` and `-z_s + 2nH`. The implementation duplicates members of the reflected family and omits translated-source terms. `N_IMAGES` is fixed at three. Existing bounded and fast-multipole tests call the same bounded kernel as their oracle. | Reflecting-wall toxin concentration and geometry are not established. This directly compromises VADI wall confinement, core/halo reach, and exposure near both walls. | `src/diffusion/greens_function.cpp:317-345`; `src/diffusion/greens_function.h:96-98`; `tests/test_greens_function.cpp:296-300`; `tests/test_fmm.cpp:118-127` (task:cdd046b0-61fe-459b-9725-e2b2364d6654; task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae) | Generate both families without duplication and either converge adaptively or justify a truncation bound. Test against an independently coded 30 to 100-pair analytical series and a finite-volume reference. Require both wall-normal derivatives near zero and preregistered concentration error across source depths, decay, flow, and wall proximity. |
| 2 | **P0 / Critical** | **Default `uptake_limit=none` permits biomass growth unfunded by realized carbon removal.** The lead audit counted 94 of 96 checked-in JSON configurations without an explicit `uptake_limit`; generated configurations may override this. Funded delivery exists but is rejected with GPU until parity is implemented. | Quantitative density, carrying capacity, washout, oxygen demand, lysis, and mass-balance claims can be driven by substrate not removed from the field. The documentation is candid, but the scientific default remains unsafe. | `src/fixes/fix_metabolism.h:70-76`; `AGENTS.md:127-130`; `src/gpu/chemistry_pipeline.cpp:185-191`; `src/io/input_parser.h:88-92`; GPU rejection at `src/io/input_parser.cpp:638-664` (task:cdd046b0-61fe-459b-9725-e2b2364d6654; task:180ea40c-de07-4423-ada9-8675f07c6e84) | Make funded uptake mandatory in scientific configurations, with explicit legacy opt-in. In a closed one-agent domain, require cumulative biomass carbon demand to match initial minus final field carbon within numerical tolerance. Require reaction clipping below a declared threshold. Do not pool funded CPU and unfunded GPU arms. |
| 3 | **P0 / Critical** | **Configuration ingestion fails open.** A missing or unreadable file returns defaults; malformed JSON may fall back to the legacy parser; unknown keys warn unless an optional environment variable is set. | A campaign can complete successfully while running a scientifically different experiment. This is a control-plane failure, not merely input ergonomics. | `src/io/input_parser.cpp:1794-1815`; `src/main.cpp:53-70`; `src/io/config_json.cpp:545-553,651-726`; `src/io/input_parser.cpp:1770-1779` (task:3969f802-a51c-4b43-8a0d-d2ec122df88a) | Require a readable strict JSON config for production. Reject malformed tokens, trailing content, duplicate/unknown keys, invalid numeric types, and non-finite values. Put legacy parsing and demo defaults behind explicit flags. CI must fail on a missing file, `"6O"`, `NaN`, a misspelled mechanism key, and trailing garbage. |
| 4 | **High** | **GPU partial per-species diffusion can suppress CPU fallback.** `apply_diffusion()` returns whether any species ran, although eligibility differs by boundary and line length. At `nz=1025`, a Dirichlet species can execute while a Robin or Flux species is skipped; the aggregate success then suppresses host diffusion. | An accepted mixed-boundary GPU configuration can leave a species undiffused for that step, preventing valid CPU/GPU pooling. | `src/gpu/chemical_field_gpu.cpp:298-330`; `src/gpu/diffusion_gpu.cpp:161-167`; `src/gpu/chemistry_pipeline.cpp:199-217` (task:180ea40c-de07-4423-ada9-8675f07c6e84) | Preflight all species and use all-or-nothing fallback, or return a per-species completion mask. Add the exact `nz=1025` mixed Dirichlet plus Robin/Flux reproducer and require field and ledger parity with CPU. |
| 5 | **High** | **GPU mechanics omits active persistent CDI corpses while dispatch still selects GPU.** CPU mechanics includes recent CDI corpses in repulsion; the kernel drops state 3 from both sides of a pair. | Cell positions and contact structure diverge by backend whenever active corpse persistence is present. | `src/fixes/fix_mechanics.cpp:19-24,88-101,164-195`; `src/gpu/mechanics_kernel.cu:94-128`; `src/gpu/README.md:100,140` (task:180ea40c-de07-4423-ada9-8675f07c6e84) | Immediately route mechanics to CPU when an active persistent corpse is possible. Later add corpse state and death time to the device model. Test one active corpse and neighboring live cells for force/displacement parity. |
| 6 | **High** | **Checkpointing restores population and chemistry, but not random-number-generator state.** Current restart is a reseeded fork, not an exact continuation. | Split runs cannot be interpreted as the same stochastic trajectory. Exact resume and experimental branching are conflated operationally. | `src/core/simulation.cpp:530-532,771-794,874-944`; `src/io/hdf5_reader.h:19-99`; `docs/BRANCHING_FROM_CHECKPOINTS.md:53-55` (task:3969f802-a51c-4b43-8a0d-d2ec122df88a) | Expose `resume_exact` and `fork_reseeded`. Serialize both RNG streams and any distribution state. Require uninterrupted versus split stochastic equivalence in serial and supported parallel modes where exactness is promised. |
| 7 | **High** | **“1% per division” is applied to both parent and daughter.** Both are marked `just_divided`; each receives an exponential SOS draw. At 0.01 this is `1 - exp(-0.01)^2 = 1.98%` before basal and other hazards. | Lysis and toxin exposure are nearly doubled relative to the labeled per-division estimand, changing calibration and selection. | `src/fixes/fix_metabolism.cpp:893-897`; `src/fixes/fix_bacteriocin.cpp:86-107`; `docs/PARAMETERS.md:680-686` (task:cdd046b0-61fe-459b-9725-e2b2364d6654) | Decide whether the parameter means per division event or per resulting cell, rename if needed, and implement once. With all other hazards off, simulate at least one million divisions across several `bio_dt` values and require the declared estimand. |
| 8 | **High** | **Comet-tail analysis fabricates geometry.** It flattens toxin arrays and assigns an arbitrary 0 to 1 mm one-dimensional x coordinate with y=z=0. | The reported anisotropy follows storage order rather than the actual three-dimensional grid, flow vector, source, boundaries, or spacing. | `python/gut_ibm_tools/validation.py:49-67`; `python/gut_ibm_tools/analysis.py:201-257` (task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae) | Read grid shape, spacing, origin, and boundary metadata. Compute source-registered three-dimensional anisotropy along the known velocity vector with periodic handling. Test on synthetic plumes rotated across axes and on a simulator-generated grid. |
| 9 | **High** | **Resident retention does not measure resident-strain descendant retention.** It intersects unique lineage identifiers at the first and last saved steps and labels final noninitial lineages transient. | The 70 to 80% comparison does not represent resident strains, descendants, immigrant cohorts, censoring, or metagenomic detection. | `python/gut_ibm_tools/validation.py:109-156`; `python/gut_ibm_tools/validation_regression.py:44-55` (task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae) | Define residents and immigrants explicitly, follow descendants, specify observation intervals and detection/censoring, and validate the estimator on known pedigrees. |
| 10 | **High** | **The flagship immigration scenario has no immigration and disables grid output.** Type 2 is present at time zero, and no immigration block exists. | The scenario cannot test periodic immigrant establishment as described, and it cannot audit the advertised toxin comet-tail field. | `examples/diversity_paradox/README.md:11-22`; `examples/diversity_paradox/input.json:38-63` (task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae) | Add explicit immigration timing/cohorts and enable required toxin-grid output. Test arrival counts, labels, and availability of every preregistered observable. |
| 11 | **Medium** | **Taylor-Aris toggle is dead.** The formula exists but has no caller while documentation says it is enabled and captures shear enhancement. | A named VADI transport mechanism has no effect. | `src/fields/advection.cpp:80-97`; `docs/PARAMETERS.md:124-139` (task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae) | Wire the effective diffusivity into the intended kernel and test toggle sensitivity, or remove the toggle and narrow the claim. |
| 12 | **Medium** | **VBF drag and carrying capacity are stored but unused.** | The VBF provides nutrient coupling, but does not presently deliver the advertised mechanical drag or local carrying-capacity control. | `src/fields/vbf.cpp:176-180,248-254`; `src/fields/vbf.h:85`; `docs/MECHANISMS.md:891-908` (task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae) | Connect each feature to mechanics/density with directional tests and outputs, or remove the executable claim. |
| 13 | **Medium** | **Requested HDF5 output can fail open.** | A long run can complete without its requested scientific record. | `src/io/hdf5_writer.cpp:471-485,534-541` (task:3969f802-a51c-4b43-8a0d-d2ec122df88a) | Make `hdf5.enabled=true` fatal on creation/write failure by default; permit best effort only explicitly and record it in exit status and manifests. |
| 14 | **Medium** | **Source-archive provenance may be `unknown-git-unavailable`; compiler contract is incomplete.** GCC 13+ is effectively required for `<format>`, but the primary README does not say so. | Results from source archives may lack an authoritative code identity, and documented generic builds can fail after configuration. | `cmake/generate_git_sha.cmake:5-34`; `CMakeLists.txt:12-14`; `src/io/hdf5_writer.cpp:563-567`; `README.md:113-133`; `deploy/aws/Dockerfile:15-18,25-26` (task:3969f802-a51c-4b43-8a0d-d2ec122df88a) | Consume `.AUDIT_SOURCE.json` or a release identity file, fail scientific releases with unknown identity, record source digest, document GCC 13+, and probe `<format>` at configure time. |
| 15 | **Medium** | **Python manifest writes are non-atomic.** | A crash during rewrite can destroy resumability and campaign provenance. | `python/gut_ibm_tools/batch_manifest.py:107-155` (task:3969f802-a51c-4b43-8a0d-d2ec122df88a) | Write, flush, `fsync`, and atomically replace a same-directory temporary file; retain one validated prior generation. |
| 16 | **Medium** | **Ethanolamine absolute units are off by 1000.** Both concentration and `K_m` preserve their ratio, but not their documented physical values. | The current Monod penalty may be numerically unchanged, but absolute labels and future cross-mechanism coupling are wrong. | `src/io/input_parser.cpp:194-196`; `src/fixes/fix_metabolism.h:66-68`; `docs/MECHANISMS.md:340-343`; `docs/UNITS_AUDIT.md:111` (task:cdd046b0-61fe-459b-9725-e2b2364d6654) | Multiply both by 1000, assert current ratio-based penalty invariance, and test every other ethanolamine-dependent output. |

### Independent Neumann residual check

The lead auditor independently recomputed a representative screened no-flow kernel using `H = 100 µm`, source `z/H = 0.3`, x offset `20 µm`, `D = 4×10^-11 m²/s`, and decay `10^-4 s^-1`. The current `N_IMAGES=3` series gave normalized wall slopes `(dC/dz)H/C = +0.0072` at the lower wall and `-0.103` at the upper wall, versus approximately zero for a correct 30-pair series. Across source depths `0.05` to `0.95 H`, the median absolute normalized slope was `0.044` and the maximum was `0.216`. These values demonstrate the defect under a representative setting. They are **not** universal concentration or flux error bounds.

One additional scale fact should not be misclassified as a new defect. At `grid_dx = 2 µm`, one agent per voxel corresponds to `1.25×10^11 cells/mL`. Earlier Damköhler analysis showed that geometric coherence, rather than density alone, controls whether a local sink matters. This is relevant when redesigning funded uptake and multiscale density experiments, but it does not replace the verified carbon-closure defect.

## 4. Scientific risks that are not proven defects

| Risk | Why it matters | Cheapest discriminating test |
|---|---|---|
| **Quasi-steady-state approximation (QSSA) timescale separation** | With default toxin `D=4×10^-11 m²/s`, retardation 10, and a 200 µm cutoff, the rough diffusion time `L²/D_eff` is about 10,000 s or 2.8 h, compared with a typical 60 s biology step. Decay relaxation is also roughly 10,000 s. A steady plume may therefore appear too quickly after release. Advection, decay, distance, and release shape prevent treating this estimate as proof of failure. | Compare one-burst QSSA with a transient finite-volume toxin solver at 1, 5, 30, and 180 minutes across retardation, decay, flow, and source-target distance. Define a dimensionless acceptance region before production use. `src/diffusion/greens_function.cpp:273-298`; `src/diffusion/qssa_solver.h:37-45` (task:cdd046b0-61fe-459b-9725-e2b2364d6654) |
| **Variable flow collapsed to source-local constant velocity** | The Green kernel evaluates velocity at the source once, although flow varies with z and peristalsis. This may misstate transport across velocity gradients. | Compare against a steady finite-volume advection-diffusion solution for sources near low- and high-velocity regions. Swap source and target heights to expose source-local asymmetry. `src/diffusion/greens_function.cpp:285-315` (task:cdd046b0-61fe-459b-9725-e2b2364d6654) |
| **Robin/flux oxygen starts at zero** | Non-Dirichlet oxygen modes begin anoxic, potentially overlapping early growth, metabolic switching, and SOS history. This may be intended for some experiments. | Compare zero initialization with the corresponding steady Robin/flux profile and quantify early respiration, fermentation, division, lysis, and required burn-in. `src/io/input_parser.cpp:675-683` (task:cdd046b0-61fe-459b-9725-e2b2364d6654) |
| **Lie reaction-diffusion splitting plus clipping** | Forward reaction, non-negativity clipping, then full-step diffusion is first-order. Strong sinks can create timestep dependence, while clipping can hide it. | Use a manufactured reaction-diffusion solution at `dt`, `dt/2`, and `dt/4`; report convergence order, clipped mass, closure, and biological outputs. Require negligible clipping for production. `src/gpu/chemistry_pipeline.cpp:171-213`; `src/io/input_parser.h:88-92` (task:cdd046b0-61fe-459b-9725-e2b2364d6654) |
| **Adaptive SOS timestep uses absolute counts** | Thresholds at 5 and 20 SOS cells can give different timestep histories for equivalent concentrations in different simulated volumes. | Replicate the same density and SOS fraction at 1×, 10×, and 100× volume. Compare timestep and event-rate distributions; prefer a fraction or hazard-density rule. `src/core/simulation.cpp:1179-1190` (task:cdd046b0-61fe-459b-9725-e2b2364d6654) |
| **Daughter placement near bounded z walls** | Division offsets daughters in all dimensions, but periodic correction does not act on bounded z. Near-wall daughters may leave the domain or be clamped. | Divide agents one radius from each z wall over many seeds; assert in-domain placement and compare immediate loss with midplane controls. `src/fixes/fix_metabolism.cpp:876-885`; `src/core/domain.cpp:78-80` (task:cdd046b0-61fe-459b-9725-e2b2364d6654) |
| **Central structural non-identifiability** | Burst, retardation, decay, receptor affinity/expression, immunity, nutrient supply, flow, refugia, HGT, and mutation can move the same patchiness and retention endpoints. Current outputs and sensitivity tests do not uniquely assign causality. | Run a preregistered factorial ablation and compare held-out predictive likelihood or scores against resource-only, toxin-only, neutral priority-effect, and imposed-threshold models. Use independently measured unit parameters. (task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae) |

## 5. CUDA assessment

### Supported operating envelope

The CUDA backend is genuine. It executes device kernels for metabolism, receptor hazard probabilities, reaction integration, directional diffusion, oxygen/VBF chemistry, near-field Green deposition, and soft-sphere mechanics. Host execution for advection, mutation, conjugation, division, migration, HDF5, random receptor-death draws, and parts of fast-multipole and MPI work is not itself a defect. The architecture is host-authoritative and intentionally hybrid. (task:180ea40c-de07-4423-ada9-8675f07c6e84)

Until the two high-severity defects are fixed, safe experimental use should be declared narrowly: full supported species configuration, grid-mode toxin evaluation, line lengths where every diffusing species is eligible, no active persistent CDI corpses, and no claim that funded delivery CPU runs are interchangeable with GPU runs that reject that mode. The partial-diffusion and corpse-mechanics cases require guards before general backend pooling.

### Correctness defects and required guards

1. **Per-species diffusion eligibility:** make the operation all-or-nothing or return a completion mask. The acceptance case is `nz=1025` with mixed Dirichlet and Robin/Flux boundaries, requiring all species to match CPU or the full step to fall back.
2. **CDI corpse mechanics:** dispatch to CPU whenever active persistent corpse semantics are required, until the device representation includes them. Require parity for live-cell displacement around an active corpse.

These are semantic mismatches. By contrast, host-side random draws, host FMM tree construction, host migration, and explicit unsupported-mode fallbacks are design boundaries, not failures, if the dispatch is correct and documented.

### Performance architecture after correctness

The largest likely gains are architectural rather than kernel arithmetic:

- `AgentPoolGpu::resize()` reallocates arrays on each call, and full agent state is repeatedly packed and uploaded at several mechanism boundaries (`src/gpu/agent_pool_gpu.cpp:13-38`; `src/gpu/device_memory.h:44-54`).
- Grid coupling is executed on the GPU, downloaded, then recomputed on the CPU; an upload immediately after MPI migration precedes host washout and is uploaded again next step (`src/core/simulation.cpp:1633-1652,1684-1698`).
- Upload/download helpers synchronize immediately and use pageable vector staging, limiting copy-compute overlap (`src/gpu/device_memory.h:70-105`).
- Spatial hashing computes keys on device but downloads them for CPU counts, prefix sums, and sort before re-uploading CSR arrays (`src/gpu/spatial_hash_gpu.cpp:18-83`).
- Several wrappers allocate scratch buffers in hot paths and synchronize after short launch groups.
- CUDA-aware MPI is opt-in and limited to replicated reaction grids. Slab x diffusion, migration, source exchange, and other reductions use host paths. No NCCL implementation exists. These are scaling limits, not correctness defects. (task:180ea40c-de07-4423-ada9-8675f07c6e84)

Correctness must precede optimization. Then preserve capacity in device buffers, eliminate redundant transfers, move compressed sparse row (CSR) construction to device, use pinned staging and asynchronous streams where measurement supports them, and profile allocator calls, bytes, barriers, and end-to-end throughput.

### Physical-GPU evidence

The ordinary GPU parity workflow can report success after a no-device skip, so it is compile and dispatch evidence rather than guaranteed physical-device evidence. The dedicated AWS Batch T4 workflow is materially stronger: it requires `nvidia-smi`, requires non-skipped GPU CTests, and has multiple recent successful pull-request runs. Its path filters are narrow and it did not execute at exact audited merge SHA `7bfa7c6`. Latest `main` reports successful CI, but exact-SHA device execution for the audited snapshot remains an evidence gap, not a demonstrated failure. Widen the device-gate path set to shared core, field, diffusion, fix, CMake, and test files that can alter CUDA behavior. (task:180ea40c-de07-4423-ada9-8675f07c6e84; task:3969f802-a51c-4b43-8a0d-d2ec122df88a)

## 6. Claim-evidence verdict

EARI/VADI components are substantially implemented:

- toxin release routes and finite burst inventory;
- charge/retardation-dependent effective diffusion;
- receptor-specific toxin/ligand occupancy and kill hazard;
- receptor-linked metabolic penalties and evolvable resistance;
- physical advection and boundary loss;
- nutrient and chemistry coupling;
- HGT, mutation, event provenance, and spatial outputs.

The full causal claim remains underidentified. Multiple free mechanisms can produce similar abundance, patchiness, exclusion, and retention outcomes. The current validation layer does not resolve that ambiguity:

- EARI/VADI continuous integration invokes only the FISH target check, not the full target suite (`scripts/validate_eari_vadi.sh:23-34`).
- Spatial and retention thresholds are encoded from internal targets and exercised on synthetic fixtures (`python/tests/test_eari_vadi_regression.py:30-45`; `python/gut_ibm_tools/validation_regression.py:18-55`). These are regression checks, not external biological validation.
- The short validation README says one hour, but `total_time: 900` is 15 minutes (`examples/eari_vadi_validation/README.md:3-5`; `examples/eari_vadi_validation/input.json:8-10`).
- There is no held-out biological dataset, fitted uncertainty model, or comparison with competing simpler causal models.
- The reflecting-wall kernel, comet-tail geometry, resident-retention estimand, and flagship immigration scenario are currently defective.

The FISH observation model is useful as an arithmetic and feasibility pipeline, but its copy number, gain, brightness, noise, and threshold assumptions are internal. CI confirms behavior under those assumptions rather than validating a microscope or a biological dataset. The honest paper-level claim is therefore: **GutIBM implements a falsifiable, coupled EARI/VADI hypothesis and produces mechanism-resolved observables.** The unsupported claim is: **EARI/VADI has been validated as the causal explanation for limited Enterobacteriaceae diversity.** (task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae)

## 7. Ranked repair plan for Devin

### P0: correctness and scientific control

1. **Correct the Neumann slab kernel.** Implement both image families and remove duplicates. Add an independent analytical oracle, wall-normal derivative tests at both walls, image-count convergence, near-wall and midplane cases, decay/no-decay cases, and a finite-volume comparison. Do not use `concentration_bounded()` as the reference.
2. **Fail closed on configuration.** Make strict JSON mandatory for production. Reject missing/unreadable input, unknown or duplicate keys, malformed numerics, trailing data, invalid integer forms, and non-finite values. Legacy/demo behavior must require explicit command-line selection.
3. **Require funded uptake for scientific configurations.** Make `delivery` or a validated `sherwood` mode explicit. Gate campaigns on carbon closure and negligible clipping. Keep legacy `none` only as a named compatibility mode. Until funded GPU parity exists, prohibit CPU/GPU pooling for these campaigns.
4. **Add both GPU correctness guards.** Use all-or-nothing or masked diffusion fallback, with the `nz=1025` mixed-boundary regression. Route active CDI corpse mechanics to CPU, with a one-corpse parity regression.
5. **Resolve SOS semantics.** Decide per division event versus per resulting cell, rename/document the parameter, and validate the frequency across timesteps with all other hazards disabled.

**P0 release gate:** all five acceptance suites pass; no self-oracle in the boundary test; malformed configurations exit nonzero; closed-domain carbon closure is within a declared tolerance; the two GPU reproducers either match CPU or demonstrably fall back; SOS frequency matches the named estimand.

### P1: validation semantics and numerical evidence

1. **Correct observables and scenario:** reconstruct toxin geometry from actual HDF5 coordinates and flow; define resident-strain descendant retention with immigrant cohorts and censoring; add real immigration and toxin-grid output to the flagship scenario; correct the 15-minute label.
2. **Separate exact resume from fork:** serialize RNG state for exact continuation and retain reseeded branching as an explicit mode. Test stochastic uninterrupted versus split runs where exactness is promised.
3. **Make scientific output and provenance enforceable:** requested HDF5 should fail fatally by default; enforce schema compatibility; consume archive source identity; record compiler, standard library, dependencies, flags, device/driver, MPI/OpenMP environment, executable or container digest, and output policy.
4. **Validate QSSA against a transient solver:** report exposure and field error over burst age, flow, retardation, decay, and distance. State a valid operating region rather than a universal QSSA claim.
5. **Run convergence ladders:** biology and chemistry timestep; grid spacing; domain and wall distance; cutoff and image count; FMM order/opening; source heterogeneity; and stochastic replicate count. Include mass, wall flux, event ledgers, exposure, lineage establishment, and spatial metrics.

**P1 release gate:** every named validation metric has a valid estimand and synthetic truth test; one simulator-generated reference artifact is SHA-tagged; exact resume and fork differ by explicit policy; failed HDF5 requests cannot silently complete; QSSA and numerical convergence bounds accompany every reported operating envelope.

### P2: performance, hardening, and external validity

1. **Optimize CUDA only after parity:** capacity-preserving buffers, fewer full-state uploads, removal of redundant coupling/post-migration transfers, pinned asynchronous copies, device CSR construction, reused scratch storage, and broader profiling.
2. **Harden native testing:** AddressSanitizer plus UndefinedBehaviorSanitizer, targeted ThreadSanitizer, static analysis, minimum-compiler matrix, full-state/thread-count parity, and atomic Python manifests.
3. **Widen physical-GPU coverage:** trigger the AWS T4 gate on shared core, fixes, fields, diffusion, build, and test paths. Report compiled, skipped, and physically executed counts separately. Add the two edge regressions and sanitizer runs where supported.
4. **Build observation models and external evidence:** calibrate FISH imaging on controls; use domain-aware three-dimensional spatial nulls; collect unit-level toxin, receptor, growth, flow, and HGT measurements; reserve held-out biological data; compare against simpler models.

## 8. Minimum evidence for a defensible methods paper

1. **Independent partial differential equation verification.** Verify infinite-domain advection-diffusion-decay, the corrected two-wall Neumann solution, wall-normal flux, mass/source balance, image convergence, finite-volume agreement, and FMM error across Péclet and Damköhler ranges. Error criteria must be set before examining the final biological endpoints.
2. **Unit-level empirical parameterization separated from validation.** Measure toxin diffusion/retardation and decay in mucus; receptor expression/binding/kill under measured ligand levels; growth penalties for receptor states across carbon, iron, corrinoid, acetate, and oxygen conditions; flow-cell clearance; and HGT versus shear. Estimate uncertainty distributions. Hold out at least one condition in each relevant assay class.
3. **Factorial causal ablation.** At minimum cross toxin off/on, receptor penalty off/on, and flow off/on or mode. Include VBF competition, HGT/mutation, refugia/crypts, and resistance genotype as secondary factors. The central evidence is the preregistered **toxin × receptor penalty × flow interaction**, not reproduction of one endpoint.
4. **Competing-model comparison.** Evaluate EARI/VADI against resource-only, toxin-only, neutral priority-effect, and imposed-threshold alternatives using held-out likelihood or predictive scores. Fit unit parameters to unit assays, then predict colony outcomes without endpoint retuning.
5. **Correct observation models.** Use real three-dimensional coordinates and source registration for toxin anisotropy; use domain- and periodicity-aware spatial nulls; define resident strains, descendants, immigrants, detection, and censoring; calibrate FISH noise and detection externally.
6. **Stochastic replicate distributions.** Report distributions, extinction/establishment probabilities, Monte Carlo uncertainty, early termination, and ledger closure. Use 30 to 100 replicates per key contrast only if justified by power or Monte Carlo error, not as an automatic number.
7. **One frozen reproducibility package.** Freeze source SHA, source-tree digest, compiler/container, exact configurations, seeds, raw outputs, calibration and holdout data, analysis environment, preregistered metrics, and scripts that produce every result. Include one CPU reference and backend parity evidence before pooling accelerated results. (task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae)

## 9. Final go/no-go verdict

### Okay now

- **Go:** exploratory mechanism development and hypothesis generation.
- **Go:** software architecture, accounting, and test development.
- **Conditional go:** relative comparisons inside a declared operating envelope, using identical validated semantics, strict configurations, clear uptake mode, no affected GPU edge case, and no interpretation that exceeds the tested numerical regime.

### Not okay yet

- **No-go:** quantitative carrying-capacity, density, washout-threshold, or mass-balance claims under default unfunded uptake.
- **No-go:** reflecting-wall toxin footprints or VADI boundary-confinement claims from the current image series.
- **No-go:** pooling CPU/GPU results for mixed-eligibility diffusion, active persistent CDI corpse mechanics, or funded/unfunded uptake differences.
- **No-go:** causal or validated EARI/VADI explanation of Enterobacteriaceae diversity.
- **No-go:** empirical resident-retention or comet-tail claims from the present validation estimands.

**Bottom line:** preserve the architecture and ledgers, fix correctness before performance, and narrow claims until independent physics, valid observation models, external holdouts, and competing-model tests are complete. The repository is a strong platform for that work. It is not yet the validation result.

## Appendix A. Evidence index

The task citations identify the four independent source audits. File references are to the audited snapshot unless noted.

| Topic | Exact repository evidence | Audit source |
|---|---|---|
| Snapshot provenance | `.AUDIT_SOURCE.json:2-6` | All four source audits |
| Neumann image construction and shared oracle | `src/diffusion/greens_function.cpp:317-345`; `src/diffusion/greens_function.h:96-98`; `tests/test_greens_function.cpp:296-300`; `tests/test_fmm.cpp:118-127` | task:cdd046b0-61fe-459b-9725-e2b2364d6654; task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae |
| Uptake default and closure | `src/fixes/fix_metabolism.h:70-76`; `AGENTS.md:127-130`; `src/gpu/chemistry_pipeline.cpp:185-191`; `src/io/input_parser.h:88-92`; `src/io/input_parser.cpp:638-664` | task:cdd046b0-61fe-459b-9725-e2b2364d6654; task:180ea40c-de07-4423-ada9-8675f07c6e84 |
| Fail-open configuration | `src/io/input_parser.cpp:1794-1815`; `src/io/config_json.cpp:545-553,651-726`; `src/main.cpp:53-70` | task:3969f802-a51c-4b43-8a0d-d2ec122df88a |
| GPU diffusion aggregate return | `src/gpu/chemical_field_gpu.cpp:298-330`; `src/gpu/diffusion_gpu.cpp:161-167`; `src/gpu/chemistry_pipeline.cpp:199-217` | task:180ea40c-de07-4423-ada9-8675f07c6e84 |
| CDI corpse mechanics | `src/fixes/fix_mechanics.cpp:19-24,88-101,164-195`; `src/gpu/mechanics_kernel.cu:94-128` | task:180ea40c-de07-4423-ada9-8675f07c6e84 |
| Checkpoint and RNG semantics | `src/core/simulation.cpp:530-532,771-794,874-944`; `src/io/hdf5_reader.h:19-99`; `docs/BRANCHING_FROM_CHECKPOINTS.md:53-55` | task:3969f802-a51c-4b43-8a0d-d2ec122df88a |
| SOS per-division semantics | `src/fixes/fix_metabolism.cpp:893-897`; `src/fixes/fix_bacteriocin.cpp:86-107`; `docs/PARAMETERS.md:680-686` | task:cdd046b0-61fe-459b-9725-e2b2364d6654 |
| Invalid comet-tail geometry | `python/gut_ibm_tools/validation.py:49-67`; `python/gut_ibm_tools/analysis.py:201-257` | task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae |
| Retention estimand | `python/gut_ibm_tools/validation.py:109-156`; `python/gut_ibm_tools/validation_regression.py:44-55` | task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae |
| Flagship scenario | `examples/diversity_paradox/README.md:11-22`; `examples/diversity_paradox/input.json:38-63` | task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae |
| Taylor-Aris dead path | `src/fields/advection.cpp:80-97`; `docs/PARAMETERS.md:124-139` | task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae |
| VBF unused drag/capacity | `src/fields/vbf.cpp:176-180,248-254`; `src/fields/vbf.h:85`; `docs/MECHANISMS.md:891-908` | task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae |
| HDF5 fail-open | `src/io/hdf5_writer.cpp:471-485,534-541` | task:3969f802-a51c-4b43-8a0d-d2ec122df88a |
| Archive provenance and compiler | `cmake/generate_git_sha.cmake:5-34`; `CMakeLists.txt:8-14`; `README.md:113-133`; `deploy/aws/Dockerfile:15-18,25-26`; `src/io/hdf5_writer.cpp:563-567` | task:3969f802-a51c-4b43-8a0d-d2ec122df88a |
| Atomic manifests | `python/gut_ibm_tools/batch_manifest.py:107-155` | task:3969f802-a51c-4b43-8a0d-d2ec122df88a |
| Ethanolamine units | `src/io/input_parser.cpp:194-196`; `src/fixes/fix_metabolism.h:66-68`; `docs/MECHANISMS.md:340-343`; `docs/UNITS_AUDIT.md:111` | task:cdd046b0-61fe-459b-9725-e2b2364d6654 |
| QSSA and variable flow | `src/diffusion/greens_function.cpp:273-315`; `src/diffusion/qssa_solver.h:37-45` | task:cdd046b0-61fe-459b-9725-e2b2364d6654 |
| Oxygen initialization and splitting | `src/io/input_parser.cpp:675-683`; `src/gpu/chemistry_pipeline.cpp:171-213` | task:cdd046b0-61fe-459b-9725-e2b2364d6654 |
| Adaptive SOS and daughter placement | `src/core/simulation.cpp:1179-1190`; `src/fixes/fix_metabolism.cpp:876-885`; `src/core/domain.cpp:78-80` | task:cdd046b0-61fe-459b-9725-e2b2364d6654 |
| CUDA architecture/performance | `src/gpu/agent_pool_gpu.cpp:13-38`; `src/gpu/device_memory.h:44-105`; `src/gpu/spatial_hash_gpu.cpp:18-83`; `src/core/simulation.cpp:1633-1652,1684-1698` | task:180ea40c-de07-4423-ada9-8675f07c6e84 |
| GPU CI and device gate | `.github/workflows/ci.yml:120-162`; `.github/workflows/gpu-device-tests.yml:3-53`; `scripts/compare_gpu_parity.sh:34-53`; `scripts/run_gpu_device_tests.sh:218-241,341-357` | task:180ea40c-de07-4423-ada9-8675f07c6e84; task:3969f802-a51c-4b43-8a0d-d2ec122df88a |
| EARI/VADI validation limits | `scripts/validate_eari_vadi.sh:23-34`; `python/tests/test_eari_vadi_regression.py:30-45`; `examples/eari_vadi_validation/README.md:3-19`; `examples/eari_vadi_validation/input.json:8-10` | task:cdd046b0-61fe-459b-9725-e2b2364d6654; task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae |
| Mechanism and ledger strengths | `src/core/simulation.cpp:1645-1682,1713-1768`; `src/fixes/fix_receptor.cpp:114-210`; `src/io/hdf5_writer.cpp:834-966,1006-1051`; `tests/CMakeLists.txt:34-46` | task:cdd046b0-61fe-459b-9725-e2b2364d6654; task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae |

## Appendix B. Audit-source task index

- **Scientific model validity:** task:cdd046b0-61fe-459b-9725-e2b2364d6654
- **CUDA/GPU architecture:** task:180ea40c-de07-4423-ada9-8675f07c6e84
- **Software engineering and reproducibility:** task:3969f802-a51c-4b43-8a0d-d2ec122df88a
- **Claim-evidence and identifiability:** task:cb11a8fa-2165-44ed-a7cd-25a87f33c0ae
