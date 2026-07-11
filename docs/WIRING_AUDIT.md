# Mechanism Wiring Audit

**Question this answers:** *Are the mechanisms from all the specs (1–7, EARI,
VADI) wired together in a way that makes logical sense, and how do we build
tests to be sure?*

This document is a systematic pass over every coupled quantity in the model.
For each species/mechanism it records the **source(s)**, **sink(s)**, and
**feedback** wiring points (with `file:line`), whether the loop is *closed*
(bounded steady state reachable) or *open* (source-only or sink-only), and how
the wiring is now guarded by tests.

The guiding principle (from `.agents/skills/ci-test-design`): a mechanism is
"wired" only if we can show both **(a)** a golden/bounded outcome for a fixed
config and **(b)** directional sensitivity — turning the coupling on/up moves
the downstream quantity in the biologically correct direction. A config key
that parses but never changes an outcome is *dead wiring*.

---

## 1. Chemical species mass-balance map

| Species | Source(s) | Sink(s) | Loop closed? | Notes |
|---------|-----------|---------|--------------|-------|
| **carbon** | VBF mucin liberation `vbf.cpp apply_carbon_source`; z=0 boundary | **Agent uptake `fix_metabolism.cpp grow_agent` (canonical, Spec 6)**; **VBF carbon sink** `vbf.cpp apply_carbon_sink` (Spec 5 §1) | Yes | Implicit diffusion transports local source/sink changes with `D=5e-10 m²/s`; the metabolism Fix is the single per-agent uptake site. |
| **b12 / corrinoid** | z=0 boundary | **none — not consumed (Spec 6 §3)** | N/A (constant pool) | The ~1 µM bioavailable pool remains spatially uniform under implicit diffusion (`D=5e-10 m²/s`). |
| **iron** | z=0 boundary; siderophore liberation | VBF first-order sink `apply_iron_sink`; **agent uptake `grow_agent` (canonical, Spec 6)** | Yes | Implicit diffusion (`D=7e-10 m²/s`) couples local depletion to the epithelial supply. |
| **oxygen** | Epithelial Dirichlet boundary `chemical_field.cpp apply_boundaries` | **Agent Pirt respiration** `qssa_solver.cpp solve_nutrient_depletion`; first-order VBF background sink `apply_oxygen_sink` | Yes | Implicit diffusion (`D=2.1e-9 m²/s`) replaces checkerboard-prone local-only updates and preserves a smooth boundary-fed gradient. |
| **acetate** | VBF fermentation `apply_acetate_coupling`; agent overflow | VBF cross-feeding; MetE uptake | Yes | Closed. |
| **mucin** | Goblet secretion `apply_mucin_secretion` | VBF degradation → carbon | Yes | Closed. |
| **bacteriocins** (btuB/fepA/cirA/fhuA) | Producer burst `fix_bacteriocin` → QSSA deposition | First-order decay; diffusion out | Yes | QSSA analytic field, recomputed each step. |

---

## 2. Cross-module coupling map (are the specs talking to each other?)

| Coupling | Producer module | Consumer module | Wired at | Guarded by |
|----------|-----------------|-----------------|----------|------------|
| Agent growth → carbon/iron depletion | metabolism Fix | ChemicalField reac | `grow_agent` (yield-based, canonical) | `test_mechanism_wiring`, `smoke` |
| O₂ field → aerobic growth boost | ChemicalField | metabolism Fix | `fix_metabolism` O2 Monod boost | `test_O2_growth_boost` |
| Agents → O₂ depletion | QSSA | ChemicalField | `solve_nutrient_depletion` | `test_mechanism_wiring::test_o2_consumption_wired`, `test_qssa_stoichiometry` |
| VBF continuum → carbon/iron/O₂/acetate/mucin | VBF | ChemicalField | `apply_nutrient_coupling` | `test_mucin_liberation`, `test_mechanism_wiring` |
| Toxin field → receptor killing | QSSA bacteriocin | receptor Fix | `fix_receptor` | `test_receptor`, `smoke::test_receptor_killing` |
| Fur (iron) → receptor expression → toxin susceptibility | cell-bio Fur | receptor Fix | `fix_fur`, `fix_receptor` | `test_fur` |
| μ_realized < γ_flow → washout (VADI) | metabolism + advection | Simulation washout | `simulation.cpp check_washout` | `smoke::test_metabolic_washout_*` | `washout_trap` (issue #160) |
| Density blow-up → halt (Spec 5 §4) | Simulation | run loop | `simulation.cpp run()` dysbiosis check | `test_mechanism_wiring::test_dysbiosis_halt` |
| Local nutrient reactions → neighboring grid cells (Spec 7) | metabolism/QSSA/VBF | ChemicalField | `sum_reactions_across_ranks` → `apply_diffusion` | `test_nutrient_diffusion`, `test_mpi_multi_rank`, `gpu_smoke` |

---

## 3. Findings

### 3.1 Spec 5 gaps addressed (updated by Spec 6)
- **§1 Carbon sink** — Monod-saturating VBF sink (`vbf_carbon_sink_vmax`,
  `vbf_carbon_sink_km`). Spec 6 **activates it by default** (`5.5e-5`, just
  above mucin liberation) → carbon reaches a bounded ~1 mM equilibrium.
- **§3 B12 source** — **removed by Spec 6**. The corrinoid field is now a
  constant ~1 µM pool (neither produced nor consumed), so a source term is no
  longer meaningful; `vbf_b12_production` was deleted.
- **§4 Dysbiosis threshold** — density-based run halt (`dysbiosis_threshold`,
  cells/mL). Off by default.

Spec 6 re-baselined the `eari-vadi` golden references (see
`python/tests/fixtures/eari_vadi_ci*_golden.json`) because activating the carbon
sink and correcting the corrinoid pool shift the reference metrics.

### 3.2 O₂ consumption — wired in Spec 1, made density-tracking post-Spec 6
Agent O₂ consumption lives in `src/diffusion/qssa_solver.cpp`
`solve_nutrient_depletion` and after Spec 6 is the **only** per-agent term there.
Two defects made it fail to track agent *density* (present-but-non-growing cells
still respire), which surfaced once the population began collapsing (6→3→1
cells) yet O₂ stayed pinned:

1. **Growth-only respiration.** The term was `o2_use = q_consumption · max(μ,0)`,
   so a non-growing cell (μ→0, e.g. washing out) consumed ~zero O₂ and the field
   stopped tracking cell count. Fixed with a **Pirt** model: `o2_use =
   q_consumption · max(μ,0) + q_maintenance`, where `q_maintenance` is a basal,
   density-coupled term applied per living cell regardless of growth.

2. **Zero-order VBF background sink.** `apply_oxygen_sink` removed a constant
   `vbf_sink` (mol/m³/s) independent of `[O₂]`. At the default scale that
   removes more O₂ than is present in one 60 s bio step, hard-zeroing the entire
   interior and pinning the mean at the epithelial Dirichlet layer — so the
   background sink alone dictated the field and per-agent respiration was masked
   entirely. Fixed by making it **first-order** (`reac −= vbf_sink · [O₂]`,
   `vbf_sink` now a 1/s rate constant), mirroring the iron sink; the sink is now
   self-limiting so a smooth gradient survives and agent respiration is visible
   on top of it.

`test_o2_consumption_wired` proves agent respiration is live with the background
sink off; `test_o2_tracks_density_over_background` proves the mean O₂ still falls
with agent count with the (first-order) background sink **on**;
`test_qssa_stoichiometry::test_o2_maintenance_tracks_density` proves a
non-growing cell respires and N cells draw N× the maintenance flux.

### 3.3 Spec 7 nutrient transport — implicit and rank-consistent

The Spec 7 explicit-stencil proposal used an incorrect CFL estimate. At the model defaults, O₂ has `D·dt/dx² = 5040`, not `0.005`; a stable 3-D explicit solve would require roughly 30,000 substeps per 60 s biological step. Nutrient transport therefore uses backward-Euler directional splitting: periodic x/y cyclic tridiagonal solves, epithelial Dirichlet z=0, and a luminal zero-flux z face. This is L-stable and O(cells) per species; concentrations are clamped nonnegative after the solve. Configured exponential z-gradients remain prescribed background profiles while diffusion smooths local departures from them.

Every rank stores the full chemical grid but contributes reactions from only its local agents. `ChemicalField::sum_reactions_across_ranks` performs an `MPI_Allreduce` before the identical global VBF terms are added, so all ranks diffuse the same global field without double-counting VBF. CUDA keeps metabolism on device, transfers reaction grids to the host for rank summation and diffusion, then synchronizes concentrations back.

### 3.4 ✅ Resolved by Spec 6 — the nutrient double-count
Carbon, iron and B12 were previously consumed in **two** places every CPU step:

1. `src/fixes/fix_metabolism.cpp` `grow_agent()` — `rate = d_biomass · yield_x / (cell_vol · dt)`.
2. `src/diffusion/qssa_solver.cpp` `solve_nutrient_depletion()` — `Δ = μ·biomass · x_stoichiometry`
   (no `cell_vol` division, no `dt`).

Because the two paths used **inconsistent units** (the metabolism term is
`≈ dt/cell_vol ≈ 5e17×` larger), the metabolism Fix term dominated and the QSSA
nutrient terms were numerically negligible — the same uptake applied twice.

**Spec 6 resolution:** the **metabolism Fix is the single canonical per-agent
uptake site** for carbon and iron. The carbon/iron/B12 terms in
`solve_nutrient_depletion` (and the GPU `nutrient_depletion_kernel`) were
removed; only the dimensionally-clean O₂ respiration term remains there. B12 is
no longer consumed anywhere (constant corrinoid pool, §3.1). The former
`iron/b12/carbon_stoichiometry` knobs on `QSSAConfig` were deleted. See
`docs/PARAMETERS.md` → "Nutrient Cycle (Spec 6)".

---

## 4. How the wiring is tested (`tests/test_mechanism_wiring.cpp`)

Two invariants, applied to the couplings above:

1. **Mass balance / bounded steady state.**
   - `test_carbon_sink_bounds_accumulation`: with a constant carbon source and
     no sink, carbon grows to `source·dt·N` (unbounded); with the sink it
     reaches the analytic Monod steady state `km·S/(vmax−S)`. Catches the
     "source without sink" class directly.
   - `test_all_species_bounded_steady_state`: full `Simulation` run with O₂,
     the carbon sink, and the constant corrinoid pool active — asserts no
     tracked species goes NaN, negative, or explodes.

2. **Directional coupling sensitivity.**
   - `test_carbon_sink_sensitivity`: larger Vmax removes more carbon and the
     term is exactly zero when the key is zero (dead-wiring guard).
   - `test_corrinoid_field_constant`: agent metabolism leaves the prescribed
     ~1 µM corrinoid pool unchanged.
   - `test_o2_consumption_wired`: more agents ⇒ lower mean O₂ field (agent
     metabolism actually feeds back into chemistry through the full step loop).
   - `test_dysbiosis_halt`: the threshold halts the run early and is inert when
     `0` (compared against a control run).

3. **Transport regression and execution parity.**
   - `test_nutrient_diffusion` checks a quantitative point-source golden profile, uniform-field and configured-gradient invariance, boundary gradients, nonnegativity at `dt=300 s`, and enable/coefficient sensitivity.
   - `test_mpi_multi_rank` checks that rank-local reactions sum to one identical diffused field on both ranks.
   - `gpu_smoke` includes chemical concentrations in the CPU/GPU fingerprint, guarding GPU diffusion wiring in `module_chemistry` (CPU fallback when line lengths exceed the PCR limit).

Every new config key also has an ingestion probe in
`tests/test_config_ingestion.cpp`, whose completeness guard fails CI if any
parsed key lacks a probe.

### Extending this audit
When adding a new species or coupling: add a row to §1/§2, then add (a) a
bounded-steady-state assertion and (b) a directional-sensitivity assertion in
`test_mechanism_wiring.cpp`. A coupling without both is considered unverified.
