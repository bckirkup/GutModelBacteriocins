# Mechanism Wiring Audit

**Question this answers:** *Are the mechanisms from all the specs (1–5, EARI,
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
| **carbon** | VBF mucin liberation `vbf.cpp apply_carbon_source`; z-gradient boundary | **Agent uptake `fix_metabolism.cpp grow_agent` (canonical, Spec 6)**; **VBF carbon sink** `vbf.cpp apply_carbon_sink` (Spec 5 §1) | Yes | Spec 6 makes the metabolism Fix the single per-agent uptake site (the duplicate `solve_nutrient_depletion` term was removed). Monod VBF sink now **active by default** (`vbf_carbon_sink_vmax=5.5e-5`) → bounded ~1 mM equilibrium. |
| **b12 / corrinoid** | none (constant field) | **none — not consumed (Spec 6 §3)** | N/A (constant pool) | Spec 6: the B12 field represents the total bioavailable corrinoid pool (~1 µM), produced by the anaerobic majority far faster than E. coli demand. Neither produced nor depleted; pinned at `1e-6`. Replaces the Spec 5 `vbf_b12_production` source (removed). |
| **iron** | z-gradient boundary; siderophore liberation | VBF first-order sink `apply_iron_sink`; **agent uptake `grow_agent` (canonical, Spec 6)** | Yes | Spec 6: uptake consolidated to the metabolism Fix (duplicate `solve_nutrient_depletion` term removed). |
| **oxygen** | Epithelial Dirichlet boundary `chemical_field.cpp apply_boundaries`; z-gradient | **Agent consumption** `qssa_solver.cpp solve_nutrient_depletion` (Spec 1; the sole remaining per-agent term there after Spec 6); VBF background sink `apply_oxygen_sink` | Yes | Already wired in Spec 1 — see §3 below. |
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
| μ_realized < γ_flow → washout (VADI) | metabolism + advection | Simulation washout | `simulation.cpp check_washout` | `smoke::test_metabolic_washout_*` |
| Density blow-up → halt (Spec 5 §4) | Simulation | run loop | `simulation.cpp run()` dysbiosis check | `test_mechanism_wiring::test_dysbiosis_halt` |

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

### 3.2 Spec 5 §2 (O₂ consumption) was already wired
Spec 5 describes O₂ as "read for the growth boost but never consumed." That is
**stale** — agent O₂ consumption was added in Spec 1 and lives in
`src/diffusion/qssa_solver.cpp` `solve_nutrient_depletion` (`o2_use =
q_consumption · max(μ,0)`, deposited as `-o2_use/cell_vol`). After Spec 6 this
is the **only** per-agent term in that function. `test_o2_consumption_wired`
proves it is live (more agents ⇒ lower O₂ field).

### 3.3 ✅ Resolved by Spec 6 — the nutrient double-count
Carbon, iron and B12 were previously consumed in **two** places every CPU step:

1. `src/fixes/fix_metabolism.cpp` `grow_agent()` — `Δ = d_biomass · yield_x / cell_vol`
   (per-volume, includes `dt` through `d_biomass`).
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
     carbon sink and B12 source active — asserts no tracked species goes NaN,
     negative, or explodes.

2. **Directional coupling sensitivity.**
   - `test_carbon_sink_sensitivity`, `test_b12_source_directional`: larger
     Vmax removes more carbon; larger production adds more B12; the term is
     exactly zero when the key is zero (dead-wiring guard).
   - `test_o2_consumption_wired`: more agents ⇒ lower mean O₂ field (agent
     metabolism actually feeds back into chemistry through the full step loop).
   - `test_dysbiosis_halt`: the threshold halts the run early and is inert when
     `0` (compared against a control run).

Every new config key also has an ingestion probe in
`tests/test_config_ingestion.cpp`, whose completeness guard fails CI if any
parsed key lacks a probe.

### Extending this audit
When adding a new species or coupling: add a row to §1/§2, then add (a) a
bounded-steady-state assertion and (b) a directional-sensitivity assertion in
`test_mechanism_wiring.cpp`. A coupling without both is considered unverified.
