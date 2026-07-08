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
| **carbon** | VBF mucin liberation `vbf.cpp apply_carbon_source`; z-gradient boundary | Agent uptake `fix_metabolism.cpp grow_agent` + `qssa_solver.cpp solve_nutrient_depletion`; **VBF carbon sink** `vbf.cpp apply_carbon_sink` (Spec 5 §1, added) | **Now yes** (was source-only) | Before Spec 5 §1 the only sink was the sparse agent population, so carbon accumulated without bound. The Monod VBF sink closes it. Disabled by default (`vbf_carbon_sink_vmax=0`). |
| **b12** | **VBF B12 production** `vbf.cpp apply_b12_source` (Spec 5 §3, added) | Agent uptake `grow_agent` + `solve_nutrient_depletion` | **Now yes** (was sink-only) | Before Spec 5 §3, B12 had an initial pool but no source → drained to zero → all cells forced onto MetE. Constant source restores homeostasis. Disabled by default (`vbf_b12_production=0`). |
| **iron** | z-gradient boundary; siderophore liberation | VBF first-order sink `apply_iron_sink`; agent uptake | Yes | Closed pre-Spec-5. |
| **oxygen** | Epithelial Dirichlet boundary `chemical_field.cpp apply_boundaries`; z-gradient | **Agent consumption** `qssa_solver.cpp solve_nutrient_depletion:411` (Spec 1); VBF background sink `apply_oxygen_sink` | Yes | **Spec 5 §2 is already implemented** — see §3 below. |
| **acetate** | VBF fermentation `apply_acetate_coupling`; agent overflow | VBF cross-feeding; MetE uptake | Yes | Closed. |
| **mucin** | Goblet secretion `apply_mucin_secretion` | VBF degradation → carbon | Yes | Closed. |
| **bacteriocins** (btuB/fepA/cirA/fhuA) | Producer burst `fix_bacteriocin` → QSSA deposition | First-order decay; diffusion out | Yes | QSSA analytic field, recomputed each step. |

---

## 2. Cross-module coupling map (are the specs talking to each other?)

| Coupling | Producer module | Consumer module | Wired at | Guarded by |
|----------|-----------------|-----------------|----------|------------|
| Agent growth → nutrient depletion | metabolism Fix | ChemicalField reac | `grow_agent`, `solve_nutrient_depletion` | `test_mechanism_wiring` (O2), `smoke` |
| O₂ field → aerobic growth boost | ChemicalField | metabolism Fix | `fix_metabolism` O2 Monod boost | `test_O2_growth_boost` |
| Agents → O₂ depletion | metabolism/QSSA | ChemicalField | `solve_nutrient_depletion:411` | `test_mechanism_wiring::test_o2_consumption_wired` |
| VBF continuum → carbon/iron/O₂/acetate/mucin/B12 | VBF | ChemicalField | `apply_nutrient_coupling` | `test_mucin_liberation`, `test_mechanism_wiring` |
| Toxin field → receptor killing | QSSA bacteriocin | receptor Fix | `fix_receptor` | `test_receptor`, `smoke::test_receptor_killing` |
| Fur (iron) → receptor expression → toxin susceptibility | cell-bio Fur | receptor Fix | `fix_fur`, `fix_receptor` | `test_fur` |
| μ_realized < γ_flow → washout (VADI) | metabolism + advection | Simulation washout | `simulation.cpp check_washout` | `smoke::test_metabolic_washout_*` |
| Density blow-up → halt (Spec 5 §4) | Simulation | run loop | `simulation.cpp run()` dysbiosis check | `test_mechanism_wiring::test_dysbiosis_halt` |

---

## 3. Findings

### 3.1 Spec 5 gaps addressed
- **§1 Carbon sink** — implemented as a Monod-saturating VBF sink
  (`vbf_carbon_sink_vmax`, `vbf_carbon_sink_km`). Off by default.
- **§3 B12 source** — implemented as a constant VBF production term
  (`vbf_b12_production`). Off by default.
- **§4 Dysbiosis threshold** — implemented as a density-based run halt
  (`dysbiosis_threshold`, cells/mL). Off by default.

Defaults keep every new term **disabled** so existing golden/fingerprint tests
(`config_diversity`, `eari-vadi`) are unchanged. To run the biologically
"closed" configuration Spec 5 intends, set the three keys to positive values
(see `docs/PARAMETERS.md`). Recommendation: promote these to non-zero defaults
in a follow-up once the golden references are re-baselined.

### 3.2 Spec 5 §2 (O₂ consumption) was already wired
Spec 5 describes O₂ as "read for the growth boost but never consumed." That is
**stale** — agent O₂ consumption was added in Spec 1 and lives in
`src/diffusion/qssa_solver.cpp:411-417` (`o2_use = q_consumption · max(μ,0)`,
deposited as `-o2_use/cell_vol`). `test_o2_consumption_wired` now proves it is
live (more agents ⇒ lower O₂ field). No code change needed; flagged here so the
spec can be annotated.

### 3.3 ⚠ Open question for the maintainer — possible double-counting of nutrient uptake
Carbon, iron and B12 are each consumed in **two** places every CPU step:

1. `src/fixes/fix_metabolism.cpp` `grow_agent()` — `Δ = d_biomass · yield_x / cell_vol`
   (per-volume, includes `dt` through `d_biomass`).
2. `src/diffusion/qssa_solver.cpp` `solve_nutrient_depletion()` — `Δ = μ·biomass · x_stoichiometry`
   (no `cell_vol` division, no `dt`).

The default `yield_*` and `*_stoichiometry` values are equal (carbon 0.5, iron
1e-6, B12 1e-9), so this looks like the **same** uptake applied twice, with
**inconsistent units** between the two paths (only the O₂ term in
`solve_nutrient_depletion` divides by `cell_vol`). This is not obviously wrong —
they could model distinct processes — but it is the kind of coupling that
"doesn't make logical sense" on its face. **I did not change this**, because
either fix (removing one path, or reconciling units) materially changes core
dynamics and every golden fingerprint. It needs a maintainer decision.

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
