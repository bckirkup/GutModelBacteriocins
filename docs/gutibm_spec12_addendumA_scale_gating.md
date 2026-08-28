# GutIBM Spec 12 — Addendum A: Scale-Gating of Agent–Field Coupling

**Status**: Analysis complete, implementation proposed. Supersedes the Spec 12
characterization of `agent_carbon_coupling` as a "no-op".
**Date**: 2026-08-28
**Depends on**: `gutibm_spec12_density_limitation.md`
**Related**: `gutibm_spec13_multiscale.md` (the patch grid is the scale this analysis implicates)

---

## 1. Why this addendum exists

Spec 12 Change 1 added an agent-density-coupled VBF carbon sink. Testing found it
had effectively no influence on the carbon field: diffusion refills a voxel faster
than the agents in it can drain it. We recorded this as "`agent_carbon_coupling`
is a no-op."

That conclusion is correct as an observation and wrong as an explanation. The
coupling is not broken and it is not mis-scaled. It is **scale-gated**: the
physics that makes it negligible at our current resolution is the same physics
that makes it dominant a factor of ten larger in radius. This addendum quantifies
the gate, and — more importantly — identifies that our previous reasoning about
*why* the coupling failed was wrong in a way that changes what we should do about it.

The practical output is (a) a corrected statement for the paper, (b) a runtime
diagnostic so the model self-reports when the coupling matters instead of us
guessing, and (c) a specific reason not to spend effort "fixing" the coupling.

---

## 2. The controlling dimensionless group

For a substrate consumed by a population of cells and supplied by diffusion, the
ratio of reaction rate to diffusive supply over a length scale `R` is the
Damköhler number:

```
Da = rho * q_cell * R^2 / (D * C_inf)
```

| symbol  | meaning                              | units        |
|---------|--------------------------------------|--------------|
| `rho`   | **local** cell number density        | cells/m^3    |
| `q_cell`| per-cell uptake rate (saturated)     | mol/cell/s   |
| `R`     | radius of the populated region       | m            |
| `D`     | substrate diffusivity in mucus       | m^2/s        |
| `C_inf` | substrate concentration outside      | mol/m^3      |

Setting `Da = 1` gives a characteristic radius:

```
R_50 = sqrt( D * C_inf / (rho * q_cell) )
```

**The `R^2` is the entire story.** `Da` is only linear in density and in per-cell
uptake, but quadratic in the size of the populated region. Any argument about
whether agent coupling "works" that does not specify a length scale is not a
well-posed argument.

### 2.1 What `Da = 1` actually means

We verified the closed form against a nonlinear Michaelis–Menten
reaction–diffusion solve in spherical geometry with proper exterior matching
(`c'(1) = 1 - c(1)`), damped Newton iteration, and grid-refinement checks
(worst relative difference between the two finest runs: 2e-10).
task:2e313c14-d016-46db-a920-5c2933d935bf

At `Da = 1`:
- centre depletion `1 - C(0)/C_inf` = **48%**
- effectiveness factor `eta` = **0.985**

So at the nominal "crossover" the aggregate is still consuming 98.5% of what it
would consume in a perfectly mixed medium. **`R_50` is a half-depletion radius,
not a threshold for anything.** An anoxic core appears at roughly `1.8 * R_50`.
Do not write `Da = 1` in the paper as though it were an onset condition.

Zero-order kinetics are an excellent approximation here (error −1.3% at R10,
−2.1% at R50, −5.3% at R90) because `K_m << C_inf` through most of the
oxygenated mucus. The Monod correction only becomes decisive below ~0.2 µM.

**Rule of thumb for the spec text**:

> `R_50 ~ sqrt(D*C_inf/(rho*q_cell))` is the radius at which a uniformly
> populated aggregate depletes its own centre by about half, valid when
> `K_m << C_inf`. It is not an anoxic-core threshold; the core goes anoxic
> at roughly `1.8x` that radius.

---

## 3. The corrected diagnosis: it is geometry, not density

This is the part that changes our interpretation, and it is worth stating plainly
because it contradicts what we assumed in Spec 12.

We had been reasoning as though the coupling was negligible because our
populations are too *dilute* — healthy culture density is 1e4–1e5 CFU/mL, which is
less than one cell per patch, so of course a single agent cannot move the field.

**That reasoning is wrong.** Compute the local density a single agent implies at
our grid resolution:

| `grid_dx` | voxel volume | local density of ONE agent |
|-----------|--------------|----------------------------|
| 2 µm      | 8.0e-18 m^3  | **1.25e11 cells/mL**       |
| 5 µm      | 1.25e-16 m^3 | **8.0e9 cells/mL**         |

A single agent in a 2 µm voxel *already* represents the dense-microcolony packing
regime (1e11 cells/mL). Density is not the problem — at voxel resolution we are
**already** at the densities where the literature reports oxygen limitation.

The problem is purely that `R = dx/2 = 1 µm` is tiny, and `Da` goes as `R^2`:

```
Da(single agent, dx = 2 µm) = 0.003
Da(single agent, dx = 5 µm) = 0.001
```

At that same local density, `R_50 = 19.6 µm` — about **10 voxels in radius**,
requiring roughly **3,900 agents** in a coherent contiguous blob.

**So the correct statement is:** the sink per agent is the right magnitude; what
is missing is *spatially coherent occupancy over ~10 voxels*. The gate is a
geometric coherence requirement, not a density requirement. This matters because
the two diagnoses imply opposite fixes:

- "too dilute" → raise `agent_carbon_coupling`, raise carrying capacity, or
  increase per-agent uptake. **All of these are wrong** and would produce
  unphysical per-cell metabolism to compensate for a geometry deficit.
- "not spatially coherent over enough voxels" → the coupling coefficient is
  fine; either resolve larger aggregates (Spec 13) or accept the continuum VBF
  sink as the correct sub-crossover representation. **Change no parameters.**

Devin: this is the reason Addendum A adds a diagnostic and *not* a parameter
retune. Please do not "fix" the coupling by scaling `agent_carbon_coupling` up.

---

## 4. Numerical anchors

Literature-anchored central values. task:1114ae84-1eb2-4f68-88f3-378a90ee6d96

### 4.1 Recommended parameters

| Parameter | Central | Range | Notes |
|---|---|---|---|
| `D_O2` in mucus | 2.0e-9 m^2/s | 1.0–2.1e-9 | O2 is small and neutral; mucus does **not** meaningfully retard it. Do not impose a 10–100x retardation factor. |
| `D_glucose` in mucus | 0.6e-9 m^2/s | — | Glucose crosses 5 wt% mucin gels; prefer this to a large viscosity penalty. |
| `C_O2` at epithelial surface | 51 µM (40 mmHg) | — | Mouse cecal tissue surface. |
| `C_O2` in lumen | 1.3 µM (1 mmHg) | 0.64–3.8 µM | Ascending colon ~14 µM; distal colon 0.64–3.8 µM. |
| O2 decay length in mucus | 100 µm | 50–300 µm | No universal measured value; consistent with 10–125 µm bacterial zones. |
| `K_m` O2, cytochrome *bo* | 0.20 µM | — | |
| `K_m` O2, cytochrome *bd* | 0.025 µM | — | High-affinity; expression shifts to *bd* microaerobically. |
| `q_O2` aerobic | 4.75 mmol/gDW/h | 4.75–20 | 20 is max respiratory capacity, not typical. |
| `m_DW` per cell | 0.30 pg | — | |
| `q_cell` O2 | 8.33e-19 mol/cell/s | ~1e-19–2e-18 | **Dominant uncertainty.** |
| `C_glucose` free monosaccharide | 0.10 mM | 0.01–1 mM | **Provisional and poorly constrained.** Mucin is polymeric; its total carbohydrate cannot be used as dissolved glucose. |
| `K_m` glucose (whole-cell) | scan 1–100 µM | — | Transporter composition adapts under limitation. |

### 4.2 Crossover radii

Nonlinear Michaelis–Menten solve, spherical geometry:

| Local density | O2 R10 / R50 / R90 (µm) | O2 dead-core (µm) | Carbon R50 (µm) |
|---:|---:|---:|---:|
| 1e9 cells/mL  | 99 / 224 / 311 | 393 | 276 |
| 1e10 cells/mL | 31 / 71 / 98   | 124 | 87  |
| 1e11 cells/mL | 9.9 / 22.4 / 31.1 | 39.3 | 27.6 |
| 1e12 cells/mL | 3.1 / 7.1 / 9.8 | 12.4 | 8.7 |

Independent literature-anchored recomputation gives O2 `R_50` = 169 µm / 16.9 µm /
5.34 µm at 1e9 / 1e11 / 1e12 cells/mL — within ~25% of the solver values. The
difference is driven by choice of `D`, `C`, dry mass, and OUR, **not** by `K_m`.

Planar (slab) geometry, which matches our mucus layer with the epithelium as a
face source — O2 penetration depth:

| Density | Penetration (µm) |
|---:|---:|
| 1e9  | 431 |
| 1e10 | 136 |
| 1e11 | 43  |
| 1e12 | 14  |

The 1e11 value of 43 µm falls inside the FISH-anchored 20–60 µm founder placement
band, which is a useful independent consistency check on the model geometry.

### 4.3 Two counter-intuitive results worth remembering

**Carbon is an *easier* sink to satisfy than oxygen, not harder.** Carbon `R_50`
runs ~23% *larger* than O2 at every density (276 vs 224 µm at 1e9 cells/mL). The
lower diffusivity of glucose is more than offset by its higher ambient
concentration relative to per-cell demand. So the original carbon no-op is *not*
explained by carbon being a weak sink — both species share essentially the same
geometric gate. This strengthens the scale-gating account rather than weakening it.

**`K_m` is nearly irrelevant; `q_cell` dominates.** One-factor-at-a-time spans in
`R_50`: `q_cell` 10.0x, `C_inf` 3.09x, `D` 1.45x, `K_m` 1.19x. Sensitivity effort
belongs on per-cell physiology (is the cell growing, carbon-limited, or in
maintenance?) and on the local boundary concentration. Refining mucus diffusivity
is not worth the effort.

---

## 5. Implementation: a diagnostic, not a parameter change

The goal is for the model to **report** its own local `Da` so we stop arguing
about whether coupling matters and instead read it off the output.

### 5.1 Change A1 — per-voxel Damköhler diagnostic

**What**: For each field with agent coupling enabled, compute and output a local
`Da` estimate per voxel, using a *coherence-aware* length scale rather than `dx`.

The naive choice `R = dx/2` is what makes the number uninformative — it always
reports `Da << 1` by construction. Instead, estimate the local coherent aggregate
radius from the occupancy of a neighbourhood:

```
// For voxel i, over a cubic neighbourhood of half-width W voxels (W = 5 default,
// so 11^3; chosen because R_50 ~ 10 voxels at 2 um for our densities):
//   N_nb    = total agents in neighbourhood
//   V_nb    = (2W+1)^3 * dx^3
//   rho_loc = N_nb / V_nb                      [cells/m^3]
// Equivalent-sphere radius of the OCCUPIED volume, not the neighbourhood volume:
//   V_occ   = N_occupied_voxels_in_nb * dx^3
//   R_coh   = cbrt( 3 * V_occ / (4 * pi) )     [m]
// Then:
//   Da_loc  = rho_loc * q_cell_eff * R_coh^2 / (D_eff * C_local)
```

`q_cell_eff` must be the *realized* per-agent uptake actually applied by
`FixMetabolism` / the O2 respiration term in `qssa_solver.cpp`, averaged over the
neighbourhood — not the config maximum. If we report `Da` using `q_max` while the
model applies a Monod-limited rate, the diagnostic will overstate coupling
precisely where it matters. Please read the applied value.

For oxygen, `q_cell_eff` should be taken as
`oxygen.q_consumption * max(mu_realized, 0) + oxygen.q_maintenance`, matching
`src/diffusion/qssa_solver.cpp:496-497`.

**Files to modify**:
- `src/diffusion/qssa_solver.cpp` — after the existing agent-respiration loop
  (around line 480–503), accumulate per-voxel agent counts and summed realized
  `q`. The loop already walks agents and has `a.grid_cell`; extend it rather than
  adding a second pass.
- `src/diffusion/qssa_solver.h` — expose `const std::vector<Real>& damkohler_o2()`
  and `damkohler_carbon()`.
- `src/io/hdf5_writer.cpp` — write `damkohler_o2` / `damkohler_carbon` fields
  alongside existing per-voxel field output, and per-interval scalar summaries:
  `da_max`, `da_p99`, `da_mean_occupied`, and
  `frac_voxels_da_gt_0p1` (fraction of *occupied* voxels with `Da_loc > 0.1`).
- `src/io/chem_environment_config.h` — add a small config struct:
  ```
  struct DamkohlerDiagnosticConfig {
    bool enabled = false;          // OFF by default
    Int  neighbourhood_half_width = 5;   // W
    bool output_field = false;     // per-voxel field is large; scalars by default
  };
  ```
- `src/io/input_parser.cpp` — parse `damkohler.enabled`,
  `damkohler.neighbourhood_half_width`, `damkohler.output_field`.

**Cost control**: the per-voxel field is the same size as a chemical field. Default
`output_field = false` and emit only the scalar summaries. The neighbourhood sum is
a separable box filter — please implement it as three 1-D passes
(O(N) not O(N*W^3)); a naive triple loop at W=5 is 1,331 reads per voxel and will
dominate step time.

**Test**:
- Uniform dilute population, 1 agent per 2 µm voxel over a 4-voxel-radius blob:
  `da_max < 0.05`.
- Same density over a 15-voxel-radius blob: `da_max > 0.5`.
- **This pair is the actual regression test for the whole addendum** — same
  density, different extent, order-of-magnitude different `Da`. If both report
  similar `Da`, the coherence length is not being computed correctly and the
  diagnostic is worthless.
- With `damkohler.enabled = false`, output must be byte-identical to current `main`.

### 5.2 Change A2 — aggregate-radius sweep harness

**What**: A campaign that initializes a spherical clump of agents of specified
radius and density, runs briefly with fields on and growth/motility/division OFF,
and records centre-vs-edge field concentration.

This measures the model's *actual* crossover and validates it against the
analytical prediction. It is the experiment that converts "we think coupling is
scale-gated" into "the model's coupling turns on at radius X, as predicted."

**Config**: reuse the existing single-colony example structure
(`examples/single_colony/`). Add `init.clump_radius` and `init.clump_density`
if not already expressible; otherwise generate the agent list in Python.

**Sweep**: `clump_radius` in {2, 5, 10, 20, 40, 80} µm at
`clump_density` in {1e10, 1e11, 1e12} cells/mL, with
`agent_carbon_coupling` and O2 respiration ON, `vbf_sink` at its calibrated value.

**Metric**: `C(centre)/C(edge)` for O2 and carbon vs `clump_radius`. Expect a
sigmoid crossing 0.5 near the predicted `R_50` for each density.

**Expected outcome**: crossing near 22 µm at 1e11 cells/mL and near 7 µm at 1e12.
If the model's crossing is far from prediction, that is a real finding about the
field solver (likely the ordering/implicitness issue we already hit with the VBF
O2 sink) and should be reported rather than tuned away.

**Important**: growth and division must be OFF. If agents divide during the run,
the clump geometry changes and the measurement is meaningless.

### 5.3 What NOT to do

- Do **not** increase `agent_carbon_coupling` to make the effect visible.
- Do **not** increase per-cell uptake rates beyond the literature range in
  section 4.1 to force depletion at small radii.
- Do **not** remove the continuum VBF sink. Below the crossover it is the
  *correct* physics: the anaerobic majority is genuinely a distributed sink and
  representing it as a continuum is right, not a compromise.

---

## 6. Revised language for Spec 12 and for the paper

Replace the "no-op" characterization with:

> Agent-resolved coupling into the chemical fields is scale-gated by the
> Damköhler number, which grows as the square of the coherent aggregate radius.
> At the model's voxel resolution a single agent already represents local packing
> of order 1e11 cells/mL, so the coupling is not limited by cell density; it is
> limited by spatial coherence, requiring sustained occupancy over roughly ten
> voxels (~20 µm, ~4,000 agents) before the local sink measurably perturbs the
> field. Below that scale the continuum background sink is the appropriate
> representation, and explicit per-agent coupling is expected to be inert. Above
> it — the regime of mucosal microcolonies and of the Spec 13 patch grid — the
> coupling governs local niche construction.

This is a stronger claim than "the coupling doesn't work" and it is falsifiable by
the Change A2 sweep.

---

## 7. Honest limitations

- The analysis is **steady-state**. Real mucus has peristaltic mixing and
  turnover on timescales that may not permit a steady gradient to establish.
  Contraction cycles in Spec 13 Layer 2 could reset gradients faster than they form.
- **Uniform density within the aggregate.** Real microcolonies are heterogeneous.
- **Single limiting substrate at a time.** Cross-coupling (O2 limitation changing
  carbon yield via the metabolic switch) is exactly what Spec 12 Change 2
  introduces, and is not in this analysis.
- The 2026 microfluidic study (Scheidweiler et al.) supports the Damköhler
  mechanism **qualitatively only**. It does not report local packing density, a
  calibrated per-cell uptake rate, or a clonal aggregate radius, so it cannot be
  used to validate our specific crossover. Do not cite it as quantitative
  confirmation.
- `C_glucose` is the weakest input in the whole analysis. The carbon crossover
  numbers should be treated as provisional until free monosaccharide
  concentrations in mucus are pinned down.

---

## 8. Priority

| Change | Effort | Value | Priority |
|---|---|---|---|
| A1 diagnostic (scalars only) | Low | High — ends the argument, self-reporting | **1** |
| A2 sweep harness | Medium | High — validates the crossover | **2** |
| A1 per-voxel field output | Low | Medium — useful for figures | 3 |

Both changes are additive and default-off. No existing behaviour changes.
