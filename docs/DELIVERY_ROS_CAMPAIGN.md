# Delivery-limited uptake and ROS mortality: consolidated findings

Status: current as of `bee9819` (merge of #339). This document consolidates the
measurement campaigns behind PRs #330–#339 so the numbers, the retired claims,
and the remaining open items survive outside individual PR descriptions.

Campaign scripts and per-arm metrics are committed under
`experiments/delivery_ros/`. The raw HDF5 outputs are not committed; every
number below is reproducible from the scripts plus the recorded commit.

## 1. What was wrong, and what each fix established

Four independent defects were interacting. Each was diagnosed from the ledgers
rather than inferred, and each changed the interpretation of every run that
preceded it.

| PR | Defect | Mechanism |
|---|---|---|
| #331 | Delivery uptake was a bookkeeping quantity | The per-agent sink rate was rescaled to imitate the analytic ceiling. At 2 µm the implicit solve is saturated (`k·dt ≈ 1e5`), so realized removal was nearly independent of the rate: cutting it 4.5× moved removal 28%. Replaced by prescribing the analytic mass inside the solve, so funding equals realized removal identically. |
| #332 | ROS-driven SOS was priced off ambient oxygen | Hazard was `k_ROS × [O₂]_ambient × mu_realized` — oxygen the cell never acquired. With a Robin epithelium the mucus sits at 40–50 µM, so the shipped `k_ROS = 1e2` gives ≈2.5% lysis per cell per 60 s step. Added a funded driver and a four-component SOS hazard ledger. |
| #334 | The Sherwood ceiling read a drained voxel | `4πDrC` is defined against the far field, but `C` came from the agent's own voxel, double-counting the `1/r` depletion the derivation already contains. Added a far-field averaging radius. Per-agent invariance improved 99% → 0.25%, and **population-scale behaviour did not move at all** — see §4. |
| #335, #336 | The prescribed mass was removed from one voxel | The prescribed mass contains no `dx`, but it was deposited into a voxel holding `C·dx³`: ~3–6e-17 mol per cell-step against 2.8e-19 mol in a 2 µm voxel, so negativity was structurally guaranteed and the retry path halved the draw. Deposition was regularized over a 10 µm physical support (default since #336). |
| #338 | Rationing gave up instead of enforcing positivity | Prescribed mass was removed unconditionally and, after two retries, a negative solve was accepted. Mass was conserved, so this was a positivity failure, not an accounting one — the field borrowed against undelivered oxygen. Replaced by local-first reduction with a global bisection guarantee, and a reported per-species rationing factor. |
| #339 | The funded hazard could not consume the calibration | The funded driver normalized per unit biomass (kg/mol); the lysis calibration is per generation (mol⁻¹). Added the absolute funded-flux form. |

## 2. Claims retired by these measurements

Do not reuse the following, in this repository's history or in prior analysis:

- **"Population collapse is metabolic."** It was not. Every agent loss in the
  pre-#332 series was `mortality_lysis` driven by the uncited ambient ROS term.
  Identical configs differing only in `oxygen.k_ROS`: 4 founders die at 3540 s
  with the shipped value and survive 6 h with it at zero; at 80 founders,
  21 versus 76 final cells at 2 µm and 68 versus 171 at 6 µm. ROS accounted for
  **99.5%** of cumulative SOS hazard in the controls.
- **"The maintenance factor sets survival"** (the 1.0 / 4.1 / 8.3 / 15 ladder).
  All four arms were indistinguishable because carbon maintenance shortfall was
  exactly zero in every arm: the brake had nothing to act on.
- **"Fermentative death"** and any fermentation fraction measured before #331.
  Mean realized fermentation fraction on identical configs moved 0.85 → 0.25–0.29
  once funding was honest.
- **"Carrying capacity is 76–171 cells, and grid-dependent."** That range was a
  discretisation floor. With the regularized support it is 201 / 207 / 204 at
  2 / 4 / 6 µm — above every radius-off arm, including the coarsest.
- **Any carbon carrying capacity from #314.** It was measured under
  single-voxel deposition and moves with #335/#336.

## 3. Current measured quantities

All from merged main with the shipped 10 µm support, `oxygen.k_ROS = 0`, no
fitting.

| Quantity | Value | Conditions |
|---|---|---|
| `Q_O2_per_gen` | **1.63e-14 mol** (≈16 fmol) | 2 founders; grid-invariant to 0.4% across 2/4/6 µm; rationing factor 1.0; zero clips; field strictly non-negative |
| Funded O₂ fraction | 0.988 at every resolution | 80 founders, 6 h |
| Carrying capacity | 201 / 207 / 204 cells | 80 founders, 6 h, 2/4/6 µm, 180×180×108 µm patch |
| `T_gen` | ≈2.3 h uncrowded, ≈3.3 h at 80 founders | from completed divisions |
| Growth/maintenance split of funded O₂ | ≈50 / 50 | 2–20 founders |
| Bloom-density funding | 1.5% of the analytic draw, all of it maintenance | 500 agents; growth respiration exactly zero |
| `k_ROS_funded` | 6.2e11 / 1.24e12 / 3.1e12 mol⁻¹ | 1% / 2% / 5% lysis per generation against the measured Q |

`Q_O2_per_gen` is only defined uncrowded: it rises to 2.31e-14 by 20 founders,
and at 80 the delivery rationing factor halves. At bloom density the answer
degenerates to epithelial supply divided by the crowd, so Q must be measured at
healthy mucus density (a handful of cells per patch), not at bloom.

### The patch is not self-sustaining

With ROS off, **every** loss is `outflow_boundary` — 32 of 32 at 2 founders,
462 at 80 — with zero mortality of any kind. `T_gen ≈ 2.3 h` against a 5400 s
radial mucus turnover means a single patch cannot persist without reseeding, so
no patch-scale carrying capacity should be quoted as a persistence result. This
is the washout axis Spec 13 exists to close, and it is now the dominant channel.

### Why `mu_realized` declines during a run

It is transport down the carbon gradient, not a hidden throttle. The often-quoted
"7× decay" was inflated by the step-0 sample, which is the `mu_max` placeholder;
the real decline is **2.0×** (1.62e-4 → 8.2e-5 /s over 12 h at 2–16 cells).

Carbon holds a stationary axial gradient (4.95e-3 mol/m³ at the epithelium to
1.86e-3 in the lumen) against `km_carbon = 5e-3`, so `monod_carbon` runs
0.50 → 0.27 across the domain, and the constant `maintenance_rate = 1e-5`
subtraction sharpens the drop. Radial mucus turnover advects every cell from
z ≈ 24 µm to the outflow boundary, so each cell slides down that gradient.
Within a *single* snapshot, `mu_realized` falls monotonically with depth
(1.45e-4 /s at 48 µm to 9.3e-5 /s at 168 µm), which is what makes it spatial
rather than temporal. Oxygen is flat at 3.8e-2 mol/m³ with the rationing factor
at 1.0, so it is not oxygen.

**Carbon is funded at 99.5% of demand with zero maintenance shortfall** at these
densities. Carbon-limited growth here is therefore set by the field's
*concentration* being comparable to `Km`, not by the delivery path — so a carbon
capacity ladder should sweep the terms that set the gradient (mucin liberation,
the VBF sink) rather than the uptake mode.

## 4. Methodological lessons that cost time

- **Per-agent invariance does not imply population-scale invariance.** #334
  achieved 0.25% funded-fraction agreement across resolutions and moved final
  population by 1 cell (124 vs 123) in a real run, because it corrected the
  concentration the ceiling *reads* and left untouched the voxel the mass is
  *removed from*. Acceptance gates for grid artifacts must be population-scale.
- **Assert both halves of an invariance claim.** An invariance test passes
  trivially against an inert code path; the contrast (radius off reproducing the
  old numbers exactly) is what proves the path is live.
- **A conserved-mass check does not detect a positivity failure.** #338's field
  was mass-exact and still negative over two-thirds of cells under stress.
- **The retry/rationing ledger is the discriminator.** Both #335 and #338 were
  diagnosed from withheld mass and retry counts, not from concentrations.
  Retry-count evidence does not generalize across campaigns, however: outside
  the 80-founder campaign, retries saturate at the per-solve ceiling and are
  identical at both radii.
- **Coefficients calibrated against a grid-dependent divisor bake in the
  artifact.** `Q_O2_per_gen` was 3.3× grid-dependent before #335/#336/#338;
  pinning the lysis coefficient against it then would have been the same class
  of error as pricing ROS off ambient oxygen.

## 5. Open items

- **Shipped defaults do not reflect these findings.** `oxygen.ros_driver`
  defaults to `ambient` with the uncited `k_ROS = 1e2`, i.e. a fresh run using
  defaults still reproduces the retired mortality; and
  `metabolism.uptake_limit` defaults to `none`, so the whole funded-delivery
  path — including the 10 µm support — is off unless a config opts in. Both are
  scientific decisions rather than code defects, and both are unresolved.
- **`k_ROS_funded` has no chosen value.** The 1/2/5% per-generation targets are
  in-vitro priors (1–8% colicin-expressing fraction; ~0.9% spontaneous
  SOS-positive) applied out of domain: there is no direct measurement of
  spontaneous colicin-mediated lysis in gut mucus. The intended validation is
  in-vivo colicin-mediated displacement kinetics at Spec 13 scale, so the
  coefficient should be an *output* of the metapopulation work, not an input.
  Respiratory ROS is also only one induction route — bile salts, nutrient
  stress, and SOS-independent operon regulation are unrepresented — so a fitted
  coefficient stands in for several routes, and real lysis is expected to be
  patchy.
- **`Q_O2_per_gen` averages over a non-stationary `T_gen`** (§3), since cells
  advect down the carbon gradient during the measurement window.
- **Re-run the #314 carbon ladder** on the current code, sweeping the
  gradient-setting terms.
- **CUDA parity** for delivery-limited uptake is covered by an explicit refusal
  rather than an implementation.
