# Oxygen respiration is refused in full: the O₂ sink never got the delivery fix

Measurement, not calibration. No new runs: all fifteen supply-ladder arms
(PR #314, 24 h, seed 1001, `48x48x150` at `grid_dx = 2 µm`) already ran with
`oxygen.enabled = true` and `oxygen.metabolic_switch_enabled = true`, and
per-species reaction clips are accounted, so the answer was already on disk.
Reproduce with:

```bash
python experiments/density_limitation/oxygen_audit/audit_oxygen_clip.py <ladder_root>
```

## Result

`clip_O2` is cumulative oxygen refused by the `conc >= 0` clamp after the
explicit reaction update; `demand_O2` is an independent trajectory integral of
`q_maintenance·N + q_consumption·µ·(1−ferm)·N` from the summary series;
`bdry_O2` is all oxygen ever delivered across the epithelial boundary.

| arm | N | clip_O2 (mol) | demand_O2 (mol) | clip/demand | clip/boundary | clip_C |
|---|---:|---:|---:|---:|---:|---:|
| s1_v1 | 28 | 5.724e-12 | 4.074e-12 | 1.41 | 5.2e+03 | 0 |
| s2_v1 | 321 | 3.750e-11 | 2.840e-11 | 1.32 | 3.3e+04 | 0 |
| s4_v1 | 1155 | 1.316e-10 | 1.036e-10 | 1.27 | 1.1e+05 | 0 |
| s8_v1 | 2684 | 3.052e-10 | 2.448e-10 | 1.25 | 2.4e+05 | 0 |
| s16_v1 | 5608 | 7.099e-10 | 5.707e-10 | 1.24 | 4.7e+05 | 0 |

The remaining ten arms (flora 0.5x and 0x) give the same ratios; see the script
output for the full table.

**Essentially the whole respiratory sink is refused.** `clip/demand` sits at
1.24–1.41 across a 200x population range. The ratio exceeding 1 is expected: the
demand integral uses interval-mean `N` and a population-weighted mean µ, so it
under-resolves within-interval growth. The measurement is that the two
quantities are the same size, i.e. the fraction of respiration that reaches the
field is not small but *nil*. For contrast, carbon — which runs in
`uptake_limit = delivery` since #308/#310 — clips exactly zero in the same runs,
so this is oxygen-specific and not a general closure failure.

**The refused amount is 5e3–5e5 times all oxygen the domain ever received**, so
no boundary condition can fund the sink as written.

## The oxygen field is static

Mean O₂ moves from `4.5821e-06` to `4.5864e-06 mol/m³` at N=28 and to
`4.5842e-06` at N=5608 — a 0.04% spread across the whole ladder, and *upward* in
every arm: the field is still equilibrating toward its boundary value while
consumption is invisible. Agents cannot deplete oxygen at all. There is
therefore no O₂ gradient, no anoxic microcolony core, and no density feedback
through oxygen anywhere in the results this project has produced.

## Consequence for Spec 12's phase plane

With a static field, `f_O2 = C/(Km + C) = 4.586e-6/(1e-6 + 4.586e-6) = 0.821`
everywhere and for all time, so the only live term in
`resp_capacity = f_O2 · mu_crit / max(µ, mu_crit)` is µ. The fermentation
fractions reported in #313 and #314 (0.33 at N=28 rising to 0.79 at N=5608) are
**growth-rate driven, not oxygen driven**; describing them as O₂-driven was
wrong and is withdrawn here. `mu_crit`, the aerobic-overflow axis and
`anaerobic_mu_factor` remain untested by any run to date.

## Why this is not a units bug alone

`oxygen.epithelial_conc = 55.0e-6` is annotated `~42 mmHg`, but 42 mmHg of
dissolved O₂ is ~55 µM = `5.5e-2 mol/m³`, so the value and its annotation differ
by 10³ (`Km = 1.0e-6 mol/m³` = 1 nM carries the same offset, which is why the
Monod *ratio* still looks sane while the absolute inventory does not). Fixing
the units does not fix the sink: at `5.5e-2 mol/m³` an 8 fL voxel holds
`4.4e-19` mol while one agent demands `q_maintenance · 60 s = 6e-17` mol per
step, still ~100x the voxel inventory. An explicit pre-diffusion point sink
cannot be funded by within-step resupply — the exact argument that produced the
delivery-limited carbon path in #308 and the total-vs-perturbation correction in
#310.

## What this does and does not license

- It licenses: oxygen needs an `uptake_limit`-style implicit delivery path, and
  the `epithelial_conc`/`Km` units need reconciling against the annotation.
- It does not license: adopting any apical O₂ value (Spec 13's 2–10 mmHg prior
  or the tissue-side 30/39 mmHg electrode values). The boundary is not the
  binding problem, and choosing a number now would look like a fix.
- It does not invalidate the carbon results. The supply ladder is a carbon
  measurement, carbon clip is zero throughout, and the aerobic/anaerobic split
  it reports is a function of µ, which is measured.
