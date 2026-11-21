# Spec 13 Phase 3 — Layer-3 colonic gradient readout

Three-region chain (cecum/ascending → transverse → descending/sigmoid →
stool), each region a `MucusSegment`, coupled by a count-based luminal
compartment with per-lineage BI-genotype resolution and anaerobic growth.
Deliverable gates: `layer3.enabled` (off by default), `layer3.uniform_profile`
(gradient-collapsed sanity arm), `layer3.checkpoint_file` (luminal counts +
genotypes + per-region state round-trip).

Sources of parameterization, per the design contract:

- **Within-region**: Phase-1 admissible corner
  (`analysis/spec13_occupancy/REPORT.md`) — crypt fraction 0.10 / exposed
  0.45+0.45, disruption 0.07/0.5/0.5, single-cell loss 0.02/0.10/0.20 h⁻¹,
  supply 2.0/1.0/0.5, contraction ≤ 0.6 min⁻¹, agent loss ≤ 0.3,
  transit half-life 4 h, reattach 0.05, establish 0.01, ratio 10.
- **Axial profiles** (S7 — sourced endpoints, fixed interpolation):
  transit 12 / 8 / 10 h (table midpoints of 8–15 / 5–10 / 5–15 h),
  pH 6.25 / 6.60 / 6.80 (midpoints of 6.0–6.5 / 6.4–6.8 / 6.6–7.0),
  lumen carbon 8 / 5 / 2.5 mol m⁻³ ×10⁻³, lumen volume 0.40 / 0.35 / 0.30 L,
  contraction k0 0.5 / 0.4 / 0.3 min⁻¹, outflow survival 0.92 / 0.90 / 0.88.
- **Mucosa→lumen mechanisms**: mucus turnover k₀ per region (0.3–0.5 min⁻¹
  range of the Layer-2 table), contraction + HAPC detachment
  (α_HAPC = 0.2 min⁻¹ when covered; ~6 events/day, 95 % antegrade, 80 % in
  the 08:00–20:00 day window), growth-edge shedding f_edge = 0.2.
- **Luminal chemistry**: shared `metabolic_mode` helpers — the same
  `yield_carbon`, acid pKa/Ki/max and anaerobic factors as Layer 1
  (`cfg.fixes.metabolism`, `cfg.chem_env.oxygen`); implemented once in
  `src/fixes/metabolic_mode.h`, used by both layers.
- **Stool normalization**: 128 g/day (median of the 100–200 g range).
- **Flatness bound** (declared, pre-run): max/min regional
  Enterobacteriaceae fraction ≤ 4.0. The same community-density denominator
  applies to every region in these arms, so its value cancels in the ratio.

Runs: `run_readout.py` — {gradient, uniform} × seeds {1001,1002,1003},
72 h at dt = 60 s, mucosal purge (90 % of luminal cells removed) at 48 h.
All runs terminated `HorizonReached`; the chain ledger (initial + mucosal
births + luminal births + reattached deposits reconciled against stocks +
stool + deaths + purged) **closed exactly at every step of every arm**
(`ledger.closed = true`, stocks = expected to printed precision).

## Verdicts (tail 24 h, post-purge recovery window)

*Correction 2026-10-21: the realized tail spans 47h50–71h50 (last sampled
row t = 258,600 s), not 48–72 h; a 48–72 h sensitivity preserves the
stool verdict. The (c) retention verdict rests on configured transit
sums (30 h gradient / 24 h uniform) — "effective ~29 h" was not computed
by the analyzer and measured microbial retention time is not
established by these runs.*

| Observable | Target | Gradient arm | Uniform arm | Met? |
|---|---|---|---|---|
| (a) Stool CFU/g | 1e6–1e8 | 8.28e6 (all seeds) | 8.12e6 (all seeds) | **yes** |
| (b) Mucosal flatness | ≤ 4.0 | 2.07–2.96 | 1.15–1.38 | **yes** |
| (c) Retention / shedding | transit 30–70 h (lit.) | effective ~29 h in-chain; distal residence ~5.7 h | same | **yes, marginal** |

### Mucosal fraction: historical versus corrected denominator (GAP-04)

The six archived runs used a community-density proxy of **1e7 cells/mL**,
as recorded in each `out/*_prov.json`. Those records describe the original
runs; their final-sample fractions are not 24 h means. The table compares
the historical tail fractions in `readout_summary.json` at `28242fd` with
the corrected summary: the same regional mean CFU/mL divided by the
declared **1e8 cells/mL** FISH-scale total community density.

Regional means over 47h50–71h50 inclusive (145 samples per region).
Each entry is **historical % → reinterpreted %**:

| Arm | Seed | Cecum | Transverse | Descending |
|---|---:|---:|---:|---:|
| Gradient | 1001 | 13.722 → 1.372 | 19.393 → 1.939 | 40.581 → 4.058 |
| Gradient | 1002 | 13.767 → 1.377 | 20.701 → 2.070 | 30.156 → 3.016 |
| Gradient | 1003 | 13.475 → 1.347 | 19.336 → 1.934 | 27.827 → 2.783 |
| Uniform | 1001 | 19.744 → 1.974 | 16.723 → 1.672 | 17.712 → 1.771 |
| Uniform | 1002 | 19.025 → 1.903 | 21.499 → 2.150 | 15.602 → 1.560 |
| Uniform | 1003 | 17.232 → 1.723 | 19.599 → 1.960 | 17.045 → 1.704 |

`run_readout.py --analyze-only` (Python 3.12) recomputed all 18 fractions
from the archived series; each is one tenth of the historical fraction
within floating-point rounding. No trajectories were rerun. Regional mean
CFU/mL and stool CFU/g are unchanged, as are the six closed ledgers and
all non-fraction summary fields. Both arms still pass stool and flatness;
`verdicts.csv` is byte-identical to the historical file.

### Post-purge window (`> 48 h`) vs `71h50m` sample vs `72 h` endpoint (GAP-05 / #467)

Because the `48 h` (`172,800 s`) luminal purge executes after step `2880` is emitted (`purged_cum = 0` at both `t = 172,200 s` and `t = 172,800 s`, jumping to `5.68e9` at `t = 173,400 s`), `run_readout.py` now reports three distinct temporal summaries alongside the preserved 145-sample legacy window (`47h50m–71h50m`):

1. **Strictly post-purge window (`t > 48 h`, `173,400–258,600 s`, 143 samples/region):** Stool output is `8.12e6 CFU/g` (`gradient`) and `8.00e6 CFU/g` (`uniform`), with mucosal flatness `2.05–2.94` (`gradient`) and `1.15–1.38` (`uniform`) — preserving every Phase-3 gate verdict.
2. **Last time-series CSV sample (`71h50m`, `258,600 s`) vs exact `72 h` (`259,200 s`) `prov.json` endpoint:** Terminal contraction/HAPC shedding lowers instantaneous cecal binary occupancy to `90.25%–96.00%` at `71h50m` and `89.75%–96.25%` at exact `72 h` (`gradient`), distinct from the `99.45%–99.64%` post-purge window mean.
3. **Microcolony density & count-thresholded reference occupancy (`N >= 1, 10, 50` cells/patch):** Occupied patches carry `16.5–19.0` (`cecum`), `22.4–28.4` (`transverse`), and `36.1–57.3` (`descending`) cells per occupied patch in the `gradient` arm (`22.0–29.6` in `uniform`), confirming multi-cell microcolonies rather than single transient immigrants, while high-threshold (`N >= 50` Poisson reference) occupancy separates descending mucosa (`2.2%–76.3%`) from cecum/transverse (`<= 0.11%`). A dedicated `2 x 3` `reattach_p0` (`0.1x, 0.3x, 1.0x`) x (`gradient`, `uniform`) sensitivity generator is staged in `experiments/spec13_occupancy_sensitivity/prepare.py`.

| Arm | Seed | Region | Post-48h `N>=1` (%) | `71h50m` (%) | `72h` prov (%) | Cells / occ patch | Poisson ref `N>=10` / `N>=50` (%) |
|---|---:|---|---:|---:|---:|---:|---:|
| Gradient | 1001 | Cecum | 99.45 | 90.25 | 89.75 | 18.44 | 97.16 / 0.00 |
| Gradient | 1001 | Transverse | 99.85 | 99.75 | 100.00 | 26.11 | 99.97 / 0.02 |
| Gradient | 1001 | Descending | 99.81 | 100.00 | 99.75 | 57.31 | 100.00 / 76.31 |
| Gradient | 1002 | Cecum | 99.64 | 96.00 | 96.25 | 16.47 | 93.12 / 0.00 |
| Gradient | 1002 | Transverse | 99.84 | 100.00 | 100.00 | 28.42 | 99.99 / 0.11 |
| Gradient | 1002 | Descending | 99.82 | 100.00 | 99.75 | 43.29 | 100.00 / 19.36 |
| Gradient | 1003 | Cecum | 99.50 | 90.50 | 91.25 | 18.95 | 96.69 / 0.00 |
| Gradient | 1003 | Transverse | 99.85 | 100.00 | 100.00 | 22.36 | 99.81 / 0.00 |
| Gradient | 1003 | Descending | 99.82 | 100.00 | 100.00 | 36.09 | 100.00 / 2.21 |
| Uniform | 1001 | Cecum | 99.84 | 100.00 | 100.00 | 26.22 | 99.98 / 0.01 |
| Uniform | 1001 | Transverse | 99.81 | 99.75 | 100.00 | 22.54 | 99.74 / 0.00 |
| Uniform | 1001 | Descending | 99.72 | 100.00 | 100.00 | 24.83 | 99.95 / 0.00 |
| Uniform | 1002 | Cecum | 99.85 | 100.00 | 100.00 | 21.98 | 99.73 / 0.00 |
| Uniform | 1002 | Transverse | 99.82 | 100.00 | 100.00 | 29.58 | 100.00 / 0.09 |
| Uniform | 1002 | Descending | 99.71 | 99.75 | 100.00 | 22.16 | 99.81 / 0.00 |
| Uniform | 1003 | Cecum | 99.82 | 100.00 | 100.00 | 24.12 | 99.89 / 0.00 |
| Uniform | 1003 | Transverse | 99.83 | 100.00 | 99.75 | 22.61 | 99.85 / 0.00 |
| Uniform | 1003 | Descending | 99.71 | 99.75 | 100.00 | 22.26 | 99.56 / 0.01 |

## S9 observable set — what each pins

| Observable | Result | Parameter(s) it pins |
|---|---|---|
| Stool CFU/g | 8.1–8.3e6 | distal outflow × `descending_outflow_survival` × `stool_g_per_day` |
| Flatness | 2.1–3.0 (grad) / 1.2–1.4 (unif) | regional `patch_supply_scale`, `flora_density_cfu_ml` proxy |
| Effective chain transit | 30 h configured (gradient) / 24 h (uniform); realized retention not measured by the analyzer | `transit_h` profile |
| Diurnal envelope of stool exports | day/night = 1.06 (grad), 1.42 (unif) | `hapc_*` vs basal outflow — HAPC modulation is visible but second-order to basal washout |
| Growth-in-transit amplification | luminal births / stool exports = 3.28 (grad), 2.96 (unif) | `lumen_carbon_mol_m3` profile, `anaerobic_mu_factor`, yield/acid params |
| Purge recovery | 12 h post-purge stool rate = 12–13 % of pre-purge | `reattach_p0`, `lumen_seed_cells`-independent regrowth, `purge_fraction` |

## Findings beyond the gate

1. **The mucosa saturates under luminal feedback.** Occupancy passes the
   Phase-1 0.1–1 % band early and saturates high by 24 h.
   *Corrected numbers 2026-10-21: initial patches are Bernoulli-seeded at
   an expected 1 %, so the first strict exceedance of 1 % (0.167–0.667 h
   across series) is not an unambiguous growth crossing; at exactly 12 h
   occupancy is 71.5–87.25 % (gradient) and 48–82 % (uniform); at 24 h
   it is 100 % (gradient) and 96.75–99.75 % (uniform) — not "~78–84 % by
   12 h / ~100 % by 24 h in both arms".* The Phase-1 "stationary" window was transient: in standalone
   Layer 2 the pool drains (4 h half-life) and cannot sustain reseeding;
   the Layer-3 lumen *grows*, so `p_attach = p0·(1−occ)` feedback reaches
   its attracting equilibrium at occ → 1. Ledger stays closed throughout —
   this is coupling physics, not a bookkeeping defect.
2. **Fraction level is not the flatness gate.** The denominator correction
   changes the estimated level, not the relative spread across regions.
   The declared density interval (1e8–2e8 cells/mL) is an observation
   assumption, not a measurement from these runs. For the empirical
   context see `docs/gutibm_density_reconciliation.md` and
   `literature_mucosal_density_consensus.md`.
3. **Luminal growth dominates stool production.** Luminal births
   (≈1.77e10) exceed mucosal births (≈5e5) by ~10⁴×; stool composition is
   set by the distal compartment's growth + outflow, as S10 requires.
4. **Diurnal contrast is weak** (1.06–1.42 vs the 4× that 80 %-daytime
   HAPC placement would produce if HAPC detachment dominated exports).
   Basal outflow at rate dt/τ dominates the stool stream; HAPC shows in
   the mucosal detachment ledger (`contraction_cum`) more than in stool.
5. **Seed robustness**: all observables agree across seeds to <0.1 % —
   population-scale counts are deterministic at these sizes; per-lineage
   composition (not scored here) carries the stochasticity.

## Non-goals honored

No antibiotics, no live-patch promotion, no GPU patches, no per-segment
free parameters (profiles are fixed axial values), no fit to the Elliott
or Swidsinski figures, no default-mode behavior change
(`layer3.enabled` defaults false; Layer-2 keys unchanged when unset).
S11 stands: the luminal compartment is anaerobic-only; no O₂-dependent
result is reported.

## Reproduce

```
python3 analysis/spec13_gradient/run_readout.py --analyze-only
```

Per-arm configs are generated under `configs/`; raw rows under
`out/*_timeseries.csv`, run records under `out/*_prov.json`, aggregate in
`readout_summary.json`, per-arm verdicts in `verdicts.csv`.
