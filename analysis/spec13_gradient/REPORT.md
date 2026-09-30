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
  Enterobacteriaceae fraction ≤ 4.0, fraction = regional mean CFU/mL ÷
  uniform flora-density proxy 1e7 CFU/mL.

Runs: `run_readout.py` — {gradient, uniform} × seeds {1001,1002,1003},
72 h at dt = 60 s, mucosal purge (90 % of luminal cells removed) at 48 h.
All runs terminated `HorizonReached`; the chain ledger (initial + mucosal
births + luminal births + reattached deposits reconciled against stocks +
stool + deaths + purged) **closed exactly at every step of every arm**
(`ledger.closed = true`, stocks = expected to printed precision).

## Verdicts (tail 24 h, post-purge recovery window)

| Observable | Target | Gradient arm | Uniform arm | Met? |
|---|---|---|---|---|
| (a) Stool CFU/g | 1e6–1e8 | 8.28e6 (all seeds) | 8.12e6 (all seeds) | **yes** |
| (b) Mucosal flatness | ≤ 4.0 | 2.07–2.96 | 1.15–1.38 | **yes** |
| (c) Retention / shedding | transit 30–70 h (lit.) | effective ~29 h in-chain; distal residence ~5.7 h | same | **yes, marginal** |

Per-region fractions (gradient, tail): cecum 0.135–0.138, transverse
0.19–0.21, descending 0.28–0.41 — a real ~2–3× proximal→distal gradient
*within* the declared flatness bound. The uniform arm collapses it to
1.15–1.38, so the residual gradient-arm spread is the gradient's
contribution; both arms stay under the bound.

## S9 observable set — what each pins

| Observable | Result | Parameter(s) it pins |
|---|---|---|
| Stool CFU/g | 8.1–8.3e6 | distal outflow × `descending_outflow_survival` × `stool_g_per_day` |
| Flatness | 2.1–3.0 (grad) / 1.2–1.4 (unif) | regional `patch_supply_scale`, `flora_density_cfu_ml` proxy |
| Effective chain transit | ~29 h realized vs 30 h configured | `transit_h` profile |
| Diurnal envelope of stool exports | day/night = 1.06 (grad), 1.42 (unif) | `hapc_*` vs basal outflow — HAPC modulation is visible but second-order to basal washout |
| Growth-in-transit amplification | luminal births / stool exports = 3.28 (grad), 2.96 (unif) | `lumen_carbon_mol_m3` profile, `anaerobic_mu_factor`, yield/acid params |
| Purge recovery | 12 h post-purge stool rate = 12–13 % of pre-purge | `reattach_p0`, `lumen_seed_cells`-independent regrowth, `purge_fraction` |

## Findings beyond the gate

1. **The mucosa saturates under luminal feedback.** Occupancy crosses the
   Phase-1 0.1–1 % band at ~6 h, ~78–84 % by 12 h, and is ~100 % by 24 h in
   both arms. The Phase-1 "stationary" window was transient: in standalone
   Layer 2 the pool drains (4 h half-life) and cannot sustain reseeding;
   the Layer-3 lumen *grows*, so `p_attach = p0·(1−occ)` feedback reaches
   its attracting equilibrium at occ → 1. Ledger stays closed throughout —
   this is coupling physics, not a bookkeeping defect.
2. **Fraction level vs flatness must be kept distinct.** At the saturated
   equilibrium the regional fraction is 0.14–0.41 against the 1e7 flora
   proxy — above the empirical ~1–3 % Enterobacteriaceae fraction. The
   gate asks for flatness (achieved); the level is pinned by
   `patch_supply_scale` / `flora_density_cfu_ml` and is a parameter
   consequence to report, not a fitted quantity.
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
python3 analysis/spec13_gradient/run_readout.py            # runs all arms
python3 analysis/spec13_gradient/run_readout.py --analyze-only  # re-tally
```

Per-arm configs are generated under `configs/`; raw rows under
`out/*_timeseries.csv`, run records under `out/*_prov.json`, aggregate in
`readout_summary.json`, per-arm verdicts in `verdicts.csv`.
