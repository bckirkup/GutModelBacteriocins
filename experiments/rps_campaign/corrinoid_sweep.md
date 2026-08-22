# Corrinoid / BtuB competition sweep

Post-#299 code (no ColE1 SOS cascade). Flux-0.3x epithelial delivery arm,
VBF sink 1x, gradient off — otherwise identical to `vbf_f030_s100_*.json`.
`kd_b12_btuB` (alias `kd_corrinoid_btuB`) swept over the range
`docs/PARAMETERS.md` recommends; 3 seeds per arm (1001/2027/3313); 24 h
horizon; `kill_rate_colicin` left at its default 1e-3 /s.

The ambient corrinoid field is pinned at 1e-3 mol/m^3, so

```
apparent_Kd = kd_colicinE_btuB * (1 + [B12] / kd_b12_btuB)
```

gives competitive factors 1001 / 101 / 11 / 2 across the four arms.

Configs: `corr_{kd1e6,kd1e5,kd1e4,kd1e3}_s{1001,2027,3313}.json`
Image: `sha256:49eb3713…` (main + the env-gated receptor trace, off by default)
Analysis: `corrinoid_analysis.py`, `corrinoid_analysis.json`

## Results, all arms compared at matched simulated time 61 800 s (17.2 h)

Seed means. `E[kills]` is the summed per-step Bernoulli kill probability over
all live type-2 agents — a far lower-variance estimator of colicin pressure
than the handful of realized deaths. Occupancy percentiles are over type-2
(susceptible) agent-observations only, i.e. sampled at the targets rather than
taken from the global field maximum.

| `kd_b12_btuB` | factor | apparent_Kd | occ p99 | occ max | E[kills] | kills | lysis | SOS | divisions | type 1 | type 2 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1e-6 (current) | 1001 | 5.005e-4 | 5.42e-3 | 8.18e-3 | 0.66 | 0.0 | 12.3 | 12 | 426 | 168 | 73 |
| 1e-5 | 101 | 5.050e-5 | 2.58e-2 | 7.03e-2 | 2.98 | 3.0 | 12.3 | 12 | 416 | 161 | 68 |
| 1e-4 | 11 | 5.500e-6 | 2.71e-1 | 4.63e-1 | 16.42 | 14.3 | 10.7 | 10 | 387 | 157 | 39 |
| 1e-3 | 2 | 1.000e-6 | 5.01e-1 | 8.02e-1 | 19.31 | 17.7 | 13.0 | 13 | 373 | 185 | 3 |

Per-seed detail is in the script output; the seed spread is wide (E[kills]
0.28–1.01 at 1e-6, 7.1–26.7 at 1e-4) but the arm separation is far wider.

## What this establishes

1. **Occupancy, not the kill rate, was the limiting quantity.** Realized kills
   track E[kills] closely in every arm (2.98 vs 3.0; 16.42 vs 14.3; 19.31 vs
   17.7), so the counter is honest and the deficit at the default parameter is
   exposure. With occupancy freed, the *existing* `kill_rate_colicin` = 1e-3 /s
   eradicates the susceptible population — type 2 falls from 73 survivors at
   1001-fold competition to 3 at 2-fold, having started at 20 in both. The
   proposed 50x increase to 0.05 /s is therefore not needed and would be
   fitted to a suppression it does not address.
2. **The sweep is graded, not binary.** Unlike `vbf_carbon_sink_vmax` — which
   saturated by 0.5x — every step of this parameter moves occupancy, expected
   kills and the surviving target population monotonically. It is the strongest
   graded lever on colicin efficacy found so far.
3. **Producer cost is now small in every arm** (10–13 lysis events over 17 h,
   against 104–148 pre-#299), and producers grow in all four arms, so the
   cost/benefit no longer inverts on an artifact.

## What it does not establish

The sweep does not say which value is right. The model spans "colicin is
irrelevant" to "colicin eradicates the susceptible population" across the
documented uncertainty of a single parameter, which means the campaign cannot
inherit 1e-6 silently: either the free-corrinoid concentration and BtuB
affinity get pinned against literature — noting that the pinned 1e-3 mol/m^3
pool is total corrinoid, most of which in the colon is non-cobalamin analogue
and may not compete at BtuB with the same affinity — or `kd_corrinoid_btuB`
becomes an explicit campaign factor.

Caveats: `uptake_limit = none` in these runs, matching the earlier delivery and
sink arms rather than the new `sherwood` funding, so they are comparable to
those arms and not to a Sherwood-funded campaign. Eight of twelve runs were
guard-halted before 24 h, which is why every comparison above is at 17.2 h.
