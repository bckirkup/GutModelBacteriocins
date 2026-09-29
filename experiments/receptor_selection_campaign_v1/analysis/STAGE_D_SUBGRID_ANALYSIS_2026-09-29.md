# Stage D sub-0.25% grid combined analysis — 98 runs (2026-09-29)

Combined readout of four arrays on the same execution source `29fd563c`,
image digest `fafc77ee`, job definition `gutibm-cuda-campaign:14`:

- original Stage D `47460bfa` (5 arms x seeds 20260911/13/17, indices 0-14)
- seed extension `734ab0fb` (5 arms x seeds 20260919-20261007, indices 15-64)
- grid refinement `adfd5af1` (7 arms x seeds 20261009/1011/1013, indices 65-85)
- sub-0.25% fill-in `1fc9a916` (4 arms x same seeds, indices 86-97)

The sub-grid probes the 0-0.25% realm at 0.05/0.1/0.15/0.20% nominal lysis
on the same three seeds as the refinement array, so the existing same-seed
null, 0% carrier, and 0.25% arms complete the six-point local grid with
same-seed pairing.

## Integrity

- 98/98 outputs readable; `execution_source_sha` matches `29fd563c` in all;
  `chemistry_placement=device_delivery` in all.
- All 12 new runs succeeded first attempt (exit 0) and terminated on the
  **dysbiosis guard** at t=20000 s, matching every prior Stage D arm.
  Paired contrasts use the shared calendar window ending at `t_common`
  per seed.
- n per arm: null 16, target 0 16, targets 0.01/0.02/0.05 13, each
  refine/subgrid target (0.0005-0.03 except 0.01/0.02/0.05) 3. Screening
  replicates — medians and ranges, no p-values.

## Sub-0.25% per-arm medians

| arm / target | n | divisions | lysis | realized lysis / producer division | kills / lysis | own-end slope (dec/h) |
|---|---|---|---|---|---|---|
| producer 0.0005 | 3 | 18,970 | 3 | 0.00039 | 72.7 | 0.0117 |
| producer 0.0010 | 3 | 18,731 | 5 | 0.00064 | 46.2 | 0.0148 |
| producer 0.0015 | 3 | 27,646 | 16 | 0.00140 | 63.9 | 0.0518 |
| producer 0.0020 | 3 | 26,195 | 25 | 0.00216 | 59.0 | 0.0707 |

## Paired contrast vs same-seed plasmid-free null (shared `t_common` window)

| target | n | Δ slope med (dec/h) | Δ slope range | Δ tail median (dec) |
|---|---|---|---|---|
| 0.000 | 16 | +0.0150 | [+0.0090, +0.0231] | -0.0003 |
| 0.0005 | 3 | +0.0251 | [+0.0242, +0.0267] | +0.0062 |
| 0.0010 | 3 | +0.0353 | [+0.0283, +0.0366] | +0.0083 |
| 0.0015 | 3 | +0.0451 | [+0.0417, +0.0667] | +0.0238 |
| 0.0020 | 3 | +0.0739 | [+0.0603, +0.0893] | +0.0327 |
| 0.0025 | 3 | +0.0734 | [+0.0685, +0.0818] | +0.0358 |

Per-seed Δ-slope series across 0 -> 0.0025 are monotone increasing through
0.0020 at all three seeds; the 0.0025 point sits level with 0.0020
(+0.0739 vs +0.0734 medians, overlapping ranges).

## Readout

- **No threshold; the response is linear through the origin on
  0.05-0.15%.** Δ slope over the inert carrier rises +0.010, +0.020,
  +0.030 dec/h at 0.05/0.1/0.15% — a straight line with slope ~20 dec/h
  per unit target whose extrapolated intercept at 0 is +0.0007 dec/h,
  inside carrier noise. Selection response is proportional to lysis
  probability at the lowest levels resolvable.
- **The concavity starts right at 0.20-0.25%.** Δ slope flattens between
  the 0.20% and 0.25% arms (+0.0739 vs +0.0734 medians; every seed's
  ranges overlap), before the response resumes rising at 0.5% (+0.102).
  The knee the coarse grid implied near ~0.25% is now localized at
  0.20-0.25% nominal.
- **Kill efficiency is highest at the lowest lysis.** kills/lysis runs
  ~46-73 across the sub-0.25% arms versus 54 at 0.25% declining to 3.3 at
  5% — each rare lysis event lands in an undepleted susceptible pool,
  consistent with the linear regime; the knee appears as the susceptible
  pool begins to deplete.
- **Realized lysis tracks nominal at the low end.** Realized
  lysis/producer division medians run ~0.64-1.08x nominal on the
  sub-grid arms (individual seeds 0.36-1.12x; small event counts of
  3-25 lysis events per run drive the spread).

## Caveats

- n=3 per sub-grid arm, screening replicates; the 0.20-vs-0.25% ordering
  is within seed noise at 2/3 seeds.
- All arms halt on the dysbiosis guard at t=20000 s; contrasts are on the
  shared window ending at `t_common`, so they measure relative
  repression rate, not endpoint density.

Metrics: `analysis/stage_D_subgrid/` (98-run `run_metrics`,
`paired_metrics`, `gate_status`, `missing_outputs`).
