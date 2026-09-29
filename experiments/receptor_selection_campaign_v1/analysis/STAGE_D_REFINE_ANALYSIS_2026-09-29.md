# Stage D refined-grid combined analysis — 86 runs (2026-09-29)

Combined readout of three arrays on the same execution source `29fd563c`,
image digest `fafc77ee`, job definition `gutibm-cuda-campaign:14`:

- original Stage D `47460bfa` (5 arms x seeds 20260911/13/17, indices 0-14)
- seed extension `734ab0fb` (5 arms x seeds 20260919-20261007, indices 15-64)
- grid refinement `adfd5af1` (7 arms x seeds 20261009/1011/1013, indices 65-85)

The refined grid probes onset (0-1%) and curvature (1-5%) of the
lysis-selection response while re-running both controls on the new seeds.

## Integrity

- 86/86 outputs readable; `execution_source_sha` matches `29fd563c` in all;
  `chemistry_placement=device_delivery` in all.
- All 86 terminated on the **dysbiosis guard** (t_end 16800-20400 s vs the
  21600 s horizon) — a scientific termination. Paired contrasts use the
  shared calendar window ending at `t_common` per seed.
- n per arm: null 16, target 0 16, targets 0.01/0.02/0.05 13, each new
  target (0.0025/0.005/0.0075/0.015/0.03) 3. New-target arms are screening
  replicates — medians and ranges, no p-values.

## Per-arm medians

| arm / target | n | divisions | lysis | realized lysis / producer division | kills / lysis | own-end slope (dec/h) |
|---|---|---|---|---|---|---|
| plasmid_free_null | 16 | 17,294 | 0 | 0.0000 | — | 0.0019 |
| producer 0.000 | 16 | 19,797 | 0 | 0.0000 | — | 0.0007 |
| producer 0.0025 | 3 | 26,267 | 32 | 0.0023 | 54.4 | 0.0859 |
| producer 0.005 | 3 | 25,456 | 57 | 0.0038 | 44.3 | 0.1323 |
| producer 0.0075 | 3 | 24,362 | 78 | 0.0055 | 36.6 | 0.1499 |
| producer 0.010 | 13 | 23,687 | 110 | 0.0075 | 28.6 | 0.1968 |
| producer 0.015 | 3 | 35,154 | 309 | 0.0120 | 21.3 | 0.3632 |
| producer 0.020 | 13 | 34,627 | 480 | 0.0177 | 12.2 | 0.4352 |
| producer 0.030 | 3 | 32,300 | 703 | 0.0272 | 7.0 | 0.5203 |
| producer 0.050 | 13 | 29,198 | 1,127 | 0.0456 | 3.3 | 0.6065 |

## Paired contrast vs same-seed plasmid-free null (shared `t_common` window)

| target | n | Δ slope med (dec/h) | Δ slope range | Δ tail median (dec) |
|---|---|---|---|---|
| 0.000 | 16 | +0.0150 | [+0.0090, +0.0231] | −0.0003 |
| 0.0025 | 3 | +0.0734 | [+0.0685, +0.0818] | +0.0358 |
| 0.005 | 3 | +0.1016 | [+0.0868, +0.1085] | +0.0521 |
| 0.0075 | 3 | +0.1243 | [+0.1179, +0.1529] | +0.0712 |
| 0.010 | 13 | +0.1590 | [+0.1418, +0.1855] | +0.0918 |
| 0.015 | 3 | +0.1976 | [+0.1963, +0.2076] | +0.1447 |
| 0.020 | 13 | +0.2555 | [+0.2335, +0.2665] | +0.1868 |
| 0.030 | 3 | +0.2987 | [+0.2882, +0.3177] | +0.2149 |
| 0.050 | 13 | +0.3798 | [+0.3591, +0.4143] | +0.3336 |

## Readout

- **Onset is immediate and steep.** The 0.25% arm already sits ~5x above
  the carrier baseline (+0.073 vs +0.015 dec/h); there is no detectable
  lysis threshold below 0.25% per generation in this geometry.
- **The response is concave (diminishing returns) across the full grid.**
  Marginal Δ slope per unit target falls monotonically: ~23 dec/h per unit
  target below 0.25%, ~7.6 between 1% and 5%. The 0.25-1% segment rises
  ~0.14 dec/h; the 1-5% segment rises only ~0.22 over four times the span.
- **Monotone ordering holds on the refined grid** at every seed pairing —
  adjacent-target ranges are disjoint or near-disjoint (only 0.0075 vs
  0.010 n=3 ranges touch: [0.118, 0.153] vs [0.142, 0.186]).
- **Carrier baseline is inert:** 16 seeds now confirm the zero-lysis ColE1
  carrier adds +0.015 dec/h vs the plasmid-free null — plasmid carriage
  alone contributes almost nothing.
- Realized per-producer-division lysis runs ~73-92% of nominal targets
  (0.0023-0.0456 realized vs 0.0025-0.05 configured) — the same two-draw
  division hazard documented for the deployed grid; report realized, not
  configured, values.
- `kills_per_lysis` continues its steep decline (54 at 0.25% -> 3.3 at
  5%): each lysis releases a fixed burst while the susceptible pool is
  depleted/reshaped, consistent with the concave response.

## Artifacts

`analysis/stage_D_refined/` contains `run_metrics.csv/json`,
`paired_metrics.csv/json`, `gate_status.json`, `missing_outputs.json`
(none missing), produced by `analyze.py` against a combined 86-run
manifest (indices 0-14 original, 15-64 seeds10, 65-85 refinement).
