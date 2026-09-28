# Stage D combined analysis — 13 seeds per arm (2026-09-28)

Combined readout of the original 15-job Stage D array
(`47460bfa-1558-47b6-9a48-47b1bd9a7b0d`, seeds 20260911/13/17) and the
50-job seed-extension array (`734ab0fb-7df4-4023-80d7-e703075d8b23`, seeds
20260919–20261007 odd days). Same execution SHA `29fd563c`, image digest
`fafc77ee`, job definition `gutibm-cuda-campaign:14` throughout.

## Integrity

- 65/65 outputs readable; `execution_source_sha` matches `29fd563c` in all 65;
  `chemistry_placement=device_delivery` in all 65.
- All 65 runs terminated on the **dysbiosis guard** (~16800–20400 s simulated
  vs the 21600 s horizon) — a scientific termination, not a Batch failure.
  Producer and null arms stop at different steps, so all paired contrasts use
  the shared window ending at `t_common` per seed (as analyzer-enforced).

## Per-arm medians (13 seeds each)

| arm | divisions | lysis | realized lysis / producer division | kills / lysis | own-end slope (dec/h) |
|---|---|---|---|---|---|
| plasmid_free_null | 17,498 | 0 | 0.0000 | — | 0.0019 |
| producer target 0.000 | 19,807 | 0 | 0.0000 | — | 0.0010 |
| producer target 0.010 | 23,687 | 110 | 0.0075 | 28.6 | 0.1968 |
| producer target 0.020 | 34,627 | 480 | 0.0177 | 12.2 | 0.4352 |
| producer target 0.050 | 29,198 | 1,127 | 0.0456 | 3.3 | 0.6065 |

## Paired contrast vs same-seed plasmid-free null (shared `t_common` window)

| lysis target | n | Δ slope med (dec/h) | Δ slope IQR | Δ tail median (dec) |
|---|---|---|---|---|
| 0.000 | 13 | +0.0151 | [0.0105, 0.0178] | −0.0007 |
| 0.010 | 13 | +0.1590 | [0.1495, 0.1700] | +0.0918 |
| 0.020 | 13 | +0.2555 | [0.2386, 0.2625] | +0.1868 |
| 0.050 | 13 | +0.3798 | [0.3686, 0.3907] | +0.3336 |

## Readout

- **Monotone producer-fraction (lysis) response confirmed at 13 seeds:**
  Δ slope increases strictly with the lysis target in every seed pairing
  (IQRs are disjoint across adjacent targets).
- The P=0 carrier control is nearly inert vs the plasmid-free null
  (+0.015 dec/h) — plasmid carriage alone contributes almost nothing; the
  response is driven by realized lysis.
- Realized per-division lysis runs at ~75–91% of nominal targets
  (0.0075/0.0177/0.0456 vs 0.01/0.02/0.05) — the two-draw division hazard
  documented in AGENTS.md applies on top of the configured probability;
  report realized, not configured, values.
- `kills_per_lysis` falls steeply with lysis rate (28.6 → 3.3): each lysed
  producer releases a fixed burst, so marginal kills per lysis decline as
  the susceptible pool is depleted/reshaped.

## Artifacts

`analysis/stage_D_combined/` contains `run_metrics.csv/json`,
`paired_metrics.csv/json`, `gate_status.json`, `missing_outputs.json`
(none missing) produced by `analyze.py` against a combined 65-run manifest
(stage_D results 0–14 + stage_D2 results remapped to indices 15–64).
