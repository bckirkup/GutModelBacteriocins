# Stage D dysbiosis-guard diagnostic — 1e13 cells/mL (2026-09-29)

Guard diagnostic following the 13-seed Stage D combined readout, in which every
run terminated on the dysbiosis guard at ~16800–20400 s against the 21600 s
horizon. Question: does the arm ordering survive if runs actually reach the
intended six hours? Eight runs: 4 arms x 2 seeds with **only**
`dysbiosis_threshold` changed (1e10 -> 1e13 cells/mL). This is a diagnostic,
not pooled replicates — the threshold changes the stopping regime, so these
outputs are not merged into the 65-run Stage D set.

## What was submitted

- **Array job:** `de21b0b0-26c0-4c39-9959-d088572718a7`
  (`gutibm-receptor-v1-d-guard-1e13`, array size 8)
- **Queue / job def:** `gutibm-gpu-campaign` / `gutibm-cuda-campaign:14` —
  same revision and image digest
  `sha256:fafc77ee802882bee6e47321515057aed20c0326d622bfb287862ef15304f748`
  as both prior Stage D arrays; timeout 7200 s
- **Execution source SHA:** `29fd563ce73b28059cb439694e92466f451ddb4a`
  (inputs generated in a detached worktree at that commit)
- **Arms:** `D_target{0.000,0.010,0.050}_producer` + `D_plasmid_free_null`
  (the 0.020 arm is dropped per the diagnostic spec)
- **Seeds:** 20260911, 20260913 — original contract seeds, so each run pairs
  against an already-analyzed same-seed Stage D output
- **Inputs:** `s3://gutibm-inputs-994254241749/receptor-selection-v1/29fd563c/stage-d-guard-1e13/inputs/<i>/input.json`
  plus `manifest.json` at the prefix root
- **Outputs:** `s3://gutibm-outputs-994254241749/receptor-selection-v1/29fd563c/stage-d-guard-1e13/outputs/<i>/output.h5.gz`
- **Generator:** `prepare_stage_D_guard_1e13.py` (record of derivation)

## Constancy verification

All 8 generated inputs were compared field-by-field against the deployed
`stage-d/inputs` objects for the same (arm, seed) pairs. Every input is
byte-identical except `dysbiosis_threshold` (1e10 -> 1e13); run_id, seed,
paired_seed, strains, and all physics/config keys are unchanged.

## Horizon check — PASS

All 8 runs report `termination_cause=horizon_reached` at `t_end_s=21600`
(`/run_provenance`: `completed_total_time=1`, `halt_reason_code=0`),
`execution_source_sha_match` true, `chemistry_placement=device_delivery`.
Previously every run guard-halted at ~16800–20400 s; at 1e13 the density never
re-arms the guard and all arms run the full six hours.

| arm | seed 20260911 t_end | seed 20260913 t_end | cause |
|---|---|---|---|
| plasmid_free_null | 21600 | 21600 | horizon_reached |
| producer target 0.000 | 21600 | 21600 | horizon_reached |
| producer target 0.010 | 21600 | 21600 | horizon_reached |
| producer target 0.050 | 21600 | 21600 | horizon_reached |

## Ordering check — PASS (retained, strengthened)

Paired Δ slope vs same-seed plasmid-free null, shared `t_common` window
(7200 s tail window as in the combined analysis):

| lysis target | prior Δ slope @16800 (s11 / s13) | new Δ slope @21600 (s11 / s13) |
|---|---|---|
| 0.000 | +0.020 / +0.016 | −0.012 / −0.023 |
| 0.010 | +0.171 / +0.155 | +0.381 / +0.395 |
| 0.050 | +0.389 / +0.376 | +0.706 / +0.702 |

- Monotone ordering holds in **both** seeds: P=0 carrier ~inert, then strictly
  increasing with the lysis target; the 0.010 and 0.050 deltas are far outside
  each other's error for a two-seed diagnostic.
- The response **strengthens** at the full horizon — selection accumulates over
  the extra ~1.3–1.3 h of simulated time the guard used to cut off.
- The P=0 carrier shifts to slightly negative (−0.01 to −0.02 dec/h) vs the
  prior +0.015–0.020 — still inert relative to the lysis arms; no reversal.
- Realized per-producer-division lysis now reads ≈ nominal (0.0096–0.0105 for
  the 1% arm, 0.048–0.051 for the 5% arm); kills/lysis 31.8–33.4 at 1% and
  3.0–3.4 at 5%, consistent with the prior trend.

## Readout

**Guard diagnostic PASSED.** All arms reach 21,600 s at threshold 1e13, the
monotone producer-fraction ordering is retained, and the claim of behavior
through the intended six hours is supported — at higher, not lower, contrast
than the guard-truncated readout. Not flagged; the 13-seed conclusion stands
with the caveat that the prior window understates the response.

## Operational note

Spot capacity churned mid-array: children 1 and 2 wrote partial status
heartbeats (~91%) before their instance was reclaimed, then re-ran from step 0
on new capacity (no checkpointing configured). All children ultimately
SUCCEEDED; outputs are complete and provenance-clean.

## Artifacts

`analysis/stage_D_guard_1e13/` — `run_metrics.csv/json`,
`paired_metrics.csv/json`, `gate_status.json`, `missing_outputs.json`
produced by `analyze.py` against the 8-run manifest.
