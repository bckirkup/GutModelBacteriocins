# Stage D grid-refinement submission (2026-09-29)

## What was submitted

A refined Stage D lysis grid — 7 arms x 3 new seeds = 21 Batch array
children — probing onset and curvature between the deployed 0/1/2/5%
targets while preserving both controls:

- `D_plasmid_free_null` (control, re-run on new seeds)
- `D_target0.0000_producer` (ColE1 zero-lysis carrier control, re-run on new seeds)
- `D_target{0.0025,0.0050,0.0075,0.0150,0.0300}_producer` (new arms)

- **Array job:** `adfd5af1-3fd1-4ea8-b3d2-6b7d40609d66`
  (`gutibm-receptor-selection-v1-stage-d-refine`, array size 21)
- **Queue:** `gutibm-gpu-campaign`
- **Job definition:** `gutibm-cuda-campaign:14` — same revision as the
  prior Stage D arrays; image digest
  `sha256:fafc77ee802882bee6e47321515057aed20c0326d622bfb287862ef15304f748`,
  timeout 7200 s
- **Execution source SHA:** `29fd563ce73b28059cb439694e92466f451ddb4a`
  (unchanged — inputs were generated in a detached worktree at that commit)
- **Inputs:** `s3://gutibm-inputs-994254241749/receptor-selection-v1/29fd563c/stage-d-refine/inputs/<i>/input.json`
  plus `manifest.json` at the prefix root
- **Outputs:** `s3://gutibm-outputs-994254241749/receptor-selection-v1/29fd563c/stage-d-refine/outputs/<i>/output.h5.gz`
- **New seeds:** 20261009, 20261011, 20261013 (disjoint from the 13 seeds
  already run: 20260911/13/17 plus the 20260919–20261007 extension block)
- **Approval record:** `approval_stage_D_refine.json`
- **Generator:** `prepare_stage_D_refine.py` (records exactly how inputs
  were derived; ran in a `29fd563c` worktree)

## Lysis-target mapping

Each division draws the lysis hazard on both mother and daughter, so the
realized per-generation lysis is `1-(1-p)^2`; `sos_lysis_prob = 1 - sqrt(1-t)`
with `sos_basal_rate = 0` reproduces nominal target `t`:

| nominal target | sos_lysis_prob |
|---|---|
| 0 (carrier) | 0 |
| 0.0025 | 0.0012507822280910519 |
| 0.005 | 0.0025031328369998773 |
| 0.0075 | 0.0037570577414361983 |
| 0.015 | 0.007528337936039575 |
| 0.03 | 0.015114219820389518 |

## Constancy verification

Before submission, all 21 generated inputs were compared field-by-field
against the 15 deployed stage-d inputs downloaded from S3. After removing
`seed`, `_campaign.run_id`, and `_campaign.paired_seed`, the two control
arms are byte-identical to their deployed counterparts, and the five new
producer arms differ from the deployed producer template only in
`sos_lysis_prob` and the `_campaign.axes` lysis fields (0 mismatches).
`selected_amplitude` is float 15.0, matching deployed axes metadata.

## Expected analysis

Combined with the deployed grid, the refined arms fill in onset
(0–1%) and curvature (1–5%) of the lysis–selection response: 16 seeds at
target 0 and the null, 3 seeds at each new target. Paired contrasts use
the shared calendar window `t_common = min(t_end_producer, t_end_null)`
per seed, as before.
