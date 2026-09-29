# Stage D sub-0.25% fill-in submission (2026-09-29)

## What was submitted

Four new producer targets filling the 0-0.25% lysis interval on the same
seeds as the grid-refinement array — 4 arms x 3 seeds = 12 Batch array
children. Together with the stage-d-refine runs on seeds 20261009/1011/1013
(plasmid-free null, zero-lysis carrier, target 0.25%), this completes the
six-point sub-grid 0 / 0.05 / 0.10 / 0.15 / 0.20 / 0.25% with same-seed
paired controls.

- **Array job:** `1fc9a916-fb4a-4655-a664-47241ab083d5`
  (`gutibm-receptor-selection-v1-stage-d-subgrid`, array size 12)
- **Queue:** `gutibm-gpu-campaign`
- **Job definition:** `gutibm-cuda-campaign:14` — image digest
  `sha256:fafc77ee802882bee6e47321515057aed20c0326d622bfb287862ef15304f748`,
  timeout 7200 s
- **Execution source SHA:** `29fd563ce73b28059cb439694e92466f451ddb4a`
  (unchanged — inputs generated in a detached worktree at that commit)
- **Inputs:** `s3://gutibm-inputs-994254241749/receptor-selection-v1/29fd563c/stage-d-subgrid/inputs/<i>/input.json`
  plus `manifest.json` at the prefix root
- **Outputs:** `s3://gutibm-outputs-994254241749/receptor-selection-v1/29fd563c/stage-d-subgrid/outputs/<i>/output.h5.gz`
- **Seeds:** 20261009, 20261011, 20261013 (same as stage-d-refine — no new
  seeds; paired contrasts reuse the existing same-seed nulls)
- **Approval record:** `approval_stage_D_subgrid.json`
- **Generator:** `prepare_stage_D_subgrid.py` (ran in a `29fd563c` worktree)

## Lysis-target mapping

`sos_lysis_prob = 1 - sqrt(1 - t)` with `sos_basal_rate = 0` (two-draw
division hazard), same construction as the deployed grid:

| nominal target | sos_lysis_prob |
|---|---|
| 0.0005 | 0.00025003125781486446 |
| 0.001 | 0.0005001250625390474 |
| 0.0015 | 0.0007502814611354269 |
| 0.002 | 0.0010005005006258338 |

## Constancy verification

All 12 generated inputs were compared field-by-field against the deployed
stage-d producer inputs downloaded from S3. After removing `seed`,
`_campaign.run_id`, `_campaign.paired_seed`, `sos_lysis_prob`, and the
`_campaign.axes` lysis fields, every input is byte-identical (0 mismatches).

## Expected analysis

Extends the combined 86-run manifest to 98 runs (indices 86-97). The
sub-grid tests whether the lysis-selection response has a low-dose
threshold or remains linear through the origin on the 0-0.25% interval;
paired contrasts reuse each run's same-seed plasmid-free null and the
shared `t_common` window.
