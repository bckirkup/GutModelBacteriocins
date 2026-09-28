# Stage D seed-extension submission (2026-09-28)

## What was submitted

The same 5 Stage D arms — `D_target{0.000,0.010,0.020,0.050}_producer` and
`D_plasmid_free_null` — rerun against 10 new seeds (50 Batch array children),
everything else held constant.

- **Array job:** `734ab0fb-7df4-4023-80d7-e703075d8b23`
  (`gutibm-receptor-selection-v1-stage-d-seeds10`, array size 50)
- **Queue:** `gutibm-gpu-campaign`
- **Job definition:** `gutibm-cuda-campaign:14` — same revision as the
  original Stage D array `47460bfa-1558-47b6-9a48-47b1bd9a7b0d`; image digest
  `sha256:fafc77ee802882bee6e47321515057aed20c0326d622bfb287862ef15304f748`,
  timeout 7200 s
- **Execution source SHA:** `29fd563ce73b28059cb439694e92466f451ddb4a`
  (unchanged — inputs were generated in a detached worktree at that commit)
- **Inputs:** `s3://gutibm-inputs-994254241749/receptor-selection-v1/29fd563c/stage-d-seeds10/inputs/<i>/input.json`
  plus `manifest.json` at the prefix root
- **Outputs:** `s3://gutibm-outputs-994254241749/receptor-selection-v1/29fd563c/stage-d-seeds10/outputs/<i>/output.h5.gz`
- **New seeds:** 20260919, 20260921, 20260923, 20260925, 20260927,
  20260929, 20261001, 20261003, 20261005, 20261007
  (disjoint from the original 20260911 / 20260913 / 20260917)
- **Approval record:** `approval_stage_D_seeds10.json`
- **Generator:** `prepare_stage_D_seeds10.py` (records exactly how inputs were
  derived; ran in a `29fd563c` worktree via `stage_d_runs()`)

## Constancy verification

Before submission, all 50 generated inputs were compared field-by-field
against the 15 deployed stage-d inputs downloaded from S3. After removing
`seed`, `_campaign.run_id`, and `_campaign.paired_seed`, every new input is
byte-identical to the corresponding deployed arm input (0 mismatches).
`selected_amplitude` was cast to float (15.0) to match the deployed axes
metadata exactly.

## Expected analysis

After completion, combine with the original 15 runs for 13 seeds per arm
(5 arms x 13 seeds total). The dysbiosis guard halted the original runs at
~20000 s against the 21600 s horizon — paired contrasts must again use the
shared calendar window `t_common = min(t_end_producer, t_end_null)` per seed.
