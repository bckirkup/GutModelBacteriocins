# Stage D retardation-sensitivity submission (2026-09-29)

## What was submitted

Two arrays testing whether the Stage D lysis gradient moves along the
mucin-charge transport axis. The deployed gradient was measured entirely
at the selected amplitude 15, where `retardation_from_pI(9.0)` resolves
ColE1 retardation to ~13.5. Amplitude 0 and amplitude 60 give the two
SPEC13 endpoints — retardation 1.2 (literature D_eff ~3.3e-11 m2/s) and
~50.2 (the pI-derived shipped value).

- **Array 1:** `f0e1c478-0f49-49af-9a2b-b0b96524f248`
  (`gutibm-receptor-selection-v1-stage-d-retardation`, array size 18) —
  amplitude {0, 60} x targets {0.25%, 1%, 5%} x seeds {20261009, 20261011,
  20261013}
- **Array 2:** `11eb5ef8-3fac-4f39-a0a8-2aaf640eb255`
  (`gutibm-receptor-selection-v1-stage-d-retardation-mid`, array size 6) —
  amplitude 15 x targets {1%, 5%} x the same seeds; same-seed amp-15
  anchors, since the deployed amp-15 grid ran 1%/5% on other seeds
- **Queue:** `gutibm-gpu-campaign`
- **Job definition:** `gutibm-cuda-campaign:14` — image digest
  `sha256:fafc77ee802882bee6e47321515057aed20c0326d622bfb287862ef15304f748`,
  timeout 7200 s
- **Execution source SHA:** `29fd563ce73b28059cb439694e92466f451ddb4a`
  (unchanged — inputs generated in a detached worktree at that commit)
- **Inputs:**
  `s3://gutibm-inputs-994254241749/receptor-selection-v1/29fd563c/stage-d-retardation/inputs/<i>/input.json`
  and `.../stage-d-retardation-mid/inputs/<i>/input.json`, plus
  `manifest.json` at each prefix root
- **Outputs:**
  `s3://gutibm-outputs-994254241749/receptor-selection-v1/29fd563c/stage-d-retardation/outputs/<i>/output.h5.gz`
  and `.../stage-d-retardation-mid/outputs/<i>/output.h5.gz`
- **Approval record:** `approval_stage_D_retardation.json`
- **Generators:** `prepare_stage_D_retardation.py`,
  `prepare_stage_D_retardation_mid.py` (ran in a `29fd563c` worktree)

## Design

Every cell of the 3-amplitude x 3-target grid pairs within seed against
the refine array's same-seed plasmid-free null. The existing refine runs
supply the amp-15 row at 0.25% (indices 3-5), the zero-lysis carriers
(0-2), and the nulls (18-20). Nulls and carriers need no amplitude
re-runs: a plasmid-free null releases no colicin, and a P=0 carrier emits
no toxin, so the mucin-charge amplitude does not act on either.

If the gradient moves between amplitudes 0 and 60, the Stage D response
measures the transport parameter, not the lysis prior — the documented
confound this probe exists to resolve.

## Constancy verification

All 24 generated inputs were compared field-by-field against the deployed
stage-d producer inputs downloaded from S3. After removing `seed`,
`_campaign.run_id`, `_campaign.paired_seed`, `sos_lysis_prob`,
`bacteriocin.mucin_charge.amplitude`, and the `_campaign.axes` fields,
every input is byte-identical (0 mismatches outside `_campaign`).

## Expected analysis

Scratch combined manifest of 33 runs (9 refine controls/anchors + 18
endpoint arms + 6 mid fills); paired contrasts reuse each run's same-seed
plasmid-free null and the shared `t_common` window. Readout: per-target
per-seed delta slope at amplitudes 0 / 15 / 60.
