# Gate C adjudication — 2026-09-27 (issue #420)

Verdict: **intrinsic assay PASS; ecological refinement READY_FOR_APPROVAL.**
Stage D inputs are prepared at the shared identity and remain unsubmitted —
they require a stage-specific operator approval binding SHA + digest +
amplitude.

## Identity under adjudication (shared by both streams)

- execution_source_sha: `29fd563ce73b28059cb439694e92466f451ddb4a`
  (merge of PR #427, main at submission)
- container_image_digest:
  `994254241749.dkr.ecr.us-east-1.amazonaws.com/gutibm@sha256:fafc77ee802882bee6e47321515057aed20c0326d622bfb287862ef15304f748`
- Job definitions: `gutibm-cuda-campaign:12` (3600 s, assay) and
  `gutibm-cuda-campaign:14` (7200 s, refinement), both pinning that digest.

## Execution evidence

- Assay array `cd7fc9b3-e6b6-43e6-ab23-a98ef1b1e4ee`: 18/18 children
  SUCCEEDED; every `outputs/<i>/output.h5.gz` present under
  `s3://gutibm-outputs-994254241749/single-source-assay-v1/29fd563c/`.
- Refinement array `2355877a-a3a0-4765-8863-b84346868e0d`: 12/12 children
  SUCCEEDED; every `outputs/<i>/output.h5.gz` present under
  `s3://gutibm-outputs-994254241749/adaptive-ecology-refinement-v1/29fd563c/`.
- Deployment regeneration of both package manifests reproduces the submitted
  S3 manifests exactly (all `input_sha256` identical), so the committed
  generators are the source of what ran.
- Per-run `/run_provenance` was read on all 30 outputs:
  - assay: `git_sha=29fd563c…`, `chemistry_placement=device`,
    `mpi_rank_count=1`, `termination_cause=horizon_reached` at t=3600 s;
    no metabolism, no `host_forced_delivery`.
  - refinement: `git_sha=29fd563c…`, `chemistry_placement=device_delivery`,
    `termination_cause=dysbiosis_guard` — producers halt at t=20400 s,
    plasmid-free nulls at t=17100 s, both before the 21600 s horizon.
    Paired contrasts therefore use `t_common = min(t_end_producer,
    t_end_null)` per seed, as declared; the guard halt is a scientific
    termination condition, not a Batch failure.

## Intrinsic verdict — `single_source_transport_gate: true` (PASS)

`experiments/single_source_transport_assay_v1/analysis/assay_gate.json`:

- all 18 outputs returned and authenticated (SHA, placement `device`, rank,
  seed/amplitude, Kd 1e-4, B12 1e-3, burst tau 300 s, ColE1 identity,
  schedules);
- exactly one strain-1 lysis event per producer run; source position and
  event time identical across all five amplitudes of each seed;
- 36 shared paired snapshot times across the three seeds (12 per seed); at
  every paired time and seed: r50 ordered 0>15>60 by
  >=0.5 um, r90 ordered 0>15>60 by >=0.1 um, near-field mean strictly
  increasing 0<15<60, and interpolation monotonicity across 0/15/20/30/60
  within numeric tolerance;
- the three toxin-free nulls hold |bacteriocin_BtuB| <= 1e-15 everywhere.

## Ecological verdict — `C_transport_gate: READY_FOR_APPROVAL`

`analysis/refinement_gate.json` + `analysis/amplitude_selection.json`:

- 12/12 outputs readable, zero authentication violations;
- per-seed producer-null delta slopes positive at every amplitude;
- **selected amplitude: 15.0, basis "qualified"** — median delta slope
  0.300 log10/h, median kills/lysis 10.04, median final susceptible 1317
  (all criteria pass); amplitude 20 also qualifies (9.04 kills/lysis, inside
  the 10% tie band) but loses the kills-per-lysis tie-break; amplitude 30 is
  disqualified by the susceptible floor (median 636, seed minimum 574);
- `C_transport_gate` remains `false` in the JSON by design: the analyzer
  never auto-promotes — promotion requires the operator approval file binding
  `authorize_c_transport_gate=true` + SHA + digest + selected amplitude.

## Tooling gap found and fixed

`analyze_assay.py` wrote `execution_source_sha`/`container_image_digest` into
`assay_metrics.json` but not into `assay_gate.json`, while
`analyze_refinement.py`'s `intrinsic_gate()` requires exactly those fields on
the gate file to bind the two streams to one identity — so the refinement
reported `BLOCKED_MISSING_INTRINSIC_GATE` despite clean evidence. The assay
analyzer now stamps both fields onto the gate object. No gate criterion,
tolerance, or selection rule changed; the pre-patch verdict was already PASS
with zero blockers (the stamp is metadata binding, not a criterion).

## PR #418 caveats — verified before trusting Stage D

- `preflight.py` does **not** bind `promotion_inputs`/`--selected-*` to the
  decision record: `check_decision_record()` validates the record's contents,
  and `check_stage_d()` only requires the 15 D inputs to share one
  (Kd, B12, amplitude) tuple — not that it equals (1e-4, 1e-3, 15).
  Mitigation applied here: every generated Stage D input was read back and
  confirmed `kd_corrinoid_btuB=1e-4`, `b12_initial_conc=1e-3`,
  `bacteriocin.mucin_charge.amplitude=15.0`; the stage_D manifest's
  `promotion_inputs` records the same (status `USER_RECORDED_PROMOTION`).
- `analyze.py` authenticates per-output `run_provenance/git_sha` against the
  per-input expected SHA (`execution_source_sha_match`), i.e. the manifest
  chain, not `campaign_manifest.json` alone. It does **not** authenticate the
  image digest per output — `run_provenance` does not record a container
  digest, so the digest binding lives in the digest-pinned job definition
  plus the binary-stamped git_sha; the digest on paired rows is copied from
  the campaign manifest for reference only. This is an inherent artifact
  limit, not a bypass: no output field exists to compare against.

## Stage D preparation (pending approval — NOT submitted)

- Regenerated `experiments/receptor_selection_campaign_v1` in deployment mode
  on the same identity (`--deployment --execution-source-sha 29fd563c…
  --image-digest sha256:fafc77ee… --selected-kd 1e-4 --selected-b12 1e-3
  --selected-amplitude 15.0`): 15 jobs under `generated/stage_D/jobs/`
  (12 ColE1 carriers: targets {0, 1%, 2%, 5%} x 3 seeds; 3 same-image
  plasmid-free nulls).
- `preflight.py --deployment --execution-source-sha 29fd563c…` → PASS,
  65 runs checked, zero errors.
- `aws_commands.py --stage D` remains refused without an approval JSON
  carrying `C_transport_gate=true` + matching baseline/SHA/digest — the
  submission path is fail-closed as designed.
- Per package convention, deployment-regenerated files under `generated/`
  are left uncommitted.

## Deliberately not done (non-goals honored)

- No Stage D submission (needs the stage-specific operator approval).
- No new image, no HEAD advance, no Gate C re-runs.
- `approval.json` and the campaign decision record untouched — the record
  still shows `single_source_transport_assay_v1: PENDING_ADAPTIVE_RUNS` and
  `D_release: BLOCKED`, which is what the operator approval step updates.
