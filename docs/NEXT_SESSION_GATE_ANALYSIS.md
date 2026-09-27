# Gate C analysis and Stage D unblock — follow-up session prompt

You are working in the bckirkup/GutModelBacteriocins repository on the Gate C
analysis phase of the receptor-selection campaign (issue 420). Both Gate C
evidence streams have now been executed on AWS Batch at one shared identity.
Your job is to analyze the outputs, adjudicate the two gates, and — only if
both pass — prepare Stage D. Do NOT submit Stage D jobs without an explicit
stage-specific operator approval.

## Settled inputs (do not re-derive)

- Execution source SHA (both campaigns): `29fd563ce73b28059cb439694e92466f451ddb4a`
  (merge commit of PR #427, `main` at time of submission).
- Container image: `994254241749.dkr.ecr.us-east-1.amazonaws.com/gutibm@sha256:fafc77ee802882bee6e47321515057aed20c0326d622bfb287862ef15304f748`
  (built from that SHA, binary embeds the SHA without `-dirty`).
- Job definitions: `gutibm-cuda-campaign:12` (3600 s) for the assay,
  `gutibm-cuda-campaign:14` (7200 s) for the refinement; both pin the digest above.
- Intrinsic assay: 18-job array `cd7fc9b3-e6b6-43e6-ab23-a98ef1b1e4ee`
  (`gutibm-single-source-assay-v1-adaptive`), inputs/outputs under
  `s3://gutibm-outputs-994254241749/single-source-assay-v1/29fd563c/outputs/<i>/output.h5.gz`
  (manifest at `s3://gutibm-inputs-994254241749/single-source-assay-v1/29fd563c/manifest.json`).
- Ecological refinement: 12-job array `2355877a-a3a0-4765-8863-b84346868e0d`
  (`gutibm-adaptive-ecology-refinement-v1`), outputs under
  `s3://gutibm-outputs-994254241749/adaptive-ecology-refinement-v1/29fd563c/outputs/<i>/output.h5.gz`.
- The 2026-09-07 assay attempt (SHA `d7b16c30...`, digest `sha256:78657698...`)
  is adjudicated scientifically supportive but cannot authorize the gate —
  it lacks amplitudes 20 and 30. It is historical evidence only; see
  `experiments/single_source_transport_assay_v1/analysis/ASSAY_RUN_2026-09-07_ADJUDICATION.md`.
- Placement contracts: the intrinsic assay must authenticate
  `chemistry_placement=device` (it runs bacteriocin + receptor only — do NOT
  add metabolism to change the label; `host`/`host_forced_delivery` are
  forbidden). The ecological refinement must authenticate
  `chemistry_placement=device_delivery` (it legitimately runs
  `metabolism.uptake_limit=delivery`).
- Every ecological/Stage D run halts on the dysbiosis guard near t=20000 s
  against the 21600 s horizon; paired contrasts use the shared calendar window
  ending at `t_common = min(t_end_producer, t_end_null)`.

## Deliverable (exactly one)

A gate adjudication for Gate C: the `single_source_transport_gate` verdict from
the intrinsic assay and the `C_transport_gate` verdict from the ecological
refinement, committed as analysis artifacts under each package's `generated/`
or `analysis/` tree, plus Stage D preparation only if both gates pass.

## Steps

1. Confirm all 30 array children completed and every
   `<i>/output.h5.gz` exists under both output prefixes; record per-run
   `run_provenance/git_sha` and `chemistry_placement` during analysis.
2. Download assay outputs into
   `experiments/single_source_transport_assay_v1/generated/results/<i>/`
   and run
   `python3 analyze_assay.py --results-root generated/results --output-dir analysis`.
   It writes `analysis/assay_gate.json` carrying `single_source_transport_gate`.
   Treat non-PASS as blocking; report the blocker list verbatim.
3. Download refinement outputs into
   `experiments/adaptive_ecology_refinement_v1/generated/results/<i>/`
   and run
   `python3 analyze_refinement.py --results-root generated/results --results-dir analysis --intrinsic-gate-file <assay analysis/assay_gate.json>`.
   Expected statuses: PASS / READY_FOR_APPROVAL / BLOCKED_*
   (BLOCKED_MISSING_OUTPUTS, BLOCKED_AUTHENTICATION,
   BLOCKED_MISSING_INTRINSIC_GATE, BLOCKED_NO_QUALIFYING_AMPLITUDE).
   The analyzer never writes an approval file, and a fallback-retained
   amplitude never satisfies the gate — only the predeclared selection rule
   resolving to an amplitude qualifies. Exit 1 (READY_FOR_APPROVAL) means the
   science checks passed and the gate awaits the recorded operator approval
   naming this exact SHA, digest, and selected amplitude; that is a working
   fail-closed guard, not an analysis failure.
4. If both gates are satisfied at the shared SHA/digest above, record the
   selected Stage D amplitude and proceed to Stage D prep in
   `experiments/receptor_selection_campaign_v1` (15-job inputs in
   `generated/stage_D/`; authoritative selections: `kd_corrinoid_btuB=1e-4`,
   `b12_initial_conc=1e-3`, and the amplitude selected by the refinement gate).
   Regenerate Stage D in deployment mode on the SAME execution SHA and image
   digest — do not build a new image and do not advance HEAD for this.
5. Stage D submission itself requires a stage-specific approval file that
   binds the SHA, digest, and selected amplitude. Prepare it; do not submit.

## PR #418 review caveats to carry forward (verify before trusting Stage D)

- `experiments/receptor_selection_campaign_v1/preflight.py`: deployment
  preflight must validate generated `promotion_inputs` and every Stage D input
  against the authoritative decision record selections (Kd 1e-4, B12 1e-3,
  selected amplitude) — any `--selected-*` CLI value otherwise passes the
  static decision checks and could launch the wrong campaign. Verify the
  current preflight performs this match before running it for real.
- `experiments/receptor_selection_campaign_v1/analyze.py`: same-image paired
  metrics must authenticate the container image per output — read
  `run_provenance/container_image_digest` (or the passed pinned digest) from
  both treatment and control outputs and require equality with the deployment
  manifest; do not populate same-image metadata solely from
  `campaign_manifest.json`.

## Report immediately if

- Any output authenticates to a different SHA, a different image digest,
  `host`, `host_forced_delivery`, or (assay only) `device_delivery`.
- Any array child is missing its `output.h5.gz` or its provenance.
- The refinement analyzer reports a BLOCKED_* status.

## Non-goals

- No Stage D submission (prep only).
- No new image build, no new execution SHA, no re-runs of Gate C.
- No changes to `approval.json`; historical approvals are not authorization.

## Stop when

Both gate verdicts are committed to the repo (analysis artifacts), the
selected amplitude is recorded, and — if the gates pass — the Stage D
deployment-prepared inputs exist pending explicit operator approval.
