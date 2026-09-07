# Single-source transport assay v1

This intrinsic assay isolates mucin-charge-dependent colicin transport around one controlled release source. It does **not** enable metabolism, select a mucin amplitude, authorize the separate ecological refinement, or release Stage D.

## Correct physical contract

Each run uses a 100 µm cube, 2 µm grid, 1 h simulated duration, one GPU, one Message Passing Interface (MPI) rank, and only the `bacteriocin` and `receptor` fixes. One source cell and three bystanders all have `mu_max=0` and `BtuB=0`; producer arms give only the source cell `ColE1`, while nulls remove it. `kd_corrinoid_btuB=1e-4`, `b12_initial_conc=1e-3`, and `burst_release_tau=300 s`. Grids are written every 120 s and provenance every step.

`/run_provenance/chemistry_placement = device` is the required physical GPU placement. This assay has no metabolism, so it must not set `metabolism.uptake_limit=delivery` and must not require `device_delivery`. `host` and `host_forced_delivery` are not accepted.

## Adaptive design (issue 420)

The deployment package has **18 explicit jobs**:

- producer amplitudes `0, 15, 20, 30, 60` × three paired seeds = 15;
- three toxin-free nulls, one per seed, at amplitude 60.

All IDs begin `SS_amp`. Source position, exact lysis time, and snapshot times are paired across all five producer amplitudes. Missing any arm or paired time blocks the assay.

The original endpoint gate is preserved without reinterpretation: at every paired time and seed, `r50` must decrease across `0 > 15 > 60` by at least 0.5 µm per endpoint interval, `r90` by at least 0.1 µm, and the 0–10 µm near-field mean must strictly increase. The added interpolation criterion requires `r50/r90` to be nonincreasing and near field nondecreasing across `0,15,20,30,60`, using only a `1e-9` numeric tolerance. The endpoint margins are **not** imposed on adjacent interpolation pairs.

Shell averaging, 40 µm support, periodic x/y minimum-image distance, nonperiodic z, r50/r90 definitions, and the exact `event_time_s < t <= event_time_s + 1500 s` window are fixed in `assay_contract.json`.

## Historical attempt and corrected adjudication

The original 12-job `d7b16c30…` attempt and its committed `analysis/assay_gate.json`, `assay_metrics.json`, and CSV remain unchanged as historical evidence. Its only recorded blocker was the inherited, scientifically irrelevant expectation `device_delivery`; the outputs actually recorded the correct physical placement `device`. All scientific criteria then declared passed. See `analysis/ASSAY_RUN_2026-09-07_ADJUDICATION.md` for the separate adjudication.

Those outputs support the original endpoint transport conclusion, but they do not contain amplitudes 20/30. Therefore the 18-job adaptive assay remains formally pending. Stage D remains blocked and its amplitude selection is unchanged.

## Separate ecological refinement

Issue 420 also defines a distinct **12-job ecological refinement**: amplitudes 15/20/30 × three seeds plus three same-image nulls. It retains ecological `C_amp` identities and is not generated, preflighted, or submitted by this intrinsic `SS_amp` package. Do not use this package to submit it.

## Planning and deployment

Planning generation is safe in a codeload tree and uses explicit identity placeholders:

```bash
python3 prepare_assay.py --clean
python3 preflight_assay.py
```

Planning preflight passes with warnings because it cannot prove a deployable commit/image. Deployment requires a normal clean git checkout whose exact HEAD descends from PR417 event-time commit `f9a908e0eb5323cfaeea82501a8ca16b1a95910c`, plus an immutable image URI:

```bash
EXECUTION_SOURCE_SHA=$(git rev-parse HEAD)
IMAGE_URI='<repository>@sha256:<64 hex>'
python3 prepare_assay.py --clean --deployment   --execution-source-sha "$EXECUTION_SOURCE_SHA" --image-digest "$IMAGE_URI"
python3 preflight_assay.py --deployment --execution-source-sha "$EXECUTION_SOURCE_SHA"
python3 aws_commands_assay.py   --execution-source-sha "$EXECUTION_SOURCE_SHA" --image-uri "$IMAGE_URI"   --input-prefix 's3://<bucket>/single-source-assay-v1/inputs'   --output-prefix 's3://<bucket>/single-source-assay-v1/outputs'   --job-queue '<queue>' --job-definition '<definition>'
```

`aws_commands_assay.py` only prints an upload command and an 18-element AWS Batch array command compatible with `deploy/aws/entry.sh`; it never calls AWS. It fixes `REQUIRE_GPU=1`, `MPI_RANKS=1`, and a 3600 s attempt timeout and rejects ecological paths/IDs. An optional authorization JSON may set `authorize_single_source_assay_submission=true` and match the exact SHA, image URI, and `jobs=18`.

Returned outputs belong only at `generated/results/<array-index>/output.h5.gz` (uncompressed `.h5` is accepted for analysis). Run `python3 analyze_assay.py`; exit codes are 0 pass, 1 scientific fail, and 2 blocked.
