# Adaptive ecology refinement v1

This package deploys the **ecological** Gate C evidence stream of issue 420: a same-image producer-versus-plasmid-free-null contrast at mucin-charge amplitudes 15, 20, and 30 that selects the Stage D amplitude. It does **not** supply the intrinsic transport stream (that is `single_source_transport_assay_v1`, `SS_amp` identities) and it never reads or modifies the historical `receptor_selection_campaign_v1` Stage C record.

## The 12 jobs

- Nine producer runs: amplitudes `15, 20, 30` × three seeds (`20260911`, `20260913`, `20260917`), IDs `C_amp<A>_producer_s<seed>`.
- Three plasmid-free nulls, one per seed, pinned to amplitude 15 (a null releases no colicin, so one null amplitude is a valid shared control and is byte-identical across the axis), IDs `C_amp15_plasmid_free_null_s<seed>`.

## Fixed settings

The ecological Stage C configuration is reproduced exactly: 100 µm cube, 2 µm grid, `total_time=21600 s`, `bio_dt=60 s`, `metabolism.uptake_limit=delivery`, `carbon.boundary_conc=0.05`, `kd_corrinoid_btuB=1e-4`, `kd_colicinE_btuB=5e-7`, `b12_initial_conc=1e-3`, fixes `metabolism, bacteriocin, receptor, mechanics`, two strains of 60 cells at `mu_max=5.5e-4`, one GPU, one MPI rank, 7200 s attempt timeout, and the contract HDF5 schedule (summary 1, agents 10, grid 60, lineage/genome 0, provenance 10, `grid_species=["bacteriocin_BtuB"]`).

`chemistry_placement` must be **`device_delivery`**: this stream funds uptake by delivered nutrient flux, so `device` (the intrinsic assay's physical placement), `host`, and `host_forced_delivery` are all refused.

## Command sequence

Planning (safe in a codeload tree, explicit identity placeholders):

```bash
python3 prepare_refinement.py --clean
python3 preflight_refinement.py
```

Deployment, only after the image exists on the merge commit:

```bash
EXECUTION_SOURCE_SHA=$(git rev-parse HEAD)
IMAGE_URI='<repository>@sha256:<64 hex>'
python3 prepare_refinement.py --clean --deployment \
  --execution-source-sha "$EXECUTION_SOURCE_SHA" --image-digest "$IMAGE_URI"
python3 preflight_refinement.py --deployment --execution-source-sha "$EXECUTION_SOURCE_SHA"
python3 aws_commands_refinement.py \
  --execution-source-sha "$EXECUTION_SOURCE_SHA" --image-uri "$IMAGE_URI" \
  --input-prefix 's3://<bucket>/adaptive-ecology-refinement-v1/inputs' \
  --output-prefix 's3://<bucket>/adaptive-ecology-refinement-v1/outputs' \
  --job-queue '<queue>' --job-definition '<definition>'
```

`aws_commands_refinement.py` **only prints** the upload and 12-element AWS Batch array commands for review — it never executes any AWS command. It rejects intrinsic (`SS_amp`/`single-source`) and stale Stage C paths and IDs, and requires the deployment preflight to pass first.

After outputs land at `generated/results/<array-index>/output.h5.gz` (plain `.h5` also accepted):

```bash
python3 analyze_refinement.py \
  [--intrinsic-gate-file intrinsic_gate.json] [--approval-file approval.json] \
  [--results-dir analysis]
```

The gate reports `BLOCKED_MISSING_OUTPUTS`, `BLOCKED_AUTHENTICATION`, `BLOCKED_MISSING_INTRINSIC_GATE` (until the intrinsic assay's `single_source_transport_gate` is shown PASS at the same execution SHA and image digest), `BLOCKED_NO_QUALIFYING_AMPLITUDE`, `READY_FOR_APPROVAL`, or `PASS`. **Stage D stays BLOCKED until this gate is PASS**; a fallback-retained amplitude never satisfies it, and the analyzer never writes an approval file. Exit codes: 0 PASS, 1 READY_FOR_APPROVAL, 2 blocked.

## Relationship to other packages

- `single_source_transport_assay_v1`: intrinsic transport stream (`SS_amp`, `chemistry_placement=device`). This package refuses those identities and paths.
- `receptor_selection_campaign_v1`: historical staged campaign at amplitudes 0/15/60 on an older image. **Its generated files are an immutable historical record and are left untouched.**
