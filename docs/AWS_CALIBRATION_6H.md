# Six-hour Spot calibration run

This run is a calibration, not a campaign result. It measures the real Stage 3
throughput, exercises checkpoint resume end to end, and gives a first structural
look to compare with the finite-colony toxin thresholds in
[T2_COLONY_CHALLENGE.md](T2_COLONY_CHALLENGE.md).

Do not change the domain, grid spacing, strains, FMM settings, or GPU settings.
The input is otherwise the `3a_burnin.json` campaign configuration:

```text
experiments/diversity_campaign/stage3_campaign/calibration_6h.json
```

The six-hour wall-clock limit is intentional. `total_time` remains 172800 s
(two simulated days), so the run should be stopped by the Batch timeout or a
deliberate termination rather than by a shortened biological horizon.

Provenance is dumped every 360 biological steps, matching the closed-restart
cadence. The writer accumulates kill events between dumps, so this cadence does
not discard deaths. The resulting throughput measurement includes the
provenance work that an equivalent production burn-in will perform and is
therefore directly transferable to that configuration. It is not a
provenance-free benchmark.

## Storage before submission

The Stage 3 domain is:

```text
(0.002 / 2e-6) × (0.002 / 2e-6) × (0.0001 / 2e-6)
= 1000 × 1000 × 50
= 50,000,000 grid cells
```

The closed restart writer includes all 11 chemical species. The grid payload
alone is therefore approximately:

```text
50,000,000 × 11 × 8 bytes = 4.40 GB
```

before the writer's metadata overhead and gzip-4 compression. With
`restart.interval_steps = 360`, the unchanged two-day `total_time` contains:

```text
172800 / 60 = 2880 steps
2880 / 360 = 8 closed restarts
```

The immutable S3 checkpoint prefix can therefore retain up to approximately
**35.2 GB uncompressed grid payload** for a complete two-day run, before
metadata and compression. The six-hour wall limit will usually produce fewer
than eight artifacts; after the measured throughput is known, use:

```text
checkpoint_count = min(8, floor(21600 × steps_per_second / 360))
```

The entrypoint prunes older uploaded local files and keeps the newest local
checkpoint, so local scratch is not an eight-checkpoint accumulation. S3 keeps
all immutable steps. Actual compressed sizes must be recorded from S3; they
cannot be inferred reliably from the uncompressed estimate.

## Submit one calibration job

Run these commands from a clone with AWS credentials. They are instructions
only; this document was prepared without submitting or describing any AWS job.

The campaign stack uses:

```text
Region:         us-east-1
Queue:          gutibm-gpu-campaign
Job definition: gutibm-cuda-campaign
```

The retry budget must be the ten-attempt policy from the AWS retry-budget
change before using this as a resume calibration.

```bash
source deploy/aws/env.sh

export CALIBRATION_NAME="calibration_6h_seed1001"
export INPUT_URI="s3://${INPUT_BUCKET}/campaign/calibration_6h/input.json"
export OUTPUT_URI="s3://${OUTPUT_BUCKET}/campaign/calibration_6h/${CALIBRATION_NAME}/output.h5.gz"
export CHECKPOINT_PREFIX="s3://${OUTPUT_BUCKET}/campaign/calibration_6h/${CALIBRATION_NAME}/ckpt/"

aws s3 cp \
  experiments/diversity_campaign/stage3_campaign/calibration_6h.json \
  "${INPUT_URI}"

OVERRIDES="$(mktemp)"
cat > "${OVERRIDES}" <<EOF
{
  "environment": [
    {"name": "INPUT_S3_URI", "value": "${INPUT_URI}"},
    {"name": "OUTPUT_S3_URI", "value": "${OUTPUT_URI}"},
    {"name": "CHECKPOINT_S3_PREFIX", "value": "${CHECKPOINT_PREFIX}"},
    {"name": "CHECKPOINT_INTERVAL_SECONDS", "value": "300"},
    {"name": "MPI_RANKS", "value": "1"},
    {"name": "REQUIRE_GPU", "value": "1"}
  ]
}
EOF

JOB_ID="$(aws batch submit-job \
  --job-name "${CALIBRATION_NAME}" \
  --job-queue "${JOB_QUEUE_CAMPAIGN}" \
  --job-definition "${JOB_DEFINITION_CAMPAIGN}" \
  --timeout attemptDurationSeconds=21600 \
  --container-overrides "file://${OVERRIDES}" \
  --query jobId --output text)"
rm -f "${OVERRIDES}"
echo "Submitted ${JOB_ID}"
```

The input, output, checkpoint, and status locations are:

```text
Input:      ${INPUT_URI}
Output:     ${OUTPUT_URI}
Checkpoints: ${CHECKPOINT_PREFIX}step_NNNNNN.h5
Pointer:    ${CHECKPOINT_PREFIX}latest.json
Status:     ${CHECKPOINT_PREFIX}status.json
```

## Watch and deliberately interrupt

Watch the job and its S3 heartbeat:

```bash
CHECKPOINT_S3_PREFIX="${CHECKPOINT_PREFIX}" \
  bash deploy/aws/04_watch_job.sh "${JOB_ID}"
```

Before interrupting, wait until `latest.json` exists and record its step and
URI:

```bash
aws s3 cp "${CHECKPOINT_PREFIX}latest.json" - \
  | tee calibration_latest_before.json
```

Terminate the running job deliberately, before the six-hour timeout:

```bash
aws batch terminate-job \
  --job-id "${JOB_ID}" \
  --reason "calibration resume validation"
```

Batch should send the container through the same graceful-stop path used for a
Spot notice: GutIBM drains the current step, the entrypoint uploads the newest
closed restart, and the nonzero exit permits a retry. Continue watching the
same Batch job ID. Evidence of a real resume is:

1. the retry attempt's log contains `Resuming from restart artifact`;
2. the URI is the immutable step named by `latest.json`, not a newly invented
   step-zero file;
3. the final `status.json` has `"resume_from_checkpoint": true`;
4. the resumed progress begins at or after the recorded checkpoint step;
5. the final `latest.json` advances beyond the pre-interruption step.

Resume is a population-state continuation, not a bit-identical continuation:
the RNG is reseeded from the configured job seed.

## Completion checks

Record these from the final status, Batch attempts, logs, and S3:

- completed simulation step and simulated seconds;
- wall-clock elapsed seconds and `rate` (simulated seconds per wall second);
- mean wall seconds per biological step (`wall_seconds / completed_steps`);
- number of Batch attempts and whether the interruption was classified as
  `Host EC2*` or a deliberate termination;
- seconds of simulated time between the recorded pre-interruption step and the
  first resumed progress step;
- every immutable checkpoint key and its S3 `ContentLength`;
- the final checkpoint's step, time, SHA-256, and `latest.json` pointer;
- provenance groups under `provenance/step_*` in the live output HDF5 if an
  attempt reaches successful completion. Closed Tier-2 checkpoints intentionally
  contain agents and grid, not the live provenance event buffer. The 360-step
  dump is lossless because kill events accumulate until the writer clears them;
  see `src/io/hdf5_writer.cpp:567-624`.
- interpret the throughput as the rate for the production-equivalent
  provenance cadence, not as a provenance-disabled upper bound. For scale,
  #220 measured a 15.3% receptor-path diagnostic cost when provenance arrays
  were built every step; this calibration deliberately avoids that
  every-step-dump overhead.
- population count and producer count at the final checkpoint;
- whether any observed colony reaches the 113, 527, or 1361 producer
  thresholds from the T2 surface.

For S3 checkpoint sizes:

```bash
aws s3 ls --summarize --human-readable --recursive "${CHECKPOINT_PREFIX}"
aws s3api head-object \
  --bucket "${OUTPUT_BUCKET}" \
  --key "campaign/calibration_6h/${CALIBRATION_NAME}/ckpt/step_XXXXXX.h5"
```

The container's local restart directory is ephemeral and older uploaded files
are deliberately pruned. Capture the newest local file size while the attempt
is still running if container-level access is available; after termination,
S3 is the durable size record.

Build the colony catalog from the final immutable checkpoint. The checkpoint
contains the full agents layer even when the wall-clock timeout prevented the
live output file from being uploaded:

```bash
aws s3 cp "${CHECKPOINT_PREFIX}latest.json" calibration_6h_latest.json
FINAL_CHECKPOINT_URI="$(jq -r '.uri' calibration_6h_latest.json)"
aws s3 cp "${FINAL_CHECKPOINT_URI}" calibration_6h_checkpoint.h5
python -m gut_ibm_tools.colony calibration_6h_checkpoint.h5 \
  --agent-output calibration_6h_agents.csv \
  --colony-output calibration_6h_colonies.csv \
  --diagnostics-output calibration_6h_colony_diagnostics.csv
```

If the job reaches a successful completion and uploads `OUTPUT_URI`, inspect
that file separately for the enabled per-step provenance layer:

```bash
aws s3 cp "${OUTPUT_URI}" calibration_6h_output.h5.gz
gunzip -f calibration_6h_output.h5.gz
h5ls -r calibration_6h_output.h5 | grep '^.*provenance/'
```

Interpret the catalog as a structural reconnaissance only. The T2 thresholds
are mechanism thresholds; this six-hour calibration does not establish that a
full campaign reaches them.
