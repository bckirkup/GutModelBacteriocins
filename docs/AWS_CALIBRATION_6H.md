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

Provenance is dumped every 360 biological steps. The writer accumulates kill
events between dumps, so this cadence does not discard deaths. Closed Tier-2
restarts use a tighter cadence so a Spot reclaim can be followed by a durable
resume. The resulting throughput measurement includes the provenance work that
an equivalent production burn-in will perform and is therefore directly
transferable to that configuration. It is not a provenance-free benchmark.

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

before the writer's metadata overhead and gzip-4 compression. The measured
steady-state rate is approximately 25.3 wall seconds per biological step.
The previous 360-step cadence therefore took approximately 2.53 wall hours
between local checkpoint writes, which was longer than the observed Spot host
lifetime. The calibration now uses:

```text
restart.interval_steps = 30
30 × 25.3 s = 759 s = 12.65 wall minutes between local writes
759 s + 300 s sync poll = 1059 s = 17.65 wall minutes maximum
durable-checkpoint exposure
```

The 300-second term matters because the entrypoint's S3 sync loop can take up
to one poll period to upload a newly written local checkpoint. The observed
host lifetime was approximately 37 minutes, so the maximum exposure is less
than half that lifetime rather than exceeding it.

With `restart.interval_steps = 30`, the unchanged two-day `total_time` contains:

```text
172800 / 60 = 2880 steps
2880 / 30 = 96 closed restarts
```

The entrypoint separates resume retention from coarse archive retention:

```text
CHECKPOINT_RETAIN_NEWEST_K = 2
CHECKPOINT_ARCHIVE_INTERVAL_STEPS = 360
```

The newest two checkpoints preserve a resume point if the newest upload is
interrupted. Every 360th step remains as a permanent fork origin. Older
non-archive checkpoints are deleted from S3 only after the new checkpoint and
its `latest.json` pointer have uploaded successfully. The retention path
protects the object named by `latest.json`, and a denied or failed delete is
logged while leaving the simulation alive.

For a complete two-day run, the eight archive steps are
`360, 720, ..., 2880`. The final step overlaps the newest-two set, so nine
objects are retained:

```text
8 archive checkpoints + 2 newest checkpoints - 1 overlap = 9 objects
9 × 4.40 GB = 39.6 GB uncompressed grid payload upper bound
```

At the measured rate, six wall hours reaches approximately 854 steps and
therefore produces up to 28 checkpoints. The full archive set is not yet
reached; steps 360 and 720 plus the newest two are retained, for at
most four objects or **17.6 GB uncompressed grid payload**. These figures are
upper bounds based on the 4.40 GB uncompressed grid estimate per checkpoint;
no compressed checkpoint size has yet been measured. The retained-count
calculation is:

```text
archive_count = floor(checkpoint_steps / 360)
retained_count = min(2, checkpoint_count) + archive_count
                 - archives_overlap_newest_two
```

The entrypoint prunes older uploaded local files and keeps the newest local
checkpoint, so local scratch is not a 96-checkpoint accumulation. Actual
compressed sizes must be recorded from S3; they cannot be inferred reliably
from the uncompressed estimate.

The Batch job role must have `s3:DeleteObject` on the output bucket for S3
retention to reduce the archive. The setup scripts already include that
permission in the `gutibm-s3-access` inline policy for the input and output
buckets. The scoped monitoring role cannot inspect the live job role policy,
so Benjamin should verify or re-run the setup script with administrator
credentials before relying on pruning. If deletion is denied, the entrypoint
retains the objects and continues the simulation rather than failing the run.

## Measured calibration result

The first calibration attempt ran on the GPU path and CloudWatch confirmed:

```text
GPU: ON (device 0)
```

Two progress events measured:

```text
Step 1  → Step 60: 1492.104 wall seconds / 59 steps
             = 25.3 wall seconds per biological step
             = 2.37 simulated seconds per wall second
```

The resulting projections are:

```text
6 wall hours:  approximately 854 steps = 14.23 simulated hours
2 simulated days (2880 steps): approximately 20.2 wall hours
7 simulated days (10080 steps): approximately 70.8 wall hours
```

The first attempt was reclaimed after approximately 37 wall minutes:

```text
Host EC2 (instance i-0103a811de8d956c4) terminated.
```

Batch classified the reclaim correctly and automatically started a retry
without resubmission. This proves the submit-time ten-attempt retry policy and
real `Host EC2*` classification. No checkpoint existed before that reclaim, so
checkpoint-backed resume remains unproven: there was no
`Resuming from restart artifact` line, no `latest.json`, and
`resume_from_checkpoint` remained false.

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
change before using this as a resume calibration. Note that editing
`deploy/aws/retry_strategy.json` does **not** change a live job definition:
nothing re-registers it, and `gutibm-cuda-campaign:1` still carries
`attempts: 2`. Pass the policy at submit time instead, which needs no IAM
change:

```bash
aws batch submit-job ... --retry-strategy file://deploy/aws/retry_strategy.json
```

A submit-time strategy overrides the job definition for that job only.

## Check the image before submitting

The job definition resolves the floating `gutibm:cuda` tag. A floating tag
silently going stale is a real failure mode here: it was last rebuilt on
2026-07-31 while the toxin, units, aliasing, siderophore, FMM, immigration and
restart fixes all landed afterwards, so a run submitted against it would have
simulated the pre-repair model and produced plausible output. Confirm what the
tag currently points at, and rebuild if it predates the science you intend to
measure:

```bash
aws ecr describe-images --repository-name gutibm \
  --query 'sort_by(imageDetails,&imagePushedAt)[].[imagePushedAt,imageTags]' --output json
```

Push an immutable `cuda-<commit>` tag alongside `cuda` so the image that ran a
given result stays identifiable after the floating tag moves.

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
  --retry-strategy file://deploy/aws/retry_strategy.json \
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

## Watch, and how to prove resume

Watch the job and its S3 heartbeat:

```bash
CHECKPOINT_S3_PREFIX="${CHECKPOINT_PREFIX}" \
  bash deploy/aws/04_watch_job.sh "${JOB_ID}"
```

Wait until `latest.json` exists and record its step and URI, so there is a
recorded pre-interruption position to compare against:

```bash
aws s3 cp "${CHECKPOINT_PREFIX}latest.json" - \
  | tee calibration_latest_before.json
```

Do **not** use `aws batch terminate-job` to test resume. A user-initiated
termination is not a host reclaim, so the `onReason: "*" -> EXIT` rule applies
and Batch does not retry: the run simply ends and nothing is proven. Prefer, in
order:

1. a **real Spot reclaim** during the run — the genuine article, and the only
   thing that exercises Batch's reclaim classification and the retry budget;
2. failing that, a **short second job pointed at the same
   `CHECKPOINT_S3_PREFIX`**, which exercises `resolve_resume_uri` and the
   restore path but *not* the reclaim classification. Report it as such; it is
   not an end-to-end reclaim proof. Give it a small
   `attemptDurationSeconds` so it cannot become another long run.

Evidence of a real resume is:

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
