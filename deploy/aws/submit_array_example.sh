#!/usr/bin/env bash
# Example: submit an AWS Batch array job (docs/AWS_BATCH.md Phase 3).
# Fill in account-specific names; does not create infrastructure.
#
# Pre-upload per-index inputs:
#   s3://.../jobs/0/input.json ... s3://.../jobs/$((ARRAY_SIZE-1))/input.json
# entry.sh resolves PREFIX + AWS_BATCH_JOB_ARRAY_INDEX at runtime.

set -euo pipefail

JOB_NAME="${JOB_NAME:-gutibm-stage3-kd}"
JOB_QUEUE="${JOB_QUEUE:-gutibm-gpu-spot}"
JOB_DEFINITION="${JOB_DEFINITION:-gutibm-cuda:1}"
ARRAY_SIZE="${ARRAY_SIZE:-12}"
INPUT_PREFIX="${INPUT_PREFIX:-s3://gutibm-inputs/stage3_kd_sweep/jobs}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-s3://gutibm-outputs/stage3_kd_sweep/jobs}"
CKPT_PREFIX="${CKPT_PREFIX:-s3://gutibm-outputs/stage3_kd_sweep/checkpoints}"

aws batch submit-job \
  --job-name "${JOB_NAME}" \
  --job-queue "${JOB_QUEUE}" \
  --job-definition "${JOB_DEFINITION}" \
  --array-properties "size=${ARRAY_SIZE}" \
  --container-overrides "{
    \"environment\": [
      {\"name\": \"INPUT_S3_PREFIX\", \"value\": \"${INPUT_PREFIX}\"},
      {\"name\": \"OUTPUT_S3_PREFIX\", \"value\": \"${OUTPUT_PREFIX}\"},
      {\"name\": \"CHECKPOINT_S3_PREFIX\", \"value\": \"${CKPT_PREFIX}\"},
      {\"name\": \"MPI_RANKS\", \"value\": \"1\"}
    ]
  }"

echo "Submitted ${JOB_NAME} array size=${ARRAY_SIZE}"
