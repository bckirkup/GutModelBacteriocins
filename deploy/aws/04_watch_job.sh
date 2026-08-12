#!/usr/bin/env bash
# Poll one Batch job until it finishes (or you Ctrl-C).
# Prefer gut-ibm-aws-status when status.json / prefixes are known.
#
#   bash deploy/aws/04_watch_job.sh <jobId>
#   CHECKPOINT_S3_PREFIX=s3://bucket/ckpt/0/ bash deploy/aws/04_watch_job.sh <jobId>
#   STATUS_S3_URI=s3://bucket/ckpt/0/status.json bash deploy/aws/04_watch_job.sh <jobId>

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -f "${SCRIPT_DIR}/env.sh" ]]; then
  echo "ERROR: env.sh is missing beside this script. Clone the repository and run it from the clone." >&2
  exit 1
fi
# shellcheck source=env.sh
source "${SCRIPT_DIR}/env.sh"

JOB_ID="${1:-}"
if [[ -z "${JOB_ID}" ]]; then
  echo "Usage: bash deploy/aws/04_watch_job.sh <jobId>" >&2
  exit 2
fi

have_status_cli=0
if command -v gut-ibm-aws-status >/dev/null 2>&1; then
  have_status_cli=1
elif python3 -c "import gut_ibm_tools.aws_batch_status" >/dev/null 2>&1; then
  have_status_cli=1
fi

print_snapshot() {
  if [[ "${have_status_cli}" -eq 1 ]]; then
    local args=(--region "${AWS_REGION}")
    if [[ -n "${STATUS_S3_URI:-}" ]]; then
      args+=(--status-uri "${STATUS_S3_URI}")
    fi
    if [[ -n "${CHECKPOINT_S3_PREFIX:-}" ]]; then
      args+=(--checkpoint-prefix "${CHECKPOINT_S3_PREFIX}")
    fi
    if [[ -n "${OUTPUT_S3_PREFIX:-}" ]]; then
      args+=(--output-prefix "${OUTPUT_S3_PREFIX}")
    fi
    if [[ -n "${AWS_BATCH_JOB_ARRAY_INDEX:-}" ]]; then
      args+=(--array-index "${AWS_BATCH_JOB_ARRAY_INDEX}")
    fi
    if command -v gut-ibm-aws-status >/dev/null 2>&1; then
      gut-ibm-aws-status "${JOB_ID}" "${args[@]}" || true
    else
      python3 -m gut_ibm_tools.aws_batch_status "${JOB_ID}" "${args[@]}" || true
    fi
    return 0
  fi
  aws batch describe-jobs --jobs "${JOB_ID}" --region "${AWS_REGION}" \
    --query 'jobs[0].{status:status,reason:statusReason,log:container.logStreamName}' \
    --output text
}

echo "Watching job ${JOB_ID} (Ctrl-C to stop polling; job keeps running)"
while true; do
  echo "---- $(date -u +%Y-%m-%dT%H:%M:%SZ) ----"
  print_snapshot
  STATUS="$(aws batch describe-jobs --jobs "${JOB_ID}" --region "${AWS_REGION}" \
    --query 'jobs[0].status' --output text)"
  case "${STATUS}" in
    SUCCEEDED)
      echo "OK: job succeeded"
      echo "Download: aws s3 cp s3://${OUTPUT_BUCKET}/practice/smoke_gpu/output.h5.gz ./output.h5.gz"
      exit 0
      ;;
    FAILED)
      echo "FAILED: status=${STATUS}" >&2
      echo "Open Batch → Jobs → this job → CloudWatch logs in the AWS Console." >&2
      echo "If status.json shows spot_interruption or Batch reason is Host EC2*, reclaim likely." >&2
      exit 1
      ;;
    SUBMITTED|PENDING|RUNNABLE|STARTING|RUNNING|*)
      sleep 20
      ;;
  esac
done
