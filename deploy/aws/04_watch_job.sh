#!/usr/bin/env bash
# Poll one Batch job until it finishes (or you Ctrl-C).
#
#   bash deploy/aws/04_watch_job.sh <jobId>

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=env.sh
source "${ROOT}/deploy/aws/env.sh"

JOB_ID="${1:-}"
if [[ -z "${JOB_ID}" ]]; then
  echo "Usage: bash deploy/aws/04_watch_job.sh <jobId>" >&2
  exit 2
fi

echo "Watching job ${JOB_ID} (Ctrl-C to stop polling; job keeps running)"
while true; do
  OUT="$(aws batch describe-jobs --jobs "${JOB_ID}" --region "${AWS_REGION}" \
    --query 'jobs[0].{status:status,reason:statusReason,log:container.logStreamName}' \
    --output text)"
  echo "$(date -u +%H:%M:%S) ${OUT}"
  STATUS="$(echo "${OUT}" | awk '{print $1}')"
  case "${STATUS}" in
    SUCCEEDED)
      echo "OK: job succeeded"
      echo "Download: aws s3 cp s3://${OUTPUT_BUCKET}/practice/smoke_gpu/output.h5.gz ./output.h5.gz"
      exit 0
      ;;
    FAILED)
      echo "FAILED: ${OUT}" >&2
      echo "Open Batch → Jobs → this job → CloudWatch logs in the AWS Console." >&2
      exit 1
      ;;
    SUBMITTED|PENDING|RUNNABLE|STARTING|RUNNING|*)
      # Keep polling for known in-progress statuses and any unexpected value.
      sleep 20
      ;;
  esac
done
