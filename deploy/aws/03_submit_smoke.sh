#!/usr/bin/env bash
# Upload experiments/smoke_gpu.json and submit one Batch job.
# From repo root (laptop or CloudShell with a clone):
#
#   bash deploy/aws/03_submit_smoke.sh
#
# Prints jobId. Then poll with:
#   bash deploy/aws/04_watch_job.sh <jobId>

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=env.sh
source "${ROOT}/deploy/aws/env.sh"

SMOKE_JSON="${ROOT}/experiments/smoke_gpu.json"
if [[ ! -f "${SMOKE_JSON}" ]]; then
  echo "ERROR: missing ${SMOKE_JSON}" >&2
  echo "Clone the repo (or upload that file) before running this script." >&2
  exit 1
fi

INPUT_URI="s3://${INPUT_BUCKET}/practice/smoke_gpu/input.json"
OUTPUT_URI="s3://${OUTPUT_BUCKET}/practice/smoke_gpu/output.h5.gz"

echo "==> Upload smoke config"
aws s3 cp "${SMOKE_JSON}" "${INPUT_URI}"

OVERRIDES="$(mktemp)"
cat > "${OVERRIDES}" <<EOF
{
  "environment": [
    {"name": "INPUT_S3_URI", "value": "${INPUT_URI}"},
    {"name": "OUTPUT_S3_URI", "value": "${OUTPUT_URI}"},
    {"name": "MPI_RANKS", "value": "1"}
  ]
}
EOF

echo "==> Submit job"
JOB_ID="$(aws batch submit-job \
  --job-name gutibm-smoke-gpu \
  --job-queue "${JOB_QUEUE}" \
  --job-definition "${JOB_DEFINITION}" \
  --container-overrides "file://${OVERRIDES}" \
  --region "${AWS_REGION}" \
  --query jobId --output text)"
rm -f "${OVERRIDES}"

echo "OK: submitted jobId=${JOB_ID}"
echo "Watch:  bash deploy/aws/04_watch_job.sh ${JOB_ID}"
echo "Output: ${OUTPUT_URI}"
