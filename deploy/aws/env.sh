#!/usr/bin/env bash
# Shared env for GutIBM AWS practice scripts. Source this; do not run it alone.
#   source deploy/aws/env.sh

export AWS_REGION="${AWS_REGION:-us-east-1}"
export AWS_DEFAULT_REGION="${AWS_REGION}"

if ! command -v aws >/dev/null 2>&1; then
  echo "ERROR: aws CLI not found. Install AWS CLI v2, or use CloudShell." >&2
  return 1 2>/dev/null || exit 1
fi

ACCOUNT="$(aws sts get-caller-identity --query Account --output text 2>/dev/null || true)"
if [[ -z "${ACCOUNT}" || "${ACCOUNT}" == "None" ]]; then
  echo "ERROR: cannot resolve AWS account. Run: aws configure   (or open CloudShell)" >&2
  return 1 2>/dev/null || exit 1
fi
export ACCOUNT

export REPO="${ACCOUNT}.dkr.ecr.${AWS_REGION}.amazonaws.com/gutibm"
export IMAGE_URI="${REPO}:cuda"
export INPUT_BUCKET="gutibm-inputs-${ACCOUNT}"
export OUTPUT_BUCKET="gutibm-outputs-${ACCOUNT}"
export JOB_QUEUE="${JOB_QUEUE:-gutibm-gpu-practice}"
export JOB_DEFINITION="${JOB_DEFINITION:-gutibm-cuda}"
export COMPUTE_ENV="${COMPUTE_ENV:-gutibm-gpu-practice-od}"

echo "AWS_REGION=${AWS_REGION}"
echo "ACCOUNT=${ACCOUNT}"
echo "IMAGE_URI=${IMAGE_URI}"
echo "INPUT_BUCKET=s3://${INPUT_BUCKET}"
echo "OUTPUT_BUCKET=s3://${OUTPUT_BUCKET}"
