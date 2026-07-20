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

# S3 buckets are overridable so you can point at an existing bucket you own.
#   - Set BUCKET to use one shared bucket for both inputs and outputs
#     (differentiate with key prefixes, e.g. s3://my-bucket/gutibm/jobs).
#   - Or set INPUT_BUCKET / OUTPUT_BUCKET individually.
#   - Unset: derive the practice defaults (gutibm-inputs-/outputs-<account>),
#     which 02_setup_practice_stack.sh creates for you.
# Existing external buckets are never created (head-bucket guard in 02 skips
# `aws s3 mb`); only ensure the job role can read/write whatever these resolve to.
if [[ -n "${BUCKET:-}" ]]; then
  export INPUT_BUCKET="${INPUT_BUCKET:-${BUCKET}}"
  export OUTPUT_BUCKET="${OUTPUT_BUCKET:-${BUCKET}}"
else
  export INPUT_BUCKET="${INPUT_BUCKET:-gutibm-inputs-${ACCOUNT}}"
  export OUTPUT_BUCKET="${OUTPUT_BUCKET:-gutibm-outputs-${ACCOUNT}}"
fi

# Practice stack (02): On-Demand g4dn.xlarge.
export JOB_QUEUE="${JOB_QUEUE:-gutibm-gpu-practice}"
export JOB_DEFINITION="${JOB_DEFINITION:-gutibm-cuda}"
export COMPUTE_ENV="${COMPUTE_ENV:-gutibm-gpu-practice-od}"

# Campaign stack (05): Spot GPU (g5.2xlarge / g4dn.2xlarge), one GPU per run.
export COMPUTE_ENV_CAMPAIGN="${COMPUTE_ENV_CAMPAIGN:-gutibm-gpu-campaign-spot}"
export COMPUTE_ENV_CAMPAIGN_OD="${COMPUTE_ENV_CAMPAIGN_OD:-gutibm-gpu-campaign-od}"
export JOB_QUEUE_CAMPAIGN="${JOB_QUEUE_CAMPAIGN:-gutibm-gpu-campaign}"
export JOB_DEFINITION_CAMPAIGN="${JOB_DEFINITION_CAMPAIGN:-gutibm-cuda-campaign}"
# Max vCPUs for the campaign CE. Default 96 = 12 concurrent g5.2xlarge (8 vCPU
# each) so the 12-run Kd sweep can run fully in parallel, one GPU per run.
export CAMPAIGN_MAX_VCPUS="${CAMPAIGN_MAX_VCPUS:-96}"
# Set to 1 to also create an On-Demand fallback CE and order it after Spot.
export CAMPAIGN_ONDEMAND_FALLBACK="${CAMPAIGN_ONDEMAND_FALLBACK:-0}"

echo "AWS_REGION=${AWS_REGION}"
echo "ACCOUNT=${ACCOUNT}"
echo "IMAGE_URI=${IMAGE_URI}"
echo "INPUT_BUCKET=s3://${INPUT_BUCKET}"
echo "OUTPUT_BUCKET=s3://${OUTPUT_BUCKET}"
