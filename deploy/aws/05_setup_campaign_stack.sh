#!/usr/bin/env bash
# Create the campaign Batch stack: Spot GPU compute environment + queue + job
# definition sized for one GPU instance per Stage 3 run (g5.2xlarge / g4dn.2xlarge).
# Safe to re-run. Region pinned to us-east-1 via env.sh.
#
# Prerequisite: run 02_setup_practice_stack.sh first. This script REUSES the IAM
# roles, default VPC/subnet, and security group that 02 creates; it only adds the
# campaign-specific Spot compute environment, job queue, and job definition.
#
# From repo root:
#   bash deploy/aws/02_setup_practice_stack.sh   # once (roles/bucket/VPC/SG)
#   bash deploy/aws/05_setup_campaign_stack.sh
#
# Overridable env vars (see env.sh):
#   CAMPAIGN_MAX_VCPUS          max vCPUs for the Spot CE (default 96 = 12x g5.2xlarge)
#   CAMPAIGN_ONDEMAND_FALLBACK  set to 1 to also create an On-Demand fallback CE
#   COMPUTE_ENV_CAMPAIGN / COMPUTE_ENV_CAMPAIGN_OD / JOB_QUEUE_CAMPAIGN / JOB_DEFINITION_CAMPAIGN

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=env.sh
source "${ROOT}/deploy/aws/env.sh"

require_role() {
  local role="$1"
  if ! aws iam get-role --role-name "${role}" >/dev/null 2>&1; then
    echo "ERROR: IAM role '${role}' not found. Run 02_setup_practice_stack.sh first." >&2
    exit 1
  fi
}

wait_ce_valid() {
  local ce="$1"
  local status=""
  for _ in $(seq 1 60); do
    status="$(aws batch describe-compute-environments \
      --compute-environments "${ce}" \
      --query 'computeEnvironments[0].status' --output text --region "${AWS_REGION}")"
    echo "  status=${status}"
    if [[ "${status}" == "VALID" ]]; then
      return 0
    fi
    if [[ "${status}" == "INVALID" ]]; then
      aws batch describe-compute-environments \
        --compute-environments "${ce}" \
        --query 'computeEnvironments[0].statusReason' --output text --region "${AWS_REGION}"
      exit 1
    fi
    sleep 5
  done
  echo "ERROR: compute environment ${ce} not VALID yet (last=${status})" >&2
  exit 1
}

echo "==> Preflight: IAM roles from 02"
require_role AWSBatchServiceRole
require_role ecsInstanceRole
require_role AmazonEC2SpotFleetTaggingRole
require_role gutibm-batch-job-role
require_role gutibm-batch-execution-role
if ! aws iam get-instance-profile --instance-profile-name ecsInstanceRole >/dev/null 2>&1; then
  echo "ERROR: instance profile 'ecsInstanceRole' not found. Run 02_setup_practice_stack.sh first." >&2
  exit 1
fi

echo "==> Default VPC + subnet + security group (created by 02)"
VPC_ID="$(aws ec2 describe-vpcs \
  --filters Name=isDefault,Values=true \
  --query 'Vpcs[0].VpcId' --output text --region "${AWS_REGION}")"
if [[ -z "${VPC_ID}" || "${VPC_ID}" == "None" ]]; then
  echo "ERROR: no default VPC in ${AWS_REGION}. Run 02_setup_practice_stack.sh first." >&2
  exit 1
fi
SUBNET_ID="$(aws ec2 describe-subnets \
  --filters "Name=vpc-id,Values=${VPC_ID}" \
  --query 'Subnets[0].SubnetId' --output text --region "${AWS_REGION}")"
SG_ID="$(aws ec2 describe-security-groups \
  --filters "Name=group-name,Values=gutibm-batch-sg" "Name=vpc-id,Values=${VPC_ID}" \
  --query 'SecurityGroups[0].GroupId' --output text --region "${AWS_REGION}" 2>/dev/null || true)"
if [[ -z "${SG_ID}" || "${SG_ID}" == "None" ]]; then
  echo "ERROR: security group 'gutibm-batch-sg' not found. Run 02_setup_practice_stack.sh first." >&2
  exit 1
fi
echo "  VPC_ID=${VPC_ID}"
echo "  SUBNET_ID=${SUBNET_ID}"
echo "  SG_ID=${SG_ID}"

echo "==> Spot GPU compute environment ${COMPUTE_ENV_CAMPAIGN} (maxvCpus=${CAMPAIGN_MAX_VCPUS})"
CE_STATUS="$(aws batch describe-compute-environments \
  --compute-environments "${COMPUTE_ENV_CAMPAIGN}" \
  --query 'computeEnvironments[0].status' --output text --region "${AWS_REGION}" 2>/dev/null || true)"
if [[ -z "${CE_STATUS}" || "${CE_STATUS}" == "None" ]]; then
  SPOT_JSON="$(mktemp)"
  # type=SPOT is required for the SPOT_CAPACITY_OPTIMIZED allocation strategy.
  cat > "${SPOT_JSON}" <<EOF
{
  "type": "SPOT",
  "allocationStrategy": "SPOT_CAPACITY_OPTIMIZED",
  "minvCpus": 0,
  "maxvCpus": ${CAMPAIGN_MAX_VCPUS},
  "desiredvCpus": 0,
  "instanceTypes": ["g5.2xlarge", "g4dn.2xlarge"],
  "subnets": ["${SUBNET_ID}"],
  "securityGroupIds": ["${SG_ID}"],
  "instanceRole": "arn:aws:iam::${ACCOUNT}:instance-profile/ecsInstanceRole",
  "spotIamFleetRole": "arn:aws:iam::${ACCOUNT}:role/AmazonEC2SpotFleetTaggingRole",
  "ec2Configuration": [{"imageType": "ECS_AL2023_NVIDIA"}]
}
EOF
  aws batch create-compute-environment \
    --compute-environment-name "${COMPUTE_ENV_CAMPAIGN}" \
    --type MANAGED \
    --state ENABLED \
    --service-role "arn:aws:iam::${ACCOUNT}:role/AWSBatchServiceRole" \
    --compute-resources "file://${SPOT_JSON}" \
    --region "${AWS_REGION}" >/dev/null
  rm -f "${SPOT_JSON}"
else
  echo "  compute env already present (status=${CE_STATUS})"
fi
wait_ce_valid "${COMPUTE_ENV_CAMPAIGN}"

CE_ORDER="order=1,computeEnvironment=${COMPUTE_ENV_CAMPAIGN}"

if [[ "${CAMPAIGN_ONDEMAND_FALLBACK}" == "1" ]]; then
  echo "==> On-Demand fallback compute environment ${COMPUTE_ENV_CAMPAIGN_OD}"
  OD_STATUS="$(aws batch describe-compute-environments \
    --compute-environments "${COMPUTE_ENV_CAMPAIGN_OD}" \
    --query 'computeEnvironments[0].status' --output text --region "${AWS_REGION}" 2>/dev/null || true)"
  if [[ -z "${OD_STATUS}" || "${OD_STATUS}" == "None" ]]; then
    OD_JSON="$(mktemp)"
    cat > "${OD_JSON}" <<EOF
{
  "type": "EC2",
  "allocationStrategy": "BEST_FIT_PROGRESSIVE",
  "minvCpus": 0,
  "maxvCpus": ${CAMPAIGN_MAX_VCPUS},
  "desiredvCpus": 0,
  "instanceTypes": ["g5.2xlarge", "g4dn.2xlarge"],
  "subnets": ["${SUBNET_ID}"],
  "securityGroupIds": ["${SG_ID}"],
  "instanceRole": "arn:aws:iam::${ACCOUNT}:instance-profile/ecsInstanceRole",
  "ec2Configuration": [{"imageType": "ECS_AL2023_NVIDIA"}]
}
EOF
    aws batch create-compute-environment \
      --compute-environment-name "${COMPUTE_ENV_CAMPAIGN_OD}" \
      --type MANAGED \
      --state ENABLED \
      --service-role "arn:aws:iam::${ACCOUNT}:role/AWSBatchServiceRole" \
      --compute-resources "file://${OD_JSON}" \
      --region "${AWS_REGION}" >/dev/null
    rm -f "${OD_JSON}"
  else
    echo "  compute env already present (status=${OD_STATUS})"
  fi
  wait_ce_valid "${COMPUTE_ENV_CAMPAIGN_OD}"
  CE_ORDER="${CE_ORDER} order=2,computeEnvironment=${COMPUTE_ENV_CAMPAIGN_OD}"
fi

echo "==> Campaign job queue ${JOB_QUEUE_CAMPAIGN}"
JQ_STATE="$(aws batch describe-job-queues \
  --job-queues "${JOB_QUEUE_CAMPAIGN}" \
  --query 'jobQueues[0].state' --output text --region "${AWS_REGION}" 2>/dev/null || true)"
if [[ -z "${JQ_STATE}" || "${JQ_STATE}" == "None" ]]; then
  # shellcheck disable=SC2086
  aws batch create-job-queue \
    --job-queue-name "${JOB_QUEUE_CAMPAIGN}" \
    --state ENABLED \
    --priority 1 \
    --compute-environment-order ${CE_ORDER} \
    --region "${AWS_REGION}" >/dev/null
else
  echo "  job queue exists (state=${JQ_STATE})"
fi

echo "==> Campaign job definition ${JOB_DEFINITION_CAMPAIGN}"
JD_JSON="$(mktemp)"
# vcpus/memory sized to a full g5.2xlarge (8 vCPU, ~32 GB) so GPU=1 => one job
# per instance. Stage 3 needs ~8 GB host + A10G VRAM; 28000 MiB leaves headroom.
cat > "${JD_JSON}" <<EOF
{
  "image": "${IMAGE_URI}",
  "vcpus": 8,
  "memory": 28000,
  "resourceRequirements": [{"type": "GPU", "value": "1"}],
  "jobRoleArn": "arn:aws:iam::${ACCOUNT}:role/gutibm-batch-job-role",
  "executionRoleArn": "arn:aws:iam::${ACCOUNT}:role/gutibm-batch-execution-role",
  "privileged": true,
  "linuxParameters": {"sharedMemorySize": 2048}
}
EOF
RETRY_JSON="$(mktemp)"
# Retry Spot reclaims (Batch surfaces them as "Host EC2*" status reasons); exit
# immediately on application failures so bad configs do not burn two attempts.
cat > "${RETRY_JSON}" <<EOF
{
  "attempts": 2,
  "evaluateOnExit": [
    {"onStatusReason": "Host EC2*", "action": "RETRY"},
    {"onReason": "*", "action": "EXIT"}
  ]
}
EOF
aws batch register-job-definition \
  --job-definition-name "${JOB_DEFINITION_CAMPAIGN}" \
  --type container \
  --platform-capabilities EC2 \
  --container-properties "file://${JD_JSON}" \
  --retry-strategy "file://${RETRY_JSON}" \
  --region "${AWS_REGION}" >/dev/null
rm -f "${JD_JSON}" "${RETRY_JSON}"

echo "OK: campaign stack ready"
echo "  computeEnv=${COMPUTE_ENV_CAMPAIGN} (Spot, maxvCpus=${CAMPAIGN_MAX_VCPUS})"
if [[ "${CAMPAIGN_ONDEMAND_FALLBACK}" == "1" ]]; then
  echo "  fallbackCE=${COMPUTE_ENV_CAMPAIGN_OD} (On-Demand)"
fi
echo "  queue=${JOB_QUEUE_CAMPAIGN}"
echo "  jobDefinition=${JOB_DEFINITION_CAMPAIGN}"
echo "  image=${IMAGE_URI}"
echo "Next: python -m gut_ibm_tools.aws_batch_export \\"
echo "        experiments/diversity_campaign/stage3_campaign/batch_kd_sweep.json --dry-run"
