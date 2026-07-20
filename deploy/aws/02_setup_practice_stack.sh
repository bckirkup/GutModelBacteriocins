#!/usr/bin/env bash
# Create S3 buckets + IAM roles + Batch practice stack (On-Demand g4dn.xlarge).
# Safe to re-run. Prefer CloudShell OR a laptop with AWS CLI.
#
# From repo root:
#   bash deploy/aws/02_setup_practice_stack.sh
#
# If you only have CloudShell and no git clone:
#   1) Upload the whole deploy/aws/ folder (Actions → Upload file), or
#   2) git clone your fork into CloudShell, then run this script.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
POLICY_DIR="${ROOT}/deploy/aws/policies"
# shellcheck source=env.sh
source "${ROOT}/deploy/aws/env.sh"

ensure_role() {
  local role="$1"
  local trust_file="$2"
  if aws iam get-role --role-name "${role}" >/dev/null 2>&1; then
    echo "  role exists: ${role}"
  else
    echo "  creating role: ${role}"
    aws iam create-role \
      --role-name "${role}" \
      --assume-role-policy-document "file://${trust_file}" >/dev/null
  fi
}

echo "==> S3 buckets"
aws s3api head-bucket --bucket "${INPUT_BUCKET}" 2>/dev/null \
  || aws s3 mb "s3://${INPUT_BUCKET}" --region "${AWS_REGION}"
aws s3api head-bucket --bucket "${OUTPUT_BUCKET}" 2>/dev/null \
  || aws s3 mb "s3://${OUTPUT_BUCKET}" --region "${AWS_REGION}"

echo "==> Default VPC + subnet + security group"
VPC_ID="$(aws ec2 describe-vpcs \
  --filters Name=isDefault,Values=true \
  --query 'Vpcs[0].VpcId' --output text --region "${AWS_REGION}")"
if [[ -z "${VPC_ID}" || "${VPC_ID}" == "None" ]]; then
  echo "ERROR: no default VPC in ${AWS_REGION}. Create a VPC or pick subnets manually." >&2
  exit 1
fi
SUBNET_ID="$(aws ec2 describe-subnets \
  --filters "Name=vpc-id,Values=${VPC_ID}" \
  --query 'Subnets[0].SubnetId' --output text --region "${AWS_REGION}")"
echo "  VPC_ID=${VPC_ID}"
echo "  SUBNET_ID=${SUBNET_ID}"

SG_ID="$(aws ec2 describe-security-groups \
  --filters "Name=group-name,Values=gutibm-batch-sg" "Name=vpc-id,Values=${VPC_ID}" \
  --query 'SecurityGroups[0].GroupId' --output text --region "${AWS_REGION}" 2>/dev/null || true)"
if [[ -z "${SG_ID}" || "${SG_ID}" == "None" ]]; then
  SG_ID="$(aws ec2 create-security-group \
    --group-name gutibm-batch-sg \
    --description "GutIBM Batch GPU egress" \
    --vpc-id "${VPC_ID}" \
    --region "${AWS_REGION}" \
    --query GroupId --output text)"
  echo "  created SG ${SG_ID}"
else
  echo "  SG exists: ${SG_ID}"
fi
# Default SG already allows all egress; extra authorize is optional/idempotent-fail-ok
aws ec2 authorize-security-group-egress \
  --group-id "${SG_ID}" \
  --ip-permissions 'IpProtocol=-1,IpRanges=[{CidrIp=0.0.0.0/0}]' \
  --region "${AWS_REGION}" >/dev/null 2>&1 || true

echo "==> IAM roles"
ensure_role AWSBatchServiceRole "${POLICY_DIR}/batch-service-trust.json"
aws iam attach-role-policy \
  --role-name AWSBatchServiceRole \
  --policy-arn arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole >/dev/null 2>&1 || true

ensure_role ecsInstanceRole "${POLICY_DIR}/ec2-trust.json"
aws iam attach-role-policy \
  --role-name ecsInstanceRole \
  --policy-arn arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role >/dev/null 2>&1 || true
if ! aws iam get-instance-profile --instance-profile-name ecsInstanceRole >/dev/null 2>&1; then
  aws iam create-instance-profile --instance-profile-name ecsInstanceRole >/dev/null
fi
aws iam add-role-to-instance-profile \
  --instance-profile-name ecsInstanceRole \
  --role-name ecsInstanceRole >/dev/null 2>&1 || true

ensure_role AmazonEC2SpotFleetTaggingRole "${POLICY_DIR}/spot-fleet-trust.json"
aws iam attach-role-policy \
  --role-name AmazonEC2SpotFleetTaggingRole \
  --policy-arn arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetTaggingRole >/dev/null 2>&1 || true

ensure_role gutibm-batch-job-role "${POLICY_DIR}/ecs-tasks-trust.json"
S3_POLICY="$(mktemp)"
cat > "${S3_POLICY}" <<EOF
{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Action": ["s3:GetObject", "s3:PutObject", "s3:ListBucket", "s3:DeleteObject"],
    "Resource": [
      "arn:aws:s3:::${INPUT_BUCKET}",
      "arn:aws:s3:::${INPUT_BUCKET}/*",
      "arn:aws:s3:::${OUTPUT_BUCKET}",
      "arn:aws:s3:::${OUTPUT_BUCKET}/*"
    ]
  }]
}
EOF
aws iam put-role-policy \
  --role-name gutibm-batch-job-role \
  --policy-name gutibm-s3-access \
  --policy-document "file://${S3_POLICY}"
rm -f "${S3_POLICY}"

ensure_role gutibm-batch-execution-role "${POLICY_DIR}/ecs-tasks-trust.json"
aws iam attach-role-policy \
  --role-name gutibm-batch-execution-role \
  --policy-arn arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy >/dev/null 2>&1 || true

echo "==> Wait for instance profile propagation"
sleep 15

echo "==> Batch compute environment ${COMPUTE_ENV}"
CE_STATUS="$(aws batch describe-compute-environments \
  --compute-environments "${COMPUTE_ENV}" \
  --query 'computeEnvironments[0].status' --output text --region "${AWS_REGION}" 2>/dev/null || true)"
if [[ -z "${CE_STATUS}" || "${CE_STATUS}" == "None" ]]; then
  CE_JSON="$(mktemp)"
  cat > "${CE_JSON}" <<EOF
{
  "type": "EC2",
  "allocationStrategy": "BEST_FIT_PROGRESSIVE",
  "minvCpus": 0,
  "maxvCpus": 4,
  "desiredvCpus": 0,
  "instanceTypes": ["g4dn.xlarge"],
  "subnets": ["${SUBNET_ID}"],
  "securityGroupIds": ["${SG_ID}"],
  "instanceRole": "arn:aws:iam::${ACCOUNT}:instance-profile/ecsInstanceRole",
  "ec2Configuration": [{"imageType": "ECS_AL2023_NVIDIA"}]
}
EOF
  aws batch create-compute-environment \
    --compute-environment-name "${COMPUTE_ENV}" \
    --type MANAGED \
    --state ENABLED \
    --service-role "arn:aws:iam::${ACCOUNT}:role/AWSBatchServiceRole" \
    --compute-resources "file://${CE_JSON}" \
    --region "${AWS_REGION}" >/dev/null
  rm -f "${CE_JSON}"
else
  echo "  compute env already present (status=${CE_STATUS})"
fi

echo "==> Wait until compute environment is VALID"
for _ in $(seq 1 60); do
  CE_STATUS="$(aws batch describe-compute-environments \
    --compute-environments "${COMPUTE_ENV}" \
    --query 'computeEnvironments[0].status' --output text --region "${AWS_REGION}")"
  echo "  status=${CE_STATUS}"
  if [[ "${CE_STATUS}" == "VALID" ]]; then
    break
  fi
  if [[ "${CE_STATUS}" == "INVALID" ]]; then
    aws batch describe-compute-environments \
      --compute-environments "${COMPUTE_ENV}" \
      --query 'computeEnvironments[0].statusReason' --output text --region "${AWS_REGION}"
    exit 1
  fi
  sleep 5
done
if [[ "${CE_STATUS}" != "VALID" ]]; then
  echo "ERROR: compute environment not VALID yet (last=${CE_STATUS})" >&2
  exit 1
fi

echo "==> Job queue ${JOB_QUEUE}"
JQ_STATE="$(aws batch describe-job-queues \
  --job-queues "${JOB_QUEUE}" \
  --query 'jobQueues[0].state' --output text --region "${AWS_REGION}" 2>/dev/null || true)"
if [[ -z "${JQ_STATE}" || "${JQ_STATE}" == "None" ]]; then
  aws batch create-job-queue \
    --job-queue-name "${JOB_QUEUE}" \
    --state ENABLED \
    --priority 1 \
    --compute-environment-order "order=1,computeEnvironment=${COMPUTE_ENV}" \
    --region "${AWS_REGION}" >/dev/null
else
  echo "  job queue exists (state=${JQ_STATE})"
fi

echo "==> Job definition ${JOB_DEFINITION}"
JD_JSON="$(mktemp)"
cat > "${JD_JSON}" <<EOF
{
  "image": "${IMAGE_URI}",
  "vcpus": 4,
  "memory": 14000,
  "resourceRequirements": [{"type": "GPU", "value": "1"}],
  "jobRoleArn": "arn:aws:iam::${ACCOUNT}:role/gutibm-batch-job-role",
  "executionRoleArn": "arn:aws:iam::${ACCOUNT}:role/gutibm-batch-execution-role",
  "privileged": true,
  "linuxParameters": {"sharedMemorySize": 2048}
}
EOF
aws batch register-job-definition \
  --job-definition-name "${JOB_DEFINITION}" \
  --type container \
  --platform-capabilities EC2 \
  --container-properties "file://${JD_JSON}" \
  --region "${AWS_REGION}" >/dev/null
rm -f "${JD_JSON}"

echo "OK: practice stack ready"
echo "  queue=${JOB_QUEUE}"
echo "  jobDefinition=${JOB_DEFINITION}"
echo "  image=${IMAGE_URI}"
echo "Next: bash deploy/aws/03_submit_smoke.sh"
