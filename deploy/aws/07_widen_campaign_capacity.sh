#!/usr/bin/env bash
# Widen the existing campaign Spot compute environment across every available
# subnet in its VPC and across single-GPU instance sizes that fit Stage 3.
#
# ADMIN ONLY — this changes an existing Batch compute environment. It does not
# change the campaign queue or job definition. Run from any working directory:
#
#   bash deploy/aws/07_widen_campaign_capacity.sh
#
# The VPC is derived from the CE's existing subnet; no default-VPC assumption is
# made. The IAM principal must be allowed to describe EC2 resources and update
# the Batch compute environment.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -f "${SCRIPT_DIR}/env.sh" ]]; then
  echo "ERROR: env.sh is missing beside this script. Clone the repository and run it from the clone." >&2
  exit 1
fi
# shellcheck source=env.sh
source "${SCRIPT_DIR}/env.sh"

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq is required to build the CE update payload." >&2
  exit 1
fi

CE="${COMPUTE_ENV_CAMPAIGN}"
TMP_RESOURCES="$(mktemp)"
RESTORE_ENABLED=0
cleanup() {
  if [[ "${RESTORE_ENABLED}" == "1" ]]; then
    echo "WARNING: attempting to leave ${CE} ENABLED after an interrupted update" >&2
    aws batch update-compute-environment \
      --compute-environment "${CE}" \
      --state ENABLED \
      --region "${AWS_REGION}" >/dev/null 2>&1 || \
      echo "ERROR: could not re-enable ${CE}; inspect it with describe-compute-environments." >&2
  fi
  rm -f "${TMP_RESOURCES}"
}
trap cleanup EXIT

echo "==> Inspect campaign compute environment ${CE}"
CE_JSON="$(aws batch describe-compute-environments \
  --compute-environments "${CE}" \
  --region "${AWS_REGION}" \
  --output json)"
if [[ "$(jq '.computeEnvironments | length' <<<"${CE_JSON}")" != "1" ]]; then
  echo "ERROR: campaign compute environment '${CE}' was not found uniquely." >&2
  exit 1
fi

CE_STATE="$(jq -r '.computeEnvironments[0].state' <<<"${CE_JSON}")"
CE_STATUS="$(jq -r '.computeEnvironments[0].status' <<<"${CE_JSON}")"
EXISTING_SUBNETS_JSON="$(jq -c '.computeEnvironments[0].computeResources.subnets // []' \
  <<<"${CE_JSON}")"
if [[ "$(jq 'length' <<<"${EXISTING_SUBNETS_JSON}")" == "0" ]]; then
  echo "ERROR: ${CE} has no configured subnet; refusing to guess a VPC." >&2
  exit 1
fi

mapfile -t EXISTING_SUBNETS < <(jq -r '.[]' <<<"${EXISTING_SUBNETS_JSON}")
SOURCE_SUBNET="${EXISTING_SUBNETS[0]}"
VPC_ID="$(aws ec2 describe-subnets \
  --subnet-ids "${SOURCE_SUBNET}" \
  --region "${AWS_REGION}" \
  --query 'Subnets[0].VpcId' \
  --output text)"
if [[ -z "${VPC_ID}" || "${VPC_ID}" == "None" ]]; then
  echo "ERROR: could not derive a VPC from CE subnet ${SOURCE_SUBNET}." >&2
  exit 1
fi

EXISTING_VPCS_JSON="$(aws ec2 describe-subnets \
  --subnet-ids "${EXISTING_SUBNETS[@]}" \
  --region "${AWS_REGION}" \
  --query 'Subnets[].VpcId' \
  --output json)"
if ! jq -e --arg vpc "${VPC_ID}" \
  'length > 0 and all(.[]; . == $vpc)' <<<"${EXISTING_VPCS_JSON}" >/dev/null; then
  echo "ERROR: existing CE subnets do not all belong to VPC ${VPC_ID}." >&2
  exit 1
fi

ALL_SUBNETS_JSON="$(aws ec2 describe-subnets \
  --filters "Name=vpc-id,Values=${VPC_ID}" \
  --region "${AWS_REGION}" \
  --query 'Subnets[?State==`available`].SubnetId' \
  --output json)"
if [[ "$(jq 'length' <<<"${ALL_SUBNETS_JSON}")" == "0" ]]; then
  echo "ERROR: VPC ${VPC_ID} has no available subnets." >&2
  exit 1
fi

# Every listed type has one GPU, at least 8 vCPUs, and at least 32 GiB host
# memory. Multi-GPU sizes are intentionally excluded: this job requests one
# GPU, so using them would pay for idle accelerators without improving fit.
INSTANCE_TYPES_JSON='[
  "g5.2xlarge", "g5.4xlarge", "g5.8xlarge", "g5.16xlarge",
  "g4dn.2xlarge", "g4dn.4xlarge", "g4dn.8xlarge", "g4dn.16xlarge"
]'

echo "  state=${CE_STATE} status=${CE_STATUS}"
echo "  source_subnet=${SOURCE_SUBNET}"
echo "  vpc=${VPC_ID}"
echo "  current_subnets=$(jq -c '.' <<<"${EXISTING_SUBNETS_JSON}")"
echo "  target_subnets=$(jq -c '.' <<<"${ALL_SUBNETS_JSON}")"
echo "  current_instance_types=$(jq -c '.computeEnvironments[0].computeResources.instanceTypes' <<<"${CE_JSON}")"
echo "  target_instance_types=$(jq -c '.' <<<"${INSTANCE_TYPES_JSON}")"

# UpdateComputeEnvironment accepts a ComputeResourceUpdate object, not the
# describe response. Omitted fields remain unchanged, so use an explicit
# allow-list rather than round-tripping output-only or create-only fields such
# as type, spotIamFleetRole, launchTemplate, and image settings. Preserve the
# capacity bounds, but omit desiredvCpus: Batch owns that scheduler-controlled
# value and resending it can fight scaling.
if ! jq -e \
  '.computeEnvironments[0].computeResources
   | (.minvCpus | type == "number")
   and (.maxvCpus | type == "number")' <<<"${CE_JSON}" >/dev/null; then
  echo "ERROR: CE is missing numeric minvCpus or maxvCpus." >&2
  exit 1
fi
jq --argjson subnets "${ALL_SUBNETS_JSON}" \
  --argjson instance_types "${INSTANCE_TYPES_JSON}" \
  '.computeEnvironments[0].computeResources
   | {
       minvCpus: .minvCpus,
       maxvCpus: .maxvCpus,
       subnets: $subnets,
       instanceTypes: $instance_types
     }' <<<"${CE_JSON}" >"${TMP_RESOURCES}"

wait_for_state() {
  local expected="$1"
  local state=""
  local status=""
  for _ in $(seq 1 60); do
    state="$(aws batch describe-compute-environments \
      --compute-environments "${CE}" \
      --region "${AWS_REGION}" \
      --query 'computeEnvironments[0].state' \
      --output text)"
    status="$(aws batch describe-compute-environments \
      --compute-environments "${CE}" \
      --region "${AWS_REGION}" \
      --query 'computeEnvironments[0].status' \
      --output text)"
    echo "  state=${state} status=${status}"
    if [[ "${state}" == "${expected}" && "${status}" != "INVALID" ]]; then
      return 0
    fi
    if [[ "${status}" == "INVALID" ]]; then
      aws batch describe-compute-environments \
        --compute-environments "${CE}" \
        --region "${AWS_REGION}" \
        --query 'computeEnvironments[0].statusReason' \
        --output text >&2
      return 1
    fi
    sleep 5
  done
  echo "ERROR: ${CE} did not reach state ${expected}." >&2
  return 1
}

update_resources() {
  local output=""
  if output="$(aws batch update-compute-environment \
    --compute-environment "${CE}" \
    --compute-resources "file://${TMP_RESOURCES}" \
    --region "${AWS_REGION}" 2>&1)"; then
    echo "${output}"
    return 0
  fi

  echo "${output}" >&2
  local api_error_code
  api_error_code="$(sed -n 's/.*An error occurred (\([^)]*\)).*/\1/p' <<<"${output}")"
  if [[ "${api_error_code}" != "ClientException" ]]; then
    echo "ERROR: CE update failed with API error code '${api_error_code:-unknown}'." >&2
    return 1
  fi

  RESTORE_ENABLED=1
  if [[ "${CE_STATE}" != "DISABLED" ]]; then
    echo "==> Retrying the infrastructure update with the CE disabled"
    aws batch update-compute-environment \
      --compute-environment "${CE}" \
      --state DISABLED \
      --region "${AWS_REGION}" >/dev/null
  fi
  wait_for_state DISABLED
  aws batch update-compute-environment \
    --compute-environment "${CE}" \
    --compute-resources "file://${TMP_RESOURCES}" \
    --region "${AWS_REGION}" >/dev/null
  aws batch update-compute-environment \
    --compute-environment "${CE}" \
    --state ENABLED \
    --region "${AWS_REGION}" >/dev/null
  RESTORE_ENABLED=0
}

echo "==> Update ${CE} without touching queue or job definitions"
update_resources
CURRENT_STATE="$(aws batch describe-compute-environments \
  --compute-environments "${CE}" \
  --region "${AWS_REGION}" \
  --query 'computeEnvironments[0].state' \
  --output text)"
if [[ "${CURRENT_STATE}" != "ENABLED" ]]; then
  echo "==> Re-enabling ${CE}"
  aws batch update-compute-environment \
    --compute-environment "${CE}" \
    --state ENABLED \
    --region "${AWS_REGION}" >/dev/null
fi
wait_for_state ENABLED
RESTORE_ENABLED=0

AFTER_JSON="$(aws batch describe-compute-environments \
  --compute-environments "${CE}" \
  --region "${AWS_REGION}" \
  --output json)"
echo "==> Campaign compute environment updated"
echo "  state=$(jq -r '.computeEnvironments[0].state' <<<"${AFTER_JSON}")"
echo "  status=$(jq -r '.computeEnvironments[0].status' <<<"${AFTER_JSON}")"
echo "  subnets=$(jq -c '.computeEnvironments[0].computeResources.subnets' <<<"${AFTER_JSON}")"
echo "  instance_types=$(jq -c '.computeEnvironments[0].computeResources.instanceTypes' <<<"${AFTER_JSON}")"
echo "OK: ${CE} is ENABLED; queue and job definitions were not modified"
