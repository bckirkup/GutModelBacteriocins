#!/usr/bin/env bash
# Replace the legacy campaign Spot compute environment with a service-linked
# role CE spanning every available subnet and the approved single-GPU sizes.
#
# ADMIN ONLY — this changes a Batch compute environment and campaign queue.
# It never changes a job definition. Run from any working directory in a clone:
#
#   bash deploy/aws/07_widen_campaign_capacity.sh
#
# A legacy-service-role CE cannot update its subnets or instance types. This
# script creates a new CE, repoints the queue, and leaves the old CE enabled.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  if [[ ! -f "${SCRIPT_DIR}/env.sh" ]]; then
    echo "ERROR: env.sh is missing beside this script. Clone the repository and run it from the clone." >&2
    exit 1
  fi
  # shellcheck source=env.sh
  source "${SCRIPT_DIR}/env.sh"
fi

cleanup() {
  local rc=$?
  trap - EXIT
  if [[ "${RESTORE_OLD}" == "1" ]]; then
    echo "WARNING: restoring legacy CE ${CE} to ENABLED before exit" >&2
    restore_old_ce || \
      echo "ERROR: could not verify ${CE} is ENABLED; inspect it manually." >&2
  fi
  rm -f "${TMP_RESOURCES}"
  exit "${rc}"
}

describe_ce() {
  aws batch describe-compute-environments \
    --compute-environments "$1" \
    --region "${AWS_REGION}" \
    --output json
}

describe_queue() {
  aws batch describe-job-queues \
    --job-queues "${QUEUE}" \
    --region "${AWS_REGION}" \
    --output json
}

wait_for_ce() {
  local name="$1"
  local data=""
  for _ in $(seq 1 60); do
    data="$(describe_ce "${name}")"
    if jq -e \
      '.computeEnvironments | length == 1 and
       .[0].status == "VALID" and .[0].state == "ENABLED"' \
      <<<"${data}" >/dev/null; then
      return 0
    fi
    if jq -e '.computeEnvironments[0].status == "INVALID"' \
      <<<"${data}" >/dev/null; then
      jq -r '.computeEnvironments[0].statusReason' <<<"${data}" >&2
      return 1
    fi
    echo "  ${name}: state=$(jq -r '.computeEnvironments[0].state' <<<"${data}") status=$(jq -r '.computeEnvironments[0].status' <<<"${data}")"
    sleep 5
  done
  echo "ERROR: ${name} did not become VALID and ENABLED." >&2
  return 1
}

restore_old_ce() {
  local data=""
  local state=""
  local status=""
  for _ in $(seq 1 60); do
    data="$(describe_ce "${CE}")"
    state="$(jq -r '.computeEnvironments[0].state' <<<"${data}")"
    status="$(jq -r '.computeEnvironments[0].status' <<<"${data}")"
    if [[ "${state}" == "ENABLED" && "${status}" == "VALID" ]]; then
      return 0
    fi
    if [[ "${status}" == "INVALID" ]]; then
      jq -r '.computeEnvironments[0].statusReason' <<<"${data}" >&2
      return 1
    fi
    if [[ "${state}" != "ENABLED" && "${status}" != "UPDATING" ]]; then
      aws batch update-compute-environment \
        --compute-environment "${CE}" \
        --state ENABLED \
        --region "${AWS_REGION}" >/dev/null 2>&1 || true
    fi
    echo "  ${CE}: state=${state} status=${status}; waiting to verify ENABLED"
    sleep 5
  done
  return 1
}

validate_copy_fields() {
  if [[ ! "${MIN_VCPUS}" =~ ^[0-9]+$ || ! "${MAX_VCPUS}" =~ ^[0-9]+$ ||
        "$(jq 'length' <<<"${SECURITY_GROUPS_JSON}")" == "0" ||
        -z "${INSTANCE_ROLE}" || "${INSTANCE_ROLE}" != *:instance-profile/* ||
        -z "${SPOT_FLEET_ROLE}" ||
        "$(jq 'length' <<<"${EC2_CONFIGURATION_JSON}")" == "0" ]]; then
    echo "ERROR: legacy CE lacks required Spot/EC2 resource settings to copy." >&2
    return 1
  fi
  if ! jq -e 'all(.[]; (.imageType | type) == "string")' \
    <<<"${EC2_CONFIGURATION_JSON}" >/dev/null; then
    echo "ERROR: legacy CE has an invalid ec2Configuration image type." >&2
    return 1
  fi
}

build_resources_payload() {
  jq -n \
    --argjson subnets "${ALL_SUBNETS_JSON}" \
    --argjson instance_types "${INSTANCE_TYPES_JSON}" \
    --argjson security_groups "${SECURITY_GROUPS_JSON}" \
    --argjson ec2_configuration "${EC2_CONFIGURATION_JSON}" \
    --arg min_vcpus "${MIN_VCPUS}" \
    --arg max_vcpus "${MAX_VCPUS}" \
    --arg instance_role "${INSTANCE_ROLE}" \
    --arg spot_fleet_role "${SPOT_FLEET_ROLE}" \
    '{
       type: "SPOT",
       allocationStrategy: "SPOT_CAPACITY_OPTIMIZED",
       minvCpus: ($min_vcpus | tonumber),
       maxvCpus: ($max_vcpus | tonumber),
       desiredvCpus: 0,
       instanceTypes: $instance_types,
       subnets: $subnets,
       securityGroupIds: $security_groups,
       instanceRole: $instance_role,
       spotIamFleetRole: $spot_fleet_role,
       ec2Configuration: $ec2_configuration
     }'
}

validate_replacement_json() {
  local replacement_json="$1"
  jq -e \
    --argjson subnets "${ALL_SUBNETS_JSON}" \
    --argjson instance_types "${INSTANCE_TYPES_JSON}" \
    '.computeEnvironments[0]
     | (.serviceRole | contains("AWSServiceRoleForBatch"))
       and (.state != "INVALID")
       and (.computeResources.subnets | sort == ($subnets | sort))
       and (.computeResources.instanceTypes | sort == ($instance_types | sort))' \
    <<<"${replacement_json}" >/dev/null
}

main() {
  local arg
  DRY_RUN=0
  for arg in "$@"; do
    case "${arg}" in
      --dry-run)
        DRY_RUN=1
        ;;
      --help|-h)
        echo "Usage: bash deploy/aws/07_widen_campaign_capacity.sh [--dry-run]"
        return 0
        ;;
      *)
        echo "ERROR: unknown argument '${arg}'. Use --help for usage." >&2
        return 2
        ;;
    esac
  done
  if ! command -v jq >/dev/null 2>&1; then
    echo "ERROR: jq is required to build and validate the CE payloads." >&2
    return 1
  fi

  CE="${COMPUTE_ENV_CAMPAIGN}"
  NEW_CE="${COMPUTE_ENV_CAMPAIGN}-v2"
  QUEUE="${JOB_QUEUE_CAMPAIGN}"
  TMP_RESOURCES="$(mktemp)"
  RESTORE_OLD=1
  trap cleanup EXIT

echo "==> Inspect legacy campaign compute environment ${CE}"
CE_JSON="$(describe_ce "${CE}")"
if [[ "$(jq '.computeEnvironments | length' <<<"${CE_JSON}")" != "1" ]]; then
  echo "ERROR: campaign compute environment '${CE}' was not found uniquely." >&2
  exit 1
fi

CE_RESOURCES="$(jq -c '.computeEnvironments[0].computeResources' <<<"${CE_JSON}")"
EXISTING_SUBNETS_JSON="$(jq -c '.subnets // []' <<<"${CE_RESOURCES}")"
if [[ "$(jq 'length' <<<"${EXISTING_SUBNETS_JSON}")" == "0" ]]; then
  echo "ERROR: ${CE} has no configured subnet; refusing to guess a VPC." >&2
  exit 1
fi
mapfile -t EXISTING_SUBNETS < <(jq -r '.[]' <<<"${EXISTING_SUBNETS_JSON}")
VPC_ID="$(aws ec2 describe-subnets \
  --subnet-ids "${EXISTING_SUBNETS[0]}" \
  --region "${AWS_REGION}" \
  --query 'Subnets[0].VpcId' \
  --output text)"
if [[ -z "${VPC_ID}" || "${VPC_ID}" == "None" ]]; then
  echo "ERROR: could not derive a VPC from CE subnet ${EXISTING_SUBNETS[0]}." >&2
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

INSTANCE_TYPES_JSON='[
  "g5.2xlarge", "g5.4xlarge", "g5.8xlarge", "g5.16xlarge",
  "g4dn.2xlarge", "g4dn.4xlarge", "g4dn.8xlarge", "g4dn.16xlarge"
]'

MIN_VCPUS="$(jq -r '.minvCpus' <<<"${CE_RESOURCES}")"
MAX_VCPUS="$(jq -r '.maxvCpus' <<<"${CE_RESOURCES}")"
SECURITY_GROUPS_JSON="$(jq -c '.securityGroupIds // []' <<<"${CE_RESOURCES}")"
INSTANCE_ROLE="$(jq -r '.instanceRole // empty' <<<"${CE_RESOURCES}")"
SPOT_FLEET_ROLE="$(jq -r '.spotIamFleetRole // empty' <<<"${CE_RESOURCES}")"
EC2_CONFIGURATION_JSON="$(jq -c '[.ec2Configuration[]? | select(.imageType != null) | {imageType: .imageType}]' <<<"${CE_RESOURCES}")"
# The AWS ComputeResource API requires spotIamFleetRole for SPOT resources.
# The service-linked role replaces only the legacy Batch serviceRole.
validate_copy_fields

echo "  old_state=$(jq -r '.computeEnvironments[0].state' <<<"${CE_JSON}")"
echo "  old_status=$(jq -r '.computeEnvironments[0].status' <<<"${CE_JSON}")"
echo "  vpc=${VPC_ID}"
echo "  target_subnets=$(jq -c '.' <<<"${ALL_SUBNETS_JSON}")"
echo "  target_instance_types=$(jq -c '.' <<<"${INSTANCE_TYPES_JSON}")"

echo "==> Ensure the old CE is ENABLED before changing the queue"
if ((DRY_RUN)); then
  echo "  dry-run: would verify/re-enable ${CE}; no update call will be made"
else
  restore_old_ce
fi

build_resources_payload >"${TMP_RESOURCES}"
if ((DRY_RUN)); then
  echo "==> Dry-run create-compute-environment payload"
  jq -n \
    --arg name "${NEW_CE}" \
    --argjson resources "$(cat "${TMP_RESOURCES}")" \
    '{
       computeEnvironmentName: $name,
       type: "MANAGED",
       state: "ENABLED",
       computeResources: $resources
     }'
fi

NEW_JSON="$(describe_ce "${NEW_CE}" 2>/dev/null || true)"
if [[ -z "${NEW_JSON}" || "$(jq '.computeEnvironments | length' <<<"${NEW_JSON}")" == "0" ]]; then
  echo "==> Create service-linked-role campaign CE ${NEW_CE}"
  echo "  no --service-role is passed; Batch uses AWSServiceRoleForBatch"
  if ((DRY_RUN)); then
    echo "  dry-run: would run aws batch create-compute-environment"
  else
    aws batch create-compute-environment \
      --compute-environment-name "${NEW_CE}" \
      --type MANAGED \
      --state ENABLED \
      --compute-resources "file://${TMP_RESOURCES}" \
      --region "${AWS_REGION}" >/dev/null
  fi
else
  if ! validate_replacement_json "${NEW_JSON}"; then
    echo "ERROR: existing ${NEW_CE} does not match the required service-linked-role capacity shape." >&2
    exit 1
  fi
  echo "  ${NEW_CE} already exists; reusing it"
fi
if ((DRY_RUN)); then
  echo "  dry-run: would wait for ${NEW_CE} to become VALID/ENABLED"
else
  wait_for_ce "${NEW_CE}"
fi

QUEUE_JSON="$(describe_queue)"
if [[ "$(jq '.jobQueues | length' <<<"${QUEUE_JSON}")" != "1" ]]; then
  echo "ERROR: campaign queue '${QUEUE}' was not found uniquely." >&2
  exit 1
fi
QUEUE_PRIORITY="$(jq -r '.jobQueues[0].priority' <<<"${QUEUE_JSON}")"
if [[ ! "${QUEUE_PRIORITY}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: campaign queue has no numeric priority." >&2
  exit 1
fi
OLD_ORDER="$(jq -r --arg old_ce "${CE}" \
  '.jobQueues[0].computeEnvironmentOrder[]
   | select(
       .computeEnvironment == $old_ce or
       ((.computeEnvironment | split("/") | .[-1]) == $old_ce)
     ) | .order' <<<"${QUEUE_JSON}" | head -n 1)"
if [[ -z "${OLD_ORDER}" ]]; then
  OLD_ORDER=1
fi
QUEUE_ORDER_ARGS=("order=${OLD_ORDER},computeEnvironment=${NEW_CE}")
while IFS='|' read -r order environment; do
  [[ -z "${environment}" ]] && continue
  QUEUE_ORDER_ARGS+=("order=${order},computeEnvironment=${environment}")
done < <(jq -r --arg old_ce "${CE}" \
  '.jobQueues[0].computeEnvironmentOrder[]
   | select(
       (.computeEnvironment != $old_ce) and
       ((.computeEnvironment | split("/") | .[-1]) != $old_ce)
     )
   | "\(.order)|\(.computeEnvironment)"' <<<"${QUEUE_JSON}")

# Drop the old CE rather than leaving it at its old order: its legacy service
# role makes its single-AZ capacity permanently unfixable. Preserve any
# unrelated fallback CEs at their existing orders. Jobs are associated with
# the queue, so queued jobs remain eligible after this repoint.
echo "==> Repoint queue ${QUEUE} to ${NEW_CE}; old CE is intentionally dropped"
if ((DRY_RUN)); then
  printf '  dry-run: would run aws batch update-job-queue --compute-environment-order %q\n' \
    "${QUEUE_ORDER_ARGS[@]}"
else
  aws batch update-job-queue \
    --job-queue "${QUEUE}" \
    --state ENABLED \
    --priority "${QUEUE_PRIORITY}" \
    --compute-environment-order "${QUEUE_ORDER_ARGS[@]}" \
    --region "${AWS_REGION}" >/dev/null
fi

for _ in $(seq 1 60); do
  ((DRY_RUN)) && break
  QUEUE_JSON="$(describe_queue)"
  if jq -e --arg new_ce "${NEW_CE}" --arg old_ce "${CE}" \
    '.jobQueues | length == 1 and .[0].status == "VALID" and
     .[0].state == "ENABLED" and
     (.[0].computeEnvironmentOrder | any(
       .computeEnvironment == $new_ce or
       ((.computeEnvironment | split("/") | .[-1]) == $new_ce)
     )) and
     (.[0].computeEnvironmentOrder | all(
       (.computeEnvironment != $old_ce) and
       ((.computeEnvironment | split("/") | .[-1]) != $old_ce)
     ))' \
    <<<"${QUEUE_JSON}" >/dev/null; then
    break
  fi
  echo "  queue ${QUEUE}: waiting for VALID/ENABLED order"
  sleep 5
done
if (( ! DRY_RUN )) && ! jq -e --arg new_ce "${NEW_CE}" --arg old_ce "${CE}" \
  '.jobQueues | length == 1 and .[0].status == "VALID" and
   .[0].state == "ENABLED" and
   (.[0].computeEnvironmentOrder | any(
     .computeEnvironment == $new_ce or
     ((.computeEnvironment | split("/") | .[-1]) == $new_ce)
   )) and
   (.[0].computeEnvironmentOrder | all(
     (.computeEnvironment != $old_ce) and
     ((.computeEnvironment | split("/") | .[-1]) != $old_ce)
   ))' \
  <<<"${QUEUE_JSON}" >/dev/null; then
  echo "ERROR: queue ${QUEUE} did not become VALID/ENABLED on ${NEW_CE}." >&2
  exit 1
fi

if ((DRY_RUN)); then
  RESTORE_OLD=0
  echo "OK: dry-run completed without AWS mutations"
  return 0
fi
restore_old_ce
RESTORE_OLD=0
FINAL_CE_JSON="$(describe_ce "${NEW_CE}")"
FINAL_QUEUE_JSON="$(describe_queue)"
echo "==> Campaign capacity replacement complete"
echo "  new_ce=${NEW_CE}"
echo "  new_ce_state=$(jq -r '.computeEnvironments[0].state' <<<"${FINAL_CE_JSON}")"
echo "  new_ce_status=$(jq -r '.computeEnvironments[0].status' <<<"${FINAL_CE_JSON}")"
echo "  subnets=$(jq -c '.computeEnvironments[0].computeResources.subnets' <<<"${FINAL_CE_JSON}")"
echo "  instance_types=$(jq -c '.computeEnvironments[0].computeResources.instanceTypes' <<<"${FINAL_CE_JSON}")"
echo "  queue_state=$(jq -r '.jobQueues[0].state' <<<"${FINAL_QUEUE_JSON}")"
echo "  queue_status=$(jq -r '.jobQueues[0].status' <<<"${FINAL_QUEUE_JSON}")"
echo "  queue_compute_environments=$(jq -c '.jobQueues[0].computeEnvironmentOrder' <<<"${FINAL_QUEUE_JSON}")"
echo "OK: new CE and queue are VALID/ENABLED; old CE remains ENABLED; job definitions were not modified"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  main "$@"
fi
