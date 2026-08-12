#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$ROOT/deploy/aws/07_widen_campaign_capacity.sh"
FIXTURE="$ROOT/tests/fixtures/campaign_compute_environment.json"

# Source only the script's reusable validation and payload helpers. Its AWS
# mutating main is guarded by BASH_SOURCE and is not executed here.
# shellcheck source=../deploy/aws/07_widen_campaign_capacity.sh
source "$SCRIPT"

fail() {
  echo "test_aws_capacity_script: FAILED: $*" >&2
  exit 1
}

CE_JSON="$(cat "$FIXTURE")"
CE_RESOURCES="$(jq -c '.computeEnvironments[0].computeResources' <<<"$CE_JSON")"
ALL_SUBNETS_JSON='[
  "subnet-0dfd15e585b9b7e84",
  "subnet-013e6340c89680943",
  "subnet-05c53fe19c7c50781",
  "subnet-0cdb5a8652bb8f0c6",
  "subnet-019f9335170c600d9",
  "subnet-0e6e3008d9fa8cf5d"
]'
INSTANCE_TYPES_JSON='[
  "g5.2xlarge", "g5.4xlarge", "g5.8xlarge", "g5.16xlarge",
  "g4dn.2xlarge", "g4dn.4xlarge", "g4dn.8xlarge", "g4dn.16xlarge"
]'
MIN_VCPUS="$(jq -r '.minvCpus' <<<"$CE_RESOURCES")"
MAX_VCPUS="$(jq -r '.maxvCpus' <<<"$CE_RESOURCES")"
SECURITY_GROUPS_JSON="$(jq -c '.securityGroupIds' <<<"$CE_RESOURCES")"
INSTANCE_ROLE="$(jq -r '.instanceRole' <<<"$CE_RESOURCES")"
SPOT_FLEET_ROLE="$(jq -r '.spotIamFleetRole' <<<"$CE_RESOURCES")"
EC2_CONFIGURATION_JSON="$(jq -c '[.ec2Configuration[] | {imageType: .imageType}]' <<<"$CE_RESOURCES")"

validate_copy_fields || fail "captured CE failed copy validation"
CREATE_PAYLOAD="$(
  jq -n \
    --arg name "gutibm-gpu-campaign-spot-v2" \
    --argjson resources "$(build_resources_payload)" \
    '{
       computeEnvironmentName: $name,
       type: "MANAGED",
       state: "ENABLED",
       computeResources: $resources
     }'
)"

jq -e '
  .computeResources.subnets | length == 6
  and (. | sort == [
    "subnet-013e6340c89680943",
    "subnet-019f9335170c600d9",
    "subnet-05c53fe19c7c50781",
    "subnet-0cdb5a8652bb8f0c6",
    "subnet-0dfd15e585b9b7e84",
    "subnet-0e6e3008d9fa8cf5d"
  ])
' <<<"$CREATE_PAYLOAD" >/dev/null ||
  fail "payload did not contain the six target subnets"
jq -e '.computeResources.instanceTypes | length == 8' <<<"$CREATE_PAYLOAD" >/dev/null ||
  fail "payload did not contain eight instance types"
jq -e '.serviceRole == null' <<<"$CREATE_PAYLOAD" >/dev/null ||
  fail "payload unexpectedly contained a legacy serviceRole"
jq -e '
  .computeResources.type == "SPOT"
  and .computeResources.allocationStrategy == "SPOT_CAPACITY_OPTIMIZED"
  and .computeResources.minvCpus == 0
  and .computeResources.maxvCpus == 96
  and .computeResources.securityGroupIds == ["sg-0842264e2d23e6b5d"]
  and .computeResources.instanceRole == "arn:aws:iam::994254241749:instance-profile/ecsInstanceRole"
  and .computeResources.spotIamFleetRole == "arn:aws:iam::994254241749:role/AmazonEC2SpotFleetTaggingRole"
  and .computeResources.ec2Configuration == [{"imageType": "ECS_AL2023_NVIDIA"}]
' <<<"$CREATE_PAYLOAD" >/dev/null ||
  fail "payload did not preserve the copied Spot/EC2 fields"

REUSE_JSON="$(
  jq --argjson resources "$(jq '.computeResources' <<<"$CREATE_PAYLOAD")" \
    '.computeEnvironments = [{
       serviceRole: "arn:aws:iam::994254241749:role/aws-service-role/batch.amazonaws.com/AWSServiceRoleForBatch",
       state: "ENABLED",
       computeResources: $resources
     }]' <<<"$CE_JSON"
)"
validate_replacement_json "$REUSE_JSON" ||
  fail "valid service-linked-role replacement failed reuse validation"

echo "AWS capacity validation and payload self-test passed."
