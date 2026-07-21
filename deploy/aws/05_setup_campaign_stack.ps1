# Create the campaign Batch stack: Spot GPU compute environment + queue + job
# definition sized for one GPU instance per Stage 3 run (g5.2xlarge / g4dn.2xlarge).
# PowerShell port of 05_setup_campaign_stack.sh. Safe to re-run.
#
# Prerequisite: run 02_setup_practice_stack.ps1 first (it creates the IAM roles,
# VPC/subnet, and security group this script reuses).
#
# From repo root:
#   .\deploy\aws\02_setup_practice_stack.ps1   # once
#   .\deploy\aws\05_setup_campaign_stack.ps1
#
# Overridable via $env: CAMPAIGN_MAX_VCPUS, CAMPAIGN_ONDEMAND_FALLBACK (=1),
# COMPUTE_ENV_CAMPAIGN / COMPUTE_ENV_CAMPAIGN_OD / JOB_QUEUE_CAMPAIGN / JOB_DEFINITION_CAMPAIGN.

$ErrorActionPreference = "Stop"

. (Join-Path $PSScriptRoot "env.ps1")

function Test-AwsText {
  param([string]$Value)
  if ([string]::IsNullOrWhiteSpace($Value) -or $Value.Trim() -eq "None") { return "" }
  return $Value.Trim()
}

function Require-Role {
  param([string]$Role)
  aws iam get-role --role-name $Role *> $null
  if ($LASTEXITCODE -ne 0) { throw "IAM role '$Role' not found. Run 02_setup_practice_stack.ps1 first." }
}

function Wait-CeValid {
  param([string]$Ce)
  $status = ""
  for ($i = 0; $i -lt 60; $i++) {
    $status = (aws batch describe-compute-environments --compute-environments $Ce --query 'computeEnvironments[0].status' --output text --region $env:AWS_REGION).Trim()
    Write-Host "  status=$status"
    if ($status -eq "VALID") { return }
    if ($status -eq "INVALID") {
      aws batch describe-compute-environments --compute-environments $Ce --query 'computeEnvironments[0].statusReason' --output text --region $env:AWS_REGION
      throw "compute environment $Ce INVALID"
    }
    Start-Sleep -Seconds 5
  }
  throw "compute environment $Ce not VALID yet (last=$status)"
}

Write-Host "==> Preflight: IAM roles from 02"
Require-Role "AWSBatchServiceRole"
Require-Role "ecsInstanceRole"
Require-Role "AmazonEC2SpotFleetTaggingRole"
Require-Role "gutibm-batch-job-role"
Require-Role "gutibm-batch-execution-role"
aws iam get-instance-profile --instance-profile-name ecsInstanceRole *> $null
if ($LASTEXITCODE -ne 0) { throw "instance profile 'ecsInstanceRole' not found. Run 02_setup_practice_stack.ps1 first." }

Write-Host "==> Default VPC + subnet + security group (created by 02)"
$VPC_ID = Test-AwsText (aws ec2 describe-vpcs --filters Name=isDefault,Values=true --query 'Vpcs[0].VpcId' --output text --region $env:AWS_REGION)
if (-not $VPC_ID) { throw "no default VPC in $($env:AWS_REGION). Run 02_setup_practice_stack.ps1 first." }
$SUBNET_ID = Test-AwsText (aws ec2 describe-subnets --filters "Name=vpc-id,Values=$VPC_ID" --query 'Subnets[0].SubnetId' --output text --region $env:AWS_REGION)
$SG_ID = Test-AwsText (aws ec2 describe-security-groups --filters "Name=group-name,Values=gutibm-batch-sg" "Name=vpc-id,Values=$VPC_ID" --query 'SecurityGroups[0].GroupId' --output text --region $env:AWS_REGION 2>$null)
if (-not $SG_ID) { throw "security group 'gutibm-batch-sg' not found. Run 02_setup_practice_stack.ps1 first." }
Write-Host "  VPC_ID=$VPC_ID"
Write-Host "  SUBNET_ID=$SUBNET_ID"
Write-Host "  SG_ID=$SG_ID"

Write-Host "==> Spot GPU compute environment $($env:COMPUTE_ENV_CAMPAIGN) (maxvCpus=$($env:CAMPAIGN_MAX_VCPUS))"
$CE_STATUS = Test-AwsText (aws batch describe-compute-environments --compute-environments $env:COMPUTE_ENV_CAMPAIGN --query 'computeEnvironments[0].status' --output text --region $env:AWS_REGION 2>$null)
if (-not $CE_STATUS) {
  $SpotJson = New-TemporaryFile
  # type=SPOT is required for the SPOT_CAPACITY_OPTIMIZED allocation strategy.
  @"
{
  "type": "SPOT",
  "allocationStrategy": "SPOT_CAPACITY_OPTIMIZED",
  "minvCpus": 0,
  "maxvCpus": $($env:CAMPAIGN_MAX_VCPUS),
  "desiredvCpus": 0,
  "instanceTypes": ["g5.2xlarge", "g4dn.2xlarge"],
  "subnets": ["$SUBNET_ID"],
  "securityGroupIds": ["$SG_ID"],
  "instanceRole": "arn:aws:iam::${ACCOUNT}:instance-profile/ecsInstanceRole",
  "spotIamFleetRole": "arn:aws:iam::${ACCOUNT}:role/AmazonEC2SpotFleetTaggingRole",
  "ec2Configuration": [{"imageType": "ECS_AL2023_NVIDIA"}]
}
"@ | Set-Content -Encoding ascii $SpotJson.FullName
  aws batch create-compute-environment --compute-environment-name $env:COMPUTE_ENV_CAMPAIGN --type MANAGED --state ENABLED --service-role "arn:aws:iam::${ACCOUNT}:role/AWSBatchServiceRole" --compute-resources "file://$($SpotJson.FullName)" --region $env:AWS_REGION | Out-Null
  Remove-Item $SpotJson.FullName -Force
} else {
  Write-Host "  compute env already present (status=$CE_STATUS)"
}
Wait-CeValid $env:COMPUTE_ENV_CAMPAIGN

$CeOrder = @("order=1,computeEnvironment=$($env:COMPUTE_ENV_CAMPAIGN)")

if ($env:CAMPAIGN_ONDEMAND_FALLBACK -eq "1") {
  Write-Host "==> On-Demand fallback compute environment $($env:COMPUTE_ENV_CAMPAIGN_OD)"
  $OD_STATUS = Test-AwsText (aws batch describe-compute-environments --compute-environments $env:COMPUTE_ENV_CAMPAIGN_OD --query 'computeEnvironments[0].status' --output text --region $env:AWS_REGION 2>$null)
  if (-not $OD_STATUS) {
    $OdJson = New-TemporaryFile
    @"
{
  "type": "EC2",
  "allocationStrategy": "BEST_FIT_PROGRESSIVE",
  "minvCpus": 0,
  "maxvCpus": $($env:CAMPAIGN_MAX_VCPUS),
  "desiredvCpus": 0,
  "instanceTypes": ["g5.2xlarge", "g4dn.2xlarge"],
  "subnets": ["$SUBNET_ID"],
  "securityGroupIds": ["$SG_ID"],
  "instanceRole": "arn:aws:iam::${ACCOUNT}:instance-profile/ecsInstanceRole",
  "ec2Configuration": [{"imageType": "ECS_AL2023_NVIDIA"}]
}
"@ | Set-Content -Encoding ascii $OdJson.FullName
    aws batch create-compute-environment --compute-environment-name $env:COMPUTE_ENV_CAMPAIGN_OD --type MANAGED --state ENABLED --service-role "arn:aws:iam::${ACCOUNT}:role/AWSBatchServiceRole" --compute-resources "file://$($OdJson.FullName)" --region $env:AWS_REGION | Out-Null
    Remove-Item $OdJson.FullName -Force
  } else {
    Write-Host "  compute env already present (status=$OD_STATUS)"
  }
  Wait-CeValid $env:COMPUTE_ENV_CAMPAIGN_OD
  $CeOrder += "order=2,computeEnvironment=$($env:COMPUTE_ENV_CAMPAIGN_OD)"
}

Write-Host "==> Campaign job queue $($env:JOB_QUEUE_CAMPAIGN)"
$JQ_STATE = Test-AwsText (aws batch describe-job-queues --job-queues $env:JOB_QUEUE_CAMPAIGN --query 'jobQueues[0].state' --output text --region $env:AWS_REGION 2>$null)
if (-not $JQ_STATE) {
  aws batch create-job-queue --job-queue-name $env:JOB_QUEUE_CAMPAIGN --state ENABLED --priority 1 --compute-environment-order @CeOrder --region $env:AWS_REGION | Out-Null
} else {
  Write-Host "  job queue exists (state=$JQ_STATE)"
}

Write-Host "==> Campaign job definition $($env:JOB_DEFINITION_CAMPAIGN)"
# vcpus/memory sized to a full g5.2xlarge (8 vCPU, ~32 GB) so GPU=1 => one job per instance.
$JDJson = New-TemporaryFile
@"
{
  "image": "$IMAGE_URI",
  "vcpus": 8,
  "memory": 28000,
  "resourceRequirements": [{"type": "GPU", "value": "1"}],
  "jobRoleArn": "arn:aws:iam::${ACCOUNT}:role/gutibm-batch-job-role",
  "executionRoleArn": "arn:aws:iam::${ACCOUNT}:role/gutibm-batch-execution-role",
  "privileged": true,
  "linuxParameters": {"sharedMemorySize": 2048}
}
"@ | Set-Content -Encoding ascii $JDJson.FullName
$RetryJson = New-TemporaryFile
# Retry Spot reclaims (surfaced as "Host EC2*"); exit on application failures.
@"
{
  "attempts": 2,
  "evaluateOnExit": [
    {"onStatusReason": "Host EC2*", "action": "RETRY"},
    {"onReason": "*", "action": "EXIT"}
  ]
}
"@ | Set-Content -Encoding ascii $RetryJson.FullName
aws batch register-job-definition --job-definition-name $env:JOB_DEFINITION_CAMPAIGN --type container --platform-capabilities EC2 --container-properties "file://$($JDJson.FullName)" --retry-strategy "file://$($RetryJson.FullName)" --region $env:AWS_REGION | Out-Null
Remove-Item $JDJson.FullName, $RetryJson.FullName -Force

Write-Host "OK: campaign stack ready"
Write-Host "  computeEnv=$($env:COMPUTE_ENV_CAMPAIGN) (Spot, maxvCpus=$($env:CAMPAIGN_MAX_VCPUS))"
if ($env:CAMPAIGN_ONDEMAND_FALLBACK -eq "1") { Write-Host "  fallbackCE=$($env:COMPUTE_ENV_CAMPAIGN_OD) (On-Demand)" }
Write-Host "  queue=$($env:JOB_QUEUE_CAMPAIGN)"
Write-Host "  jobDefinition=$($env:JOB_DEFINITION_CAMPAIGN)"
Write-Host "  image=$IMAGE_URI"
