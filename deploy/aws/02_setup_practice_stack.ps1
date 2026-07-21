# Create S3 buckets + IAM roles + Batch practice stack (On-Demand g4dn.xlarge).
# PowerShell port of 02_setup_practice_stack.sh. Safe to re-run.
#
# From repo root:
#   $env:BUCKET = "your-existing-bucket"   # optional; else derives defaults
#   .\deploy\aws\02_setup_practice_stack.ps1

$ErrorActionPreference = "Stop"

$Root      = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$PolicyDir = Join-Path $PSScriptRoot "policies"
. (Join-Path $PSScriptRoot "env.ps1")

function Test-AwsText {
  # Returns trimmed stdout, or "" when the value is empty/None.
  param([string]$Value)
  if ([string]::IsNullOrWhiteSpace($Value) -or $Value.Trim() -eq "None") { return "" }
  return $Value.Trim()
}

function Ensure-Role {
  param([string]$Role, [string]$TrustFile)
  aws iam get-role --role-name $Role *> $null
  if ($LASTEXITCODE -eq 0) {
    Write-Host "  role exists: $Role"
  } else {
    Write-Host "  creating role: $Role"
    aws iam create-role --role-name $Role --assume-role-policy-document "file://$TrustFile" | Out-Null
  }
}

Write-Host "==> S3 buckets"
aws s3api head-bucket --bucket $INPUT_BUCKET *> $null
if ($LASTEXITCODE -ne 0) { aws s3 mb "s3://$INPUT_BUCKET" --region $env:AWS_REGION }
aws s3api head-bucket --bucket $OUTPUT_BUCKET *> $null
if ($LASTEXITCODE -ne 0) { aws s3 mb "s3://$OUTPUT_BUCKET" --region $env:AWS_REGION }

Write-Host "==> Default VPC + subnet + security group"
$VPC_ID = Test-AwsText (aws ec2 describe-vpcs --filters Name=isDefault,Values=true --query 'Vpcs[0].VpcId' --output text --region $env:AWS_REGION)
if (-not $VPC_ID) { throw "no default VPC in $($env:AWS_REGION). Create a VPC or pick subnets manually." }
$SUBNET_ID = Test-AwsText (aws ec2 describe-subnets --filters "Name=vpc-id,Values=$VPC_ID" --query 'Subnets[0].SubnetId' --output text --region $env:AWS_REGION)
Write-Host "  VPC_ID=$VPC_ID"
Write-Host "  SUBNET_ID=$SUBNET_ID"

$SG_ID = Test-AwsText (aws ec2 describe-security-groups --filters "Name=group-name,Values=gutibm-batch-sg" "Name=vpc-id,Values=$VPC_ID" --query 'SecurityGroups[0].GroupId' --output text --region $env:AWS_REGION 2>$null)
if (-not $SG_ID) {
  $SG_ID = (aws ec2 create-security-group --group-name gutibm-batch-sg --description "GutIBM Batch GPU egress" --vpc-id $VPC_ID --region $env:AWS_REGION --query GroupId --output text).Trim()
  Write-Host "  created SG $SG_ID"
} else {
  Write-Host "  SG exists: $SG_ID"
}
# Default SG already allows all egress; this extra rule is optional and may fail harmlessly.
aws ec2 authorize-security-group-egress --group-id $SG_ID --ip-permissions 'IpProtocol=-1,IpRanges=[{CidrIp=0.0.0.0/0}]' --region $env:AWS_REGION *> $null

Write-Host "==> IAM roles"
Ensure-Role "AWSBatchServiceRole" (Join-Path $PolicyDir "batch-service-trust.json")
aws iam attach-role-policy --role-name AWSBatchServiceRole --policy-arn arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole *> $null

Ensure-Role "ecsInstanceRole" (Join-Path $PolicyDir "ec2-trust.json")
aws iam attach-role-policy --role-name ecsInstanceRole --policy-arn arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role *> $null
aws iam get-instance-profile --instance-profile-name ecsInstanceRole *> $null
if ($LASTEXITCODE -ne 0) { aws iam create-instance-profile --instance-profile-name ecsInstanceRole | Out-Null }
aws iam add-role-to-instance-profile --instance-profile-name ecsInstanceRole --role-name ecsInstanceRole *> $null

Ensure-Role "AmazonEC2SpotFleetTaggingRole" (Join-Path $PolicyDir "spot-fleet-trust.json")
aws iam attach-role-policy --role-name AmazonEC2SpotFleetTaggingRole --policy-arn arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetTaggingRole *> $null

Ensure-Role "gutibm-batch-job-role" (Join-Path $PolicyDir "ecs-tasks-trust.json")
$S3Policy = New-TemporaryFile
@"
{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Action": ["s3:GetObject", "s3:PutObject", "s3:ListBucket", "s3:DeleteObject"],
    "Resource": [
      "arn:aws:s3:::$INPUT_BUCKET",
      "arn:aws:s3:::$INPUT_BUCKET/*",
      "arn:aws:s3:::$OUTPUT_BUCKET",
      "arn:aws:s3:::$OUTPUT_BUCKET/*"
    ]
  }]
}
"@ | Set-Content -Encoding ascii $S3Policy.FullName
aws iam put-role-policy --role-name gutibm-batch-job-role --policy-name gutibm-s3-access --policy-document "file://$($S3Policy.FullName)"
Remove-Item $S3Policy.FullName -Force

Ensure-Role "gutibm-batch-execution-role" (Join-Path $PolicyDir "ecs-tasks-trust.json")
aws iam attach-role-policy --role-name gutibm-batch-execution-role --policy-arn arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy *> $null

Write-Host "==> Wait for instance profile propagation"
Start-Sleep -Seconds 15

Write-Host "==> Batch compute environment $($env:COMPUTE_ENV)"
$CE_STATUS = Test-AwsText (aws batch describe-compute-environments --compute-environments $env:COMPUTE_ENV --query 'computeEnvironments[0].status' --output text --region $env:AWS_REGION 2>$null)
if (-not $CE_STATUS) {
  $CEJson = New-TemporaryFile
  @"
{
  "type": "EC2",
  "allocationStrategy": "BEST_FIT_PROGRESSIVE",
  "minvCpus": 0,
  "maxvCpus": 4,
  "desiredvCpus": 0,
  "instanceTypes": ["g4dn.xlarge"],
  "subnets": ["$SUBNET_ID"],
  "securityGroupIds": ["$SG_ID"],
  "instanceRole": "arn:aws:iam::${ACCOUNT}:instance-profile/ecsInstanceRole",
  "ec2Configuration": [{"imageType": "ECS_AL2023_NVIDIA"}]
}
"@ | Set-Content -Encoding ascii $CEJson.FullName
  aws batch create-compute-environment --compute-environment-name $env:COMPUTE_ENV --type MANAGED --state ENABLED --service-role "arn:aws:iam::${ACCOUNT}:role/AWSBatchServiceRole" --compute-resources "file://$($CEJson.FullName)" --region $env:AWS_REGION | Out-Null
  Remove-Item $CEJson.FullName -Force
} else {
  Write-Host "  compute env already present (status=$CE_STATUS)"
}

Write-Host "==> Wait until compute environment is VALID"
for ($i = 0; $i -lt 60; $i++) {
  $CE_STATUS = (aws batch describe-compute-environments --compute-environments $env:COMPUTE_ENV --query 'computeEnvironments[0].status' --output text --region $env:AWS_REGION).Trim()
  Write-Host "  status=$CE_STATUS"
  if ($CE_STATUS -eq "VALID") { break }
  if ($CE_STATUS -eq "INVALID") {
    aws batch describe-compute-environments --compute-environments $env:COMPUTE_ENV --query 'computeEnvironments[0].statusReason' --output text --region $env:AWS_REGION
    throw "compute environment INVALID"
  }
  Start-Sleep -Seconds 5
}
if ($CE_STATUS -ne "VALID") { throw "compute environment not VALID yet (last=$CE_STATUS)" }

Write-Host "==> Job queue $($env:JOB_QUEUE)"
$JQ_STATE = Test-AwsText (aws batch describe-job-queues --job-queues $env:JOB_QUEUE --query 'jobQueues[0].state' --output text --region $env:AWS_REGION 2>$null)
if (-not $JQ_STATE) {
  aws batch create-job-queue --job-queue-name $env:JOB_QUEUE --state ENABLED --priority 1 --compute-environment-order "order=1,computeEnvironment=$($env:COMPUTE_ENV)" --region $env:AWS_REGION | Out-Null
} else {
  Write-Host "  job queue exists (state=$JQ_STATE)"
}

Write-Host "==> Job definition $($env:JOB_DEFINITION)"
$JDJson = New-TemporaryFile
@"
{
  "image": "$IMAGE_URI",
  "vcpus": 4,
  "memory": 14000,
  "resourceRequirements": [{"type": "GPU", "value": "1"}],
  "jobRoleArn": "arn:aws:iam::${ACCOUNT}:role/gutibm-batch-job-role",
  "executionRoleArn": "arn:aws:iam::${ACCOUNT}:role/gutibm-batch-execution-role",
  "privileged": true,
  "linuxParameters": {"sharedMemorySize": 2048}
}
"@ | Set-Content -Encoding ascii $JDJson.FullName
aws batch register-job-definition --job-definition-name $env:JOB_DEFINITION --type container --platform-capabilities EC2 --container-properties "file://$($JDJson.FullName)" --region $env:AWS_REGION | Out-Null
Remove-Item $JDJson.FullName -Force

Write-Host "OK: practice stack ready"
Write-Host "  queue=$($env:JOB_QUEUE)"
Write-Host "  jobDefinition=$($env:JOB_DEFINITION)"
Write-Host "  image=$IMAGE_URI"
Write-Host "Next: .\deploy\aws\05_setup_campaign_stack.ps1"
