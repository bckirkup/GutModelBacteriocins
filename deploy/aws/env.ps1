# Shared env for GutIBM AWS deploy scripts (PowerShell port of env.sh).
# Dot-source this; do not run it alone:
#   . .\deploy\aws\env.ps1
#
# Mirrors deploy/aws/env.sh: pins region to us-east-1, derives the account and
# image URI, and resolves the S3 bucket names (overridable via $env:BUCKET /
# $env:INPUT_BUCKET / $env:OUTPUT_BUCKET). Region is overridable via
# $env:AWS_REGION (defaults to us-east-1).

$ErrorActionPreference = "Stop"

if (-not $env:AWS_REGION) { $env:AWS_REGION = "us-east-1" }
$env:AWS_DEFAULT_REGION = $env:AWS_REGION

if (-not (Get-Command aws -ErrorAction SilentlyContinue)) {
  throw "aws CLI not found. Install AWS CLI v2 (https://aws.amazon.com/cli/)."
}

$script:ACCOUNT = (aws sts get-caller-identity --query Account --output text 2>$null)
if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($ACCOUNT) -or $ACCOUNT -eq "None") {
  throw "Cannot resolve AWS account. Run: aws configure  (or set AWS credentials)."
}
$ACCOUNT = $ACCOUNT.Trim()

$REPO      = "$ACCOUNT.dkr.ecr.$($env:AWS_REGION).amazonaws.com/gutibm"
$IMAGE_URI = "${REPO}:cuda"

# Bucket resolution: BUCKET (shared) > INPUT_BUCKET/OUTPUT_BUCKET > derived defaults.
if ($env:BUCKET) {
  if (-not $env:INPUT_BUCKET)  { $env:INPUT_BUCKET  = $env:BUCKET }
  if (-not $env:OUTPUT_BUCKET) { $env:OUTPUT_BUCKET = $env:BUCKET }
} else {
  if (-not $env:INPUT_BUCKET)  { $env:INPUT_BUCKET  = "gutibm-inputs-$ACCOUNT" }
  if (-not $env:OUTPUT_BUCKET) { $env:OUTPUT_BUCKET = "gutibm-outputs-$ACCOUNT" }
}
$INPUT_BUCKET  = $env:INPUT_BUCKET
$OUTPUT_BUCKET = $env:OUTPUT_BUCKET

# Practice stack (02): On-Demand g4dn.xlarge.
if (-not $env:JOB_QUEUE)      { $env:JOB_QUEUE      = "gutibm-gpu-practice" }
if (-not $env:JOB_DEFINITION) { $env:JOB_DEFINITION = "gutibm-cuda" }
if (-not $env:COMPUTE_ENV)    { $env:COMPUTE_ENV    = "gutibm-gpu-practice-od" }

# Campaign stack (05): Spot GPU (g5.2xlarge / g4dn.2xlarge), one GPU per run.
if (-not $env:COMPUTE_ENV_CAMPAIGN)     { $env:COMPUTE_ENV_CAMPAIGN     = "gutibm-gpu-campaign-spot" }
if (-not $env:COMPUTE_ENV_CAMPAIGN_OD)  { $env:COMPUTE_ENV_CAMPAIGN_OD  = "gutibm-gpu-campaign-od" }
if (-not $env:JOB_QUEUE_CAMPAIGN)       { $env:JOB_QUEUE_CAMPAIGN       = "gutibm-gpu-campaign" }
if (-not $env:JOB_DEFINITION_CAMPAIGN)  { $env:JOB_DEFINITION_CAMPAIGN  = "gutibm-cuda-campaign" }
if (-not $env:CAMPAIGN_MAX_VCPUS)       { $env:CAMPAIGN_MAX_VCPUS       = "96" }
if (-not $env:CAMPAIGN_ONDEMAND_FALLBACK) { $env:CAMPAIGN_ONDEMAND_FALLBACK = "0" }

Write-Host "AWS_REGION=$($env:AWS_REGION)"
Write-Host "ACCOUNT=$ACCOUNT"
Write-Host "IMAGE_URI=$IMAGE_URI"
Write-Host "INPUT_BUCKET=s3://$INPUT_BUCKET"
Write-Host "OUTPUT_BUCKET=s3://$OUTPUT_BUCKET"
