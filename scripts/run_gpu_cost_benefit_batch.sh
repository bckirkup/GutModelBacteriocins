#!/usr/bin/env bash
# Build, publish, and run the guarded GPU cost/benefit benchmark.
#
# This script is deliberately limited to the practice GPU queue. It submits
# one arm/scale job at a time and retrieves raw JSON through CloudWatch logs.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
AWS_REGION="${AWS_REGION:-us-east-1}"
export AWS_REGION AWS_DEFAULT_REGION="${AWS_REGION}"
ECR_REPOSITORY="${ECR_REPOSITORY:-}"
BATCH_QUEUE="${BATCH_QUEUE:-gutibm-gpu-practice}"
BATCH_JOB_DEFINITION="${BATCH_JOB_DEFINITION:-gutibm-gpubench}"
BATCH_LOG_GROUP="${BATCH_LOG_GROUP:-/aws/batch/job}"
BATCH_POLL_SECONDS="${BATCH_POLL_SECONDS:-20}"
BATCH_LOG_WAIT_SECONDS="${BATCH_LOG_WAIT_SECONDS:-300}"
BENCH_ARMS="${BENCH_ARMS:-A1 A2 A3 A4 A5 A6}"
BENCH_SCALES="${BENCH_SCALES:-s1 s2}"
BENCH_S3_URI="${GUTIBM_BENCH_S3_URI:-}"
OUTDIR="${BENCH_OUTDIR:-${ROOT}/bench_results/gpubench}"
IMAGE_TAG="${GUTIBM_IMAGE_TAG:-}"
DRY_RUN=0
CURRENT_JOB_ID=""
CURRENT_JOB_DEF_JSON=""
CURRENT_LOG_FILE=""
CLEANUP_RUNNING=0

die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_command() {
  local command_name="$1"
  command -v "${command_name}" >/dev/null 2>&1 \
    || die "required command not found: ${command_name}"
}

cleanup_temp_files() {
  [[ "${CLEANUP_RUNNING}" -eq 0 ]] || return 0
  CLEANUP_RUNNING=1
  rm -f "${CURRENT_JOB_DEF_JSON}" || \
    echo "WARNING: could not remove temporary Batch files" >&2
}

terminate_current_job() {
  local current_status
  [[ -n "${CURRENT_JOB_ID}" ]] || return 0
  current_status="$(aws batch describe-jobs \
    --jobs "${CURRENT_JOB_ID}" \
    --region "${AWS_REGION}" \
    --query 'jobs[0].status' \
    --output text 2>/dev/null || true)"
  case "${current_status}" in
    SUBMITTED|PENDING|RUNNABLE|STARTING|RUNNING)
      echo "Terminating non-terminal benchmark job ${CURRENT_JOB_ID}." >&2
      aws batch terminate-job \
        --job-id "${CURRENT_JOB_ID}" \
        --reason "GPU benchmark runner exited before job completion" \
        --region "${AWS_REGION}" >/dev/null 2>&1 || \
        echo "WARNING: could not terminate Batch job ${CURRENT_JOB_ID}" >&2
      ;;
  esac
}

handle_exit() {
  local exit_code=$?
  trap - EXIT INT TERM
  terminate_current_job
  cleanup_temp_files
  exit "${exit_code}"
}

handle_signal() {
  local signal="$1"
  trap - EXIT INT TERM
  terminate_current_job
  cleanup_temp_files
  if [[ "${signal}" == "INT" ]]; then
    exit 130
  fi
  exit 143
}

trap handle_exit EXIT
trap 'handle_signal INT' INT
trap 'handle_signal TERM' TERM

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --outdir)
      [[ $# -ge 2 ]] || die "--outdir requires a path"
      OUTDIR="$2"
      shift 2
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

require_command aws
require_command docker
require_command git
require_command jq
require_command python3
require_command awk
require_command grep

[[ "${GUTIBM_BENCH_AUTHORIZED:-}" == "1" ]] \
  || die "set GUTIBM_BENCH_AUTHORIZED=1 before running the benchmark"
[[ -f "${ROOT}/deploy/aws/Dockerfile.gpubench" ]] \
  || die "missing deploy/aws/Dockerfile.gpubench"
[[ -f "${ROOT}/python/gut_ibm_tools/gpu_cost_benefit.py" ]] \
  || die "missing GPU benchmark harness"

case "${BATCH_QUEUE}" in
  gutibm-gpu-practice|gutibm-gpu-practice-*) ;;
  *) die "BATCH_QUEUE must be the practice queue: ${BATCH_QUEUE}" ;;
esac
case "${BATCH_JOB_DEFINITION}" in
  gutibm-gpubench|gutibm-gpubench-*) ;;
  *) die "BATCH_JOB_DEFINITION must be a gutibm-gpubench definition: ${BATCH_JOB_DEFINITION}" ;;
esac
[[ "${BATCH_QUEUE}" != *campaign* ]] \
  || die "campaign Batch queues are forbidden: ${BATCH_QUEUE}"
[[ "${BATCH_JOB_DEFINITION}" != *campaign* ]] \
  || die "campaign Batch job definitions are forbidden: ${BATCH_JOB_DEFINITION}"
[[ "${BENCH_S3_URI}" != *campaign* ]] \
  || die "campaign S3 URIs are forbidden: ${BENCH_S3_URI}"

GIT_SHA="$(git -C "${ROOT}" rev-parse HEAD)"
GIT_SHORT_SHA="$(git -C "${ROOT}" rev-parse --short=12 HEAD)"
IMAGE_TAG="${IMAGE_TAG:-gpubench-${GIT_SHORT_SHA}}"
[[ "${IMAGE_TAG}" =~ ^gpubench-[0-9a-fA-F]{7,12}$ ]] \
  || die "GUTIBM_IMAGE_TAG must be gpubench-<short sha>: ${IMAGE_TAG}"
[[ "${IMAGE_TAG}" != *campaign* ]] \
  || die "campaign image tags are forbidden"

read -r -a ARM_ARRAY <<< "${BENCH_ARMS}"
read -r -a SCALE_ARRAY <<< "${BENCH_SCALES}"
[[ "${#ARM_ARRAY[@]}" -gt 0 && "${#SCALE_ARRAY[@]}" -gt 0 ]] \
  || die "BENCH_ARMS and BENCH_SCALES must not be empty"
for arm in "${ARM_ARRAY[@]}"; do
  [[ "${arm}" =~ ^A[1-6]$ ]] \
    || die "unsupported benchmark arm: ${arm}"
done
for scale in "${SCALE_ARRAY[@]}"; do
  [[ "${scale}" =~ ^s[012]$ ]] \
    || die "unsupported benchmark scale: ${scale}"
done

ACCOUNT="$(aws sts get-caller-identity \
  --query Account --output text --region "${AWS_REGION}")"
[[ -n "${ACCOUNT}" && "${ACCOUNT}" != "None" ]] \
  || die "could not determine the AWS account"
if [[ -z "${ECR_REPOSITORY}" ]]; then
  ECR_REPOSITORY="${ACCOUNT}.dkr.ecr.${AWS_REGION}.amazonaws.com/gutibm"
fi
ECR_REGISTRY="${ECR_REPOSITORY%/*}"
ECR_REPOSITORY_NAME="${ECR_REPOSITORY##*/}"
[[ "${ECR_REPOSITORY}" != *campaign* ]] \
  || die "campaign ECR repositories are forbidden"
[[ "${ECR_REPOSITORY_NAME}" == gutibm ]] \
  || die "ECR_REPOSITORY must name the gutibm repository: ${ECR_REPOSITORY}"
IMAGE_URI="${ECR_REPOSITORY}:${IMAGE_TAG}"

echo "Read-only AWS preflight"
echo "  account=${ACCOUNT}"
echo "  region=${AWS_REGION}"
if aws batch describe-job-queues \
  --job-queues "${BATCH_QUEUE}" \
  --region "${AWS_REGION}" \
  --query 'jobQueues[0].jobQueueName' --output text >/dev/null 2>&1; then
  echo "  practice_queue=reachable"
else
  die "practice queue is not readable: ${BATCH_QUEUE}"
fi
if aws batch describe-job-definitions \
  --job-definition-name "${BATCH_JOB_DEFINITION}" \
  --status ACTIVE --region "${AWS_REGION}" \
  --query 'jobDefinitions | length(@)' --output text >/dev/null 2>&1; then
  echo "  benchmark_job_definition=readable"
else
  echo "  benchmark_job_definition=not currently readable or registered"
fi
if aws ecr describe-repositories \
  --repository-names "${ECR_REPOSITORY_NAME}" \
  --region "${AWS_REGION}" >/dev/null 2>&1; then
  echo "  ecr_repository=reachable"
else
  die "ECR repository is not readable: ${ECR_REPOSITORY_NAME}"
fi

if [[ -n "${BENCH_S3_URI}" ]]; then
  S3_PATH="${BENCH_S3_URI#s3://}"
  S3_BUCKET="${S3_PATH%%/*}"
  [[ -n "${S3_BUCKET}" && "${S3_BUCKET}" != "${S3_PATH}" ]] \
    || die "GUTIBM_BENCH_S3_URI must include an S3 bucket and prefix"
  if aws s3api head-bucket --bucket "${S3_BUCKET}" \
    --region "${AWS_REGION}" >/dev/null 2>&1; then
    echo "  s3_bucket=reachable; write=not provable by read-only checks"
  else
    echo "  s3_bucket=not readable; optional S3 copies may fail"
  fi
else
  echo "  s3=disabled; write check not applicable"
fi

if aws ecr describe-images \
  --repository-name "${ECR_REPOSITORY_NAME}" \
  --image-ids "imageTag=${IMAGE_TAG}" \
  --region "${AWS_REGION}" \
  --query 'imageDetails[0].imageDigest' --output text 2>/dev/null \
  | grep -qv '^None$'; then
  die "refusing to overwrite existing ECR tag: ${IMAGE_TAG}"
fi

echo "Build ${IMAGE_URI}"
aws ecr get-login-password --region "${AWS_REGION}" \
  | docker login --username AWS --password-stdin "${ECR_REGISTRY}"
docker build \
  --file "${ROOT}/deploy/aws/Dockerfile.gpubench" \
  --build-arg "GUTIBM_GIT_SHA=${GIT_SHA}" \
  --tag "${IMAGE_URI}" \
  "${ROOT}"

IMAGE_REF="${IMAGE_URI}"
if [[ "${DRY_RUN}" -eq 1 ]]; then
  IMAGE_REF="${ECR_REPOSITORY}@sha256:dry-run-no-push"
else
  docker push "${IMAGE_URI}"
  IMAGE_DIGEST=""
  for _ in $(seq 1 12); do
    IMAGE_DIGEST="$(aws ecr describe-images \
      --repository-name "${ECR_REPOSITORY_NAME}" \
      --image-ids "imageTag=${IMAGE_TAG}" \
      --region "${AWS_REGION}" \
      --query 'imageDetails[0].imageDigest' \
      --output text 2>/dev/null || true)"
    [[ -n "${IMAGE_DIGEST}" && "${IMAGE_DIGEST}" != "None" ]] && break
    sleep 5
  done
  [[ -n "${IMAGE_DIGEST}" && "${IMAGE_DIGEST}" != "None" ]] \
    || die "ECR did not return a digest for ${IMAGE_TAG}"
  IMAGE_REF="${ECR_REPOSITORY}@${IMAGE_DIGEST}"
  echo "IMAGE_DIGEST=${IMAGE_DIGEST}"
fi

JOB_COUNT=$(( ${#ARM_ARRAY[@]} * ${#SCALE_ARRAY[@]} ))
ESTIMATED_HOURS_PER_JOB="${BENCH_ESTIMATED_HOURS_PER_JOB:-2}"
if ! awk -v value="${ESTIMATED_HOURS_PER_JOB}" \
  'BEGIN { exit !(value > 0) }'; then
  die "BENCH_ESTIMATED_HOURS_PER_JOB must be positive"
fi
ESTIMATED_INSTANCE_HOURS="$(awk -v jobs="${JOB_COUNT}" -v hours="${ESTIMATED_HOURS_PER_JOB}" \
  'BEGIN { printf "%.3f", jobs * hours }')"
ESTIMATED_COST="$(awk -v hours="${ESTIMATED_INSTANCE_HOURS}" \
  'BEGIN { printf "%.2f", hours * 0.526 }')"

echo "Benchmark plan (serial jobs, 4 vCPU, 1 GPU, 13000 MiB, max attempt 2h)"
for arm in "${ARM_ARRAY[@]}"; do
  for scale in "${SCALE_ARRAY[@]}"; do
    echo "  ${arm}/${scale}"
  done
done
echo "  jobs=${JOB_COUNT}"
echo "  estimated_instance_hours=${ESTIMATED_INSTANCE_HOURS}"
echo "  estimated_cost_usd_at_0.526_per_hour=${ESTIMATED_COST}"

if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo "DRY_RUN: skipping docker push, job-definition registration, and submission"
fi

CURRENT_JOB_DEF_JSON="$(mktemp)"
# The variables expand in the Batch container, not while constructing JSON.
# shellcheck disable=SC2016
BENCH_CONTAINER_SCRIPT='
set +e
job_status=0
work=/tmp/gutibm-gpubench
rm -rf "${work}"
mkdir -p "${work}/config" "${work}/results"
if [[ "${BENCH_ARM}" =~ ^A[456]$ ]]; then
  nvidia-smi || job_status=1
fi
python3 -m gut_ibm_tools.gpu_cost_benefit generate \
  --base "examples/scaling_benchmark/input_bench_${BENCH_SCALE}.json" \
  --output-dir "${work}/config" \
  --scale "${BENCH_SCALE}" \
    "examples/scaling_benchmark/input_bench_${BENCH_SCALE}.json"
if [[ "$?" -ne 0 ]]; then
  job_status=1
fi
for seed in 55 56 57; do
  result="${work}/results/${BENCH_ARM}_${BENCH_SCALE}_seed${seed}.json"
  python3 -m gut_ibm_tools.gpu_cost_benefit run \
    --manifest "${work}/config/manifest.json" \
    --arm "${BENCH_ARM}" --scale "${BENCH_SCALE}" --seed "${seed}" \
    --binary /src/build/gut_ibm --output-dir "${work}/results"
  run_status=$?
  if [[ -f "${result}" ]]; then
    echo "===BENCH_RESULT_BEGIN ${BENCH_ARM}_${BENCH_SCALE}_seed${seed}==="
    cat "${result}"
    echo "===BENCH_RESULT_END==="
  else
    echo "WARNING: missing result JSON for seed ${seed}" >&2
    job_status=1
  fi
  if [[ -n "${GUTIBM_BENCH_S3_URI:-}" ]]; then
    aws s3 cp "${result}" \
      "${GUTIBM_BENCH_S3_URI%/}/${BENCH_ARM}_${BENCH_SCALE}_seed${seed}.json" \
      || echo "WARNING: failed optional S3 copy for ${result}" >&2
    for hdf5 in "${work}/results/${BENCH_ARM}_${BENCH_SCALE}_seed${seed}."*.h5; do
      [[ -f "${hdf5}" ]] || continue
      aws s3 cp "${hdf5}" \
        "${GUTIBM_BENCH_S3_URI%/}/$(basename "${hdf5}")" \
        || echo "WARNING: failed optional S3 copy for ${hdf5}" >&2
    done
  fi
  [[ "${run_status}" -eq 0 ]] || job_status=1
done
echo "BENCH_JOB_STATUS=${job_status}"
exit "${job_status}"
'
jq -n \
  --arg name "${BATCH_JOB_DEFINITION}" \
  --arg image "${IMAGE_REF}" \
  --arg command "${BENCH_CONTAINER_SCRIPT}" \
  --arg log_group "${BATCH_LOG_GROUP}" \
  --arg s3_uri "${BENCH_S3_URI}" \
  '{
    jobDefinitionName: $name,
    type: "container",
    platformCapabilities: ["EC2"],
    containerProperties: {
      image: $image,
      command: ["-lc", $command],
      vcpus: 4,
      memory: 13000,
      resourceRequirements: [{type: "GPU", value: "1"}],
      linuxParameters: {sharedMemorySize: 2048},
      environment: (
        [
          {name: "PYTHONUNBUFFERED", value: "1"},
          {name: "BATCH_LOG_GROUP", value: $log_group}
        ]
        + (if $s3_uri == "" then [] else
          [{name: "GUTIBM_BENCH_S3_URI", value: $s3_uri}] end)
      )
    },
    timeout: {attemptDurationSeconds: 7200},
    retryStrategy: {attempts: 1}
  }' > "${CURRENT_JOB_DEF_JSON}"

if [[ "${DRY_RUN}" -eq 1 ]]; then
  exit 0
fi

JOB_DEF_ARN="$(aws batch register-job-definition \
  --cli-input-json "file://${CURRENT_JOB_DEF_JSON}" \
  --region "${AWS_REGION}" \
  --query 'jobDefinitionArn' --output text)"
[[ "${JOB_DEF_ARN}" == *":job-definition/${BATCH_JOB_DEFINITION}:"* ]] \
  || die "registered unexpected job definition: ${JOB_DEF_ARN}"
echo "JOB_DEFINITION_ARN=${JOB_DEF_ARN}"

mkdir -p "${OUTDIR}"
for arm in "${ARM_ARRAY[@]}"; do
  for scale in "${SCALE_ARRAY[@]}"; do
    CURRENT_LOG_FILE="${OUTDIR}/${arm}_${scale}.log"
    : > "${CURRENT_LOG_FILE}"
    CURRENT_JOB_ID="$(aws batch submit-job \
      --job-name "gutibm-gpubench-${arm}-${scale}-${GIT_SHORT_SHA}" \
      --job-queue "${BATCH_QUEUE}" \
      --job-definition "${JOB_DEF_ARN}" \
      --container-overrides "$(jq -cn \
        --arg arm "${arm}" --arg scale "${scale}" \
        '{environment: [{name: "BENCH_ARM", value: $arm}, {name: "BENCH_SCALE", value: $scale}]}' )" \
      --region "${AWS_REGION}" \
      --query jobId --output text)"
    [[ -n "${CURRENT_JOB_ID}" && "${CURRENT_JOB_ID}" != "None" ]] \
      || die "Batch submission returned no job ID"
    echo "JOB_ID=${CURRENT_JOB_ID} ARM=${arm} SCALE=${scale}"
    while :; do
      status="$(aws batch describe-jobs \
        --jobs "${CURRENT_JOB_ID}" --region "${AWS_REGION}" \
        --query 'jobs[0].status' --output text)"
      echo "  status=${status}"
      case "${status}" in
        SUCCEEDED|FAILED) break ;;
        SUBMITTED|PENDING|RUNNABLE|STARTING|RUNNING)
          sleep "${BATCH_POLL_SECONDS}" ;;
        *) die "unexpected Batch status: ${status}" ;;
      esac
    done
    log_stream=""
    for _ in $(seq 0 5 "${BATCH_LOG_WAIT_SECONDS}"); do
      log_stream="$(aws batch describe-jobs \
        --jobs "${CURRENT_JOB_ID}" --region "${AWS_REGION}" \
        --query 'jobs[0].container.logStreamName' --output text \
        2>/dev/null || true)"
      if [[ -n "${log_stream}" && "${log_stream}" != "None" ]]; then
        break
      fi
      sleep 5
    done
    [[ -n "${log_stream}" && "${log_stream}" != "None" ]] \
      || die "Batch job has no CloudWatch log stream: ${CURRENT_JOB_ID}"
    next_token=""
    while :; do
      if [[ -n "${next_token}" ]]; then
        events="$(aws logs get-log-events \
          --log-group-name "${BATCH_LOG_GROUP}" --log-stream-name "${log_stream}" \
          --start-from-head --next-token "${next_token}" --output json)"
      else
        events="$(aws logs get-log-events \
          --log-group-name "${BATCH_LOG_GROUP}" --log-stream-name "${log_stream}" \
          --start-from-head --output json)"
      fi
      jq -r '.events[]?.message // empty' <<< "${events}" | tee -a "${CURRENT_LOG_FILE}"
      new_token="$(jq -r '.nextForwardToken // empty' <<< "${events}")"
      [[ -z "${new_token}" || "${new_token}" == "${next_token}" ]] && break
      next_token="${new_token}"
    done
    if [[ "${arm}" =~ ^A[456]$ ]]; then
      grep -Eq 'NVIDIA-SMI [0-9]' "${CURRENT_LOG_FILE}" \
        || die "GPU benchmark log lacks NVIDIA-SMI evidence: ${CURRENT_LOG_FILE}"
    fi
    python3 - "${CURRENT_LOG_FILE}" "${OUTDIR}" "${arm}" "${scale}" <<'PY'
import json
import re
import sys
from pathlib import Path

log_path, output_dir, arm, scale = sys.argv[1:]
text = Path(log_path).read_text(encoding="utf-8")
pattern = re.compile(
    rf"===BENCH_RESULT_BEGIN ({re.escape(arm)}_{re.escape(scale)}_seed\d+)===\n"
    r"(.*?)\n===BENCH_RESULT_END===",
    re.DOTALL,
)
matches = list(pattern.finditer(text))
if len(matches) != 3:
    raise SystemExit(
        f"expected 3 benchmark result blocks for {arm}/{scale}, found {len(matches)}"
    )
for match in matches:
    identifier = match.group(1)
    payload = json.loads(match.group(2))
    Path(output_dir, f"{identifier}.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
PY
    status_line="$(grep -E 'BENCH_JOB_STATUS=[^[:space:]]+' \
      "${CURRENT_LOG_FILE}" | tail -1 || true)"
    [[ -n "${status_line}" ]] \
      || die "missing BENCH_JOB_STATUS in ${CURRENT_LOG_FILE}"
    [[ "${status}" == "SUCCEEDED" ]] \
      || die "Batch benchmark job failed: ${CURRENT_JOB_ID}"
    CURRENT_JOB_ID=""
  done
done
echo "GPU benchmark retrieval complete."
