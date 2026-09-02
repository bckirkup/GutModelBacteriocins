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
    *)
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
SUMMARY_FILE="${OUTDIR}/campaign_summary.json"
printf '%s\n' '{"jobs":[]}' > "${SUMMARY_FILE}"
FAILED_JOBS=()

reasons_to_json() {
  if [[ "$#" -eq 0 ]]; then
    printf '[]'
  else
    printf '%s\n' "$@" | jq -R . | jq -s .
  fi
}

append_campaign_summary() {
  local arm="$1"
  local scale="$2"
  local job_id="$3"
  local batch_status="$4"
  local log_file="$5"
  local extraction_file="$6"
  local shell_reasons_json="$7"
  local extraction_json='{"result_files":[],"found_block_count":0,"failure_reasons":[]}'
  local result_files
  local found_block_count
  local expected_block_count
  local combined_reasons
  if [[ -f "${extraction_file}" ]]; then
    extraction_json="$(cat "${extraction_file}")"
  fi
  result_files="$(jq -c '.result_files // []' <<< "${extraction_json}")"
  found_block_count="$(jq -c '.found_block_count // 0' <<< "${extraction_json}")"
  expected_block_count="$(jq -c '.expected_block_count // 3' <<< "${extraction_json}")"
  combined_reasons="$(jq -c \
    --argjson extraction "${extraction_json}" \
    --argjson shell "${shell_reasons_json}" \
    '($extraction.failure_reasons // []) + $shell | unique' <<< '{}')"
  jq \
    --arg arm "${arm}" \
    --arg scale "${scale}" \
    --arg job_id "${job_id}" \
    --arg batch_status "${batch_status}" \
    --arg log_file "${log_file}" \
    --argjson result_files "${result_files}" \
    --argjson expected_block_count "${expected_block_count}" \
    --argjson found_block_count "${found_block_count}" \
    --argjson failure_reasons "${combined_reasons}" \
    '
      .jobs += [{
        arm: $arm,
        scale: $scale,
        job_id: (if $job_id == "" then null else $job_id end),
        batch_status: $batch_status,
        log_file: $log_file,
        result_files: $result_files,
        expected_block_count: $expected_block_count,
        found_block_count: $found_block_count,
        failure_reasons: $failure_reasons
      }]
    ' "${SUMMARY_FILE}" > "${SUMMARY_FILE}.tmp" \
    && mv "${SUMMARY_FILE}.tmp" "${SUMMARY_FILE}"
}

run_benchmark_job() {
  local arm="$1"
  local scale="$2"
  local job_id=""
  local status="UNKNOWN"
  local log_stream=""
  local next_token=""
  local events=""
  local new_token=""
  local extraction_file="${OUTDIR}/.${arm}_${scale}.extraction.json"
  local submit_output=""
  local submit_status=0
  local describe_status=0
  local log_status=0
  local extraction_status=0
  local shell_reasons_json
  local status_line=""
  local -a reasons=()

  CURRENT_LOG_FILE="${OUTDIR}/${arm}_${scale}.log"
  : > "${CURRENT_LOG_FILE}"
  rm -f "${extraction_file}"
  submit_output="$(aws batch submit-job \
    --job-name "gutibm-gpubench-${arm}-${scale}-${GIT_SHORT_SHA}" \
    --job-queue "${BATCH_QUEUE}" \
    --job-definition "${JOB_DEF_ARN}" \
    --container-overrides "$(jq -cn \
      --arg arm "${arm}" --arg scale "${scale}" \
      '{environment: [{name: "BENCH_ARM", value: $arm}, {name: "BENCH_SCALE", value: $scale}]}' )" \
    --region "${AWS_REGION}" \
    --query jobId --output text 2>/dev/null)"
  submit_status=$?
  if [[ "${submit_status}" -ne 0 || -z "${submit_output}" ||
        "${submit_output}" == "None" ]]; then
    reasons+=("Batch submission failed")
    shell_reasons_json="$(reasons_to_json "${reasons[@]}")"
    append_campaign_summary "${arm}" "${scale}" "" "SUBMIT_FAILED" \
      "${CURRENT_LOG_FILE}" "${extraction_file}" "${shell_reasons_json}" \
      || echo "WARNING: could not update ${SUMMARY_FILE}" >&2
    return 1
  fi
  job_id="${submit_output}"
  CURRENT_JOB_ID="${job_id}"
  echo "JOB_ID=${job_id} ARM=${arm} SCALE=${scale}"
  while :; do
    status="$(aws batch describe-jobs \
      --jobs "${job_id}" --region "${AWS_REGION}" \
      --query 'jobs[0].status' --output text 2>/dev/null)"
    describe_status=$?
    if [[ "${describe_status}" -ne 0 || -z "${status}" || "${status}" == "None" ]]; then
      reasons+=("Batch status lookup failed")
      status="UNKNOWN"
      terminate_current_job
      break
    fi
    echo "  status=${status}"
    case "${status}" in
      SUCCEEDED|FAILED) break ;;
      SUBMITTED|PENDING|RUNNABLE|STARTING|RUNNING)
        sleep "${BATCH_POLL_SECONDS}" ;;
      *)
        reasons+=("unexpected Batch status: ${status}")
        terminate_current_job
        break
        ;;
    esac
  done
  [[ "${status}" == "FAILED" ]] && reasons+=("Batch job FAILED")

  for _ in $(seq 0 5 "${BATCH_LOG_WAIT_SECONDS}"); do
    log_stream="$(aws batch describe-jobs \
      --jobs "${job_id}" --region "${AWS_REGION}" \
      --query 'jobs[0].container.logStreamName' --output text \
      2>/dev/null)"
    if [[ -n "${log_stream}" && "${log_stream}" != "None" ]]; then
      break
    fi
    sleep 5
  done
  if [[ -z "${log_stream}" || "${log_stream}" == "None" ]]; then
    reasons+=("Batch job has no CloudWatch log stream")
  else
    while :; do
      if [[ -n "${next_token}" ]]; then
        events="$(aws logs get-log-events \
          --log-group-name "${BATCH_LOG_GROUP}" --log-stream-name "${log_stream}" \
          --start-from-head --next-token "${next_token}" --output json 2>/dev/null)"
      else
        events="$(aws logs get-log-events \
          --log-group-name "${BATCH_LOG_GROUP}" --log-stream-name "${log_stream}" \
          --start-from-head --output json 2>/dev/null)"
      fi
      log_status=$?
      if [[ "${log_status}" -ne 0 ]]; then
        reasons+=("CloudWatch log retrieval failed")
        break
      fi
      jq -r '.events[]?.message // empty' <<< "${events}" \
        | tee -a "${CURRENT_LOG_FILE}"
      new_token="$(jq -r '.nextForwardToken // empty' <<< "${events}")"
      [[ -z "${new_token}" || "${new_token}" == "${next_token}" ]] && break
      next_token="${new_token}"
    done
  fi
  if [[ "${arm}" =~ ^A[456]$ ]] &&
     ! grep -Eq 'NVIDIA-SMI [0-9]' "${CURRENT_LOG_FILE}"; then
    reasons+=("GPU log lacks NVIDIA-SMI evidence")
  fi
  python3 - "${CURRENT_LOG_FILE}" "${OUTDIR}" "${arm}" "${scale}" \
    "${extraction_file}" <<'PY'
import json
import re
import sys
from pathlib import Path

log_path, output_dir, arm, scale, extraction_path = sys.argv[1:]
text = Path(log_path).read_text(encoding="utf-8")
pattern = re.compile(
    rf"===BENCH_RESULT_BEGIN ({re.escape(arm)}_{re.escape(scale)}_seed\d+)===\n"
    r"(.*?)\n===BENCH_RESULT_END===",
    re.DOTALL,
)
matches = list(pattern.finditer(text))
summary = {
    "expected_block_count": 3,
    "found_block_count": 0,
    "result_files": [],
    "failure_reasons": [],
}
begin_ids = set(re.findall(
    rf"^===BENCH_RESULT_BEGIN {re.escape(arm)}_{re.escape(scale)}_seed\d+===$",
    text,
    re.MULTILINE,
))
blocks = {}
for match in matches:
    identifier = match.group(1)
    body = match.group(2)
    previous = blocks.get(identifier)
    if previous is not None:
        if previous != body:
            summary["failure_reasons"].append(
                f"conflicting duplicate result block: {identifier}"
            )
        continue
    blocks[identifier] = body
summary["found_block_count"] = len(blocks)
if len(begin_ids) > len(blocks):
    summary["failure_reasons"].append(
        f"incomplete result block(s): {len(begin_ids) - len(blocks)}"
    )
if len(blocks) < summary["expected_block_count"]:
    summary["failure_reasons"].append(
        "result block shortfall: expected "
        f"{summary['expected_block_count']}, found {len(blocks)}"
    )
for identifier, body in blocks.items():
    try:
        payload = json.loads(body)
    except json.JSONDecodeError as error:
        summary["failure_reasons"].append(
            f"unparseable JSON for {identifier}: {error.msg}"
        )
        continue
    result_path = Path(output_dir, f"{identifier}.json")
    result_path.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    summary["result_files"].append(str(result_path))
if len(summary["result_files"]) < summary["expected_block_count"]:
    summary["failure_reasons"].append(
        "valid result file shortfall: expected "
        f"{summary['expected_block_count']}, found "
        f"{len(summary['result_files'])}"
    )
Path(extraction_path).write_text(
    json.dumps(summary, indent=2) + "\n", encoding="utf-8"
)
PY
  extraction_status=$?
  [[ "${extraction_status}" -eq 0 ]] \
    || reasons+=("result extraction failed")
  status_line="$(grep -E 'BENCH_JOB_STATUS=[^[:space:]]+' \
    "${CURRENT_LOG_FILE}" | tail -1 || true)"
  if [[ -z "${status_line}" ]]; then
    reasons+=("missing BENCH_JOB_STATUS")
  elif [[ "${status_line}" != "BENCH_JOB_STATUS=0" ]]; then
    reasons+=("container reported ${status_line}")
  fi
  shell_reasons_json="$(reasons_to_json "${reasons[@]}")"
  append_campaign_summary "${arm}" "${scale}" "${job_id}" "${status}" \
    "${CURRENT_LOG_FILE}" "${extraction_file}" "${shell_reasons_json}" \
    || reasons+=("campaign summary update failed")
  CURRENT_JOB_ID=""
  CURRENT_LOG_FILE=""
  if [[ "${#reasons[@]}" -gt 0 ]]; then
    return 1
  fi
  return 0
}

for arm in "${ARM_ARRAY[@]}"; do
  for scale in "${SCALE_ARRAY[@]}"; do
    if ! run_benchmark_job "${arm}" "${scale}"; then
      FAILED_JOBS+=("${arm}/${scale}")
      echo "FAILED: ${arm}/${scale}; continuing campaign" >&2
    fi
  done
done
if [[ "${#FAILED_JOBS[@]}" -gt 0 ]]; then
  printf 'Campaign completed with failed jobs:\n' >&2
  printf '  %s\n' "${FAILED_JOBS[@]}" >&2
  exit 1
fi
echo "GPU benchmark retrieval complete."
