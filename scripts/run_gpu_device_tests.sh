#!/usr/bin/env bash
# Build, publish, and run the CUDA GPU test image on AWS Batch.
#
# The script intentionally targets only the practice GPU resources. It never
# creates, updates, submits to, or reads from the campaign queue or job
# definitions.
#
# Optional environment overrides:
#   AWS_REGION              (default: us-east-1)
#   ECR_REPOSITORY          (default: <account>.dkr.ecr.<region>.amazonaws.com/gutibm)
#   GUTIBM_IMAGE_TAG        (default: gputest-<short git sha>)
#   BATCH_QUEUE             (default: gutibm-gpu-practice)
#   BATCH_JOB_DEFINITION    (default: gutibm-gputest)
#   BATCH_LOG_GROUP         (default: /aws/batch/job)
#   BATCH_POLL_SECONDS      (default: 20)
#   BATCH_LOG_WAIT_SECONDS  (default: 300)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
AWS_REGION="${AWS_REGION:-us-east-1}"
export AWS_REGION AWS_DEFAULT_REGION="${AWS_REGION}"
ECR_REPOSITORY="${ECR_REPOSITORY:-}"
BATCH_QUEUE="${BATCH_QUEUE:-gutibm-gpu-practice}"
BATCH_JOB_DEFINITION="${BATCH_JOB_DEFINITION:-gutibm-gputest}"
BATCH_LOG_GROUP="${BATCH_LOG_GROUP:-/aws/batch/job}"
BATCH_POLL_SECONDS="${BATCH_POLL_SECONDS:-20}"
BATCH_LOG_WAIT_SECONDS="${BATCH_LOG_WAIT_SECONDS:-300}"
JOB_ID=""
STATUS=""
CLEANUP_RUNNING=0
JOB_DEF_JSON=""
LOG_FILE=""

die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_command() {
  local command_name="$1"
  command -v "${command_name}" >/dev/null 2>&1 \
    || die "required command not found: ${command_name}"
}

require_command aws
require_command docker
require_command git
require_command jq

terminate_submitted_job() {
  local current_status
  [[ -n "${JOB_ID}" && "${CLEANUP_RUNNING}" -eq 0 ]] || return 0
  CLEANUP_RUNNING=1
  current_status="$(aws batch describe-jobs \
    --jobs "${JOB_ID}" \
    --region "${AWS_REGION}" \
    --query 'jobs[0].status' \
    --output text 2>/dev/null || true)"
  case "${current_status}" in
    SUBMITTED|PENDING|RUNNABLE|STARTING|RUNNING)
      echo "Terminating non-terminal Batch job ${JOB_ID} after local exit." >&2
      aws batch terminate-job \
        --job-id "${JOB_ID}" \
        --reason "GPU device gate runner exited before job completion" \
        --region "${AWS_REGION}" >/dev/null 2>&1 || \
        echo "WARNING: could not terminate Batch job ${JOB_ID}" >&2
      ;;
    *)
      [[ -z "${current_status}" ]] || \
        echo "WARNING: unexpected status during Batch cleanup: ${current_status}" >&2
      ;;
  esac
}

cleanup_temp_files() {
  if [[ -n "${JOB_DEF_JSON}" || -n "${LOG_FILE}" ]]; then
    rm -f "${JOB_DEF_JSON}" "${LOG_FILE}" || \
      echo "WARNING: could not remove temporary GPU gate files" >&2
  fi
}

handle_exit() {
  local exit_code=$?
  trap - EXIT INT TERM
  terminate_submitted_job
  cleanup_temp_files
  exit "${exit_code}"
}

handle_signal() {
  local signal="$1"
  trap - EXIT INT TERM
  terminate_submitted_job
  cleanup_temp_files
  if [[ "${signal}" == "INT" ]]; then
    exit 130
  fi
  exit 143
}

trap handle_exit EXIT
trap 'handle_signal INT' INT
trap 'handle_signal TERM' TERM

[[ -f "${ROOT}/deploy/aws/Dockerfile.gputest" ]] \
  || die "missing deploy/aws/Dockerfile.gputest"
[[ -f "${ROOT}/.github/workflows/ci.yml" ]] \
  || die "missing .github/workflows/ci.yml"

ACCOUNT="$(aws sts get-caller-identity --query Account --output text)"
[[ -n "${ACCOUNT}" && "${ACCOUNT}" != "None" ]] \
  || die "could not determine the AWS account"

if [[ -z "${ECR_REPOSITORY}" ]]; then
  ECR_REPOSITORY="${ACCOUNT}.dkr.ecr.${AWS_REGION}.amazonaws.com/gutibm"
fi
ECR_REGISTRY="${ECR_REPOSITORY%/*}"
ECR_REPOSITORY_NAME="${ECR_REPOSITORY##*/}"

case "${ECR_REPOSITORY}" in
  */gutibm) ;;
  *) die "ECR_REPOSITORY must name the gutibm repository: ${ECR_REPOSITORY}" ;;
esac
[[ "${ECR_REPOSITORY}" != *campaign* ]] \
  || die "campaign ECR repositories are forbidden"
[[ "${BATCH_QUEUE}" != *campaign* ]] \
  || die "campaign Batch queues are forbidden: ${BATCH_QUEUE}"
[[ "${BATCH_JOB_DEFINITION}" != *campaign* ]] \
  || die "campaign Batch job definitions are forbidden: ${BATCH_JOB_DEFINITION}"
[[ "${BATCH_QUEUE}" == gutibm-gpu-practice || "${BATCH_QUEUE}" == gutibm-gpu-practice-* ]] \
  || die "BATCH_QUEUE must be the practice queue: ${BATCH_QUEUE}"
[[ "${BATCH_JOB_DEFINITION}" == gutibm-gputest || "${BATCH_JOB_DEFINITION}" == gutibm-gputest-* ]] \
  || die "BATCH_JOB_DEFINITION must be a gutibm-gputest definition: ${BATCH_JOB_DEFINITION}"

GIT_SHA="$(git -C "${ROOT}" rev-parse HEAD)"
GIT_SHORT_SHA="$(git -C "${ROOT}" rev-parse --short=12 HEAD)"
IMAGE_TAG="${GUTIBM_IMAGE_TAG:-gputest-${GIT_SHORT_SHA}}"
[[ "${IMAGE_TAG}" =~ ^[a-zA-Z0-9_][a-zA-Z0-9_.-]{0,127}$ ]] \
  || die "invalid Docker image tag: ${IMAGE_TAG}"
[[ "${IMAGE_TAG}" != *campaign* ]] || die "campaign image tags are forbidden"
IMAGE_URI="${ECR_REPOSITORY}:${IMAGE_TAG}"

GPU_TARGETS="${GPU_TARGETS:-}"
if [[ -z "${GPU_TARGETS}" ]]; then
  GPU_TARGETS="$(
    awk '
      /cmake --build build-cuda/ { capture = 1 }
      capture && /--target/ {
        for (i = 1; i <= NF; ++i) {
          if ($i == "--target") {
            for (j = i + 1; j <= NF; ++j) print $j
            exit
          }
        }
      }
    ' "${ROOT}/.github/workflows/ci.yml" | tr '\n' ' '
  )"
fi
[[ -n "${GPU_TARGETS//[[:space:]]/}" ]] || die "could not determine CUDA target list"
read -r -a GPU_TARGET_ARRAY <<< "${GPU_TARGETS}"
for target in "${GPU_TARGET_ARRAY[@]}"; do
  [[ "${target}" =~ ^[a-zA-Z0-9_.-]+$ ]] \
    || die "invalid CUDA target name: ${target}"
done

echo "GPU device gate"
echo "  image=${IMAGE_URI}"
echo "  queue=${BATCH_QUEUE}"
echo "  jobDefinition=${BATCH_JOB_DEFINITION}"
echo "  gpuTargets=${#GPU_TARGET_ARRAY[@]}"

if aws ecr describe-repositories \
  --repository-names "${ECR_REPOSITORY_NAME}" \
  --region "${AWS_REGION}" >/dev/null 2>&1; then
  :
else
  die "ECR repository gutibm does not exist in ${AWS_REGION}"
fi

if aws ecr describe-images \
  --repository-name "${ECR_REPOSITORY_NAME}" \
  --image-ids "imageTag=${IMAGE_TAG}" \
  --region "${AWS_REGION}" \
  --query 'imageDetails[0].imageDigest' \
  --output text 2>/dev/null | grep -qv '^None$'; then
  die "refusing to overwrite existing ECR tag: ${IMAGE_TAG}"
fi

BUILD_PUSH_START="$(date +%s)"
echo "==> Build ${IMAGE_URI}"
aws ecr get-login-password --region "${AWS_REGION}" \
  | docker login --username AWS --password-stdin "${ECR_REGISTRY}"
docker build \
  --file "${ROOT}/deploy/aws/Dockerfile.gputest" \
  --build-arg "GPU_TARGETS=${GPU_TARGETS}" \
  --build-arg "GUTIBM_GIT_SHA=${GIT_SHA}" \
  --tag "${IMAGE_URI}" \
  "${ROOT}"
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
BUILD_PUSH_SECONDS=$(( $(date +%s) - BUILD_PUSH_START ))
IMAGE_REF="${ECR_REPOSITORY}@${IMAGE_DIGEST}"
echo "IMAGE_DIGEST=${IMAGE_DIGEST}"
echo "IMAGE_BUILD_PUSH_SECONDS=${BUILD_PUSH_SECONDS}"

JOB_DEF_JSON="$(mktemp)"
LOG_FILE="$(mktemp)"
jq -n \
  --arg name "${BATCH_JOB_DEFINITION}" \
  --arg image "${IMAGE_REF}" \
  '{
    jobDefinitionName: $name,
    type: "container",
    platformCapabilities: ["EC2"],
    containerProperties: {
      image: $image,
      command: ["-lc", "cd /src/build && nvidia-smi && ctest -L gpu --output-on-failure; echo GPUTEST_STATUS=$?"],
      vcpus: 4,
      memory: 12000,
      resourceRequirements: [{type: "GPU", value: "1"}],
      environment: [
        {name: "REQUIRE_GPU", value: "1"},
        {name: "OMPI_ALLOW_RUN_AS_ROOT", value: "1"},
        {name: "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM", value: "1"},
        {name: "HWLOC_COMPONENTS", value: "-x86"}
      ]
    },
    retryStrategy: {attempts: 1}
  }' > "${JOB_DEF_JSON}"

echo "==> Register digest-pinned job definition"
JOB_DEF_ARN="$(aws batch register-job-definition \
  --cli-input-json "file://${JOB_DEF_JSON}" \
  --region "${AWS_REGION}" \
  --query 'jobDefinitionArn' \
  --output text)"
[[ "${JOB_DEF_ARN}" == *":job-definition/${BATCH_JOB_DEFINITION}:"* ]] \
  || die "registered unexpected job definition: ${JOB_DEF_ARN}"
echo "JOB_DEFINITION_ARN=${JOB_DEF_ARN}"

JOB_NAME="gutibm-gpu-device-${GITHUB_RUN_ID:-manual-${GIT_SHORT_SHA}}"
SUBMIT_START="$(date +%s)"
echo "==> Submit device test job"
JOB_ID="$(aws batch submit-job \
  --job-name "${JOB_NAME}" \
  --job-queue "${BATCH_QUEUE}" \
  --job-definition "${JOB_DEF_ARN}" \
  --region "${AWS_REGION}" \
  --query jobId \
  --output text)"
[[ -n "${JOB_ID}" && "${JOB_ID}" != "None" ]] || die "Batch submission returned no job ID"
echo "JOB_ID=${JOB_ID}"

STATUS=""
while :; do
  STATUS="$(aws batch describe-jobs \
    --jobs "${JOB_ID}" \
    --region "${AWS_REGION}" \
    --query 'jobs[0].status' \
    --output text)"
  echo "  status=${STATUS}"
  case "${STATUS}" in
    SUCCEEDED|FAILED) break ;;
    SUBMITTED|PENDING|RUNNABLE|STARTING|RUNNING) sleep "${BATCH_POLL_SECONDS}" ;;
    *) die "unexpected Batch status for ${JOB_ID}: ${STATUS}" ;;
  esac
done

JOB_DETAILS="$(aws batch describe-jobs --jobs "${JOB_ID}" --region "${AWS_REGION}")"
CREATED_MS="$(jq -r '.jobs[0].createdAt // empty' <<< "${JOB_DETAILS}")"
STARTED_MS="$(jq -r '.jobs[0].startedAt // empty' <<< "${JOB_DETAILS}")"
STOPPED_MS="$(jq -r '.jobs[0].stoppedAt // empty' <<< "${JOB_DETAILS}")"
if [[ "${CREATED_MS}" =~ ^[0-9]+$ && "${STARTED_MS}" =~ ^[0-9]+$ &&
  "${STOPPED_MS}" =~ ^[0-9]+$ ]]; then
  BATCH_QUEUE_WAIT_SECONDS="$(awk -v ms="$((STARTED_MS - CREATED_MS))" \
    'BEGIN { printf "%.3f", ms / 1000 }')"
  BATCH_CONTAINER_SECONDS="$(awk -v ms="$((STOPPED_MS - STARTED_MS))" \
    'BEGIN { printf "%.3f", ms / 1000 }')"
  BATCH_JOB_SECONDS="$(awk -v ms="$((STOPPED_MS - CREATED_MS))" \
    'BEGIN { printf "%.3f", ms / 1000 }')"
else
  BATCH_QUEUE_WAIT_SECONDS="unknown"
  BATCH_CONTAINER_SECONDS="unknown"
  BATCH_JOB_SECONDS=$(( $(date +%s) - SUBMIT_START ))
fi
echo "BATCH_QUEUE_WAIT_SECONDS=${BATCH_QUEUE_WAIT_SECONDS}"
echo "BATCH_CONTAINER_SECONDS=${BATCH_CONTAINER_SECONDS}"
echo "BATCH_JOB_SECONDS=${BATCH_JOB_SECONDS}"
LOG_STREAM=""
for _ in $(seq 0 5 "${BATCH_LOG_WAIT_SECONDS}"); do
  LOG_STREAM="$(aws batch describe-jobs \
    --jobs "${JOB_ID}" \
    --region "${AWS_REGION}" \
    --query 'jobs[0].container.logStreamName' \
    --output text 2>/dev/null || true)"
  if [[ -n "${LOG_STREAM}" && "${LOG_STREAM}" != "None" ]]; then
    break
  fi
  sleep 5
done
[[ -n "${LOG_STREAM}" && "${LOG_STREAM}" != "None" ]] \
  || die "Batch job has no CloudWatch log stream: ${JOB_ID}"
echo "LOG_STREAM=${LOG_STREAM}"

echo "==> CloudWatch log"
NEXT_TOKEN=""
while :; do
  if [[ -n "${NEXT_TOKEN}" ]]; then
    EVENTS="$(aws logs get-log-events \
      --log-group-name "${BATCH_LOG_GROUP}" \
      --log-stream-name "${LOG_STREAM}" \
      --start-from-head \
      --next-token "${NEXT_TOKEN}" \
      --output json)"
  else
    EVENTS="$(aws logs get-log-events \
      --log-group-name "${BATCH_LOG_GROUP}" \
      --log-stream-name "${LOG_STREAM}" \
      --start-from-head \
      --output json)"
  fi
  jq -r '.events[]?.message // empty' <<< "${EVENTS}" \
    | tee -a "${LOG_FILE}"
  NEW_TOKEN="$(jq -r '.nextForwardToken // empty' <<< "${EVENTS}")"
  [[ -z "${NEW_TOKEN}" || "${NEW_TOKEN}" == "${NEXT_TOKEN}" ]] && break
  NEXT_TOKEN="${NEW_TOKEN}"
done

STATUS_LINE="$(grep -E 'GPUTEST_STATUS=[^[:space:]]+' "${LOG_FILE}" | tail -1 || true)"
if [[ -z "${STATUS_LINE}" ]]; then
  die "CloudWatch log has no parseable GPUTEST_STATUS"
fi
STATUS_VALUE="$(sed -n 's/.*GPUTEST_STATUS=\([^[:space:]]*\).*/\1/p' <<< "${STATUS_LINE}")"
if [[ "${STATUS_VALUE}" != "0" ]]; then
  die "GPU device test status was not zero: ${STATUS_LINE}"
fi
grep -q 'NVIDIA-SMI [0-9]' "${LOG_FILE}" \
  || die "CloudWatch log does not prove a physical NVIDIA device was present"
if [[ "${STATUS}" != "SUCCEEDED" ]]; then
  die "Batch job did not SUCCEED: ${STATUS}"
fi

END_TO_END_SECONDS=$(( $(date +%s) - BUILD_PUSH_START ))
echo "END_TO_END_SECONDS=${END_TO_END_SECONDS}"
echo "GPU device gate passed."
