#!/usr/bin/env bash
# Every registered GPU job definition must carry a bounded attempt timeout,
# a GPU resource requirement, and an explicit /dev/shm allocation, and every
# GPU compute environment must scale to zero.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
AWS_DIR="${ROOT}/deploy/aws"

fail() {
  local message="$1"
  echo "${message}" >&2
  exit 1
}

# Timeout defaults exist and are positive, bounded, whole seconds.
for var in PRACTICE_JOB_TIMEOUT_SECONDS CAMPAIGN_JOB_TIMEOUT_SECONDS; do
  value="$(sed -nE "s/^export ${var}=\"\\\$\{${var}:-([0-9]+)\}\"$/\1/p" \
    "${AWS_DIR}/env.sh")"
  [[ -n "${value}" ]] || fail "env.sh does not default ${var}"
  (( value > 0 )) || fail "${var} default is not positive: ${value}"
  (( value <= 172800 )) || fail "${var} default exceeds 48 h: ${value}"
  grep -Fq "${var}" "${AWS_DIR}/env.ps1" \
    || fail "env.ps1 does not default ${var}"
done

check_stack() {
  local script="$1" timeout_var="$2"
  local path="${AWS_DIR}/${script}"

  grep -Fq -- "--timeout" "${path}" \
    || fail "${script} registers a job definition without --timeout"
  grep -Fq "attemptDurationSeconds=" "${path}" \
    || fail "${script} does not set attemptDurationSeconds"
  grep -Fq "${timeout_var}" "${path}" \
    || fail "${script} does not use ${timeout_var}"
  grep -Fq '"sharedMemorySize"' "${path}" \
    || fail "${script} does not allocate /dev/shm"
  grep -Fq '"type": "GPU"' "${path}" \
    || fail "${script} does not request a GPU"
  grep -Fq -- "--platform-capabilities EC2" "${path}" \
    || fail "${script} does not target EC2"
  grep -Fq '"minvCpus": 0' "${path}" \
    || fail "${script} compute environment does not scale to zero"
  grep -Fq '"desiredvCpus": 0' "${path}" \
    || fail "${script} compute environment does not start at zero vCPUs"
  grep -Fq 'ECS_AL2023_NVIDIA' "${path}" \
    || fail "${script} does not pin an NVIDIA-driver AMI"
}

check_stack 02_setup_practice_stack.sh PRACTICE_JOB_TIMEOUT_SECONDS
check_stack 05_setup_campaign_stack.sh CAMPAIGN_JOB_TIMEOUT_SECONDS

# PowerShell mirrors must not drift from the Bash stacks.
grep -Fq "PRACTICE_JOB_TIMEOUT_SECONDS" "${AWS_DIR}/02_setup_practice_stack.ps1" \
  || fail "02_setup_practice_stack.ps1 has no attempt timeout"
grep -Fq "CAMPAIGN_JOB_TIMEOUT_SECONDS" "${AWS_DIR}/05_setup_campaign_stack.ps1" \
  || fail "05_setup_campaign_stack.ps1 has no attempt timeout"

echo "PASS: AWS GPU job definitions are timeout-bounded and scale to zero"
