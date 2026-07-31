#!/usr/bin/env bash
# AWS Batch entrypoint (draft) — docs/AWS_BATCH.md Phase 1–4 + Observability
#
# Expected environment (set in the Batch job definition / overrides):
#   INPUT_S3_URI     s3://bucket/path/input.json   (or set INPUT_S3_PREFIX)
#   INPUT_S3_PREFIX  s3://bucket/path/jobs         → ${PREFIX}/${INDEX}/input.json
#   OUTPUT_S3_URI    s3://bucket/path/output.h5.gz (or set OUTPUT_S3_PREFIX)
#   OUTPUT_S3_PREFIX s3://bucket/path/jobs         → ${PREFIX}/${INDEX}/output.h5.gz
#   CHECKPOINT_S3_PREFIX  s3://bucket/.../ckpt/    (optional; per-index if ARRAY)
#   CHECKPOINT_S3_URI     s3://bucket/.../step_N.h5  (optional explicit restart file)
#   CHECKPOINT_INTERVAL_SECONDS  how often to scan/upload closed restarts (default 300)
#   RESTART_INTERVAL_STEPS  injected into input.json when checkpointing (default 60)
#   REQUIRE_RESTART_GRID  1 (default) = closed restarts must contain /grid/
#   STATUS_S3_URI    optional explicit status.json URI (else under CHECKPOINT prefix)
#   GUTIBM_BINARY    path to gut_ibm (default /opt/gutibm/gut_ibm)
#   MPI_RANKS        default 1
#   EXTRA_MPIRUN_ARGS  optional extra flags
#   AWS_BATCH_JOB_ARRAY_INDEX  set automatically for array jobs
#   AWS_BATCH_JOB_ID           set automatically by Batch
#   MEMORY_GUARD             1 (default) = graceful stop when free RAM/VRAM is low
#   MEMORY_MIN_AVAILABLE_MB  host/cgroup free-RAM threshold to stop (default 2048)
#   GPU_MIN_FREE_MB          nvidia free-VRAM threshold to stop (default 512; 0=skip)
#   REQUIRE_GPU              1 = fail the job if the run log has no "GPU: ON" line
#                            (guards against a silent CUDA→CPU fallback; default 0).
#                            Fresh init() and init_from_checkpoint() both print the
#                            banner, so Spot resumes are covered.

set -euo pipefail

BINARY="${GUTIBM_BINARY:-/opt/gutibm/gut_ibm}"
MPI_RANKS="${MPI_RANKS:-1}"
WORK="${GUTIBM_WORK_DIR:-/tmp/gutibm_job}"
INDEX="${AWS_BATCH_JOB_ARRAY_INDEX:-}"
JOB_ID="${AWS_BATCH_JOB_ID:-}"
CHECKPOINT_INTERVAL_SECONDS="${CHECKPOINT_INTERVAL_SECONDS:-300}"
MEMORY_GUARD="${MEMORY_GUARD:-1}"
MEMORY_MIN_AVAILABLE_MB="${MEMORY_MIN_AVAILABLE_MB:-2048}"
GPU_MIN_FREE_MB="${GPU_MIN_FREE_MB:-512}"
REQUIRE_GPU="${REQUIRE_GPU:-0}"
REQUIRE_RESTART_GRID="${REQUIRE_RESTART_GRID:-1}"
RESTART_INTERVAL_STEPS="${RESTART_INTERVAL_STEPS:-60}"
CHECKPOINT_NAME="checkpoint.h5"
LATEST_NAME="latest.json"
STATUS_NAME="status.json"
STATUS_FAILED="failed"
STATUS_MEMORY_PRESSURE="memory_pressure"
RESTART_DIR="${WORK}/restart"
SYNC_PID=""
MPIRUN_PID=""
SPOT_NOTICED=0
MEMORY_PRESSURE=0
RESUME_FROM_CHECKPOINT=0
RUN_EXIT_CODE=0
WALL_START_EPOCH="$(date -u +%s)"
MEM_AVAILABLE_MB=""
MEM_CGROUP_FREE_MB=""
MEM_EFFECTIVE_FREE_MB=""
GPU_FREE_MB=""
CHECKPOINT_KEY=""
mkdir -p "${WORK}" "${RESTART_DIR}"
RUN_LOG="${WORK}/run.log"
: >"${RUN_LOG}"

iso_now() { date -u +"%Y-%m-%dT%H:%M:%SZ"; }

resolve_status_uri() {
  if [[ -n "${STATUS_S3_URI:-}" ]]; then
    echo "${STATUS_S3_URI}"
    return 0
  fi
  if [[ -n "${CHECKPOINT_S3_PREFIX:-}" ]]; then
    echo "${CHECKPOINT_S3_PREFIX%/}/${STATUS_NAME}"
    return 0
  fi
  if [[ -n "${OUTPUT_S3_URI:-}" ]]; then
    # Same directory as the gzipped output.
    echo "${OUTPUT_S3_URI%/*}/${STATUS_NAME}"
    return 0
  fi
  echo ""
}

# True if path is a readable closed restart (agents required; grid when Tier 2).
hdf5_checkpoint_usable() {
  local f="$1"
  [[ -f "${f}" ]] || return 1
  local sz
  sz="$(stat -c '%s' "${f}" 2>/dev/null || echo 0)"
  if ! [[ "${sz}" =~ ^[0-9]+$ ]] || (( sz < 4096 )); then
    echo "Restart unusable (${f}): size=${sz} (<4096)" >&2
    return 1
  fi
  if command -v h5ls >/dev/null 2>&1; then
    local listing err rc=0
    err="$(mktemp)"
    # Path probes are more reliable than grepping `h5ls -r` (and surface open errors).
    if ! h5ls "${f}/agents" >/dev/null 2>"${err}"; then
      echo "Restart unusable (${f}): missing/unreadable /agents (size=${sz}) $(tr '\n' ' ' <"${err}" | head -c 180)" >&2
      rm -f "${err}"
      return 1
    fi
    if [[ "${REQUIRE_RESTART_GRID}" == "1" ]]; then
      if ! h5ls "${f}/grid" >/dev/null 2>"${err}"; then
        echo "Restart unusable (${f}): missing/unreadable /grid (size=${sz}) $(tr '\n' ' ' <"${err}" | head -c 180)" >&2
        rm -f "${err}"
        return 1
      fi
    fi
    # Full listing only for the agents-path regex seatbelt used by older probes.
    listing="$(h5ls -r "${f}" 2>"${err}")" || rc=$?
    if (( rc != 0 )); then
      echo "Restart unusable (${f}): h5ls -r failed rc=${rc} size=${sz} $(tr '\n' ' ' <"${err}" | head -c 180)" >&2
      rm -f "${err}"
      return 1
    fi
    rm -f "${err}"
    echo "${listing}" | grep -qE '^/agents/' || {
      echo "Restart unusable (${f}): no /agents/ child groups (size=${sz})" >&2
      return 1
    }
    if [[ "${REQUIRE_RESTART_GRID}" == "1" ]]; then
      echo "${listing}" | grep -qE '^/grid/' || {
        echo "Restart unusable (${f}): no /grid/ child groups (size=${sz})" >&2
        return 1
      }
    fi
  else
    # Weak fallback when hdf5-tools is missing: signature only.
    local magic
    magic="$(od -An -N8 -tx1 "${f}" 2>/dev/null | tr -d ' \n')"
    [[ "${magic}" == "894844460d0a1a0a" ]] || {
      echo "Restart unusable (${f}): bad HDF5 signature (no h5ls)" >&2
      return 1
    }
  fi
  return 0
}

file_sha256() {
  local f="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${f}" | awk '{print $1}'
  else
    echo ""
  fi
}

write_latest_json() {
  local step_name="$1"
  local s3_uri="$2"
  local local_path="$3"
  local size sha
  size="$(stat -c '%s' "${local_path}" 2>/dev/null || echo 0)"
  sha="$(file_sha256 "${local_path}")"
  GUTIBM_WORK_DIR="${WORK}" \
  GUTIBM_LATEST_STEP="${step_name}" \
  GUTIBM_LATEST_URI="${s3_uri}" \
  GUTIBM_LATEST_SIZE="${size}" \
  GUTIBM_LATEST_SHA="${sha}" \
  GUTIBM_JOB_ID="${JOB_ID}" \
  python3 - <<'PY'
import json, os
from pathlib import Path
work = Path(os.environ["GUTIBM_WORK_DIR"])
payload = {
    "step": os.environ.get("GUTIBM_LATEST_STEP"),
    "uri": os.environ.get("GUTIBM_LATEST_URI"),
    "size_bytes": int(os.environ.get("GUTIBM_LATEST_SIZE") or 0),
    "sha256": os.environ.get("GUTIBM_LATEST_SHA") or None,
    "job_id": os.environ.get("GUTIBM_JOB_ID") or None,
    "updated_at": __import__("datetime").datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
    "fidelity": "tier2_agents_grid",
    "rng": "reseeded_from_config_seed",
}
(work / "latest.json").write_text(json.dumps(payload, indent=2) + "\n")
PY
}

# Upload closed restart/step_*.h5 files immutably; update latest.json pointer only.
upload_checkpoint() {
  [[ -n "${CHECKPOINT_S3_PREFIX:-}" ]] || return 0
  [[ -d "${RESTART_DIR}" ]] || return 0
  local f base dest uploading step_name newest="" newest_uri=""
  shopt -s nullglob
  for f in "${RESTART_DIR}"/step_*.h5; do
    base="$(basename "${f}")"
    dest="${CHECKPOINT_S3_PREFIX}${base}"
    if [[ -f "${f}.uploaded" ]]; then
      newest="${f}"
      newest_uri="${dest}"
      continue
    fi
    if [[ -f "${f}.rejected" ]]; then
      continue
    fi
    if ! hdf5_checkpoint_usable "${f}"; then
      echo "Skipping unusable restart artifact: ${f}" >&2
      # Stop re-probing every sync interval; keep the file for the instance lifetime.
      : >"${f}.rejected"
      continue
    fi
    # Immutable: never overwrite an existing remote step object.
    if aws s3 ls "${dest}" >/dev/null 2>&1; then
      echo "Remote restart already exists (immutable): ${dest}"
      : >"${f}.uploaded"
      newest="${f}"
      newest_uri="${dest}"
      continue
    fi
    uploading="${dest}.uploading"
    if aws s3 cp "${f}" "${uploading}" >/dev/null 2>&1 \
      && aws s3 cp "${uploading}" "${dest}" >/dev/null 2>&1; then
      aws s3 rm "${uploading}" >/dev/null 2>&1 || true
      : >"${f}.uploaded"
      newest="${f}"
      newest_uri="${dest}"
      CHECKPOINT_UPLOADED_AT="$(iso_now)"
      CHECKPOINT_KEY="${base}"
      export CHECKPOINT_UPLOADED_AT CHECKPOINT_KEY
      echo "Uploaded immutable restart: ${dest}"
      # Free scratch space: drop older uploaded local steps (keep the newest).
      local old
      for old in "${RESTART_DIR}"/step_*.h5; do
        [[ -f "${old}" ]] || continue
        [[ "${old}" == "${f}" ]] && continue
        if [[ -f "${old}.uploaded" ]]; then
          rm -f "${old}" "${old}.uploaded" "${old}.rejected" || true
        fi
      done
    else
      echo "Restart upload failed for ${f} (left prior S3 objects untouched)" >&2
      aws s3 rm "${uploading}" >/dev/null 2>&1 || true
    fi
  done
  shopt -u nullglob

  if [[ -n "${newest}" && -n "${newest_uri}" ]]; then
    step_name="$(basename "${newest}" .h5)"
    write_latest_json "${step_name}" "${newest_uri}" "${newest}"
    aws s3 cp "${WORK}/latest.json" "${CHECKPOINT_S3_PREFIX}${LATEST_NAME}" >/dev/null 2>&1 || true
  fi
}

# Parse the latest gut_ibm progress line from RUN_LOG into env vars for status.json.
# Tolerant: missing fields leave defaults.
parse_progress_from_log() {
  SIM_TIME_S=""
  TOTAL_TIME_S=""
  PCT=""
  STEP=""
  GLOBAL_AGENTS=""
  MU_AVG=""
  RATE=""
  ETA_S=""
  local line
  line="$(grep -E '^Step ' "${RUN_LOG}" 2>/dev/null | tail -n 1 || true)"
  [[ -n "${line}" ]] || return 0
  # Example:
  # Step 10  t=600s  dt=60s  global_agents=50  local_agents=50  mu_avg=5e-4  pct=10  rate=1.2  eta_s=4500
  if [[ "${line}" =~ Step[[:space:]]+([0-9]+) ]]; then
    STEP="${BASH_REMATCH[1]}"
  fi
  if [[ "${line}" =~ t=([0-9.eE+-]+)s ]]; then
    SIM_TIME_S="${BASH_REMATCH[1]}"
  fi
  if [[ "${line}" =~ global_agents=([0-9]+) ]]; then
    GLOBAL_AGENTS="${BASH_REMATCH[1]}"
  fi
  if [[ "${line}" =~ mu_avg=([0-9.eE+-]+) ]]; then
    MU_AVG="${BASH_REMATCH[1]}"
  fi
  if [[ "${line}" =~ pct=([0-9.eE+-]+) ]]; then
    PCT="${BASH_REMATCH[1]}"
  fi
  if [[ "${line}" =~ rate=([0-9.eE+-]+) ]]; then
    RATE="${BASH_REMATCH[1]}"
  fi
  if [[ "${line}" =~ eta_s=([0-9.eE+-]+) ]]; then
    ETA_S="${BASH_REMATCH[1]}"
  fi
  if [[ -f "${WORK}/input.json" ]]; then
    TOTAL_TIME_S="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("total_time",""))' "${WORK}/input.json" 2>/dev/null || true)"
  fi
}

# Sample host MemAvailable, cgroup free (Batch memory limit), and GPU free VRAM.
sample_memory() {
  MEM_AVAILABLE_MB=""
  MEM_CGROUP_FREE_MB=""
  MEM_EFFECTIVE_FREE_MB=""
  GPU_FREE_MB=""
  if [[ -r /proc/meminfo ]]; then
    local avail_kb
    avail_kb="$(awk '/^MemAvailable:/ {print $2; exit}' /proc/meminfo 2>/dev/null || true)"
    if [[ -n "${avail_kb}" ]]; then
      MEM_AVAILABLE_MB=$((avail_kb / 1024))
    fi
  fi
  # cgroup v2 (ECS/Batch): memory.max / memory.current are bytes.
  local cg_max="" cg_cur=""
  if [[ -r /sys/fs/cgroup/memory.max ]]; then
    cg_max="$(cat /sys/fs/cgroup/memory.max 2>/dev/null || true)"
    cg_cur="$(cat /sys/fs/cgroup/memory.current 2>/dev/null || true)"
  elif [[ -r /sys/fs/cgroup/memory/memory.limit_in_bytes ]]; then
    cg_max="$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || true)"
    cg_cur="$(cat /sys/fs/cgroup/memory/memory.usage_in_bytes 2>/dev/null || true)"
  fi
  if [[ -n "${cg_max}" && "${cg_max}" != "max" && "${cg_max}" =~ ^[0-9]+$ \
        && -n "${cg_cur}" && "${cg_cur}" =~ ^[0-9]+$ ]] \
     && (( cg_max < 1000000000000 )); then
    # Ignore absurdly large "unlimited" legacy limits (often ~2^63).
    local free_b=$((cg_max - cg_cur))
    if (( free_b < 0 )); then free_b=0; fi
    MEM_CGROUP_FREE_MB=$((free_b / 1024 / 1024))
  fi
  if [[ -n "${MEM_AVAILABLE_MB}" && -n "${MEM_CGROUP_FREE_MB}" ]]; then
    if (( MEM_CGROUP_FREE_MB < MEM_AVAILABLE_MB )); then
      MEM_EFFECTIVE_FREE_MB="${MEM_CGROUP_FREE_MB}"
    else
      MEM_EFFECTIVE_FREE_MB="${MEM_AVAILABLE_MB}"
    fi
  elif [[ -n "${MEM_CGROUP_FREE_MB}" ]]; then
    MEM_EFFECTIVE_FREE_MB="${MEM_CGROUP_FREE_MB}"
  else
    MEM_EFFECTIVE_FREE_MB="${MEM_AVAILABLE_MB}"
  fi
  if command -v nvidia-smi >/dev/null 2>&1; then
    local gpu_raw
    gpu_raw="$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null \
      | head -n 1 | tr -d '[:space:]' || true)"
    if [[ "${gpu_raw}" =~ ^[0-9]+$ ]]; then
      GPU_FREE_MB="${gpu_raw}"
    fi
  fi
  export MEM_AVAILABLE_MB MEM_CGROUP_FREE_MB MEM_EFFECTIVE_FREE_MB GPU_FREE_MB
}

# Returns 0 when free host/cgroup RAM or GPU VRAM is below configured floors.
check_memory_pressure() {
  [[ "${MEMORY_GUARD}" == "1" ]] || return 1
  sample_memory
  if [[ -n "${MEM_EFFECTIVE_FREE_MB}" && "${MEM_EFFECTIVE_FREE_MB}" =~ ^[0-9]+$ ]] \
     && (( MEM_EFFECTIVE_FREE_MB < MEMORY_MIN_AVAILABLE_MB )); then
    echo "Memory pressure: effective free ${MEM_EFFECTIVE_FREE_MB} MiB < ${MEMORY_MIN_AVAILABLE_MB} MiB" >&2
    return 0
  fi
  if [[ "${GPU_MIN_FREE_MB}" != "0" && -n "${GPU_FREE_MB}" && "${GPU_FREE_MB}" =~ ^[0-9]+$ ]] \
     && (( GPU_FREE_MB < GPU_MIN_FREE_MB )); then
    echo "GPU memory pressure: free ${GPU_FREE_MB} MiB < ${GPU_MIN_FREE_MB} MiB" >&2
    return 0
  fi
  return 1
}

write_status_json() {
  local state="$1"
  sample_memory
  parse_progress_from_log
  local wall_now
  wall_now="$(date -u +%s)"
  local wall_elapsed=$((wall_now - WALL_START_EPOCH))
  local spot_action="${SPOT_ACTION_TIME:-}"
  GUTIBM_WORK_DIR="${WORK}" \
  GUTIBM_STATUS_STATE="${state}" \
  GUTIBM_JOB_ID="${JOB_ID}" \
  GUTIBM_ARRAY_INDEX="${INDEX}" \
  GUTIBM_SIM_TIME_S="${SIM_TIME_S}" \
  GUTIBM_TOTAL_TIME_S="${TOTAL_TIME_S}" \
  GUTIBM_PCT="${PCT}" \
  GUTIBM_STEP="${STEP}" \
  GUTIBM_GLOBAL_AGENTS="${GLOBAL_AGENTS}" \
  GUTIBM_MU_AVG="${MU_AVG}" \
  GUTIBM_RATE="${RATE}" \
  GUTIBM_ETA_S="${ETA_S}" \
  GUTIBM_WALL_ELAPSED_S="${wall_elapsed}" \
  GUTIBM_CHECKPOINT_UPLOADED_AT="${CHECKPOINT_UPLOADED_AT:-}" \
  GUTIBM_CHECKPOINT_KEY="${CHECKPOINT_KEY:-}" \
  GUTIBM_RESUME_FROM_CHECKPOINT="${RESUME_FROM_CHECKPOINT}" \
  GUTIBM_SPOT_INTERRUPTION="${SPOT_NOTICED}" \
  GUTIBM_SPOT_ACTION_TIME="${spot_action}" \
  GUTIBM_MEMORY_PRESSURE="${MEMORY_PRESSURE}" \
  GUTIBM_MEM_AVAILABLE_MB="${MEM_AVAILABLE_MB}" \
  GUTIBM_MEM_CGROUP_FREE_MB="${MEM_CGROUP_FREE_MB}" \
  GUTIBM_MEM_EFFECTIVE_FREE_MB="${MEM_EFFECTIVE_FREE_MB}" \
  GUTIBM_GPU_FREE_MB="${GPU_FREE_MB}" \
  GUTIBM_MEMORY_MIN_AVAILABLE_MB="${MEMORY_MIN_AVAILABLE_MB}" \
  GUTIBM_GPU_MIN_FREE_MB="${GPU_MIN_FREE_MB}" \
  GUTIBM_EXIT_CODE="${RUN_EXIT_CODE}" \
  python3 - <<'PY'
import json, os
from pathlib import Path

def num_or_none(raw):
    if raw is None or raw == "":
        return None
    try:
        if "." in raw or "e" in raw.lower() or "E" in raw:
            return float(raw)
        return int(raw)
    except ValueError:
        try:
            return float(raw)
        except ValueError:
            return None

work = Path(os.environ["GUTIBM_WORK_DIR"])
spot = os.environ.get("GUTIBM_SPOT_INTERRUPTION", "0") == "1"
mem_pressure = os.environ.get("GUTIBM_MEMORY_PRESSURE", "0") == "1"
payload = {
    "job_id": os.environ.get("GUTIBM_JOB_ID") or None,
    "array_index": os.environ.get("GUTIBM_ARRAY_INDEX") or None,
    "state": os.environ["GUTIBM_STATUS_STATE"],
    "sim_time_s": num_or_none(os.environ.get("GUTIBM_SIM_TIME_S", "")),
    "total_time_s": num_or_none(os.environ.get("GUTIBM_TOTAL_TIME_S", "")),
    "pct": num_or_none(os.environ.get("GUTIBM_PCT", "")),
    "step": num_or_none(os.environ.get("GUTIBM_STEP", "")),
    "global_agents": num_or_none(os.environ.get("GUTIBM_GLOBAL_AGENTS", "")),
    "mu_avg": num_or_none(os.environ.get("GUTIBM_MU_AVG", "")),
    "rate": num_or_none(os.environ.get("GUTIBM_RATE", "")),
    "wall_elapsed_s": num_or_none(os.environ.get("GUTIBM_WALL_ELAPSED_S", "")),
    "eta_s": num_or_none(os.environ.get("GUTIBM_ETA_S", "")),
    "checkpoint_uploaded_at": os.environ.get("GUTIBM_CHECKPOINT_UPLOADED_AT") or None,
    "checkpoint_key": os.environ.get("GUTIBM_CHECKPOINT_KEY") or None,
    "resume_from_checkpoint": os.environ.get("GUTIBM_RESUME_FROM_CHECKPOINT", "0") == "1",
    "spot_interruption": spot,
    "spot_action_time": os.environ.get("GUTIBM_SPOT_ACTION_TIME") or None,
    "memory_pressure": mem_pressure,
    "mem_available_mb": num_or_none(os.environ.get("GUTIBM_MEM_AVAILABLE_MB", "")),
    "mem_cgroup_free_mb": num_or_none(os.environ.get("GUTIBM_MEM_CGROUP_FREE_MB", "")),
    "mem_effective_free_mb": num_or_none(os.environ.get("GUTIBM_MEM_EFFECTIVE_FREE_MB", "")),
    "gpu_free_mb": num_or_none(os.environ.get("GUTIBM_GPU_FREE_MB", "")),
    "memory_min_available_mb": num_or_none(os.environ.get("GUTIBM_MEMORY_MIN_AVAILABLE_MB", "")),
    "gpu_min_free_mb": num_or_none(os.environ.get("GUTIBM_GPU_MIN_FREE_MB", "")),
    "exit_code": num_or_none(os.environ.get("GUTIBM_EXIT_CODE", "")),
    "updated_at": __import__("datetime").datetime.now(
        __import__("datetime").timezone.utc
    ).strftime("%Y-%m-%dT%H:%M:%SZ"),
}
(work / "status.json").write_text(json.dumps(payload, indent=2) + "\n")
PY
}

upload_status() {
  local uri
  uri="$(resolve_status_uri)"
  [[ -n "${uri}" && -f "${WORK}/${STATUS_NAME}" ]] || return 0
  aws s3 cp "${WORK}/${STATUS_NAME}" "${uri}" >/dev/null 2>&1 || true
}

# IMDSv2 Spot interruption notice (EC2). Returns 0 if a reclaim is pending.
check_spot_interruption() {
  local token
  token="$(curl -s -m 2 -X PUT "http://169.254.169.254/latest/api/token" \
    -H "X-aws-ec2-metadata-token-ttl-seconds: 60" 2>/dev/null || true)"
  [[ -n "${token}" ]] || return 1
  local action
  action="$(curl -s -m 2 -H "X-aws-ec2-metadata-token: ${token}" \
    "http://169.254.169.254/latest/meta-data/spot/instance-action" 2>/dev/null || true)"
  if [[ -n "${action}" && "${action}" != *"404"* && "${action}" != *"Not Found"* ]]; then
    SPOT_ACTION_TIME="$(iso_now)"
    export SPOT_ACTION_TIME
    echo "Spot interruption notice: ${action}" >&2
    return 0
  fi
  return 1
}

request_graceful_stop() {
  if [[ -n "${MPIRUN_PID}" ]] && kill -0 "${MPIRUN_PID}" 2>/dev/null; then
    echo "Sending SIGTERM to mpirun pid=${MPIRUN_PID}" >&2
    kill -TERM "${MPIRUN_PID}" 2>/dev/null || true
  fi
}

checkpoint_sync_loop() {
  while true; do
    sleep "${CHECKPOINT_INTERVAL_SECONDS}"
    if check_spot_interruption; then
      SPOT_NOTICED=1
      write_status_json "spot_notice"
      upload_status
      request_graceful_stop
      upload_checkpoint
      write_status_json "spot_notice"
      upload_status
      return 0
    fi
    if check_memory_pressure; then
      MEMORY_PRESSURE=1
      write_status_json "memory_pressure"
      upload_status
      request_graceful_stop
      upload_checkpoint
      write_status_json "memory_pressure"
      upload_status
      return 0
    fi
    upload_checkpoint
    write_status_json "running"
    upload_status
  done
}

stop_checkpoint_sync() {
  if [[ -n "${SYNC_PID}" ]]; then
    kill "${SYNC_PID}" 2>/dev/null || true
    wait "${SYNC_PID}" 2>/dev/null || true
    SYNC_PID=""
  fi
}

on_signal() {
  echo "Entrypoint received stop signal; flushing checkpoint + status" >&2
  request_graceful_stop
  # Give gut_ibm a moment to finish the current step and finalize HDF5.
  if [[ -n "${MPIRUN_PID}" ]]; then
    wait "${MPIRUN_PID}" 2>/dev/null || RUN_EXIT_CODE=$?
  fi
  stop_checkpoint_sync
  upload_checkpoint
  write_status_json "stopping"
  upload_status
  exit "${RUN_EXIT_CODE:-143}"
}

trap on_signal TERM INT
trap stop_checkpoint_sync EXIT

if [[ -z "${INPUT_S3_URI:-}" && -n "${INPUT_S3_PREFIX:-}" && -n "${INDEX}" ]]; then
  INPUT_S3_URI="${INPUT_S3_PREFIX%/}/${INDEX}/input.json"
fi
if [[ -z "${OUTPUT_S3_URI:-}" && -n "${OUTPUT_S3_PREFIX:-}" && -n "${INDEX}" ]]; then
  OUTPUT_S3_URI="${OUTPUT_S3_PREFIX%/}/${INDEX}/output.h5.gz"
fi
if [[ -n "${CHECKPOINT_S3_PREFIX:-}" && -n "${INDEX}" && "${CHECKPOINT_S3_PREFIX}" != */"${INDEX}"/* ]]; then
  CHECKPOINT_S3_PREFIX="${CHECKPOINT_S3_PREFIX%/}/${INDEX}/"
fi

if [[ -z "${INPUT_S3_URI:-}" || -z "${OUTPUT_S3_URI:-}" ]]; then
  echo "Need INPUT_S3_URI/OUTPUT_S3_URI, or PREFIX vars plus AWS_BATCH_JOB_ARRAY_INDEX" >&2
  exit 2
fi

echo "Downloading input: ${INPUT_S3_URI}"
aws s3 cp "${INPUT_S3_URI}" "${WORK}/input.json"

resolve_resume_uri() {
  if [[ -n "${CHECKPOINT_S3_URI:-}" ]]; then
    echo "${CHECKPOINT_S3_URI}"
    return 0
  fi
  [[ -n "${CHECKPOINT_S3_PREFIX:-}" ]] || return 1
  if aws s3 cp "${CHECKPOINT_S3_PREFIX}${LATEST_NAME}" "${WORK}/latest.json" >/dev/null 2>&1; then
    local uri
    uri="$(GUTIBM_WORK_DIR="${WORK}" python3 - <<'PY'
import json, os
from pathlib import Path
work = Path(os.environ["GUTIBM_WORK_DIR"])
data = json.loads((work / "latest.json").read_text())
print(data.get("uri") or "")
PY
)"
    if [[ -n "${uri}" ]]; then
      echo "${uri}"
      return 0
    fi
  fi
  # Legacy fallback: mutable checkpoint.h5 from older entrypoint revisions.
  if aws s3 ls "${CHECKPOINT_S3_PREFIX}${CHECKPOINT_NAME}" >/dev/null 2>&1; then
    echo "${CHECKPOINT_S3_PREFIX}${CHECKPOINT_NAME}"
    return 0
  fi
  return 1
}

RESUME_URI=""
if RESUME_URI="$(resolve_resume_uri)" && [[ -n "${RESUME_URI}" ]]; then
  echo "Resuming from restart artifact: ${RESUME_URI}"
  aws s3 cp "${RESUME_URI}" "${WORK}/checkpoint.h5"
  if ! hdf5_checkpoint_usable "${WORK}/checkpoint.h5"; then
    echo "Downloaded restart is unreadable/incomplete; refusing resume" >&2
    echo "Quarantining S3 object so retries do not loop on a corrupt file" >&2
    aws s3 mv \
      "${RESUME_URI}" \
      "${RESUME_URI}.corrupt.$(date -u +%Y%m%dT%H%M%SZ)" \
      >/dev/null 2>&1 || true
    rm -f "${WORK}/checkpoint.h5"
    write_status_json "${STATUS_FAILED}"
    upload_status
    exit 3
  fi
  RESUME_FROM_CHECKPOINT=1
  CHECKPOINT_KEY="$(basename "${RESUME_URI}")"
  export CHECKPOINT_KEY
  GUTIBM_WORK_DIR="${WORK}" python3 - <<'PY'
import json, os
from pathlib import Path
work = Path(os.environ["GUTIBM_WORK_DIR"])
cfg = json.loads((work / "input.json").read_text())
cfg["checkpoint_file"] = str(work / "checkpoint.h5")
cfg.setdefault("checkpoint_step", "")
(work / "input.json").write_text(json.dumps(cfg, indent=2) + "\n")
PY
fi

ENABLE_RESTART=0
if [[ -n "${CHECKPOINT_S3_PREFIX:-}" ]]; then
  ENABLE_RESTART=1
fi
GUTIBM_WORK_DIR="${WORK}" \
RESTART_INTERVAL_STEPS="${RESTART_INTERVAL_STEPS}" \
ENABLE_RESTART="${ENABLE_RESTART}" \
python3 - <<'PY'
import json, os
from pathlib import Path
work = Path(os.environ["GUTIBM_WORK_DIR"])
cfg = json.loads((work / "input.json").read_text())
cfg["hdf5_file"] = str(work / "output.h5")
cfg.setdefault("gpu_enabled", True)
if os.environ.get("ENABLE_RESTART") == "1":
    restart = cfg.get("restart")
    if not isinstance(restart, dict):
        restart = {}
    restart["enabled"] = True
    restart["directory"] = str(work / "restart")
    restart.setdefault(
        "interval_steps",
        int(os.environ.get("RESTART_INTERVAL_STEPS", "60")),
    )
    cfg["restart"] = restart
(work / "input.json").write_text(json.dumps(cfg, indent=2) + "\n")
PY

write_status_json "starting"
upload_status

# Fail fast if the container is already below the free-RAM floor (wrong instance /
# job-def memory for this grid) before burning GPU time.
if check_memory_pressure; then
  MEMORY_PRESSURE=1
  write_status_json "memory_pressure"
  upload_status
  echo "Refusing to start: effective free RAM/VRAM below guard thresholds" >&2
  exit 137
fi

if [[ -n "${CHECKPOINT_S3_PREFIX:-}" ]] || [[ -n "$(resolve_status_uri)" ]]; then
  echo "Status/checkpoint sync every ${CHECKPOINT_INTERVAL_SECONDS}s (memory guard=${MEMORY_GUARD})"
  checkpoint_sync_loop &
  SYNC_PID=$!
fi

echo "Running ${BINARY} (mpi_ranks=${MPI_RANKS})"
# shellcheck disable=SC2086
set +e
mpirun -np "${MPI_RANKS}" ${EXTRA_MPIRUN_ARGS:-} "${BINARY}" "${WORK}/input.json" \
  > >(tee -a "${RUN_LOG}") 2> >(tee -a "${RUN_LOG}" >&2) &
MPIRUN_PID=$!
wait "${MPIRUN_PID}"
RUN_EXIT_CODE=$?
set -e
MPIRUN_PID=""

stop_checkpoint_sync
# Final checkpoint push so a retry after a clean-but-incomplete exit can resume.
upload_checkpoint

if [[ "${SPOT_NOTICED}" -eq 1 ]]; then
  write_status_json "spot_notice"
  upload_status
  # Non-zero so Batch retryStrategy can resume from checkpoint (even if gut_ibm
  # exited 0 after a graceful SIGTERM drain).
  local_exit="${RUN_EXIT_CODE}"
  if [[ "${local_exit}" -eq 0 ]]; then
    local_exit=1
  fi
  echo "Exiting after Spot notice (exit=${local_exit}) so Batch can retry" >&2
  exit "${local_exit}"
fi

if [[ "${MEMORY_PRESSURE}" -eq 1 ]]; then
  write_status_json "memory_pressure"
  upload_status
  # Do NOT rely on Batch auto-retry here — same instance size will OOM again.
  # Checkpoint is on S3; resume manually on a larger CE / higher memory job def.
  local_exit="${RUN_EXIT_CODE}"
  if [[ "${local_exit}" -eq 0 ]]; then
    local_exit=137
  fi
  echo "Exiting after memory pressure (exit=${local_exit}); resize instance before resume" >&2
  exit "${local_exit}"
fi

if [[ "${RUN_EXIT_CODE}" -ne 0 ]]; then
  write_status_json "${STATUS_FAILED}"
  upload_status
  echo "gut_ibm exited with ${RUN_EXIT_CODE}" >&2
  exit "${RUN_EXIT_CODE}"
fi

# Guard against a silent CUDA→CPU fallback: gut_ibm exits 0 and writes HDF5 even
# when GPU init fails, so a green Batch job does not by itself prove the GPU path
# ran. When REQUIRE_GPU=1, insist on the "GPU: ON" banner (printed by both fresh
# init() and init_from_checkpoint() resume paths).
if [[ "${REQUIRE_GPU}" == "1" ]]; then
  if grep -qE 'GPU: ON' "${RUN_LOG}"; then
    echo "GPU activation confirmed (REQUIRE_GPU=1)"
  else
    write_status_json "${STATUS_FAILED}"
    upload_status
    echo "REQUIRE_GPU=1 but no 'GPU: ON' line in run log — CUDA fell back to CPU" >&2
    exit 42
  fi
fi

if [[ -f "${WORK}/output.h5" ]]; then
  gzip -c "${WORK}/output.h5" > "${WORK}/output.h5.gz"
  echo "Uploading ${OUTPUT_S3_URI}"
  aws s3 cp "${WORK}/output.h5.gz" "${OUTPUT_S3_URI}"
else
  write_status_json "${STATUS_FAILED}"
  upload_status
  echo "No output.h5 produced" >&2
  exit 1
fi

write_status_json "succeeded"
upload_status
echo "Done."
