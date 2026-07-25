#!/usr/bin/env bash
# AWS Batch entrypoint (draft) — docs/AWS_BATCH.md Phase 1–4 + Observability
#
# Expected environment (set in the Batch job definition / overrides):
#   INPUT_S3_URI     s3://bucket/path/input.json   (or set INPUT_S3_PREFIX)
#   INPUT_S3_PREFIX  s3://bucket/path/jobs         → ${PREFIX}/${INDEX}/input.json
#   OUTPUT_S3_URI    s3://bucket/path/output.h5.gz (or set OUTPUT_S3_PREFIX)
#   OUTPUT_S3_PREFIX s3://bucket/path/jobs         → ${PREFIX}/${INDEX}/output.h5.gz
#   CHECKPOINT_S3_PREFIX  s3://bucket/.../ckpt/    (optional; per-index if ARRAY)
#   CHECKPOINT_INTERVAL_SECONDS  checkpoint upload cadence (default 300)
#   STATUS_S3_URI    optional explicit status.json URI (else under CHECKPOINT prefix)
#   GUTIBM_BINARY    path to gut_ibm (default /opt/gutibm/gut_ibm)
#   MPI_RANKS        default 1
#   EXTRA_MPIRUN_ARGS  optional extra flags
#   AWS_BATCH_JOB_ARRAY_INDEX  set automatically for array jobs
#   AWS_BATCH_JOB_ID           set automatically by Batch

set -euo pipefail

BINARY="${GUTIBM_BINARY:-/opt/gutibm/gut_ibm}"
MPI_RANKS="${MPI_RANKS:-1}"
WORK="${GUTIBM_WORK_DIR:-/tmp/gutibm_job}"
INDEX="${AWS_BATCH_JOB_ARRAY_INDEX:-}"
JOB_ID="${AWS_BATCH_JOB_ID:-}"
CHECKPOINT_INTERVAL_SECONDS="${CHECKPOINT_INTERVAL_SECONDS:-300}"
CHECKPOINT_NAME="checkpoint.h5"
STATUS_NAME="status.json"
SYNC_PID=""
MPIRUN_PID=""
SPOT_NOTICED=0
RESUME_FROM_CHECKPOINT=0
RUN_EXIT_CODE=0
WALL_START_EPOCH="$(date -u +%s)"
mkdir -p "${WORK}"
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

# Upload the newest checkpoint (the working output.h5 already holds the step
# snapshot groups the C++ side writes; init_from_checkpoint restarts from the
# latest step group). Best-effort: a mid-write copy just misses the newest group.
upload_checkpoint() {
  [[ -n "${CHECKPOINT_S3_PREFIX:-}" && -f "${WORK}/output.h5" ]] || return 0
  if aws s3 cp "${WORK}/output.h5" "${CHECKPOINT_S3_PREFIX}${CHECKPOINT_NAME}" >/dev/null 2>&1; then
    CHECKPOINT_UPLOADED_AT="$(iso_now)"
    export CHECKPOINT_UPLOADED_AT
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

write_status_json() {
  local state="$1"
  local status_path="${WORK}/${STATUS_NAME}"
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
  GUTIBM_RESUME_FROM_CHECKPOINT="${RESUME_FROM_CHECKPOINT}" \
  GUTIBM_SPOT_INTERRUPTION="${SPOT_NOTICED}" \
  GUTIBM_SPOT_ACTION_TIME="${spot_action}" \
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
    "resume_from_checkpoint": os.environ.get("GUTIBM_RESUME_FROM_CHECKPOINT", "0") == "1",
    "spot_interruption": spot,
    "spot_action_time": os.environ.get("GUTIBM_SPOT_ACTION_TIME") or None,
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

if [[ -n "${CHECKPOINT_S3_PREFIX:-}" ]]; then
  if aws s3 ls "${CHECKPOINT_S3_PREFIX}${CHECKPOINT_NAME}" >/dev/null 2>&1; then
    echo "Resuming from checkpoint: ${CHECKPOINT_S3_PREFIX}${CHECKPOINT_NAME}"
    aws s3 cp "${CHECKPOINT_S3_PREFIX}${CHECKPOINT_NAME}" "${WORK}/checkpoint.h5"
    RESUME_FROM_CHECKPOINT=1
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
fi

GUTIBM_WORK_DIR="${WORK}" python3 - <<'PY'
import json, os
from pathlib import Path
work = Path(os.environ["GUTIBM_WORK_DIR"])
cfg = json.loads((work / "input.json").read_text())
cfg["hdf5_file"] = str(work / "output.h5")
cfg.setdefault("gpu_enabled", True)
(work / "input.json").write_text(json.dumps(cfg, indent=2) + "\n")
PY

write_status_json "starting"
upload_status

if [[ -n "${CHECKPOINT_S3_PREFIX:-}" ]] || [[ -n "$(resolve_status_uri)" ]]; then
  echo "Status/checkpoint sync every ${CHECKPOINT_INTERVAL_SECONDS}s"
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

if [[ "${RUN_EXIT_CODE}" -ne 0 ]]; then
  write_status_json "failed"
  upload_status
  echo "gut_ibm exited with ${RUN_EXIT_CODE}" >&2
  exit "${RUN_EXIT_CODE}"
fi

if [[ -f "${WORK}/output.h5" ]]; then
  gzip -c "${WORK}/output.h5" > "${WORK}/output.h5.gz"
  echo "Uploading ${OUTPUT_S3_URI}"
  aws s3 cp "${WORK}/output.h5.gz" "${OUTPUT_S3_URI}"
else
  write_status_json "failed"
  upload_status
  echo "No output.h5 produced" >&2
  exit 1
fi

write_status_json "succeeded"
upload_status
echo "Done."
