#!/usr/bin/env bash
# AWS Batch entrypoint (draft) — docs/AWS_BATCH.md Phase 1–4
#
# Expected environment (set in the Batch job definition / overrides):
#   INPUT_S3_URI     s3://bucket/path/input.json   (or set INPUT_S3_PREFIX)
#   INPUT_S3_PREFIX  s3://bucket/path/jobs         → ${PREFIX}/${INDEX}/input.json
#   OUTPUT_S3_URI    s3://bucket/path/output.h5.gz (or set OUTPUT_S3_PREFIX)
#   OUTPUT_S3_PREFIX s3://bucket/path/jobs         → ${PREFIX}/${INDEX}/output.h5.gz
#   CHECKPOINT_S3_PREFIX  s3://bucket/.../ckpt/    (optional; per-index if ARRAY)
#   GUTIBM_BINARY    path to gut_ibm (default /opt/gutibm/gut_ibm)
#   MPI_RANKS        default 1
#   EXTRA_MPIRUN_ARGS  optional extra flags
#   AWS_BATCH_JOB_ARRAY_INDEX  set automatically for array jobs

set -euo pipefail

BINARY="${GUTIBM_BINARY:-/opt/gutibm/gut_ibm}"
MPI_RANKS="${MPI_RANKS:-1}"
WORK="${GUTIBM_WORK_DIR:-/tmp/gutibm_job}"
INDEX="${AWS_BATCH_JOB_ARRAY_INDEX:-}"
mkdir -p "${WORK}"

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
  LATEST="$(aws s3 ls "${CHECKPOINT_S3_PREFIX}" | awk '{print $4}' | sort | tail -n 1 || true)"
  if [[ -n "${LATEST}" ]]; then
    echo "Resuming from checkpoint: ${CHECKPOINT_S3_PREFIX}${LATEST}"
    aws s3 cp "${CHECKPOINT_S3_PREFIX}${LATEST}" "${WORK}/checkpoint.h5"
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

echo "Running ${BINARY} (mpi_ranks=${MPI_RANKS})"
# shellcheck disable=SC2086
mpirun -np "${MPI_RANKS}" ${EXTRA_MPIRUN_ARGS:-} "${BINARY}" "${WORK}/input.json"

if [[ -f "${WORK}/output.h5" ]]; then
  gzip -c "${WORK}/output.h5" > "${WORK}/output.h5.gz"
  echo "Uploading ${OUTPUT_S3_URI}"
  aws s3 cp "${WORK}/output.h5.gz" "${OUTPUT_S3_URI}"
else
  echo "No output.h5 produced" >&2
  exit 1
fi

echo "Done."
