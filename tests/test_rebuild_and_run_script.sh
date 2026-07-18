#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/rebuild_and_run.sh"

fail() {
  echo "test_rebuild_and_run_script: FAILED: $*" >&2
  exit 1
}

[[ "$(json_kind "$ROOT/experiments/smoke_single.json")" == single ]] ||
  fail "smoke_single.json was not classified as single"
[[ "$(json_kind "$ROOT/experiments/smoke_batch.json")" == batch ]] ||
  fail "smoke_batch.json was not classified as batch"
[[ "$(json_kind "$ROOT/experiments/diversity_campaign/stage1_motility_validation/1a_motility_off.json")" == single ]] ||
  fail "stage1 1a was not classified as single"
[[ "$(json_kind "$ROOT/experiments/diversity_campaign/stage3_campaign/batch_kd_sweep.json")" == batch ]] ||
  fail "stage3 batch_kd_sweep was not classified as batch"

# Stage 3 singles request GPU (50M-cell grids); batches inherit via base_config.
stage3_gpu_cfg="$ROOT/experiments/diversity_campaign/stage3_campaign/3b_full_mechanisms.json"
"$PYTHON" - "$stage3_gpu_cfg" <<'PY' || fail "stage3 base missing gpu_enabled"
import json, sys
from pathlib import Path
cfg = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
assert cfg.get("gpu_enabled") is True, cfg.get("gpu_enabled")
assert cfg.get("gpu_device_id") == -1, cfg.get("gpu_device_id")
PY
kd_base="$("$PYTHON" -c "import json; from pathlib import Path; print(json.loads(Path('$ROOT/experiments/diversity_campaign/stage3_campaign/batch_kd_sweep.json').read_text())['base_config'])")"
[[ "$kd_base" == *3b_full_mechanisms.json ]] ||
  fail "batch_kd_sweep base_config is not 3b_full_mechanisms.json"

temporary_dir="$(mktemp -d)"
trap 'rm -rf "$temporary_dir"' EXIT
sensitivity_config="$temporary_dir/classification.json"
printf '%s\n' \
  '{"base_config":"experiments/smoke_single.json","runs":[]}' \
  >"$sensitivity_config"
[[ "$(json_kind "$sensitivity_config")" == batch ]] ||
  fail "runs key did not produce batch classification"
printf '%s\n' '{"base_config":"experiments/smoke_single.json"}' \
  >"$sensitivity_config"
[[ "$(json_kind "$sensitivity_config")" == single ]] ||
  fail "removing runs key did not change classification"

mapfile -t single_configs < <(collect_json_files single)
mapfile -t batch_configs < <(collect_json_files batch)
printf '%s\n' "${single_configs[@]}" |
  grep -qx "$ROOT/experiments/smoke_single.json" ||
  fail "single config discovery omitted smoke_single.json"
printf '%s\n' "${batch_configs[@]}" |
  grep -qx "$ROOT/experiments/smoke_batch.json" ||
  fail "batch discovery omitted smoke_batch.json"
printf '%s\n' "${single_configs[@]}" |
  grep -q "stage1_motility_validation/1a_motility_off.json" ||
  fail "single discovery omitted stage1 1a"
printf '%s\n' "${batch_configs[@]}" |
  grep -q "stage3_campaign/batch_kd_sweep.json" ||
  fail "batch discovery omitted stage3 kd sweep"

mapfile -t stages < <(collect_campaign_stages)
((${#stages[@]} == 3)) ||
  fail "expected 3 campaign stages, found ${#stages[@]}"
[[ "${stages[0]}" == *stage1_motility_validation ]] ||
  fail "first campaign stage is not stage1"
[[ "${stages[1]}" == *stage2_mechanism_validation ]] ||
  fail "second campaign stage is not stage2"
[[ "${stages[2]}" == *stage3_campaign ]] ||
  fail "third campaign stage is not stage3"

mapfile -t stage1_configs < <(
  collect_stage_single_configs \
    "$ROOT/experiments/diversity_campaign/stage1_motility_validation"
)
((${#stage1_configs[@]} == 6)) ||
  fail "stage1 should expose 6 single configs, found ${#stage1_configs[@]}"
[[ "$(basename "${stage1_configs[0]}")" == 1a_motility_off.json ]] ||
  fail "stage1 singles are not sorted starting at 1a"

stage_path="$(resolve_stage_path \
  experiments/diversity_campaign/stage2_mechanism_validation)"
[[ "$stage_path" == "$ROOT/experiments/diversity_campaign/stage2_mechanism_validation" ]] ||
  fail "resolve_stage_path did not resolve stage2"

[[ "$(printf '2\n' | prompt_batch_action 2>/dev/null)" == dry-run ]] ||
  fail "batch action prompt did not return dry-run"

fake_bin="$temporary_dir/bin"
fake_build="$temporary_dir/build"
mkdir -p "$fake_bin" "$fake_build"
cat >"$fake_bin/mpirun" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$*" >"$MPIRUN_LOG"
EOF
cat >"$fake_build/gut_ibm" <<'EOF'
#!/usr/bin/env bash
exit 0
EOF
cat >"$fake_bin/python" <<'EOF'
#!/usr/bin/env bash
argument_1="$1"
argument_2="${2:-}"
if [[ "$argument_1" == "-" ]]; then
  echo "build/gut_ibm"
elif [[ "$argument_1" == "-m" &&
        "$argument_2" == "gut_ibm_tools.batch_runner" ]]; then
  printf '%s\n' "$*" >"$BATCH_RUN_LOG"
  echo "Batch: 2 jobs"
else
  exit 2
fi
EOF
chmod +x "$fake_bin/mpirun" "$fake_bin/python" "$fake_build/gut_ibm"
export MPIRUN_LOG="$temporary_dir/mpirun.log"
export BATCH_RUN_LOG="$temporary_dir/batch-run.log"
original_build_dir="$BUILD_DIR"
BUILD_DIR="$fake_build"
PATH="$fake_bin:$PATH" run_single \
  "$ROOT/experiments/smoke_single.json" 2 >/dev/null
BUILD_DIR="$original_build_dir"
grep -q -- '--bind-to none' "$MPIRUN_LOG" ||
  fail "single-run MPI command omitted --bind-to none"
grep -q -- "-np 2 $fake_build/gut_ibm" "$MPIRUN_LOG" ||
  fail "single-run MPI command was not constructed correctly"

# Post-run gzip: missing/disabled HDF5 is a no-op; enabled path compresses.
fake_h5="$temporary_dir/campaign_out.h5"
printf 'HDF5-fake-%s' "$(printf 'x%.0s' {1..4000})" >"$fake_h5"
gzip_hdf5_file "$fake_h5"
[[ ! -f "$fake_h5" ]] || fail "gzip_hdf5_file left uncompressed HDF5"
[[ -f "${fake_h5}.gz" ]] || fail "gzip_hdf5_file did not create .h5.gz"
GZIP_HDF5=false
maybe_gzip_hdf5_from_config "$ROOT/experiments/smoke_single.json" ||
  fail "maybe_gzip with disabled HDF5 should no-op"
GZIP_HDF5=true
# smoke_single has hdf5.enabled=false → empty path / no-op
maybe_gzip_hdf5_from_config "$ROOT/experiments/smoke_single.json" ||
  fail "maybe_gzip with hdf5.enabled=false should no-op"

stage1_cfg="$ROOT/experiments/diversity_campaign/stage1_motility_validation/1a_motility_off.json"
hdf5_from_1a="$(hdf5_path_from_config "$stage1_cfg" "$temporary_dir")"
[[ "$hdf5_from_1a" == "$temporary_dir/1a_motility_off.h5" ]] ||
  fail "hdf5_path_from_config did not resolve stage1 filename"

original_python="$PYTHON"
PYTHON="$fake_bin/python"
GZIP_HDF5=true
run_batch "$ROOT/experiments/smoke_batch.json" dry-run |
  grep -q 'Batch: 2 jobs' ||
  fail "batch dry-run did not expand two jobs"
PYTHON="$original_python"
grep -q -- '--dry-run' "$BATCH_RUN_LOG" ||
  fail "batch dry-run did not pass --dry-run to the batch runner"

fake_venv="$temporary_dir/pipless-venv"
mkdir -p "$fake_venv/bin"
cat >"$fake_venv/bin/python" <<'EOF'
#!/usr/bin/env bash
command_line="$*"
if [[ "$command_line" == "-m pip --version" ]]; then
  [[ -f "$PIP_BOOTSTRAP_MARKER" ]]
elif [[ "$command_line" == "-m ensurepip --upgrade" ]]; then
  touch "$PIP_BOOTSTRAP_MARKER"
  echo "$command_line" >>"$PIP_BOOTSTRAP_LOG"
elif [[ "$command_line" == "-m pip install "* ]]; then
  [[ -f "$PIP_BOOTSTRAP_MARKER" ]] || exit 2
  echo "$command_line" >>"$PIP_BOOTSTRAP_LOG"
else
  exit 2
fi
EOF
chmod +x "$fake_venv/bin/python"
export PIP_BOOTSTRAP_MARKER="$temporary_dir/pip-ready"
export PIP_BOOTSTRAP_LOG="$temporary_dir/pip-bootstrap.log"
original_venv_dir="$VENV_DIR"
original_python="$PYTHON"
VENV_DIR="$fake_venv"
PYTHON=python3
ensure_python_environment >/dev/null
VENV_DIR="$original_venv_dir"
PYTHON="$original_python"
grep -q -- '-m ensurepip --upgrade' "$PIP_BOOTSTRAP_LOG" ||
  fail "pip-less venv did not invoke ensurepip"
grep -q -- '-m pip install --quiet -e' "$PIP_BOOTSTRAP_LOG" ||
  fail "package installation did not follow pip bootstrap"

"$ROOT/rebuild_and_run.sh" --help |
  grep -q -- '--gzip-hdf5' ||
  fail "help output omitted --gzip-hdf5"
"$ROOT/rebuild_and_run.sh" --help |
  grep -q -- '--mode prompt|single|batch|stage|none' ||
  fail "help output omitted run modes"

if GUTIBM_BUILD_DIR="$ROOT" "$ROOT/rebuild_and_run.sh" \
  --reuse-build --mode none >"$temporary_dir/safety.log" 2>&1; then
  fail "repository-root build path was accepted"
fi
grep -q 'build directory cannot be the repository root' \
  "$temporary_dir/safety.log" ||
  fail "unsafe build path did not produce the expected error"

echo "All rebuild-and-run script tests passed."
