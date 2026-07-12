#!/usr/bin/env bash
# Cleanly rebuild GutIBM, run CTest, then launch experiment JSON files.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPERIMENTS_DIR="$ROOT/experiments"
VENV_DIR="$ROOT/.venv"
BUILD_DIR="$(realpath -m "${GUTIBM_BUILD_DIR:-$ROOT/build}")"
BUILD_TYPE="${GUTIBM_BUILD_TYPE:-Release}"
CUDA_MODE="${GUTIBM_CUDA:-auto}"
RUN_MODE="${GUTIBM_RUN_MODE:-prompt}"
CONFIG_PATH="${GUTIBM_CONFIG:-}"
MPI_RANKS="${GUTIBM_MPI_RANKS:-}"
BATCH_ACTION="${GUTIBM_BATCH_ACTION:-run}"
REUSE_BUILD=false
PYTHON="${PYTHON:-python3}"
CUDA_ENABLED=OFF
NVIDIA_SMI=""
NVCC=""

usage() {
  cat <<'EOF'
Usage: ./rebuild_and_run.sh [options]

Default workflow:
  1. Remove build/
  2. Configure MPI + HDF5 and auto-detect CUDA
  3. Build GutIBM
  4. Run the complete CTest suite
  5. Prompt for single or batch JSON files under experiments/

Options:
  --cuda auto|on|off        CUDA build mode (default: auto)
  --mode prompt|single|batch|none
                            Run menu, one config, one batch, or stop after CTest
  --config PATH             JSON file for --mode single or batch
  --mpi-ranks N             MPI ranks for a single run
  --batch-action ACTION     run, dry-run, resume, or status
  --reuse-build             Skip clean/configure/build/CTest
  -h, --help                Show this help

Environment equivalents:
  GUTIBM_BUILD_DIR, GUTIBM_BUILD_TYPE, GUTIBM_BUILD_JOBS,
  GUTIBM_CUDA, GUTIBM_CUDA_ARCHITECTURES, GUTIBM_RUN_MODE,
  GUTIBM_CONFIG, GUTIBM_MPI_RANKS, GUTIBM_BATCH_ACTION, PYTHON
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

parse_args() {
  while (($# > 0)); do
    case "$1" in
      --cuda)
        (($# >= 2)) || die "--cuda requires a value"
        CUDA_MODE="$2"
        shift 2
        ;;
      --mode)
        (($# >= 2)) || die "--mode requires a value"
        RUN_MODE="$2"
        shift 2
        ;;
      --config)
        (($# >= 2)) || die "--config requires a path"
        CONFIG_PATH="$2"
        shift 2
        ;;
      --mpi-ranks)
        (($# >= 2)) || die "--mpi-ranks requires a value"
        MPI_RANKS="$2"
        shift 2
        ;;
      --batch-action)
        (($# >= 2)) || die "--batch-action requires a value"
        BATCH_ACTION="$2"
        shift 2
        ;;
      --reuse-build)
        REUSE_BUILD=true
        shift
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        die "unknown option: $1"
        ;;
    esac
  done
}

validate_options() {
  CUDA_MODE="${CUDA_MODE,,}"
  RUN_MODE="${RUN_MODE,,}"
  BATCH_ACTION="${BATCH_ACTION,,}"

  case "$CUDA_MODE" in
    auto|on|off) ;;
    *) die "CUDA mode must be auto, on, or off" ;;
  esac

  case "$RUN_MODE" in
    prompt|single|batch|none) ;;
    *) die "run mode must be prompt, single, batch, or none" ;;
  esac

  case "$BATCH_ACTION" in
    run|dry-run|resume|status) ;;
    *) die "batch action must be run, dry-run, resume, or status" ;;
  esac

  case "$BUILD_DIR" in
    "$ROOT")
      die "build directory cannot be the repository root"
      ;;
    "$ROOT"/*) ;;
    *)
      die "build directory must be inside $ROOT"
      ;;
  esac
}

find_cuda_tools() {
  NVCC="$(command -v nvcc || true)"
  if [[ -z "$NVCC" && -x /usr/local/cuda/bin/nvcc ]]; then
    NVCC=/usr/local/cuda/bin/nvcc
  fi

  NVIDIA_SMI="$(command -v nvidia-smi || true)"
  if [[ -z "$NVIDIA_SMI" && -x /usr/lib/wsl/lib/nvidia-smi ]]; then
    NVIDIA_SMI=/usr/lib/wsl/lib/nvidia-smi
  fi

  local device_available=false
  if [[ -n "$NVIDIA_SMI" ]] && "$NVIDIA_SMI" -L >/dev/null 2>&1; then
    device_available=true
  fi

  case "$CUDA_MODE" in
    auto)
      if [[ -n "$NVCC" && "$device_available" == true ]]; then
        CUDA_ENABLED=ON
      fi
      ;;
    on)
      [[ -n "$NVCC" ]] || die "CUDA requested, but nvcc was not found"
      [[ "$device_available" == true ]] ||
        die "CUDA requested, but nvidia-smi cannot see a GPU"
      CUDA_ENABLED=ON
      ;;
    off)
      CUDA_ENABLED=OFF
      ;;
  esac
}

load_reused_build_options() {
  local cache="$BUILD_DIR/CMakeCache.txt"
  [[ -f "$cache" ]] ||
    die "--reuse-build requested, but $cache is missing"
  CUDA_ENABLED="$(sed -n 's/^GUTIBM_USE_CUDA:BOOL=//p' "$cache" | head -1)"
  [[ "$CUDA_ENABLED" == ON || "$CUDA_ENABLED" == OFF ]] ||
    die "cannot determine CUDA mode from $cache"
}

detect_cuda_architectures() {
  if [[ -n "${GUTIBM_CUDA_ARCHITECTURES:-}" ]]; then
    printf '%s\n' "$GUTIBM_CUDA_ARCHITECTURES"
    return
  fi

  [[ -n "$NVIDIA_SMI" ]] || return
  "$NVIDIA_SMI" --query-gpu=compute_cap --format=csv,noheader 2>/dev/null |
    awk '/^[0-9]+\.[0-9]+$/ { gsub(/\./, ""); print }' |
    sort -u |
    paste -sd';' -
}

select_compilers() {
  if [[ -z "${CC:-}" ]]; then
    if command -v gcc-13 >/dev/null 2>&1; then
      CC=gcc-13
    else
      CC=gcc
    fi
  fi
  if [[ -z "${CXX:-}" ]]; then
    if command -v g++-13 >/dev/null 2>&1; then
      CXX=g++-13
    else
      CXX=g++
    fi
  fi
  export CC CXX
  require_command "$CC"
  require_command "$CXX"
}

ensure_python_environment() {
  require_command "$PYTHON"
  if [[ ! -x "$VENV_DIR/bin/python" ]]; then
    echo "=== Creating Python virtual environment ==="
    "$PYTHON" -m venv "$VENV_DIR" ||
      die "failed to create .venv; install the python3-venv package"
  fi

  PYTHON="$VENV_DIR/bin/python"
  echo "=== Installing GutIBM Python tools in .venv ==="
  "$PYTHON" -m pip install --quiet -e "$ROOT/python/.[dev]"
}

default_build_jobs() {
  local available
  available="$(nproc)"
  if ((available > 8)); then
    echo 8
  else
    echo "$available"
  fi
}

print_environment_summary() {
  echo "=== GutIBM environment ==="
  echo "Repository: $ROOT"
  echo "Build:      $BUILD_DIR ($BUILD_TYPE)"
  echo "Python:     $PYTHON"
  echo "Compiler:   $CC / $CXX"
  echo "MPI:        $(mpirun --version | head -1)"
  echo "CUDA:       $CUDA_ENABLED"
  if [[ "$CUDA_ENABLED" == ON ]]; then
    echo "nvcc:       $NVCC"
    "$NVIDIA_SMI" -L
  elif grep -qi microsoft /proc/version 2>/dev/null; then
    echo "CUDA note:  CPU build selected; see docs/WSL2_SETUP.md"
  fi
  echo
}

clean_configure_build_test() {
  local build_jobs="${GUTIBM_BUILD_JOBS:-$(default_build_jobs)}"
  [[ "$build_jobs" =~ ^[1-9][0-9]*$ ]] ||
    die "GUTIBM_BUILD_JOBS must be a positive integer"

  echo "=== Removing existing build directory ==="
  echo "$BUILD_DIR"
  rm -rf -- "$BUILD_DIR"

  local -a cmake_args=(
    -S "$ROOT"
    -B "$BUILD_DIR"
    "-DCMAKE_BUILD_TYPE=$BUILD_TYPE"
    -DGUTIBM_USE_MPI=ON
    -DGUTIBM_USE_HDF5=ON
    -DGUTIBM_USE_OPENMP=OFF
    "-DGUTIBM_USE_CUDA=$CUDA_ENABLED"
    "-DCMAKE_CXX_FLAGS=-Wall -Wextra"
  )

  if [[ "$CUDA_ENABLED" == ON ]]; then
    export CUDACXX="$NVCC"
    cmake_args+=("-DCMAKE_CUDA_COMPILER=$NVCC")
    local cuda_architectures
    cuda_architectures="$(detect_cuda_architectures)"
    if [[ -n "$cuda_architectures" ]]; then
      cmake_args+=("-DCMAKE_CUDA_ARCHITECTURES=$cuda_architectures")
    fi
  fi

  echo "=== Configuring GutIBM ==="
  cmake "${cmake_args[@]}"

  echo "=== Building GutIBM with $build_jobs jobs ==="
  cmake --build "$BUILD_DIR" --parallel "$build_jobs"

  echo "=== Running CTest ==="
  ctest --test-dir "$BUILD_DIR" --output-on-failure
}

json_kind() {
  "$PYTHON" - "$1" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
try:
    payload = json.loads(path.read_text(encoding="utf-8"))
except (OSError, json.JSONDecodeError):
    print("invalid")
    raise SystemExit

if not isinstance(payload, dict):
    print("invalid")
    raise SystemExit

is_batch = (
    "base_config" in payload
    and (("sweep" in payload) ^ ("runs" in payload))
)
print("batch" if is_batch else "single")
PY
}

collect_json_files() {
  local expected_kind="$1"
  "$PYTHON" - "$EXPERIMENTS_DIR" "$expected_kind" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
expected = sys.argv[2]
if not root.is_dir():
    raise SystemExit

for path in sorted(root.rglob("*.json")):
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"WARNING: skipping invalid JSON {path}: {exc}", file=sys.stderr)
        continue
    if not isinstance(payload, dict):
        print(f"WARNING: skipping non-object JSON {path}", file=sys.stderr)
        continue
    is_batch = (
        "base_config" in payload
        and (("sweep" in payload) ^ ("runs" in payload))
    )
    kind = "batch" if is_batch else "single"
    if kind == expected:
        print(path)
PY
}

choose_json() {
  local kind="$1"
  local -a files=()
  mapfile -t files < <(collect_json_files "$kind")
  if ((${#files[@]} == 0)); then
    echo "No $kind JSON files found under $EXPERIMENTS_DIR" >&2
    return 1
  fi

  echo "Available $kind JSON files:"
  local index
  for index in "${!files[@]}"; do
    printf '  %d) %s\n' "$((index + 1))" "${files[$index]#"$ROOT"/}"
  done

  local choice
  while true; do
    read -r -p "Select a file [1-${#files[@]}]: " choice
    if [[ "$choice" =~ ^[1-9][0-9]*$ ]] &&
       ((choice >= 1 && choice <= ${#files[@]})); then
      SELECTED_JSON="${files[$((choice - 1))]}"
      return
    fi
    echo "Invalid selection."
  done
}

resolve_config_path() {
  local candidate="$1"
  if [[ -f "$candidate" ]]; then
    candidate="$(realpath "$candidate")"
  elif [[ -f "$ROOT/$candidate" ]]; then
    candidate="$(realpath "$ROOT/$candidate")"
  elif [[ -f "$EXPERIMENTS_DIR/$candidate" ]]; then
    candidate="$(realpath "$EXPERIMENTS_DIR/$candidate")"
  else
    die "JSON config not found: $1"
  fi

  case "$candidate" in
    "$EXPERIMENTS_DIR"/*) ;;
    *) die "JSON config must be inside $EXPERIMENTS_DIR" ;;
  esac
  printf '%s\n' "$candidate"
}

validate_ranks() {
  [[ "$1" =~ ^[1-9][0-9]*$ ]] || die "MPI ranks must be a positive integer"
}

gpu_config_value() {
  "$PYTHON" - "$1" "$2" <<'PY'
import json
import sys
from pathlib import Path

payload = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
value = payload.get(sys.argv[2], "")
if isinstance(value, bool):
    print(str(value).lower())
else:
    print(value)
PY
}

run_single() {
  local config="$1"
  local ranks="$2"
  validate_ranks "$ranks"

  local gpu_enabled
  gpu_enabled="$(gpu_config_value "$config" gpu_enabled)"
  if [[ "$gpu_enabled" == true && "$CUDA_ENABLED" != ON ]]; then
    die "$(basename "$config") enables GPU runtime, but this build has CUDA disabled"
  fi

  local gpu_device_id
  gpu_device_id="$(gpu_config_value "$config" gpu_device_id)"
  if ((ranks > 1)) && [[ "$gpu_enabled" == true ]] &&
     [[ -n "$gpu_device_id" && "$gpu_device_id" != -1 ]]; then
    echo "WARNING: gpu_device_id=$gpu_device_id pins every MPI rank to one GPU." >&2
    echo "         Use gpu_device_id=-1 for rank-based device selection." >&2
  fi

  local -a mpirun_args=(mpirun --bind-to none)
  if [[ "$(id -u)" -eq 0 ]]; then
    mpirun_args+=(--allow-run-as-root)
  fi
  if ((ranks > $(nproc))); then
    mpirun_args+=(--oversubscribe)
  fi

  echo "=== Single experiment ==="
  echo "Config: ${config#"$ROOT"/}"
  echo "Ranks:  $ranks"
  (
    cd "$ROOT"
    OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
      "${mpirun_args[@]}" -np "$ranks" "$BUILD_DIR/gut_ibm" "$config"
  )
}

run_batch() {
  local config="$1"
  local action="$2"
  local manifest_binary
  manifest_binary="$("$PYTHON" - "$config" <<'PY'
import json
import sys
from pathlib import Path

payload = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
print(payload.get("binary", "build/gut_ibm"))
PY
)"
  local resolved_binary
  if [[ "$manifest_binary" == /* ]]; then
    resolved_binary="$(realpath -m "$manifest_binary")"
  else
    resolved_binary="$(realpath -m "$ROOT/$manifest_binary")"
  fi
  if [[ "$resolved_binary" != "$BUILD_DIR/gut_ibm" ]]; then
    echo "WARNING: manifest binary resolves to $resolved_binary" >&2
    echo "         current build executable is $BUILD_DIR/gut_ibm" >&2
  fi

  local -a action_args=()
  case "$action" in
    run) ;;
    dry-run) action_args+=(--dry-run) ;;
    resume) action_args+=(--resume) ;;
    status) action_args+=(--status) ;;
    *) die "unknown batch action: $action" ;;
  esac

  echo "=== Batch experiment ($action) ==="
  echo "Config: ${config#"$ROOT"/}"
  (
    cd "$ROOT"
    PYTHONPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}" \
      "$PYTHON" -m gut_ibm_tools.batch_runner "$config" "${action_args[@]}"
  )
}

prompt_ranks() {
  local default_ranks="${MPI_RANKS:-1}"
  local ranks
  read -r -p "MPI ranks [$default_ranks]: " ranks
  printf '%s\n' "${ranks:-$default_ranks}"
}

prompt_batch_action() {
  echo "Batch action:" >&2
  echo "  1) run" >&2
  echo "  2) dry-run" >&2
  echo "  3) resume" >&2
  echo "  4) status" >&2
  local choice
  while true; do
    read -r -p "Select an action [1-4]: " choice
    case "$choice" in
      1) echo run; return ;;
      2) echo dry-run; return ;;
      3) echo resume; return ;;
      4) echo status; return ;;
      *) echo "Invalid selection." >&2 ;;
    esac
  done
}

interactive_menu() {
  while true; do
    echo
    echo "Experiment menu:"
    echo "  1) Run one simulation JSON"
    echo "  2) Run a batch manifest"
    echo "  3) Exit"
    local choice
    read -r -p "Select an option [1-3]: " choice
    case "$choice" in
      1)
        if choose_json single; then
          run_single "$SELECTED_JSON" "$(prompt_ranks)"
        fi
        ;;
      2)
        if choose_json batch; then
          run_batch "$SELECTED_JSON" "$(prompt_batch_action)"
        fi
        ;;
      3)
        return
        ;;
      *)
        echo "Invalid selection."
        ;;
    esac
  done
}

main() {
  parse_args "$@"
  validate_options
  require_command realpath
  require_command cmake
  require_command ctest
  require_command mpirun
  select_compilers
  ensure_python_environment
  find_cuda_tools
  if [[ "$REUSE_BUILD" == true ]]; then
    load_reused_build_options
  fi
  print_environment_summary

  if [[ "$REUSE_BUILD" == true ]]; then
    [[ -x "$BUILD_DIR/gut_ibm" ]] ||
      die "--reuse-build requested, but $BUILD_DIR/gut_ibm is missing"
  else
    clean_configure_build_test
  fi

  case "$RUN_MODE" in
    none)
      echo "Build and CTest complete."
      ;;
    single)
      [[ -n "$CONFIG_PATH" ]] || die "--mode single requires --config"
      CONFIG_PATH="$(resolve_config_path "$CONFIG_PATH")"
      [[ "$(json_kind "$CONFIG_PATH")" == single ]] ||
        die "selected JSON is not a single-run configuration"
      run_single "$CONFIG_PATH" "${MPI_RANKS:-1}"
      ;;
    batch)
      [[ -n "$CONFIG_PATH" ]] || die "--mode batch requires --config"
      CONFIG_PATH="$(resolve_config_path "$CONFIG_PATH")"
      [[ "$(json_kind "$CONFIG_PATH")" == batch ]] ||
        die "selected JSON is not a batch manifest"
      run_batch "$CONFIG_PATH" "$BATCH_ACTION"
      ;;
    prompt)
      interactive_menu
      ;;
  esac
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  main "$@"
fi
