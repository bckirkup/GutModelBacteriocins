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
chmod +x "$fake_bin/mpirun" "$fake_build/gut_ibm"
export MPIRUN_LOG="$temporary_dir/mpirun.log"
original_build_dir="$BUILD_DIR"
BUILD_DIR="$fake_build"
PATH="$fake_bin:$PATH" run_single \
  "$ROOT/experiments/smoke_single.json" 2 >/dev/null
BUILD_DIR="$original_build_dir"
grep -q -- '--bind-to none' "$MPIRUN_LOG" ||
  fail "single-run MPI command omitted --bind-to none"
grep -q -- "-np 2 $fake_build/gut_ibm" "$MPIRUN_LOG" ||
  fail "single-run MPI command was not constructed correctly"

run_batch "$ROOT/experiments/smoke_batch.json" dry-run |
  grep -q 'Batch: 2 jobs' ||
  fail "batch dry-run did not expand two jobs"

"$ROOT/rebuild_and_run.sh" --help |
  grep -q -- '--mode prompt|single|batch|none' ||
  fail "help output omitted run modes"

if GUTIBM_BUILD_DIR="$ROOT" "$ROOT/rebuild_and_run.sh" \
  --reuse-build --mode none >"$temporary_dir/safety.log" 2>&1; then
  fail "repository-root build path was accepted"
fi
grep -q 'build directory cannot be the repository root' \
  "$temporary_dir/safety.log" ||
  fail "unsafe build path did not produce the expected error"

echo "All rebuild-and-run script tests passed."
