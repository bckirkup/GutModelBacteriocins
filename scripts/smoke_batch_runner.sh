#!/usr/bin/env bash
# Batch runner smoke: one short simulation via the Python batch CLI.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${GUTIBM_BUILD_DIR:-$ROOT/build}"
BATCH="${BATCH_CONFIG:-$ROOT/examples/batch_scan/batch_ci.json}"
PYTHON="${PYTHON:-python3}"

if [[ ! -x "$BUILD/gut_ibm" && ! -f "$BUILD/gut_ibm" ]]; then
  echo "ERROR: gut_ibm not found in $BUILD (set GUTIBM_BUILD_DIR)" >&2
  exit 1
fi

echo "=== Install Python batch runner ==="
"$PYTHON" -m pip install -q -e "$ROOT/python/.[dev]"

echo "=== Batch runner dry-run ==="
"$PYTHON" -m gut_ibm_tools.batch_runner "$BATCH" --dry-run

echo "=== Batch runner smoke (1 job) ==="
"$PYTHON" -m gut_ibm_tools.batch_runner "$BATCH"

echo "Batch runner smoke passed."
