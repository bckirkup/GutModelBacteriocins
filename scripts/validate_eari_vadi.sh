#!/usr/bin/env bash
# EARI/VADI + FISH observability invariant checks for CI (issues #56, #25).
# Build gut_ibm, run the short validation scenario, and check durable invariants.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build-eari-vadi"
EXAMPLE="$ROOT/examples/eari_vadi_validation/input.json"
OUTPUT="eari_vadi_validation.h5"
export CC="${CC:-gcc}"
export CXX="${CXX:-g++}"
PYTHON="${PYTHON:-python3}"

echo "=== Configuring GutIBM (Release) ==="
cmake -B "$BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DGUTIBM_USE_MPI=ON \
  -DGUTIBM_USE_HDF5=ON

echo "=== Building gut_ibm ==="
cmake --build "$BUILD" -j"$(nproc)" --target gut_ibm

echo "=== Running EARI/VADI validation scenario ==="
rm -f "$BUILD/$OUTPUT"
(
  cd "$BUILD"
  mpirun --allow-run-as-root -np 1 ./gut_ibm "$EXAMPLE"
)

echo "=== Validating HDF5 metrics and invariants ==="
"$PYTHON" -m pip install -q -e "$ROOT/python/.[dev]"
"$PYTHON" -m gut_ibm_tools.validation_regression \
  "$BUILD/$OUTPUT" \
  --check-fish-targets

echo "EARI/VADI + FISH observability invariant checks passed."
