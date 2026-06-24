#!/usr/bin/env bash
# EARI/VADI validation regression for CI (issue #56).
# Build gut_ibm, run the short validation scenario, compare HDF5 metrics to golden baseline.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build-eari-vadi"
EXAMPLE="$ROOT/examples/eari_vadi_validation/input.json"
GOLDEN="$ROOT/python/tests/fixtures/eari_vadi_ci_golden.json"
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

echo "=== Validating HDF5 metrics against golden baseline ==="
"$PYTHON" -m pip install -q -e "$ROOT/python/.[dev]"
"$PYTHON" -m gut_ibm_tools.validation_regression \
  "$BUILD/$OUTPUT" \
  --golden "$GOLDEN"

echo "EARI/VADI validation regression passed."
