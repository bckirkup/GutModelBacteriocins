#!/usr/bin/env bash
# Compare OpenMP OFF vs ON simulation fingerprints (issue #50).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_SERIAL="$ROOT/build-openmp-off"
BUILD_OMP="$ROOT/build-openmp-on"
export CC=gcc
export CXX=g++

common_cmake() {
  cmake -B "$1" \
    -DCMAKE_BUILD_TYPE=Release \
    -DGUTIBM_USE_MPI=ON \
    -DGUTIBM_USE_HDF5=ON \
    -DGUTIBM_USE_OPENMP="$2"
}

echo "=== Building serial (OpenMP OFF) ==="
common_cmake "$BUILD_SERIAL" OFF
cmake --build "$BUILD_SERIAL" -j"$(nproc)" --target test_openmp_parity

echo "=== Building OpenMP (OpenMP ON) ==="
common_cmake "$BUILD_OMP" ON
cmake --build "$BUILD_OMP" -j"$(nproc)" --target test_openmp_parity

FP_SERIAL="$("$BUILD_SERIAL/tests/test_openmp_parity" | awk '/^FINGERPRINT=/ {print $1}')"
FP_OMP="$("$BUILD_OMP/tests/test_openmp_parity" | awk '/^FINGERPRINT=/ {print $1}')"

echo "Serial: $FP_SERIAL"
echo "OpenMP: $FP_OMP"

if [[ "$FP_SERIAL" != "$FP_OMP" ]]; then
  echo "ERROR: OpenMP fingerprint mismatch" >&2
  exit 1
fi

echo "OpenMP parity check passed."
