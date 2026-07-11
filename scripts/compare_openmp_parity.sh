#!/usr/bin/env bash
# Compare OpenMP OFF vs ON simulation fingerprints (issues #50, #161).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_SERIAL="$ROOT/build-openmp-off"
BUILD_OMP="$ROOT/build-openmp-on"
export CC=gcc
export CXX=g++

common_cmake() {
  local build_dir="$1"
  local openmp_flag="$2"
  cmake -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release \
    -DGUTIBM_USE_MPI=ON \
    -DGUTIBM_USE_HDF5=ON \
    -DGUTIBM_USE_OPENMP="$openmp_flag"
}

echo "=== Building serial (OpenMP OFF) ==="
common_cmake "$BUILD_SERIAL" OFF
cmake --build "$BUILD_SERIAL" -j"$(nproc)" --target test_openmp_parity

echo "=== Building OpenMP (OpenMP ON) ==="
common_cmake "$BUILD_OMP" ON
cmake --build "$BUILD_OMP" -j"$(nproc)" --target test_openmp_parity

FP_SERIAL="$("$BUILD_SERIAL/tests/test_openmp_parity" | awk '/^FINGERPRINT=/ {print $1}')"
FP_OMP="$("$BUILD_OMP/tests/test_openmp_parity" | awk '/^FINGERPRINT=/ {print $1}')"
FP_STOCH_SERIAL="$("$BUILD_SERIAL/tests/test_openmp_parity" | awk '/^FINGERPRINT_STOCHASTIC=/ {print $1}')"
FP_STOCH_OMP="$("$BUILD_OMP/tests/test_openmp_parity" | awk '/^FINGERPRINT_STOCHASTIC=/ {print $1}')"

echo "Serial: $FP_SERIAL"
echo "OpenMP: $FP_OMP"
echo "Stochastic serial: $FP_STOCH_SERIAL"
echo "Stochastic OpenMP: $FP_STOCH_OMP"

if [[ "$FP_SERIAL" != "$FP_OMP" ]]; then
  echo "ERROR: OpenMP fingerprint mismatch" >&2
  exit 1
fi

if [[ "$FP_STOCH_SERIAL" != "$FP_STOCH_OMP" ]]; then
  echo "ERROR: OpenMP stochastic toxin-kill fingerprint mismatch" >&2
  exit 1
fi

echo "OpenMP parity check passed (deterministic + stochastic toxin-kill)."
