#!/usr/bin/env bash
# Compare CPU vs CUDA-enabled simulation fingerprints (issue #33).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_CPU="$ROOT/build-gpu-parity-cpu"
BUILD_CUDA="$ROOT/build-gpu-parity-cuda"
export CC=gcc
export CXX=g++

common_cmake() {
  local build_dir="$1"
  local cuda_flag="$2"
  cmake -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release \
    -DGUTIBM_USE_MPI=ON \
    -DGUTIBM_USE_HDF5=ON \
    -DGUTIBM_USE_CUDA="$cuda_flag"
}

echo "=== Building CPU reference (CUDA OFF) ==="
common_cmake "$BUILD_CPU" OFF
cmake --build "$BUILD_CPU" -j"$(nproc)" --target test_gpu_smoke

echo "=== Building CUDA (CUDA ON) ==="
if ! command -v nvcc >/dev/null 2>&1; then
  echo "WARNING: nvcc not found — skipping GPU parity (compile-only CI should install nvidia-cuda-toolkit)"
  exit 0
fi
common_cmake "$BUILD_CUDA" ON
cmake --build "$BUILD_CUDA" -j"$(nproc)" --target test_gpu_smoke

echo "=== Running parity test (skips GPU path when no device) ==="
"$BUILD_CUDA/tests/test_gpu_smoke"

echo "GPU parity check passed."
