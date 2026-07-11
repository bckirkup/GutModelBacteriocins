#!/usr/bin/env bash
# MPI scaling smoke for issues #154 / #156.
# Runs 2-rank and 4-rank MPI integration tests plus optional CUDA-aware path.
#
# Usage:
#   bash scripts/run_mpi_scaling_smoke.sh
#   BUILD=build-cuda bash scripts/run_mpi_scaling_smoke.sh
#   MPI_RANKS="2 4 8" bash scripts/run_mpi_scaling_smoke.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${BUILD:-$ROOT/build}"
TEST_BIN="${TEST_BIN:-$BUILD/tests}"
MPI_RANKS="${MPI_RANKS:-2 4}"
MPIRUN=(mpirun --allow-run-as-root)
export CC="${CC:-gcc}"
export CXX="${CXX:-g++}"

if [[ ! -x "$TEST_BIN/test_mpi_multi_rank" ]]; then
  echo "=== Configuring GutIBM (Release) ==="
  cmake -B "$BUILD" \
    -DCMAKE_BUILD_TYPE=Release \
    -DGUTIBM_USE_MPI=ON \
    -DGUTIBM_USE_HDF5=ON \
    -DGUTIBM_USE_OPENMP=OFF
  echo "=== Building MPI test targets ==="
  cmake --build "$BUILD" -j"$(nproc)" \
    --target test_mpi_multi_rank test_mpi_four_rank test_cuda_aware_mpi_reaction
fi

echo "=== MPI scaling smoke (BUILD=$BUILD) ==="

for ranks in $MPI_RANKS; do
  case "$ranks" in
    2)
      echo "--- mpirun -np 2 test_mpi_multi_rank ---"
      "${MPIRUN[@]}" -np 2 --oversubscribe "$TEST_BIN/test_mpi_multi_rank"
      if [[ -x "$TEST_BIN/test_cuda_aware_mpi_reaction" ]]; then
        echo "--- mpirun -np 2 test_cuda_aware_mpi_reaction ---"
        "${MPIRUN[@]}" -np 2 --oversubscribe "$TEST_BIN/test_cuda_aware_mpi_reaction"
      fi
      ;;
    4)
      if [[ -x "$TEST_BIN/test_mpi_four_rank" ]]; then
        echo "--- mpirun -np 4 test_mpi_four_rank ---"
        "${MPIRUN[@]}" -np 4 --oversubscribe "$TEST_BIN/test_mpi_four_rank"
      fi
      ;;
    *)
      echo "WARN: no dedicated test binary for np=$ranks; skipping"
      ;;
  esac
done

if [[ -x "$TEST_BIN/test_mpi_gpu_multi_rank" ]]; then
  echo "--- mpirun -np 2 test_mpi_gpu_multi_rank ---"
  "${MPIRUN[@]}" -np 2 --oversubscribe "$TEST_BIN/test_mpi_gpu_multi_rank"
fi

echo "MPI scaling smoke complete."
