#!/usr/bin/env bash
# GPU scaling & A1-A6 uptake-mode benchmark driver (post-residency #480).
#
# Emits per-run s/step, H2D/D2H/total PCIe GB per 10 steps, pinned staging
# status, and peak RSS across configurable uptake modes, seeds, agent counts,
# MPI ranks, and host/device backends.
#
# Usage:
#   bash scripts/run_gpu_scaling_benchmark.sh
#   AGENT_COUNTS="100000" UPTAKE_MODES="none sherwood delivery" SEEDS="55 56 57" \
#     STEPS=10 bash scripts/run_gpu_scaling_benchmark.sh
#   GRID_X=0.001 GRID_Y=0.001 GRID_Z=0.0001 bash scripts/run_gpu_scaling_benchmark.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${BUILD:-$ROOT/build-gpu-scaling-bench}"
OUT_DIR="${OUT_DIR:-$ROOT/benchmark_results}"
STEPS="${STEPS:-10}"
SEEDS="${SEEDS:-55 56 57}"
AGENT_COUNTS="${AGENT_COUNTS:-100000}"
UPTAKE_MODES="${UPTAKE_MODES:-none sherwood delivery}"
MPI_RANKS="${MPI_RANKS:-1}"
CHEMISTRY_DECOMPOSITION="${CHEMISTRY_DECOMPOSITION:-replicated}"
USE_FMM_THRESHOLD="${USE_FMM_THRESHOLD:-50000}"
GPU_MODES="${GPU_MODES:-false true}"
USE_OPENMP="${USE_OPENMP:-ON}"
GRID_X="${GRID_X:-0.001}"
GRID_Y="${GRID_Y:-0.001}"
GRID_Z="${GRID_Z:-0.0001}"
GRID_DX="${GRID_DX:-2e-6}"
export CC="${CC:-gcc}"
export CXX="${CXX:-g++}"

MPIRUN=(mpirun --allow-run-as-root)

sec_to_ms() { awk -v s="${1:-0}" 'BEGIN { printf "%.3f", s * 1000 }'; }
bytes_to_gb_per_10() {
  awk -v b="${1:-0}" -v n="${2:-10}" \
    'BEGIN { if (n > 0) printf "%.4f", (b / 1.0e9) * (10.0 / n); else print "0.0000" }'
}

arm_name() {
  local mode="$1"
  local gpu="$2"
  case "${gpu}:${mode}" in
    false:none)     echo "A1" ;;
    false:sherwood) echo "A2" ;;
    false:delivery) echo "A3" ;;
    true:none)      echo "A4" ;;
    true:sherwood)  echo "A5" ;;
    true:delivery)  echo "A6" ;;
    *)              echo "custom" ;;
  esac
}

timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
csv="$OUT_DIR/gpu_scaling_${timestamp}.csv"
mkdir -p "$OUT_DIR"

echo "=== Configuring GutIBM (Release + CUDA + OpenMP=${USE_OPENMP}) ==="
cmake -B "$BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DGUTIBM_USE_MPI=ON \
  -DGUTIBM_USE_HDF5=ON \
  -DGUTIBM_USE_CUDA=ON \
  -DGUTIBM_USE_OPENMP="${USE_OPENMP}"

echo "=== Building gut_ibm ==="
cmake --build "$BUILD" -j"$(nproc)" --target gut_ibm

write_config() {
  local agents="$1"
  local use_fmm="$2"
  local gpu="$3"
  local uptake="$4"
  local seed="$5"
  local out="$6"
  local total_time_s
  total_time_s="$(awk -v n="$STEPS" 'BEGIN { print n * 60 }')"
  cat >"$out" <<EOF
{
  "total_time": $total_time_s,
  "bio_dt": 60,
  "output_interval": 1.0e9,
  "seed": $seed,
  "profile_steps": true,
  "gpu_enabled": $gpu,
  "chemistry_decomposition": "$CHEMISTRY_DECOMPOSITION",
  "uptake_limit": "$uptake",
  "domain_x": $GRID_X,
  "domain_y": $GRID_Y,
  "domain_z": $GRID_Z,
  "grid_dx": $GRID_DX,
  "toxin_cutoff": 200e-6,
  "nutrient_cutoff": 50e-6,
  "use_fmm": $use_fmm,
  "fmm_theta": 0.5,
  "fmm_expansion_order": 2,
  "siderophore.enabled": false,
  "initial_strains": [
    {
      "type": 1,
      "count": $agents,
      "mu_max": 5.5e-4,
      "plasmids": ["ColE1", "ColB"],
      "conjugative": true
    }
  ],
  "hdf5_file": "gpu_scaling_bench.h5",
  "hdf5": {
    "enabled": false
  }
}
EOF
}

echo "timestamp,arm,seed,uptake_mode,chemistry_decomposition,agents,ranks,use_fmm,gpu_enabled,wall_s,step_s,step_ms,chemistry_ms,biology_ms,gpu_h2d_ms,gpu_d2h_ms,h2d_gb_per_10steps,d2h_gb_per_10steps,pcie_gb_per_10steps,pinned_staging_active,mpi_reaction_reduce_ms,rss_mb" >"$csv"

for agents in $AGENT_COUNTS; do
  use_fmm="false"
  if [[ "$agents" -ge "$USE_FMM_THRESHOLD" ]]; then
    use_fmm="true"
  fi

  for ranks in $MPI_RANKS; do
    for uptake in $UPTAKE_MODES; do
      for gpu_flag in $GPU_MODES; do
        arm="$(arm_name "$uptake" "$gpu_flag")"
        for seed in $SEEDS; do
          cfg="$BUILD/gpu_scaling_${arm}_${agents}_r${ranks}_s${seed}.json"
          write_config "$agents" "$use_fmm" "$gpu_flag" "$uptake" "$seed" "$cfg"
          log="$OUT_DIR/gpu_run_${arm}_${agents}_r${ranks}_s${seed}_${timestamp}.log"

          echo "=== arm=$arm uptake=$uptake seed=$seed agents=$agents ranks=$ranks gpu=$gpu_flag ==="
          rm -f "$BUILD/gpu_scaling_bench.h5"

          /usr/bin/time -f "TIME wall_s=%e maxrss_kb=%M" -o "$log.time" \
            "${MPIRUN[@]}" -np "$ranks" "$BUILD/gut_ibm" "$cfg" \
            >"$log" 2>&1 || {
              echo "WARN: run failed for arm=$arm seed=$seed agents=$agents ranks=$ranks (see $log)" >&2
              continue
            }

          wall_s="$(awk -F= '/^TIME wall_s=/ {print $2}' "$log.time" | awk '{print $1}')"
          rss_kb="$(awk '/^TIME wall_s=/ {sub(/^maxrss_kb=/, "", $NF); print $NF}' "$log.time")"
          rss_mb="$(awk -v kb="${rss_kb:-0}" 'BEGIN { if (kb == "" || kb <= 0) print -1; else printf "%.1f", kb/1024 }')"

          profile_line="$(grep '^PROFILE_CSV' "$log" | tail -1 || true)"
          step_s="$(echo "$profile_line" | sed -n 's/.* total_s=\([^ ]*\).*/\1/p')"
          chemistry_s="$(echo "$profile_line" | sed -n 's/.* chemistry_s=\([^ ]*\).*/\1/p')"
          biology_s="$(echo "$profile_line" | sed -n 's/.* biology_s=\([^ ]*\).*/\1/p')"
          gpu_h2d_s="$(echo "$profile_line" | sed -n 's/.* gpu_h2d_s=\([^ ]*\).*/\1/p')"
          gpu_d2h_s="$(echo "$profile_line" | sed -n 's/.* gpu_d2h_s=\([^ ]*\).*/\1/p')"
          gpu_h2d_bytes="$(echo "$profile_line" | sed -n 's/.* gpu_h2d_bytes_total=\([^ ]*\).*/\1/p')"
          gpu_d2h_bytes="$(echo "$profile_line" | sed -n 's/.* gpu_d2h_bytes_total=\([^ ]*\).*/\1/p')"
          pinned_active="$(echo "$profile_line" | sed -n 's/.* gpu_pinned_staging_active=\([^ ]*\).*/\1/p')"
          mpi_reduce_s="$(echo "$profile_line" | sed -n 's/.* mpi_reaction_reduce_s=\([^ ]*\).*/\1/p')"

          if [[ -n "$step_s" ]]; then
            step_ms="$(sec_to_ms "$step_s")"
            chemistry_ms="$(sec_to_ms "${chemistry_s:-0}")"
            biology_ms="$(sec_to_ms "${biology_s:-0}")"
            gpu_h2d_ms="$(sec_to_ms "${gpu_h2d_s:-0}")"
            gpu_d2h_ms="$(sec_to_ms "${gpu_d2h_s:-0}")"
            mpi_reduce_ms="$(sec_to_ms "${mpi_reduce_s:-0}")"
          else
            step_s="$(awk -v w="$wall_s" -v n="$STEPS" 'BEGIN { if (n > 0) printf "%.6f", w/n; else print 0 }')"
            step_ms="$(sec_to_ms "$step_s")"
            chemistry_ms=""
            biology_ms=""
            gpu_h2d_ms=""
            gpu_d2h_ms=""
            mpi_reduce_ms=""
          fi

          h2d_gb_10="$(bytes_to_gb_per_10 "${gpu_h2d_bytes:-0}" "$STEPS")"
          d2h_gb_10="$(bytes_to_gb_per_10 "${gpu_d2h_bytes:-0}" "$STEPS")"
          total_bytes="$(awk -v a="${gpu_h2d_bytes:-0}" -v b="${gpu_d2h_bytes:-0}" 'BEGIN { printf "%.0f", a + b }')"
          pcie_gb_10="$(bytes_to_gb_per_10 "$total_bytes" "$STEPS")"
          pinned_active="${pinned_active:-0}"

          echo "$timestamp,$arm,$seed,$uptake,$CHEMISTRY_DECOMPOSITION,$agents,$ranks,$use_fmm,$gpu_flag,$wall_s,$step_s,$step_ms,$chemistry_ms,$biology_ms,$gpu_h2d_ms,$gpu_d2h_ms,$h2d_gb_10,$d2h_gb_10,$pcie_gb_10,$pinned_active,$mpi_reduce_ms,$rss_mb" >>"$csv"
          echo "  wall_s=$wall_s step_s=$step_s pcie_gb_per_10steps=$pcie_gb_10 pinned=$pinned_active rss_mb=$rss_mb"
        done
      done
    done
  done
done

echo "GPU scaling benchmark complete. Results: $csv"
