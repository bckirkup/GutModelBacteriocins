#!/usr/bin/env bash
# GPU scaling benchmark driver — CPU vs GPU chemistry/sync comparison.
#
# Usage:
#   bash scripts/run_gpu_scaling_benchmark.sh
#   AGENT_COUNTS="10000 50000" bash scripts/run_gpu_scaling_benchmark.sh
#   STEPS=5 GPU_ENABLED=true bash scripts/run_gpu_scaling_benchmark.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${BUILD:-$ROOT/build-gpu-scaling-bench}"
OUT_DIR="${OUT_DIR:-$ROOT/benchmark_results}"
STEPS="${STEPS:-5}"
AGENT_COUNTS="${AGENT_COUNTS:-10000 50000 100000}"
MPI_RANKS="${MPI_RANKS:-1}"
USE_FMM_THRESHOLD="${USE_FMM_THRESHOLD:-50000}"
GPU_ENABLED="${GPU_ENABLED:-true}"
COMPARE_CPU="${COMPARE_CPU:-true}"
export CC="${CC:-gcc}"
export CXX="${CXX:-g++}"

MPIRUN=(mpirun --allow-run-as-root)

sec_to_ms() { awk -v s="${1:-0}" 'BEGIN { printf "%.3f", s * 1000 }'; }

timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
csv="$OUT_DIR/gpu_scaling_${timestamp}.csv"
mkdir -p "$OUT_DIR"

echo "=== Configuring GutIBM (Release + CUDA) ==="
cmake -B "$BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DGUTIBM_USE_MPI=ON \
  -DGUTIBM_USE_HDF5=ON \
  -DGUTIBM_USE_CUDA=ON \
  -DGUTIBM_USE_OPENMP=OFF

echo "=== Building gut_ibm ==="
cmake --build "$BUILD" -j"$(nproc)" --target gut_ibm

write_config() {
  local agents="$1"
  local use_fmm="$2"
  local gpu="$3"
  local out="$4"
  cat >"$out" <<EOF
{
  "total_time": $((STEPS * 60)),
  "bio_dt": 60,
  "output_interval": 1.0e9,
  "seed": 55,
  "profile_steps": true,
  "gpu_enabled": $gpu,
  "domain_x": 0.002,
  "domain_y": 0.002,
  "domain_z": 0.0001,
  "grid_dx": 2e-6,
  "toxin_cutoff": 200e-6,
  "nutrient_cutoff": 50e-6,
  "use_fmm": $use_fmm,
  "fmm_theta": 0.5,
  "fmm_expansion_order": 2,
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

echo "timestamp,agents,ranks,use_fmm,gpu_enabled,wall_s,step_ms,chemistry_ms,biology_ms,gpu_h2d_ms,gpu_d2h_ms,mpi_reaction_reduce_ms,rss_mb" >"$csv"

for agents in $AGENT_COUNTS; do
  use_fmm="false"
  if [[ "$agents" -ge "$USE_FMM_THRESHOLD" ]]; then
    use_fmm="true"
  fi

  for ranks in $MPI_RANKS; do
    for gpu_flag in $GPU_ENABLED; do
      if [[ "$COMPARE_CPU" == "true" && "$gpu_flag" == "true" ]]; then
        cfg_cpu="$BUILD/gpu_scaling_${agents}_r${ranks}_cpu.json"
        write_config "$agents" "$use_fmm" "false" "$cfg_cpu"
        log_cpu="$OUT_DIR/gpu_run_${agents}_r${ranks}_cpu_${timestamp}.log"
        echo "=== agents=$agents ranks=$ranks gpu=false ==="
        /usr/bin/time -f "TIME wall_s=%e maxrss_kb=%M" -o "$log_cpu.time" \
          "${MPIRUN[@]}" -np "$ranks" "$BUILD/gut_ibm" "$cfg_cpu" \
          >"$log_cpu" 2>&1 || echo "WARN: CPU run failed (see $log_cpu)" >&2
      fi

      cfg="$BUILD/gpu_scaling_${agents}_r${ranks}_gpu${gpu_flag}.json"
      write_config "$agents" "$use_fmm" "$gpu_flag" "$cfg"
      log="$OUT_DIR/gpu_run_${agents}_r${ranks}_gpu${gpu_flag}_${timestamp}.log"

      echo "=== agents=$agents ranks=$ranks gpu=$gpu_flag ==="
      rm -f "$BUILD/gpu_scaling_bench.h5"

      /usr/bin/time -f "TIME wall_s=%e maxrss_kb=%M" -o "$log.time" \
        "${MPIRUN[@]}" -np "$ranks" "$BUILD/gut_ibm" "$cfg" \
        >"$log" 2>&1 || {
          echo "WARN: run failed for agents=$agents ranks=$ranks gpu=$gpu_flag (see $log)" >&2
          continue
        }

      wall_s="$(awk -F= '/^TIME wall_s=/ {print $2}' "$log.time" | awk '{print $1}')"
      rss_kb="$(awk -F= '/maxrss_kb=/ {print $2}' "$log.time")"
      rss_mb="$(awk -v kb="${rss_kb:-0}" 'BEGIN { if (kb == "" || kb <= 0) print -1; else printf "%.1f", kb/1024 }')"

      profile_line="$(grep '^PROFILE_CSV' "$log" | tail -1 || true)"
      step_ms="$(echo "$profile_line" | sed -n 's/.* total_s=\([^ ]*\).*/\1/p')"
      chemistry_ms="$(echo "$profile_line" | sed -n 's/.* chemistry_s=\([^ ]*\).*/\1/p')"
      biology_ms="$(echo "$profile_line" | sed -n 's/.* biology_s=\([^ ]*\).*/\1/p')"
      gpu_h2d_ms="$(echo "$profile_line" | sed -n 's/.* gpu_h2d_s=\([^ ]*\).*/\1/p')"
      gpu_d2h_ms="$(echo "$profile_line" | sed -n 's/.* gpu_d2h_s=\([^ ]*\).*/\1/p')"
      mpi_reduce_ms="$(echo "$profile_line" | sed -n 's/.* mpi_reaction_reduce_s=\([^ ]*\).*/\1/p')"

      if [[ -n "$step_ms" ]]; then
        step_ms="$(sec_to_ms "$step_ms")"
        chemistry_ms="$(sec_to_ms "${chemistry_ms:-0}")"
        biology_ms="$(sec_to_ms "${biology_ms:-0}")"
        gpu_h2d_ms="$(sec_to_ms "${gpu_h2d_ms:-0}")"
        gpu_d2h_ms="$(sec_to_ms "${gpu_d2h_ms:-0}")"
        mpi_reduce_ms="$(sec_to_ms "${mpi_reduce_ms:-0}")"
      else
        step_ms="$(awk -v w="$wall_s" -v n="$STEPS" 'BEGIN { if (n > 0) printf "%.3f", (w/n)*1000; else print 0 }')"
        chemistry_ms=""
        biology_ms=""
        gpu_h2d_ms=""
        gpu_d2h_ms=""
        mpi_reduce_ms=""
      fi

      echo "$timestamp,$agents,$ranks,$use_fmm,$gpu_flag,$wall_s,$step_ms,$chemistry_ms,$biology_ms,$gpu_h2d_ms,$gpu_d2h_ms,$mpi_reduce_ms,$rss_mb" >>"$csv"
      echo "  wall_s=$wall_s step_ms=$step_ms chemistry_ms=$chemistry_ms gpu_d2h_ms=$gpu_d2h_ms"
    done
  done
done

echo "GPU scaling benchmark complete. Results: $csv"
