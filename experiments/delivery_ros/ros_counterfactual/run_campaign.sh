#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/ubuntu/gutibm-campaign/ros-counterfactual"
BINARY="/home/ubuntu/repos/GutModelBacteriocins/build-serial/gut_ibm"

for arm in A_ros0_res2 A_ctrl_res2 B_ros0_res2 B_ctrl_res2 B_ros0_res6 B_ctrl_res6; do
  arm_dir="$ROOT/$arm"
  echo "=== Running $arm ==="
  (
    cd "$arm_dir"
    /usr/bin/time -f "wrapper_wall_seconds=%e" \
      "$BINARY" input.json 2>&1 | tee full.log
  )
  /home/ubuntu/.venvs/gutibm/bin/python3 "$ROOT/analyze_campaign.py" \
    --check-arm "$arm"
done

/home/ubuntu/.venvs/gutibm/bin/python3 "$ROOT/analyze_campaign.py"
