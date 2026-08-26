# Delivery and ROS campaign method record

These text artifacts preserve the methods and reported measurements behind
PRs #330–#339. They are historical campaign records, not a portable harness.
The scripts use hardcoded `/home/ubuntu` paths and support the
`GUTIBM_BINARY` environment override where applicable. Raw HDF5 outputs and
copied binaries are deliberately excluded.

## Campaigns

### `q_o2_per_gen/`

Measures funded oxygen per generation and checks grid and density behavior.
The campaign ran at commit `ffab55b`, as recorded in its report.

### `regularized_support/`

Measures population-scale delivery behavior across grid resolutions with the
10 µm regularized support enabled and the radius-zero voxel-support control.
The campaign ran at commit `2a16cc7`, as recorded in its report.

### `ros_counterfactual/`

Measures ambient-ROS mortality against ROS-off controls, including SOS
component and accounting-ledger outcomes. The campaign ran at commit
`abf393750453673ebb1f4c4ab54f447797a647fd`, as recorded in its report.

The large ROS trajectory `metrics.json` is not copied because it exceeds the
approximately 1 MB artifact limit; the compact `metrics.csv` and reports are
preserved.
