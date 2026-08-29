# GPU cost/benefit benchmark harness

This directory contains the authoritative Axis A/B/C arm manifest generator,
the one-arm runner, and the cross-machine report merger.  Generate configs
without executing an arm:

```bash
PYTHONPATH=../../python python generate_configs.py generate \
  --base base_config.json --output-dir generated \
  --scale 1e5 ../scaling_benchmark/input_1e5.json \
  --scale 1e6 ../scaling_benchmark/input_1e6.json
```

The generated `manifest.json` retains every arm, including blocked A6.  Run
one record on a suitable host with `run_arm.py`, then merge raw result files
from any hosts with `merge_report.py`.  Missing `(arm, scale, seed)` records
remain explicit in the merged report; the merger never starts GutIBM.

No benchmark arm is launched by the generator or this source tree's tests.

B1 uses the pre-fix duplicated-reflection series only as a cost reference. It
is explicitly marked as physically wrong and pins `image_series_max_shells`
to the historical three-shell value; legacy mode also maps an otherwise
unchanged corrected-mode shell default to those three shells.

## Cost and attribution passes

The runner executes every runnable `(arm, scale, seed)` twice. The
`profile_steps=false` pass is the only source for wall-clock and
steps-per-second cost columns. The `profile_steps=true` pass supplies
per-phase timing and kernel/table/direct-evaluation attribution columns. Both
passes, their configs, HDF5 files, and statuses remain distinguishable in the
raw result record; the merger never substitutes a profiled wall time for the
cost pass.

This separation is necessary because kernel counting is instrumentation. A
controlled one-million-call host workload measured a median `+5.37%` wall-time
overhead with counting enabled (13 million kernel evaluations), so mixing the
profiled pass into B2/B4 cost columns would bias the comparison toward higher
cost where evaluation counts are largest.
