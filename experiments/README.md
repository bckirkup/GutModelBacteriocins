# Experiments

Place local simulation JSON files and batch manifests in this directory. The
top-level helper discovers JSON recursively:

```bash
./rebuild_and_run.sh
```

Files are classified by structure:

- **Single run:** a normal GutIBM simulation config.
- **Batch run:** an object containing `base_config` and exactly one of `sweep`
  or `runs`.

Starter files:

- `smoke_single.json` — one small, 60-second CPU simulation.
- `smoke_batch.json` — two explicit runs based on the smoke config.
- `diversity_campaign/` — full-scale baseline, mechanism, and corrinoid-Kd
  experiments with three-seed batch manifests.

To use the GPU, copy a simulation config and set:

```json
{
  "gpu_enabled": true,
  "gpu_device_id": -1
}
```

`gpu_device_id: -1` lets MPI ranks select devices by rank. A CUDA-enabled build
and visible NVIDIA GPU are still required.

Batch manifests resolve `base_config`, `binary`, and `output_dir` from the
repository root. Keep `"binary": "build/gut_ibm"` when using the helper's
default build directory.
