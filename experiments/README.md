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

The interactive menu groups configs by subdirectory (so
`diversity_campaign/stage1_*` … `stage3_*` appear in order) and can run an
entire diversity-campaign stage sequentially.

After each successful single-run or batch job, `rebuild_and_run.sh` **gzip**s
the HDF5 output to `*.h5.gz` and removes the uncompressed file (default; disable
with `--no-gzip-hdf5`). This is whole-file gzip — independent of HDF5-internal
`hdf5.compression: gzip` (which only affects grid layers). Gunzip before
analysis:

```bash
gunzip -k path/to/output.h5.gz   # keep .gz, or omit -k to replace
```

Starter files:

- `smoke_single.json` — one small, 60-second CPU simulation.
- `smoke_batch.json` — two explicit runs based on the smoke config.
- `diversity_campaign/` — three-stage plan (motility validation → mechanism
  validation → full 7-day campaign + Kd sweep). See
  `diversity_campaign/README.md`.

To use the GPU, copy a simulation config and set:

```json
{
  "gpu_enabled": true,
  "gpu_device_id": -1
}
```

`gpu_device_id: -1` lets MPI ranks select devices by rank. A CUDA-enabled build
and visible NVIDIA GPU are still required. Stage 3 campaign configs already
enable GPU; Stages 1–2 stay on CPU.

Batch manifests resolve `base_config`, `binary`, and `output_dir` from the
repository root. Keep `"binary": "build/gut_ibm"` when using the helper's
default build directory.
