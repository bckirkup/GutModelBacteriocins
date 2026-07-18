# Batch Runner

Resumable CLI for parameter scans and experiment grids. Runs multiple `gut_ibm`
simulations from a single batch JSON manifest, with per-job isolation and
incremental progress tracking.

## Quick start

```bash
# From repo root, after building gut_ibm
cd python && pip install -e ".[dev]"

python -m gut_ibm_tools.batch_runner examples/batch_scan/batch.json --dry-run
python -m gut_ibm_tools.batch_runner examples/batch_scan/batch.json
python -m gut_ibm_tools.batch_runner examples/batch_scan/batch.json --resume
python -m gut_ibm_tools.batch_runner examples/batch_scan/batch.json --status
```

Installed entry point: `gut-ibm-batch` (same interface).

## Output layout

```
batch_results/my_scan/
  batch_manifest.json          # progress tracker (rewrite after each job)
  jobs/
    seed=4092_total_time=300/
      input.json               # generated simulation config
      output.h5                # HDF5 output (absolute path in config)
      # or output.h5.gz when GUTIBM_GZIP_HDF5=true (rebuild_and_run default)
      run.log                    # mpirun + gut_ibm stdout/stderr
```

When `GUTIBM_GZIP_HDF5` is true (exported by `rebuild_and_run.sh` unless
`--no-gzip-hdf5`), each successful job’s `output.h5` is replaced with a
whole-file `output.h5.gz` after validation. Gunzip before analysis.

## Batch JSON schema

Batch configs are strict JSON (use `_comment` for notes). They are **not**
simulation input files — see [CONFIG_FORMAT.md](CONFIG_FORMAT.md) for per-run
simulation JSON.

### Required keys

| Key | Description |
|-----|-------------|
| `output_dir` | Artifact root (manifest + `jobs/` tree) |
| `base_config` | Template simulation JSON (repo-relative or absolute) |
| `sweep` **or** `runs` | Job definition (mutually exclusive) |

### Optional keys

| Key | Default | Description |
|-----|---------|-------------|
| `binary` | `build/gut_ibm` | Path to `gut_ibm` executable |
| `mpi_ranks` | `1` | MPI process count |
| `mpirun` | `mpirun` | MPI launcher binary |
| `mpirun_args` | `[]` | Extra launcher flags (e.g. `["--allow-run-as-root"]`) |
| `on_fail` | `continue` | `continue` or `abort` |
| `env` | `{}` | Extra environment variables for subprocess |
| `validate` | omitted | Post-run validation block (see below) |

### Cartesian sweep

```json
{
  "output_dir": "batch_results/eari_kd_seed",
  "base_config": "examples/eari_vadi_validation/input.json",
  "binary": "build/gut_ibm",
  "mpi_ranks": 1,
  "sweep": {
    "seed": [4092, 4093, 4094],
    "kd_colicinE_btuB": [2e-10, 8e-10]
  }
}
```

Job count = product of sweep list lengths. Auto-generated `job_id` from sorted
parameter keys (e.g. `kd_colicinE_btuB=2e-10_seed=4092`).

### Explicit run list

```json
{
  "output_dir": "batch_results/custom",
  "base_config": "examples/diversity_paradox/input.json",
  "binary": "build/gut_ibm",
  "runs": [
    { "id": "baseline", "overrides": {} },
    { "id": "short", "overrides": { "total_time": 3600 } },
    { "id": "high_kd", "overrides": { "kd_colicinE_btuB": 1e-9 } }
  ]
}
```

### Override merge rules

1. Load `base_config`; strip `_`-prefixed documentation keys.
2. Apply per-job `overrides` with shallow top-level merge.
3. Support **dot-path keys** for nested fields:
   `"initial_strains.0.count": 200`
4. Force per-job `hdf5_file` to `{job_dir}/output.h5` (absolute path).

### Optional validation block

After a successful simulation, call `validation_regression.run_validation()`:

```json
{
  "validate": {
    "golden": "python/tests/fixtures/eari_vadi_ci_golden.json",
    "fish_golden": "python/tests/fixtures/eari_vadi_ci_fish_golden.json",
    "check_targets": false,
    "check_fish_targets": true
  }
}
```

Validation failures mark the job `failed` even when `gut_ibm` exits 0.

## Resume and interrupt safety

- Manifest is flushed after **every** job completes.
- `Ctrl+C` / `SIGTERM`: current job is marked `interrupted`, manifest saved, exit `130`.
- `--resume` skips `done` jobs; retries `pending`, `failed`, and `interrupted`.
- Re-running without `--resume` on an existing `output_dir` errors (prevents accidental overwrite).

## Status output

Progress lines go to **stderr** only:

```
[3/12] kd_colicinE_btuB=8e-10 seed=4092 ... OK (41.2s)
[4/12] kd_colicinE_btuB=8e-10 seed=4093 ... FAIL (exit 1)
^C interrupted after 4/12 (resume with --resume)
```

`--status` prints manifest summary: `3/12 done, 1 failed, 0 interrupted, 8 pending`.

## CI smoke

- **python-lint**: `--dry-run` on `examples/batch_scan/batch.json`
- **integration-tests**: `scripts/smoke_batch_runner.sh` (single 60 s job via `batch_ci.json`)

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | All jobs completed successfully |
| 1 | One or more job failures |
| 2 | Config or manifest error |
| 130 | Interrupted |
