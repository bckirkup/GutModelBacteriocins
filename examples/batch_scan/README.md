# Batch Parameter Scan Examples

Starter manifests for the Python batch runner (`gut_ibm_tools.batch_runner`).

| File | Purpose |
|------|---------|
| `batch.json` | 2×2 sweep over `seed` and `total_time` on the CI validation base config |
| `batch_ci.json` | Single 60 s job used by `scripts/smoke_batch_runner.sh` in CI |

## Run locally

```bash
# Preview expanded jobs (no build required)
python -m gut_ibm_tools.batch_runner examples/batch_scan/batch.json --dry-run

# Run after building gut_ibm
mkdir -p build && cd build
cmake .. -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON && make -j$(nproc) gut_ibm
cd ..
python -m gut_ibm_tools.batch_runner examples/batch_scan/batch.json
```

Results land in `batch_results/batch_scan/` with a resumable `batch_manifest.json`.

See [docs/BATCH_RUNNER.md](../../docs/BATCH_RUNNER.md) for the full schema.
