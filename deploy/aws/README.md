# AWS Batch (CUDA Spot) — scaffolding

Planning doc: [`docs/AWS_BATCH.md`](../../docs/AWS_BATCH.md).

This directory will hold the container entrypoint and job helpers for
**Phase 1+** of that plan. Files here are intentionally thin until the region
and default instance type are pinned.

| File | Status |
|------|--------|
| `Dockerfile` | Draft multi-stage CUDA + MPI + HDF5 build |
| `entry.sh` | Draft S3 download → `gut_ibm` → upload |
| `submit_array_example.sh` | Example `aws batch submit-job` (array) |

Do not treat these as production until Phase 1 smoke on Batch succeeds.
