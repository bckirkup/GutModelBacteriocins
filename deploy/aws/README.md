# AWS Batch (CUDA Spot) — scaffolding

Planning doc: [`docs/AWS_BATCH.md`](../../docs/AWS_BATCH.md).

**Pinned:** region `us-east-1`. Infra style = laptop **Docker + AWS CLI**
(same habit as Fargate/ECR deploys). Compute backend is **Batch EC2 GPU**, not
Fargate (no NVIDIA GPUs on Fargate).

| File | Status |
|------|--------|
| `Dockerfile` | Draft multi-stage CUDA + MPI + HDF5 build |
| `entry.sh` | Draft S3 download → `gut_ibm` → upload |
| `submit_array_example.sh` | Example array submit (after smoke works) |

## Practice first

Use [`experiments/smoke_gpu.json`](../../experiments/smoke_gpu.json) on
**`g4dn.xlarge`** before any Stage 3 job. See the Practice path section in
`docs/AWS_BATCH.md`.
