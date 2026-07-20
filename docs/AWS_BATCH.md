# AWS Batch Deployment Plan (CUDA, Spot / low priority)

Plan for running production GutIBM jobs on AWS when desktop/WSL cannot finish
full campaigns (especially Stage 3: 7-day biology, 2 mm domain, GPU chemistry).

**Status:** planning + Phase 1 practice path (Jul 2026).
**Region (pinned):** `us-east-1` (co-locate ECR, S3, Batch).
**Primary workload (later):** `experiments/diversity_campaign/stage3_campaign/`.
**Practice workload (now):** `experiments/smoke_gpu.json` (+ optional
`smoke_gpu_batch.json`).

## Goals

| Goal | Detail |
|------|--------|
| Finish full runs | Stage 3 / diversity paradox scale that OOMs or crawls on a desktop |
| Low priority / cost | Prefer Spot interruptible capacity over On-Demand |
| CUDA | Production chemistry path with `gpu_enabled: true` |
| Batch semantics | Many independent sims (seeds × Kd × mechanisms), not one forever-job |
| Resume | Survive Spot reclaim via HDF5 checkpoint + S3 |
| Fit existing tools | Same local **Docker + AWS CLI → ECR** habit as other projects; reuse `batch_runner` JSON |

Non-goals for v1: multi-node MPI rings, NCCL multi-GPU, interactive login nodes,
Terraform/CDK, or replacing local Stage 1–2 validation.

## Decisions (locked / recommended)

| Topic | Choice | Notes |
|-------|--------|-------|
| Region | **`us-east-1`** | Locked |
| Infra style | **AWS CLI + Docker + checked-in JSON** (no Terraform for v1) | Matches local Crusher→ECR workflow; Batch resources created once via CLI |
| Scheduler | **AWS Batch** (EC2 GPU), not Fargate | Fargate has **no NVIDIA GPUs**; keep the Docker/ECR muscle memory, change only the compute backend |
| Practice instance | **`g4dn.xlarge`** (T4, 16 GB VRAM) | Cheapest common GPU Spot class for tiny smokes |
| Campaign instance | **`g5.2xlarge`** (A10G, 24 GB) | Stage 3 chemistry ~8 GB VRAM + headroom; allow `g4dn.2xlarge` in the same compute env if Spot is thinner on g5 |
| First jobs | **`experiments/smoke_gpu.json`** then `smoke_gpu_batch.json` | Prove CUDA path before Stage 3 |
| Image build | **Laptop `docker build` + `aws ecr get-login-password` push** | Same as existing Fargate deploys; GHA→ECR later if desired |
| First campaign (after smoke) | **`batch_baseline.json`** (3 seeds) | Smaller than the 12-run Kd sweep; validate cost/wall time once |

## Why not Fargate (even though you already use it)

Crusher_to_the_Bridge on Fargate is the right pattern for **CPU** containers:
build locally → push ECR → run managed tasks.

GutIBM CUDA needs host NVIDIA drivers + GPU device mapping. That is **Batch
managed EC2** with `ECS_AL2023_NVIDIA` (or equivalent GPU AMI), job definition
`resourceRequirements: [{type: GPU, value: 1}]`. Workflow stays familiar:

```
docker build → docker tag → aws ecr … push → aws batch submit-job
```

instead of `aws ecs run-task` / Fargate.

## Constraints from the code today

1. **Chemistry grid is per-rank full domain.** Stage 3 at `dx = 2 µm` on
   2 mm × 2 mm × 100 µm ≈ **50 M cells** → ~**8 GB host + 8 GB VRAM** chem
   alone (`experiments/diversity_campaign/README.md`).
2. **Stage 3 batches default to `mpi_ranks: 1`.** Prefer **1 GPU process per
   job** for v1.
3. **GPU FMM far-field still walks on CPU** for large trees — need host RAM,
   not a GPU-only tiny box.
4. **Checkpoint restart exists** (`checkpoint_file` / `checkpoint_step`).
5. Build CUDA arches for target GPUs: T4=`75`, A10G=`86` (image defaults to
   `75;86;89`).

## Recommended architecture (v1)

```
Laptop (AWS CLI + Docker)              AWS us-east-1
─────────────────────────              ────────────────────────────────
docker build deploy/aws/Dockerfile ──► ECR  …/gutibm:cuda
aws s3 cp smoke_gpu.json           ──► S3   gutibm-inputs-…
aws batch submit-job               ──► Batch Spot/OnDemand GPU queue
                                         │
                                    g4dn.xlarge (practice) /
                                    g5.2xlarge (campaign)
                                         │
                                    entry.sh → gut_ibm
                                         │
                                    S3 gutibm-outputs-… / CloudWatch Logs
```

| Piece | Choice | Why |
|-------|--------|-----|
| Scheduler | AWS Batch managed CE | GPU + Spot + array jobs; scales to 0 |
| Compute | EC2 Spot GPU (`SPOT_PRICE_CAPACITY_OPTIMIZED`) | Low priority / cost |
| Practice CE | Prefer `g4dn.xlarge` (optionally On-Demand for first green run) | Fast feedback, cheap fails |
| Campaign CE | `g5.2xlarge` + `g4dn.2xlarge` allowed | Stage 3 VRAM/RAM |
| Container | ECR + `deploy/aws/Dockerfile` | Same push path as Fargate apps |
| Storage | S3 prefixes in `us-east-1` | Spot-safe artifacts |
| IAM | Job role scoped to those prefixes | No keys in the image |
| IaC | None for v1 — CLI once, commit the JSON blobs under `deploy/aws/` | Avoid infra yak-shave before CUDA works |

## Instance sizing

### Practice / smoke (do this first)

| Instance | GPU | When |
|----------|-----|------|
| **`g4dn.xlarge`** | T4 16 GB | **Default for Phase 1–1b** (`smoke_gpu`, Stage 1-sized grids) |

`experiments/smoke_gpu.json` is ~20 agents, 100 µm domain, `dx = 5 µm` — chemistry
memory is negligible. Goal is “CUDA path ran on Batch,” not science.

### Full Stage 3 (later)

| Instance | GPU | Fit |
|----------|-----|-----|
| `g4dn.2xlarge` | T4 16 GB / 32 GB RAM | Usable Stage 3 candidate |
| **`g5.2xlarge`** | A10G 24 GB / 32 GB RAM | **Campaign default** |
| `g6.2xlarge` | L4 24 GB | If Spot is cheaper in `us-east-1` |

Compute environment tip: list **families or several sizes** and let Spot
allocation pick (`SPOT_PRICE_CAPACITY_OPTIMIZED`). Job definition still asks for
`GPU=1` + enough memory.

## Practice path (Phase 1) — concrete

Do **not** start with Stage 3. Order:

1. **Local CUDA binary smoke** (optional if you have a GPU laptop/WSL):

   ```bash
   ./rebuild_and_run.sh --cuda on --reuse-build --mode single \
     --config experiments/smoke_gpu.json --mpi-ranks 1
   ```

2. **Build & push image** (same muscle memory as Crusher→ECR).

   Prefer the paste-safe script (avoids fragile multi-line CLI / JSON):

   ```bash
   bash deploy/aws/01_push_image.sh
   ```

   Manual equivalent (quote `CUDA_ARCHS`; trailing spaces after `\` break pastes):

   ```bash
   export AWS_REGION=us-east-1
   ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
   REPO="${ACCOUNT}.dkr.ecr.${AWS_REGION}.amazonaws.com/gutibm"
   aws ecr create-repository --repository-name gutibm --region "$AWS_REGION" || true
   aws ecr get-login-password --region "$AWS_REGION" | docker login --username AWS --password-stdin "${ACCOUNT}.dkr.ecr.${AWS_REGION}.amazonaws.com"
   docker build -f deploy/aws/Dockerfile -t gutibm:cuda --build-arg 'CUDA_ARCHS=75;86;89' .
   docker tag gutibm:cuda "${REPO}:cuda"
   docker push "${REPO}:cuda"
   ```

3. **One-time Batch GPU stack** — prefer:

   ```bash
   bash deploy/aws/02_setup_practice_stack.sh
   ```

   That creates S3 buckets, IAM roles, On-Demand `g4dn.xlarge` CE
   (`ECS_AL2023_NVIDIA`), queue `gutibm-gpu-practice`, and job definition
   `gutibm-cuda`. Flip the CE to Spot after the first green run if desired.

4. **Upload smoke config & submit one job**:

   ```bash
   bash deploy/aws/03_submit_smoke.sh
   bash deploy/aws/04_watch_job.sh <jobId>
   ```

5. **Pass criteria for Phase 1:**
   - Job `SUCCEEDED`
   - CloudWatch / run log shows GPU init (not silent CPU fallback)
   - `output.h5.gz` lands in S3 and gunzips / opens with `gut_ibm_tools`

6. **Phase 1b — tiny array:** expand `smoke_gpu_batch.json` (2 seeds) to two
   S3 inputs; submit array size 2. Still on `g4dn.xlarge`.

Only after 1b is boring: move CE to include `g5.2xlarge` / Spot-heavy and try
one Stage 3 baseline seed.

## Job model (campaigns)

Prefer **Batch array index = one simulation** from a `batch_runner` manifest.
Do **not** run a multi-hour local `batch_runner` loop as one Spot task until
checkpoint→S3 resume is solid.

## Spot interruption strategy

1. Periodic HDF5 checkpoints → S3 prefix per job id.
2. `entry.sh` already sketches resume via `CHECKPOINT_S3_PREFIX`.
3. Batch `retryStrategy` for Spot reclaim once exit codes are stable.
4. Gzip finals before upload (existing habit).

For Phase 1 smokes, Spot reclaim is unlikely to matter (seconds–minutes). Use
**On-Demand `g4dn.xlarge` for the first green run** if Spot capacity is flaky;
flip the CE to Spot afterward.

## Container

See `deploy/aws/Dockerfile` and `entry.sh`. Multi-arch default
`CMAKE_CUDA_ARCHITECTURES=75;86;89` covers practice T4 and campaign A10G/L4.

## IAM and networking (minimal)

- Job role: S3 R/W on `gutibm-inputs-*` / `gutibm-outputs-*` prefixes only.
- Logs → CloudWatch (no SSH).
- S3 gateway endpoint in the VPC when HDF5 gets large (Stage 3).
- Everything in **`us-east-1`**.

## Cost control knobs

| Knob | Effect |
|------|--------|
| Practice on `g4dn.xlarge` first | Failures stay cheap |
| Spot after first green run | Primary savings |
| `minvCpus=0` | Scale to zero when idle |
| Array granularity | One bad seed ≠ whole sweep |
| HDF5 schedule | Agents/grid dumps dominate S3 |
| Gzip HDF5 | Storage + transfer |

## Phased delivery

### Phase 0 — Decide

- [x] Architecture: Batch + GPU EC2 + ECR + S3 (not Fargate)
- [x] Region: `us-east-1`
- [x] Infra style: CLI + Docker (no Terraform v1)
- [x] Practice instance: `g4dn.xlarge`
- [x] Campaign instance: `g5.2xlarge` (CE also allows `g4dn.2xlarge`)
- [x] Practice configs: `experiments/smoke_gpu.json`, `smoke_gpu_batch.json`

### Phase 1 — CUDA smoke on Batch (current focus)

- [x] Draft `deploy/aws/Dockerfile`, `entry.sh`
- [x] Paste-safe practice scripts (`deploy/aws/01`–`04` + README)
- [ ] Create ECR repo + push `gutibm:cuda` from laptop (`01_push_image.sh`)
- [ ] Create Batch GPU CE/queue/job definition in `us-east-1` (`02_setup_practice_stack.sh`)
- [ ] One On-Demand or Spot job with `smoke_gpu.json` (`03` + `04`)
- [ ] Confirm GPU path in logs + S3 output
- [ ] Phase 1b: array of 2 from `smoke_gpu_batch.json`

### Phase 2 — Single Stage 3 seed on `g5.2xlarge`

- [x] Campaign CE/queue/job-def sized for `g5.2xlarge` (one GPU per run) via `05_setup_campaign_stack.sh`
- [ ] One `3a_baseline` / `batch_baseline` seed (submit a single index of `batch_baseline.json`)
- [ ] Record wall time / $ in Measured section below

### Phase 3 — Array export from batch manifests

- [x] Helper: `batch_*.json` → S3 job tree + `submit-job` (`gut_ibm_tools.aws_batch_export` / `gut-ibm-aws-export`)
- [x] Parity with `batch_runner --dry-run` (shared `parse_batch_config` / `build_job_config`; unit-tested for `batch_kd_sweep.json` + `batch_baseline.json`)

### Phase 4 — Spot resilience

- [x] Checkpoint → S3 + auto-resume (`entry.sh` background sync of `output.h5` → `CHECKPOINT_S3_PREFIX`; resume half already existed)
- [x] Retry policy; Spot-only (or Spot-with-fallback) CE (`05` Spot CE `SPOT_CAPACITY_OPTIMIZED` + `retryStrategy` retrying `Host EC2*`; set `CAMPAIGN_ONDEMAND_FALLBACK=1` for an On-Demand fallback CE)

### Phase 5 — campaign deploy (concrete)

Region `us-east-1`. Assumes the image is already in ECR (`01_push_image.sh`).
Run from the repo root, one line at a time (paste-safe scripts). To use an
**existing bucket you own** instead of the derived `gutibm-inputs/outputs-<account>`
buckets, export `BUCKET` (shared) or `INPUT_BUCKET` / `OUTPUT_BUCKET` before
sourcing/running — differentiate with key prefixes, e.g.
`s3://my-bucket/gutibm/jobs`. Existing buckets are never recreated.

1. **Roles + bucket** (re-run `02`; idempotent, creates nothing you already own):

   ```bash
   bash deploy/aws/02_setup_practice_stack.sh
   ```

2. **Campaign Spot GPU stack** (CE + queue + job def, one GPU per run):

   ```bash
   bash deploy/aws/05_setup_campaign_stack.sh
   ```

   Creates `gutibm-gpu-campaign-spot` (Spot `g5.2xlarge`/`g4dn.2xlarge`,
   `SPOT_CAPACITY_OPTIMIZED`, `maxvCpus=${CAMPAIGN_MAX_VCPUS:-96}`), queue
   `gutibm-gpu-campaign`, and job def `gutibm-cuda-campaign` (8 vCPU / 28 GB /
   `GPU=1`, `retryStrategy` for Spot reclaims).

3. **Dry-run the array export** (no uploads, no submit — verify count + inputs):

   ```bash
   source deploy/aws/env.sh
   python -m gut_ibm_tools.aws_batch_export \
     experiments/diversity_campaign/stage3_campaign/batch_kd_sweep.json \
     --input-prefix "s3://${INPUT_BUCKET}/campaign/kd_sweep/jobs" \
     --output-prefix "s3://${OUTPUT_BUCKET}/campaign/kd_sweep/out" \
     --checkpoint-prefix "s3://${OUTPUT_BUCKET}/campaign/kd_sweep/ckpt" \
     --job-queue "${JOB_QUEUE_CAMPAIGN}" \
     --job-definition "${JOB_DEFINITION_CAMPAIGN}" \
     --dry-run
   ```

4. **Submit** (drop `--dry-run`): uploads per-index `input.json` and submits one
   array job of size 12 (`gpu_enabled:true`, `MPI_RANKS=1`, prefixes wired to
   `entry.sh`):

   ```bash
   python -m gut_ibm_tools.aws_batch_export \
     experiments/diversity_campaign/stage3_campaign/batch_kd_sweep.json \
     --input-prefix "s3://${INPUT_BUCKET}/campaign/kd_sweep/jobs" \
     --output-prefix "s3://${OUTPUT_BUCKET}/campaign/kd_sweep/out" \
     --checkpoint-prefix "s3://${OUTPUT_BUCKET}/campaign/kd_sweep/ckpt" \
     --job-queue "${JOB_QUEUE_CAMPAIGN}" \
     --job-definition "${JOB_DEFINITION_CAMPAIGN}"
   ```

5. **Watch** the printed array `jobId`:

   ```bash
   bash deploy/aws/04_watch_job.sh <jobId>
   ```

Outputs land at `s3://${OUTPUT_BUCKET}/campaign/kd_sweep/out/<index>/output.h5.gz`.
Start smaller with `batch_baseline.json` (3 runs) to measure cost/wall time first.

### Phase 6 — Optional later

- [ ] GitHub Actions → ECR
- [ ] Terraform/CDK if the CLI stack becomes painful to recreate
- [ ] Multi-rank MPI when chemistry is domain-decomposed

## Still open (low urgency)

1. **Checkpoint cadence** for 7-day runs (S3 PUT cost vs reclaim loss).
2. **Whether campaign CE is Spot-only** or 80/20 Spot/On-Demand after smoke.
3. **Bucket naming** (`gutibm-inputs-<account>` vs a shared research bucket).

## Relation to existing docs

| Doc | Role |
|-----|------|
| `experiments/smoke_gpu.json` | Phase 1 CUDA practice config |
| `experiments/diversity_campaign/README.md` | Memory floors, Stage 3 batches |
| `docs/BATCH_RUNNER.md` | Local sweep semantics to mirror as array jobs |
| `docs/SCALING.md` | Agent/grid scaling |
| `src/gpu/README.md` | Device kernels |
| `deploy/aws/` | Dockerfile, entrypoint, **paste-safe** `01`–`04` practice scripts |

## Measured results

_Fill in after Phase 1 / 2:_

| Config | Instance | Spot? | Wall time | $/run | Notes |
|--------|----------|-------|-----------|-------|-------|
| `smoke_gpu` | `g4dn.xlarge` | | | | Phase 1 |
| `3a_baseline` seed | `g5.2xlarge` | | | | Phase 2 |
