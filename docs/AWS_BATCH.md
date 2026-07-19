# AWS Batch Deployment Plan (CUDA, Spot / low priority)

Plan for running production GutIBM jobs on AWS when desktop/WSL cannot finish
full campaigns (especially Stage 3: 7-day biology, 2 mm domain, GPU chemistry).

**Status:** planning document (Jul 2026). No Terraform/CloudFormation yet.
**Primary workload:** `experiments/diversity_campaign/stage3_campaign/` batches
and similar parameter scans via `gut_ibm_tools.batch_runner`.

## Goals

| Goal | Detail |
|------|--------|
| Finish full runs | Stage 3 / diversity paradox scale that OOMs or crawls on a desktop |
| Low priority / cost | Prefer Spot interruptible capacity over On-Demand |
| CUDA | Production chemistry path with `gpu_enabled: true` |
| Batch semantics | Many independent sims (seeds × Kd × mechanisms), not one forever-job |
| Resume | Survive Spot reclaim via HDF5 checkpoint + S3 |
| Fit existing tools | Reuse `batch_runner`, experiment JSON, gzip HDF5 habit |

Non-goals for v1: multi-node MPI rings, NCCL multi-GPU, interactive login nodes,
or replacing local Stage 1–2 validation (those stay on laptop/CI).

## Constraints from the code today

These drive instance and job shape more than AWS product choice.

1. **Chemistry grid is per-rank full domain.** There is no slab-local chemistry
   yet. Stage 3 at `dx = 2 µm` on 2 mm × 2 mm × 100 µm ≈ **50 M cells** →
   ~**8 GB host + 8 GB VRAM** for chemistry alone (`experiments/diversity_campaign/README.md`).
2. **Stage 3 batches default to `mpi_ranks: 1`.** Raising ranks multiplies host
   RAM without splitting the grid. Prefer **1 GPU process per job** for v1.
3. **GPU FMM far-field still walks on CPU** for large trees. Expect hybrid
   CPU+GPU load; do not pick GPU-only tiny instances with starved host RAM.
4. **Checkpoint restart exists** (`checkpoint_file` / `checkpoint_step` in
   input JSON; `hdf5_checkpoint` CTest). Use it for Spot reclaim.
5. **CUDA arches in CMake default to `60;70;80`.** Rebuild image with the
   compute capability of the chosen instance family (see below).
6. **Batch runner is local-process oriented.** It launches `mpirun`/`gut_ibm`
   and writes under `output_dir`. On AWS, either wrap one job per Batch array
   index, or run `batch_runner --resume` on a long-lived worker that checkpoints
   to S3. Prefer **one sim per Batch job** for Spot.

## Recommended architecture (v1)

```
Desktop / CI                  AWS
─────────────                 ─────────────────────────────────────
git push ──► build image ──► Amazon ECR
                              │
experiment JSON + batch ──► S3://gutibm-inputs/...
manifest (optional)           │
                              ▼
                         AWS Batch
                         (Spot GPU queue, low priority)
                              │
                    Array job: one index = one sim
                              │
                    Docker: gut_ibm (+ optional mpirun -np 1)
                              │
                    checkpoints + output.h5(.gz) ──► S3://gutibm-outputs/...
                              │
                         Athena / local download
                         + python gut_ibm_tools analysis
```

| Piece | Choice | Why |
|-------|--------|-----|
| Scheduler | **AWS Batch** managed queue | Native array jobs, Spot, retries, no cluster daemon to babysit |
| Compute | **EC2 Spot GPU** via Batch compute environment | “Low priority” = interruptible capacity |
| Container | **ECR** image: CUDA runtime + OpenMPI + parallel HDF5 + `gut_ibm` | Matches CI stack; portable |
| Storage | **S3** for inputs, checkpoints, final HDF5 | Cheap, Spot-safe; EFS optional later |
| Secrets | IAM role on Batch job (no long-lived keys in image) | Least privilege to S3 prefixes |
| Orchestration glue | Thin `scripts/aws_batch_entry.sh` + optional SubmitJob helper | Keep science configs in `experiments/` |

Defer **ParallelCluster / Slurm** until you need `mpirun -np 8+` across nodes
or CUDA-aware multi-rank rings. Stage 3’s memory model does not benefit from
multi-rank today.

## Instance sizing (GPU Spot)

Chemistry VRAM floor (~8 GB) + scratch ⇒ **≥12 GB VRAM practical, 16–24 GB recommended**
for full Stage 3 (`dx = 2 µm`). Host needs **≥24 GB RAM** comfortable (16 GB is
tight once temps/HDF5 buffers appear).

| Instance (typical) | GPU | VRAM | vCPU / RAM | Fit |
|--------------------|-----|------|------------|-----|
| `g4dn.xlarge` | T4 | 16 GB | 4 / 16 GB | Marginal host RAM; OK for `dx≥4 µm` or 1 mm domain |
| `g4dn.2xlarge` | T4 | 16 GB | 8 / 32 GB | Good Stage 3 candidate (T4 = CC 7.5) |
| `g5.xlarge` | A10G | 24 GB | 4 / 16 GB | VRAM great; host RAM tight |
| `g5.2xlarge` | A10G | 24 GB | 8 / 32 GB | **Preferred default** for full Stage 3 |
| `g6.2xlarge` | L4 | 24 GB | 8 / 32 GB | Newer alternative if Spot is cheaper in-region |

**v1 default:** Batch compute environment allowing `g5.2xlarge` + `g4dn.2xlarge`
(Spot), single GPU, 1 job per instance (or 1 GPU job definition with
`resourceRequirements: GPU=1`).

Build CUDA arch:

| Family | `CMAKE_CUDA_ARCHITECTURES` |
|--------|----------------------------|
| g4dn (T4) | `75` |
| g5 (A10G) | `86` |
| g6 (L4) | `89` |

Ship a multi-arch image (`75;86;89`) if you want one ECR tag for all three.

Laptop-safe configs (`grid_dx` 4–5 µm) can use cheaper / smaller Spot if needed;
keep full-resolution jobs on the 24 GB VRAM class.

## Job model

### Prefer: AWS Batch array job = one simulation

Map each `batch_runner` job (e.g. `seed=42_kd_corrinoid_btuB=1e-09`) to one
array index:

1. Upload base config + overrides (or a generated `input.json`) to S3.
2. Submit `gutibm-stage3` array job size = N.
3. Entry script downloads index-specific JSON, runs
   `./gut_ibm /tmp/input.json` (or `mpirun -np 1`), uploads `output.h5.gz`
   and periodic checkpoints.

Pros: Spot kill loses one sim, not the whole sweep; natural Batch retries;
matches `on_fail: continue` semantics.

### Avoid for Spot: one fat EC2 running local `batch_runner` for hours

Works for On-Demand or a reserved “campaign node,” but a reclaim mid-sweep
wastes wall time unless every child sim checkpoints and the runner resumes from
S3. Use only as a v0 smoke path.

### Bridge from existing manifests

Keep authoring sweeps in `experiments/.../batch_*.json`. Add a small exporter
(future) that expands the Cartesian product to:

- `jobs/<job_id>/input.json` objects in S3, and
- a Batch `arrayProperties.size`,

without teaching AWS about GutIBM’s JSON schema.

## Spot interruption strategy

1. Enable frequent HDF5 checkpoints in campaign configs (e.g. genome/summary
   intervals already exist; add an explicit checkpoint schedule / copy to
   `s3://.../checkpoints/<job_id>/step_XXXXXX.h5` every N steps).
2. On start, entry script checks S3 for latest checkpoint; if present, set
   `checkpoint_file` + `checkpoint_step` and continue.
3. Batch job definition: `retryStrategy` with Spot reclaim as retryable
   (AWS Batch supports retry on `ExitCode` / status reason patterns — tune once
   the wrapper’s exit codes are fixed).
4. Gzip final outputs (`GUTIBM_GZIP_HDF5` / existing helper) before S3 upload to
   cut egress and storage.

Until checkpoint-to-S3 is automated, use shorter Stage 3 proxies or On-Demand
for the first end-to-end validation.

## Container sketch

Base: NVIDIA CUDA devel/runtime matching driver on Batch AMI (Batch GPU AMIs
provide the driver; image needs matching user-mode libs).

Rough contents (to be added under `deploy/aws/` in an implementation PR):

```dockerfile
# Conceptual — not production-pinned yet
FROM nvidia/cuda:12.4.1-devel-ubuntu22.04
RUN apt-get update && apt-get install -y \
    cmake g++ openmpi-bin libopenmpi-dev libhdf5-mpi-dev python3 python3-pip
COPY . /src
WORKDIR /src
RUN CC=gcc CXX=g++ cmake -B build \
      -DGUTIBM_USE_MPI=ON -DGUTIBM_USE_HDF5=ON -DGUTIBM_USE_CUDA=ON \
      -DCMAKE_CUDA_ARCHITECTURES="75;86;89" \
      -DCMAKE_BUILD_TYPE=Release \
 && cmake --build build -j$(nproc) --target gut_ibm \
 && pip3 install ./python
ENTRYPOINT ["/src/deploy/aws/entry.sh"]
```

Runtime image can multi-stage copy `build/gut_ibm` + MPI/HDF5 runtime libs onto
`nvidia/cuda:*-runtime-*` to shrink pull time.

Validate inside the image before Batch:

```bash
nvidia-smi
./build/gut_ibm experiments/smoke_single.json   # with gpu_enabled if GPU present
ctest -R 'gpu_smoke|hdf5_checkpoint' --output-on-failure   # in devel build
```

## IAM and networking (minimal)

- Job role: `s3:GetObject` / `PutObject` / `ListBucket` on
  `gutibm-*/inputs/*`, `.../outputs/*`, `.../checkpoints/*` only.
- No inbound SSH required for v1 (use Batch logs → CloudWatch).
- Egress: pull ECR, talk to S3 (gateway endpoint in VPC recommended to avoid
  NAT cost on large HDF5).
- Region: pick one with deep GPU Spot (often `us-east-1` / `us-west-2`); pin it
  in docs once chosen so ECR/S3/Batch co-locate.

## Cost control knobs

| Knob | Effect |
|------|--------|
| Spot + Batch low-priority queue | Primary savings |
| Array job granularity | Failed index ≠ failed campaign |
| `grid_dx` / domain size | Dominant chem memory & runtime |
| HDF5 schedule | Agents/grid dumps dominate I/O and S3 |
| `profile_steps` off in production | Less log noise |
| Auto-terminate / no min vCPU | Scale compute env to 0 when idle |
| Gzip HDF5 | Storage + transfer |

Budget sanity check: price a **single** Stage 3 seed on `g5.2xlarge` Spot for
the measured wall time before launching the 12-run Kd sweep.

## Phased delivery

### Phase 0 — Decide (this doc)

- [x] Architecture: Batch + Spot GPU + ECR + S3
- [ ] Pin region + default instance type
- [ ] Confirm first campaign (e.g. `batch_kd_sweep.json` vs baseline × 3)

### Phase 1 — Container that runs GPU smoke on Batch

- [ ] `deploy/aws/Dockerfile` (+ `.dockerignore`)
- [ ] `deploy/aws/entry.sh` (download config → run → upload)
- [ ] Build/push to ECR; submit one On-Demand GPU job with `experiments/smoke`
  or Stage 1 config
- [ ] CloudWatch logs show `gpu_enabled` path (no CUDA-absent fallback)

### Phase 2 — S3-backed single Stage 3 job

- [ ] Upload one Stage 3 `input.json`; run to completion or timed smoke
- [ ] Upload `output.h5.gz` to S3; download + `gut_ibm_tools` read smoke
- [ ] Document wall time / Spot price sample in this file’s “Measured” section

### Phase 3 — Array jobs from batch manifests

- [ ] Script: expand `batch_*.json` → S3 job inputs + `aws batch submit-job`
- [ ] Map `AWS_BATCH_JOB_ARRAY_INDEX` → job id
- [ ] Parity with local `batch_runner --dry-run` job list

### Phase 4 — Spot resilience

- [ ] Periodic checkpoint copy to S3
- [ ] Entry script auto-resume from latest checkpoint
- [ ] Batch retry policy for Spot reclaim
- [ ] Flip compute environment to Spot-only (or Spot with On-Demand fallback %)

### Phase 5 — Optional later

- [ ] Multi-rank MPI when chemistry is domain-decomposed
- [ ] CUDA-aware MPI AMI / custom Batch AMI
- [ ] ParallelCluster only if multi-node becomes necessary
- [ ] Cost dashboards + budget alarms

## Open decisions

1. **Region** — where GPU Spot is deepest for your account.
2. **Default instance** — `g5.2xlarge` vs `g4dn.2xlarge` (price vs VRAM headroom).
3. **Who builds the image** — GitHub Actions → ECR on tag, vs laptop `docker build`.
4. **Checkpoint frequency vs S3 PUT cost** for 7-day runs.
5. **Whether v1 requires Terraform** or click-ops + checked-in JSON job definitions.

## Relation to existing docs

| Doc | Role |
|-----|------|
| `experiments/diversity_campaign/README.md` | Memory floors, Stage 3 batches |
| `docs/BATCH_RUNNER.md` | Local sweep / resume semantics to mirror |
| `docs/SCALING.md` | Agent/grid scaling; MPI notes |
| `docs/WSL2_SETUP.md` | Local GPU; not used on Batch hosts |
| `src/gpu/README.md` | What actually runs on device |
| `AGENTS.md` | Build flags (`GUTIBM_USE_CUDA`, `CC=gcc CXX=g++`) |

## Measured results

_Fill in after Phase 2:_

| Config | Instance | Spot? | Wall time | $/run | Notes |
|--------|----------|-------|-----------|-------|-------|
| | | | | | |
