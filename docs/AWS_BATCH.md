# AWS Batch Deployment Plan (CUDA, Spot / low priority)

Plan for running production GutIBM jobs on AWS when desktop/WSL cannot finish
full campaigns (especially Stage 3: 7-day biology, 2 mm domain, GPU chemistry).

**Status:** Phase 1 CUDA smoke green on Batch (Jul 2026); campaign path next.
**Region (pinned):** `us-east-1` (co-locate ECR, S3, Batch).
**Primary workload (later):** `experiments/diversity_campaign/stage3_campaign/`.
**Practice workload (now):** `experiments/smoke_gpu.json` (+ optional
`smoke_gpu_batch.json`).

For the pre-campaign throughput and resume calibration, see
[AWS_CALIBRATION_6H.md](AWS_CALIBRATION_6H.md).

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
| Campaign instances | **`g5.2xlarge`**, `g5.4xlarge`, `g5.8xlarge`, `g5.16xlarge`, `g4dn.2xlarge`, `g4dn.4xlarge`, `g4dn.8xlarge`, `g4dn.16xlarge` | Single-GPU sizes with at least 8 vCPUs and 32 GiB host RAM; multi-GPU sizes are excluded because one requested GPU would leave paid accelerators idle |
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
| Campaign CE | Multi-AZ Spot with single-GPU campaign sizes | Stage 3 VRAM/RAM; the CE should search every available subnet in its VPC |
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
`GPU=1` + enough memory. Keep the CE in multiple AZs: a single-AZ GPU Spot
search is the narrowest possible search for the scarcest capacity class and
can present as a job queued forever with a healthy CE and no error anywhere.

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
   - CloudWatch / run log shows GPU init (not silent CPU fallback); smoke
     submits with `REQUIRE_GPU=1` so missing `GPU: ON` fails the job
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

Closed **restart artifacts** (not the live analysis `output.h5`) are the Spot/fork
primitive:

1. C++ writes `restart/step_NNNNNN.h5` on `restart.interval_steps` (Tier 2: agents +
   lineage + genome + **grid** + clock), flush+close, atomic rename.
2. `entry.sh` scans that directory every `CHECKPOINT_INTERVAL_SECONDS` (default 300),
   uploads each file **immutably** to `CHECKPOINT_S3_PREFIX/step_*.h5`, and updates
   `latest.json` (uri, size, sha256, fidelity). Never overwrites an existing step object.
3. Resume reads `latest.json` (or `CHECKPOINT_S3_URI`) → validates with `h5ls`
   (`/agents/` + `/grid/`) → injects `checkpoint_file`. Corrupt objects are quarantined
   (`.corrupt.<ts>`) and the job exits 3 — it does **not** clobber S3 with an empty stub.
4. IMDSv2 Spot notice → SIGTERM → final closed restart upload → non-zero exit so Batch retries.
5. Enable **S3 versioning** on the outputs bucket as a seatbelt.
6. Campaign job def `retryStrategy` retries `Host EC2*` reclaim reasons.

**Science contract:** Tier 2 restores population + chemistry. RNG is **reseeded** from
the job `seed` (not bit-identical continuation). Document forks as population-state
branches unless Tier 3 RNG serialize lands later.

**Cadence (Stage 3):** `restart.interval_steps: 60` (~1 h sim at `bio_dt=60`) is a
good default; entry sync every 300 s uploads whatever closed files appeared.

## Midstream burn-in and forks

Shared colonization is expensive (~$70/seed for a full week). Prefer:

1. Run [`3a_burnin.json`](../experiments/diversity_campaign/stage3_campaign/3a_burnin.json)
   (2 simulated days) with checkpoints → `s3://…/campaign/burnin_seed1001/ckpt/`.
2. Fork Kd / motility / ensemble seeds from `latest.json`:

```bash
gut-ibm-aws-branch \
  --from "s3://${OUTPUT_BUCKET}/campaign/burnin_seed1001/ckpt/latest.json" \
  --overlay experiments/diversity_campaign/stage3_campaign/3_kd_sweep_1e-7.json \
  --seed 2001 --total-time 604800 \
  --out-prefix "s3://${OUTPUT_BUCKET}/campaign/fork_kd1e-7_seed2001" \
  --job-queue gutibm-gpu-campaign \
  --job-definition gutibm-cuda-campaign
```

Keep a few full-from-scratch seeds as statistical controls so shared-burn-in bias
stays visible.

## Phase Observability — progress, usefulness, cost

**Run AWS commands from a local clone with an AWS profile.** This cloud agent
environment does not manage live Batch jobs.

For AWS CloudShell, clone the repository rather than uploading only one script.
The capacity-widening script also requires `jq`:

```bash
git clone https://github.com/bckirkup/GutModelBacteriocins.git
cd GutModelBacteriocins
jq --version
bash deploy/aws/07_widen_campaign_capacity.sh
```

If `env.sh` is not beside a numbered deployment script, it will stop and ask
you to clone the repository instead of attempting to source a misleading path.

### What you get while a job is RUNNING

| Signal | Where |
|--------|--------|
| `pct=` / `rate=` / `eta_s=` | gut_ibm stdout → CloudWatch (each `output_interval`) |
| `status.json` heartbeat | S3 beside `latest.json` / restart prefix (or next to output) |
| Batch state + reclaim reason | `describe-jobs` / `gut-ibm-aws-status` |
| Usefulness triage hints | Warnings from `gut-ibm-aws-status` (not auto-cancel) |

Progress line shape (parsed by `entry.sh`):

```text
Step 10  t=600s  dt=60s  global_agents=50  local_agents=50  mu_avg=5e-4  pct=10  rate=1.2  eta_s=4500
```

For long campaigns, consider `output_interval` 600 s so CloudWatch and heartbeats
update more often than the default 3600 s.

### Operator commands

```bash
# One-screen report (Batch + status.json + Spot / usefulness hints)
gut-ibm-aws-status <jobId> \
  --checkpoint-prefix "s3://${OUTPUT_BUCKET}/campaign/kd_sweep/ckpt" \
  --array-index 0

# Post-hoc array QA: download outputs, join seeds, assert fingerprints differ
gut-ibm-aws-qa \
  --output-prefix "s3://${OUTPUT_BUCKET}/practice/smoke_gpu_batch/out" \
  --input-prefix  "s3://${INPUT_BUCKET}/practice/smoke_gpu_batch/jobs"

# Watch loop (uses gut-ibm-aws-status when installed)
CHECKPOINT_S3_PREFIX="s3://${OUTPUT_BUCKET}/campaign/kd_sweep/ckpt/0/" \
  bash deploy/aws/04_watch_job.sh <jobId>

# Rough $/run before submit (checked-in us-east-1 table)
gut-ibm-aws-estimate --instance-type g5.2xlarge --wall-hours 24 --array-size 12
gut-ibm-aws-estimate --list-prices
# Optional interrupt sensitivity:
gut-ibm-aws-estimate --instance-type g5.2xlarge --wall-hours 24 \
  --interrupt-every-hours 6 --checkpoint-overhead-hours 0.25

# Fork from an immutable burn-in restart
gut-ibm-aws-branch \
  --from "s3://${OUTPUT_BUCKET}/campaign/burnin_seed1001/ckpt/latest.json" \
  --overlay experiments/diversity_campaign/stage3_campaign/3_kd_sweep_1e-7.json \
  --seed 2001 --total-time 604800 \
  --out-prefix "s3://${OUTPUT_BUCKET}/campaign/fork_kd1e-7_seed2001" \
  --job-queue gutibm-gpu-campaign \
  --job-definition gutibm-cuda-campaign \
  --dry-run
```

### Spot reclaim vs crash vs biology failure

| Symptom | Likely cause |
|---------|----------------|
| Batch `statusReason` starts with `Host EC2*`; `status.json` has `spot_interruption=true` | Spot reclaim (retry should resume from `latest.json`) |
| Resume attempt exits 134 / HDF5 "bad object header"; S3 object suddenly tiny | Old live-copy path — fixed by closed immutable restarts; quarantine and re-burn-in if needed |
| `status.json` has `memory_pressure=true` / state `memory_pressure` | Graceful stop before OOM; **do not** auto-retry same size — use larger instance / job-def memory, then resume from `latest.json` |
| Job FAILED, stale `status.json` (age ≫ sync interval), no Spot/memory flags | Hang / hard OOM kill without flush |
| `global_agents≤1` or population-stop in logs; job SUCCEEDED early | Biology collapse (usefulness warning) |
| `mu_avg` very low for many heartbeats | Washout-risk triage hint |
| CUDA init missing / silent CPU fallback in logs | Config/image GPU path problem |

### Memory guard (avoid hard OOM kills)

Job definitions already cap container memory (`14000` practice / `28000` campaign MiB)
and Stage 3 needs ~8 GB host + ~8 GB VRAM (`experiments/diversity_campaign/README.md`).
That sizing alone does **not** prevent mid-run growth from being OOM-killed.

`entry.sh` also:

1. Samples `MemAvailable`, cgroup free (Batch limit), and `nvidia-smi` free VRAM into `status.json`.
2. Before start and each sync interval, if effective free RAM &lt; `MEMORY_MIN_AVAILABLE_MB` (default **2048**) or GPU free &lt; `GPU_MIN_FREE_MB` (default **512**), SIGTERM → checkpoint + `memory_pressure` status → exit **without** Batch Spot-style auto-retry (same box would fail again).
3. Disable with `MEMORY_GUARD=0`; tune floors via env on the job definition / submit overrides.

`gut-ibm-aws-status` warns on `memory_pressure`, or soft `low_memory` / `low_gpu_memory` when free drops below 4 GiB / 1 GiB.

Usefulness warnings from `gut-ibm-aws-status` are **triage hints only** (v1 does
not auto-cancel jobs for biology signals).

### Filling the Measured results table

After practice / baseline smokes on a machine with credentials:

1. Note instance type, Spot vs On-Demand, wall time from Batch `startedAt`→`stoppedAt`.
2. Estimate `$` with `gut-ibm-aws-estimate --wall-hours …` (or multiply wall hours × table rate).
3. Paste into the Measured table below.
4. Use those wall hours for campaign cost planning before large arrays.

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
- [x] Campaign instances: multi-AZ CE with single-GPU sizes that fit 8 vCPU / 28 GB / GPU=1
- [x] Practice configs: `experiments/smoke_gpu.json`, `smoke_gpu_batch.json`

### Phase 1 — CUDA smoke on Batch (done Jul 2026)

- [x] Draft `deploy/aws/Dockerfile`, `entry.sh`
- [x] Paste-safe practice scripts (`deploy/aws/01`–`04` + README)
- [x] Create ECR repo + push `gutibm:cuda` from laptop (`01_push_image.sh`)
- [x] Create Batch GPU CE/queue/job definition in `us-east-1` (`02_setup_practice_stack.sh`
      previously created CE/JD/roles; PowerUser completed missing ECR/S3/queue)
- [x] One On-Demand job with `smoke_gpu.json` (`03` + watch; job
      `1fa5bae6-8860-441b-b1bb-9eeb3c8f31af`)
- [x] Confirm GPU path in logs + S3 output (`GPU: ON (device 0)`,
      `REQUIRE_GPU=1`, `output.h5.gz` opens in `gut_ibm_tools`)
- [x] Parser-fix re-smoke after `\uXXXX` support (`04cc2257-33f3-4bac-930c-205a053d9991`;
      no legacy-fallback warning; nested `initial_strains` applied)
- [x] Phase 1b: array of 2 from `smoke_gpu_batch.json`
      (`f12ed608-efea-4e69-954d-5d6a18988ce3`, both children SUCCEEDED)

### Phase 2 — Single Stage 3 seed on `g5.2xlarge`

- [x] Campaign CE/queue/job-def sized for `g5.2xlarge` (one GPU per run) via `05_setup_campaign_stack.sh`
- [ ] One `3a_baseline` / `batch_baseline` seed (submit a single index of `batch_baseline.json`)
- [ ] Record wall time / $ in Measured section below

### Phase 3 — Array export from batch manifests

- [x] Helper: `batch_*.json` → S3 job tree + `submit-job` (`gut_ibm_tools.aws_batch_export` / `gut-ibm-aws-export`)
- [x] Parity with `batch_runner --dry-run` (shared `parse_batch_config` / `build_job_config`; unit-tested for `batch_kd_sweep.json` + `batch_baseline.json`)

### Phase 4 — Spot resilience

- [x] Checkpoint → S3 + auto-resume (C++ closed `restart/step_*.h5`; `entry.sh` immutable upload + `latest.json`)
- [x] Midstream branch CLI (`gut-ibm-aws-branch`)
- [x] Retry policy; Spot-only (or Spot-with-fallback) CE (`05` Spot CE `SPOT_CAPACITY_OPTIMIZED` + `retryStrategy` retrying `Host EC2*`; set `CAMPAIGN_ONDEMAND_FALLBACK=1` for an On-Demand fallback CE)
- [x] IMDSv2 Spot notice → SIGTERM + checkpoint/status flush (`entry.sh`)

### Phase Observability

- [x] Progress lines: `pct` / `rate` / `eta_s` in `simulation.cpp`
- [x] S3 `status.json` heartbeat from `entry.sh`
- [x] `gut-ibm-aws-status` + richer `04_watch_job.sh`
- [x] `gut-ibm-aws-estimate` checked-in us-east-1 price table
- [x] Usefulness triage warnings (population / low μ / stale heartbeat / Spot)

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

   Creates `gutibm-gpu-campaign-spot` (Spot single-GPU sizes across every
   available subnet in the selected VPC, `SPOT_CAPACITY_OPTIMIZED`,
   `maxvCpus=${CAMPAIGN_MAX_VCPUS:-96}`), queue
   `gutibm-gpu-campaign`, and job def `gutibm-cuda-campaign` (8 vCPU / 28 GB /
   `GPU=1`, `retryStrategy` for Spot reclaims).

3. **Widen an existing campaign CE (admin credentials only):**

   ```bash
   bash deploy/aws/07_widen_campaign_capacity.sh
   ```

   The script derives the VPC from the CE's existing subnet, discovers every
   available subnet in that VPC, and creates or reuses a service-linked-role
   replacement CE with the single-GPU instance-type list. It waits for the new
   CE to become `VALID`/`ENABLED`, repoints the queue to it, and prints the
   resulting configuration. It never changes job definitions. The scoped
   GutIBM role intentionally lacks the EC2 discovery and Batch update
   permissions required by this admin-only operation.

   A legacy campaign CE created with `AWSBatchServiceRole` cannot update its
   subnets or instance types at all. The script therefore creates a distinct
   `gutibm-gpu-campaign-spot-v2` CE without `--service-role`, so Batch uses
   `AWSServiceRoleForBatch`, then repoints the queue to the new CE and drops
   the legacy CE from the queue order. The old CE is left enabled; it is not
   retained as a lower-order fallback because its single-AZ capacity cannot be
   repaired. The existing queued job is associated with the job queue, not a
   CE; [AWS Batch schedules jobs using the queue's current
   compute-environment order](https://docs.aws.amazon.com/batch/latest/userguide/compute_environments.html),
   so repointing the queue does not require resubmitting it.

   The practice stack's `02_setup_practice_stack.{sh,ps1}` still creates its
   separate practice CE with the legacy role. That CE is intentionally outside
   this campaign-capacity replacement; do not use it for Stage 3 campaign jobs.

4. **Dry-run the array export** (no uploads, no submit — verify count + inputs):

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

1. Whether campaign CE is Spot-only or 80/20 Spot/On-Demand after smoke.
2. Bucket naming (`gutibm-inputs-<account>` vs a shared research bucket).
3. CloudWatch custom metrics / dashboards (heartbeat JSON is enough for v1).
4. Auto-cancel on usefulness warnings (triage-only today).
5. **PowerUser IAM gap:** SSO `PowerUserAccess` cannot `iam:PutRolePolicy`. For the
   Phase 1 smoke, S3 access was granted with **bucket policies** on
   `gutibm-inputs/outputs-<account>` allowing `gutibm-batch-job-role`. An admin
   should eventually run `02_setup_practice_stack.sh` (or put the matching
   inline role policy) so identity-based access matches the scripts.
6. ~~**JSON parse warning on Batch**~~ — fixed in `ConfigJson` (`\uXXXX` /
   surrogate pairs). Re-smoke + Phase 1b array after image push showed no
   legacy-fallback warning and nested keys applied.

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

| Config | Instance | Spot? | Wall time | $/run | Notes |
|--------|----------|-------|-----------|-------|-------|
| `smoke_gpu` (first) | `g4dn.xlarge` | No (OD) | ~3.1 s container; ~4 min queue→done | ~$0.05 (0.1 h OD table) | `1fa5bae6-…`; GPU ON; had JSON legacy fallback before parser fix |
| `smoke_gpu` (parser fix) | `g4dn.xlarge` | No (OD) | ~3.3 s container | ~$0.05 | `04cc2257-…`; no JSON warning; 20 agents from nested strains; `REQUIRE_GPU=1` |
| `smoke_gpu_batch` ×2 | `g4dn.xlarge` | No (OD) | ~4–5 s / child | ~$0.05 | Array `f12ed608-…`; both SUCCEEDED; seeds 4092/4093 → final agents 13 vs 12; `gut-ibm-aws-qa` fingerprints differ |
| `3a_baseline` seed | `g5.2xlarge` | | | | Phase 2 |
