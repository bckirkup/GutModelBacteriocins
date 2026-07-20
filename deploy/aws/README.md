# AWS Batch practice path (novice / paste-safe)

Planning detail: [`docs/AWS_BATCH.md`](../../docs/AWS_BATCH.md).

**Region:** `us-east-1` only.  
**First goal:** one tiny GPU smoke job → `SUCCEEDED` + file in S3.

## Why your multi-line pastes break

Long commands with `\` line continuations fail easily when:

- a trailing space sneaks in after `\`
- CloudShell / Word / Slack wraps or “smart-quotes” the text
- nested `"..."` JSON inside a shell command gets mangled

**Fix:** do not paste giant multi-line CLI blobs. Run the numbered scripts below (one short command each). JSON lives in files under `deploy/aws/policies/` or temp files the scripts write.

## Split of work

| Step | Where | Command |
|------|--------|---------|
| Build + push CUDA image | **Laptop** (needs Docker) | `bash deploy/aws/01_push_image.sh` |
| Create S3 / IAM / Batch | Laptop **or CloudShell** | `bash deploy/aws/02_setup_practice_stack.sh` |
| Submit smoke job | Laptop **or CloudShell** | `bash deploy/aws/03_submit_smoke.sh` |
| Watch job | Laptop **or CloudShell** | `bash deploy/aws/04_watch_job.sh <jobId>` |

CloudShell is great for AWS CLI, but **cannot build** the CUDA Docker image. Do step 01 on a machine with Docker.

## One-time laptop prep

1. Install Docker Desktop and start it.
2. Install [AWS CLI v2](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).
3. Run this **one line**, then answer the prompts (region = `us-east-1`):

```bash
aws configure
```

4. Check identity with this **one line**:

```bash
aws sts get-caller-identity
```

5. Clone this repo and `cd` into it.

## Step-by-step (copy one line at a time)

From the **repo root**:

```bash
bash deploy/aws/01_push_image.sh
```

```bash
bash deploy/aws/02_setup_practice_stack.sh
```

```bash
bash deploy/aws/03_submit_smoke.sh
```

Copy the printed `jobId`, then:

```bash
bash deploy/aws/04_watch_job.sh PASTE_JOB_ID_HERE
```

When it says `SUCCEEDED`, download with this **one line** (ACCOUNT is your 12-digit account id):

```bash
aws s3 cp s3://gutibm-outputs-ACCOUNT/practice/smoke_gpu/output.h5.gz ./output.h5.gz
```

Or, after sourcing env:

```bash
source deploy/aws/env.sh
```

```bash
aws s3 cp "s3://${OUTPUT_BUCKET}/practice/smoke_gpu/output.h5.gz" ./output.h5.gz
```

## CloudShell-only for steps 02–04

If the image is already in ECR from your laptop:

1. Open CloudShell (console terminal icon), region `us-east-1`.
2. Clone the repo **or** upload the `deploy/aws/` folder plus `experiments/smoke_gpu.json`.
3. Run the same `02` / `03` / `04` one-liners above.

## Pass checklist

- Batch job status = `SUCCEEDED`
- CloudWatch logs show the run (GPU path, not a silent CPU fallback)
- `output.h5.gz` exists under `s3://gutibm-outputs-<account>/practice/...`

## Scripts in this folder

| File | Role |
|------|------|
| `env.sh` | Shared account/region/bucket names (sourced by others) |
| `01_push_image.sh` | Local Docker build → ECR |
| `02_setup_practice_stack.sh` | S3 + IAM + Batch CE/queue/job def (On-Demand `g4dn.xlarge`) |
| `03_submit_smoke.sh` | Upload `smoke_gpu.json` + submit one job |
| `04_watch_job.sh` | Poll until SUCCEEDED/FAILED |
| `05_setup_campaign_stack.sh` | Spot GPU CE + queue + job def for Stage 3 (`g5.2xlarge`, one GPU/run) |
| `entry.sh` | Container entrypoint (S3 → `gut_ibm` → S3; checkpoint sync) |
| `Dockerfile` | CUDA + MPI + HDF5 image |
| `submit_array_example.sh` | Later: array jobs (after smoke works) |
| `policies/*.json` | IAM trust documents (no paste required) |

## Use your own S3 bucket (optional)

`env.sh` derives `gutibm-inputs-<account>` / `gutibm-outputs-<account>` by
default, but you can point at an **existing bucket you own** — it is never
recreated (the `head-bucket` guard in `02` skips `aws s3 mb`), and the job role
policy is scoped to whatever these resolve to. Export **before** running the
scripts:

| Env var | Effect |
|---------|--------|
| `BUCKET` | One shared bucket for both inputs and outputs (use key prefixes, e.g. `s3://my-bucket/gutibm/jobs`) |
| `INPUT_BUCKET` / `OUTPUT_BUCKET` | Set input/output buckets individually |

```bash
export BUCKET=my-existing-bucket
bash deploy/aws/02_setup_practice_stack.sh
```

## Campaign (Stage 3) commands

After the smoke path works, run the campaign stack (details + full export
commands in [`docs/AWS_BATCH.md`](../../docs/AWS_BATCH.md) “Phase 5”). Region
`us-east-1`, one GPU instance per run, Spot + checkpoint resilience.

| Step | Command |
|------|---------|
| Roles + bucket (idempotent) | `bash deploy/aws/02_setup_practice_stack.sh` |
| Campaign Spot GPU stack | `bash deploy/aws/05_setup_campaign_stack.sh` |
| Export dry-run (parity check) | `python -m gut_ibm_tools.aws_batch_export experiments/diversity_campaign/stage3_campaign/batch_kd_sweep.json --input-prefix … --output-prefix … --job-queue "$JOB_QUEUE_CAMPAIGN" --job-definition "$JOB_DEFINITION_CAMPAIGN" --dry-run` |
| Export + submit array | same command without `--dry-run` (add `--checkpoint-prefix …` for Spot resume) |
| Watch | `bash deploy/aws/04_watch_job.sh <jobId>` |

Overridable campaign env vars: `CAMPAIGN_MAX_VCPUS` (default 96),
`CAMPAIGN_ONDEMAND_FALLBACK=1` (add On-Demand fallback CE),
`COMPUTE_ENV_CAMPAIGN` / `JOB_QUEUE_CAMPAIGN` / `JOB_DEFINITION_CAMPAIGN`.

## After smoke works

See `docs/AWS_BATCH.md` for Spot flip, array jobs, and Stage 3 on `g5.2xlarge`. Do not start there.
