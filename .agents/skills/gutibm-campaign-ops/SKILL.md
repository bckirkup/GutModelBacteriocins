---
name: gutibm-campaign-ops
description: Operate and audit GutIBM AWS campaigns safely. Covers read-only monitoring, checkpoint durability, Spot classification, capacity, throughput, storage, and fork readiness.
---

# GutIBM Campaign Operations

Use this skill for AWS Batch campaign monitoring, checkpoint inspection,
capacity diagnosis, interruption handling, and fork readiness. The default
operation is **read-only**. A campaign submission, cancellation, termination,
queue repoint, or On-Demand fallback always requires explicit user
authorization: GPU Spot runs are expensive and slow, and full runs are not the
routine unit of experimentation.

Read [`docs/AWS_BATCH.md`](../../../docs/AWS_BATCH.md) for deployment,
checkpoint-upload, resume, and Spot mechanics. Read the scripts under
[`deploy/aws/`](../../../deploy/aws/) before interpreting their behavior.
Do not duplicate their commands here; this skill supplies the operational
judgment around them.

## Operating modes and authorization

### Read-only by default

These actions need no authorization:

- inspect Batch jobs, queues, compute environments, job definitions, and
  `statusReason`;
- inspect `status.json`, CloudWatch output, S3 checkpoint objects, and
  `latest.json`;
- compare object size, SHA-256, timestamps, and checkpoint step;
- run an existing script's `--dry-run` if it is documented as mutation-free;
- calculate runtime, cadence, storage, and retry implications.

These actions require explicit user authorization immediately before execution:

- submit or resubmit a job or array;
- cancel or terminate a job;
- repoint or replace a queue or compute environment;
- widen capacity or change instance types/subnets;
- fall back to On-Demand capacity.

When authorization is absent, report the evidence and recommended next action;
do not perform the action.

## Durable progress: prove it from S3

Never claim that a checkpoint exists from a log line. A message saying that a
checkpoint was flushed has been wrong in this repository.

Count a checkpoint as durable only when both are present and consistent:

1. the immutable S3 `step_NNNNNN.h5` object; and
2. the corresponding `latest.json` update containing its URI, size, SHA-256,
   and fidelity metadata.

Use the closed restart, not the live analysis `output.h5`, as the Spot and fork
primitive. Follow the upload, validation, quarantine, and final-interruption
ordering in
[`AWS_BATCH.md#spot-interruption-strategy`](../../../docs/AWS_BATCH.md#spot-interruption-strategy).
Do not infer durability from stdout, a local file, or a `status.json` field
alone.

For a fork, record the origin step, object SHA-256, new seed, overlay revision,
and output prefix. Tier-2 restart restores population and chemistry, but RNG is
reseeded: a fork is a population-state branch, not a bit-identical continuation.

## Classify a stopped job before reacting

Inspect the Batch `statusReason` and attempt history before submitting anything:

| Evidence | Classification | Response |
|---|---|---|
| `statusReason` begins `Host EC2*` | Spot reclaim | The configured retry/resume path should use `latest.json`; inspect the newest durable checkpoint before authorizing any intervention. |
| Attempt timeout | Timeout | It does not retry under the Spot-reclaim rule; do not label it a reclaim or resubmit automatically. Job definitions register a bounded attempt budget (`CAMPAIGN_JOB_TIMEOUT_SECONDS`, default 24 h; practice 3 h), so a timeout means the run exceeded that budget — resume from the newest durable checkpoint rather than restarting the horizon. |
| Application/configuration failure | Application failure | Neither automatic Spot retry nor timeout reasoning applies; inspect logs and artifact state first. |

Do not respond symmetrically to these cases. A mistaken resubmission wastes
capacity; a missed Spot resume loses a run.

## Checkpoint cadence from measured throughput

Cadence must be shorter than the observed host lifetime in wall time. The
observed host lifetime was **37 minutes = 2,220 seconds**. A previous
**2.5-hour** checkpoint cadence therefore produced **zero durable progress**,
regardless of how many retries were configured.

Calculate the wall time between checkpoints:

```text
checkpoint_wall_seconds = interval_steps × measured_wall_seconds_per_step
```

Require that value to be comfortably below 2,220 seconds, with margin for
startup, upload, and scheduling. Do not choose a fixed interval without
checking the measured rate for the actual instance pool. At the measured
range, 60 steps costs approximately 1,008–1,518 seconds (16.8–25.3 minutes);
120 steps can cost 2,016–3,036 seconds and is not safe at the slow end.

## Throughput and storage are ranges

Measured single-GPU `g5`/`g4dn` throughput is **16.8–25.3 wall-seconds per
biological step**, a 1.5× spread. The spread came from widening the instance
pool for availability. Treat it as a planning range, not a constant.

Compression is population-dependent:

```text
step 360: 22:1
step 720: 6.09:1
step 1080: 4.89:1
```

An early lag-phase measurement extrapolates to a wildly wrong campaign figure:
that mistake was made in the dangerous direction by treating a cheap early
checkpoint as representative of the later, larger population. Plan wall time
and storage from the full measured ranges and retain enough capacity for
hundreds of GB in a seven-day all-checkpoint campaign.

## Heartbeats and apparently quiet jobs

Heartbeat quiet is not equivalent to hung. Console progress is emitted on a
simulated-time schedule, so a healthy run can show no new progress line for
tens of minutes.

Check the wall-clock heartbeat (`status.json` and its update time), Batch state,
and CloudWatch activity. Compare elapsed wall time against the expected
checkpoint/upload cadence and the measured seconds-per-step. Do not terminate
or resubmit because the last console progress line is old.

## Operator scripts

Scripts in `deploy/aws/` are run by the user from CloudShell, not by this
session. Before recommending a script:

- resolve paths relative to the script location; never assume the caller is
  inside a clone or a particular current directory;
- require a `--dry-run` path that mutates nothing;
- validate payloads against the captured fixture
  `tests/fixtures/campaign_compute_environment.json`, not an imagined account;
- test every guard against captured real values before shipping it to an
  operator.

One guard in this repository could only ever fail because it matched
`/instance-profile/`, while IAM ARNs use `:instance-profile/`. A guard is not
validated merely because its failure branch executes; prove it accepts the
captured valid value and rejects the invalid value.

Use the user's CloudShell/AWS profile workflow documented in
[`AWS_BATCH.md`](../../../docs/AWS_BATCH.md). This skill does not authorize
cloud mutations.

## Batch capacity facts

- A legacy `AWSBatchServiceRole` compute environment cannot be widened in
  place. Replace it and repoint the queue.
- Queued jobs belong to the queue, so repointing the queue causes them to use
  the replacement compute environment without resubmission.
- Single-AZ GPU Spot can starve for hours.
- Widening to six subnets and eight single-GPU instance types fixed the
  availability problem, at the cost of the measured 1.5× throughput spread.

Treat these as facts to check against the current account, not as permission
to mutate it. Capacity, quota, subnet, IAM, image, and job-definition
validation remains an authorized operator action.

## Security and evidence handling

- Never print, echo, commit, or paste credentials.
- Never print role values, external IDs, access keys, or session tokens.
- Never embed a token in a git remote.
- Prefer redacted command output containing job IDs, states, timestamps, object
  URIs, sizes, and hashes only.
- Do not place secrets in status artifacts or fork records.

## Current campaign state

Do not duplicate the measured calibration inventory or fork choices here:

- [`docs/BRANCHING_FROM_CHECKPOINTS.md`](../../../docs/BRANCHING_FROM_CHECKPOINTS.md)
  contains the checkpoint inventory and fork rules.
- [`docs/CARBON_BUDGET.md`](../../../docs/CARBON_BUDGET.md) contains the nutrient
  and epithelial-boundary audit.
- [`docs/PRE_SUBMISSION_CHECKLIST.md`](../../../docs/PRE_SUBMISSION_CHECKLIST.md)
  from PR #243 lists the real MPI, GPU, AWS, Spot, and campaign-duration checks
  that CI cannot establish.

Before authorizing a campaign action, combine the durable-artifact evidence
above with the fork record and the external checks in the pre-submission
checklist.
