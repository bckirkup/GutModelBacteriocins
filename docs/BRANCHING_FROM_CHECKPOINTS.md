# Branching from calibration checkpoints

This is the operator runbook for turning a calibration checkpoint into a
short fork or a longer campaign. It records the measured limits of the
calibration state so that a fork answers a defined question instead of silently
continuing an unchecked growth regime.

The AWS submission and resume mechanics are canonical in
[`AWS_BATCH.md`](AWS_BATCH.md), especially
[Midstream burn-in and forks](AWS_BATCH.md#midstream-burn-in-and-forks) and
[Spot interruption strategy](AWS_BATCH.md#spot-interruption-strategy). Follow
those sections for the command, checkpoint validation, immutable upload,
`latest.json`, retry, and interruption ordering. This document supplies the
choice of origin and the scientific records to make before submission.

## 1. Choose an origin explicitly

The measured calibration inventory is:

```text
s3://gutibm-outputs-994254241749/campaign/calibration_6h/calibration_6h_seed1001/ckpt/
```

| Object | Population | Use |
|---|---:|---|
| `step_000360.h5` | 654 agents | Valid origin; lag phase |
| `step_000720.h5` | 5,377 agents | Valid origin |
| `step_001080.h5` | 42,937 agents | Do not use as an origin |

Fork from **step 360 or step 720**. Step 1080 and later are deliberately
excluded: the population is in unchecked exponential growth, so measurements
from those states describe the runaway rather than established ecology.

`latest.json` currently points at step 1080. Do **not** use it to select an
origin for this calibration; name `step_000360.h5` or `step_000720.h5`
explicitly, and then follow the AWS checkpoint validation procedure.

Before submitting, record:

- the exact origin step;
- the origin object's SHA-256;
- the overlay/configuration and its git revision;
- the fork seed;
- the question and planned observables.

The seed is essential: a fork is a **population-state continuation**, not a
bit-identical continuation. Tier-2 restart restores the population and
chemistry, but RNG is reseeded from the fork job's seed.

## 2. Plan wall time before requesting capacity

Measured single-GPU Spot throughput was **16.8–25.3 wall-seconds per biological
step** on `g5` and `g4dn`. The 1.5× spread came from widening the instance pool
for availability. It is a planning range, not a fixed runtime.

| Horizon | Biological steps | Estimated wall time |
|---|---:|---:|
| 2 simulated days | 2,880 | approximately 13–20 hours |
| 7 simulated days | 10,080 | approximately 47–71 hours |

Use the shorter horizon for L2 questions such as local exclusion, immigrant
washout, or mechanism directionality. Reserve seven-day runs for L3 retention
and sweep questions, and keep scratch controls from a full-from-scratch
population so shared-burn-in bias remains visible.

Storage is also a range, not one campaign-wide ratio. Compression measured:

- 22:1 at step 360;
- 6.09:1 at step 720;
- 4.89:1 at step 1080.

Compression falls as the population grows. A seven-day campaign retaining every
checkpoint can therefore plausibly consume **hundreds of GB**; retention is
part of the operating plan, not optional cleanup.

## 3. Submit and verify the fork

1. Select step 360 or 720 and record its SHA-256.
2. Choose a new seed and record it as a population-state branch.
3. Prepare the smallest overlay that changes the parameter under study.
4. Submit through the fork path documented in
   [`AWS_BATCH.md`](AWS_BATCH.md#midstream-burn-in-and-forks); do not create a
   second submission recipe here.
5. Confirm the output prefix, checkpoint prefix, queue, job definition, and
   requested horizon before launch.
6. After launch, verify the durable checkpoint object and `latest.json`, rather
   than trusting stdout.
7. At completion or interruption, retain the origin step, origin SHA-256, fork
   seed, overlay revision, job/attempt identifiers, and final checkpoint URI.

Read event and nutrient counters from the summary artifacts using their
semantics:

- `events/*` and nutrient-flux `*_interval` fields cover only the current
  reporting window;
- `events/cumulative_*` and nutrient-flux `*_cumulative` fields are the
  continuous totals across closed windows and resumes.

Do not report a sequence of interval values as a running total.

## 4. Durability rules

These rules come from observed failures, not general caution:

- **A log line is not a checkpoint.** A 2.5-hour checkpoint cadence was paired
  with a 37-minute host lifetime and produced zero durable progress. Require the
  S3 checkpoint object **and** the corresponding `latest.json` update.
- **Closed restart artifacts are the durable primitive.** Use the immutable
  restart object, not the live analysis `output.h5`, as the branch origin. The
  upload, validation, and Spot ordering are defined in
  [`AWS_BATCH.md`](AWS_BATCH.md#spot-interruption-strategy).
- **Use cumulative counters for totals.** Nutrient-flux and event counters are
  per-window in their interval fields; the `_cumulative` fields are the values
  that survive resumes.
- **Classify interruption before interpreting retry behavior.** A Spot reclaim
  has a Batch `statusReason` beginning `Host EC2*`, triggers the configured
  retry/resume path, and should leave a closed restart. An attempt timeout is a
  different failure and does not retry under that Spot-reclaim rule.

When a density halt occurs, its closed restart is itself a usable origin. Keep
the halt metadata with the artifact:

- `halt_reason_code`;
- `halt_density_cells_per_mL`.

## 5. Define what the fork measures

State the toxin and the spatial observable before looking at results.

The **300–3,000 producer thresholds** in
[`T2_COLONY_CHALLENGE.md`](T2_COLONY_CHALLENGE.md) are **ColE1/BtuB-specific**
and use a screening length of approximately 45 µm. They do not apply to every
bacteriocin. The live calibration hazard was approximately 98.9% ColB/FepA,
with a screening length of approximately 186 µm; at 10 µm, per-target hazard
correlated **negatively** with local producer count (−0.25). A fork asking
about bacteriocin efficacy should therefore measure the diffuse ColB halo and
must name the toxin, rather than applying the ColE1/BtuB thresholds to ColB.

Prefer direct producer counts within a fixed radius of each susceptible target
over DBSCAN colony counts. DBSCAN returned 88 producers at one radius and 4,625
at twice that radius: the result is a clustering knob, not a stable biological
measurement.

For each fork, record at minimum:

- toxin and susceptible receptor;
- origin step and SHA-256;
- fork seed;
- overlay parameters;
- interval and cumulative event fields;
- spatial radius and target population used for local counts;
- whether the question is about a diffuse halo, a tight neighbourhood, growth,
  washout, retention, or diversity.

## 6. Open items: do not assume these away

- The dysbiosis threshold is **disabled by default** pending Benjamin's
  healthy-density ceiling in cells/mL. The killed run peaked near
  `1.575e8 cells/mL`, approximately three orders of magnitude below lumen
  density. The earlier `1.6e11` figure was a 1,000× arithmetic error and is
  wrong. If a fork enables the halt, retain its reason and tripped density in
  the closed restart.
- The epithelial carbon boundary is an **unmetered Dirichlet reservoir**.
  Growth rates from any fork inherit that subsidy until the boundary question
  is settled. Use [`CARBON_BUDGET.md`](CARBON_BUDGET.md) for the measured
  accounting and the open boundary decision; do not present subsidized growth
  as a closed ecological carbon budget.

## 7. Fork record

Before calling a result interpretable, put this record beside the analysis:

```text
origin_step:
origin_sha256:
fork_seed:
overlay:
git_revision:
toxin_and_receptor:
horizon_steps:
checkpoint_prefix:
question:
primary_observables:
interval_fields_used:
cumulative_fields_used:
spot_or_timeout:
halt_metadata_if_any:
```

If the origin, seed, durability evidence, toxin identity, or spatial
measurement definition is missing, treat the fork as operational output rather
than a reproducible scientific result.
