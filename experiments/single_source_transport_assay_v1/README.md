# Single-source transport assay v1

This assay isolates the mucin-charge dependence of colicin transport. It exists
because the receptor-selection Stage C transport analysis measured radii from
all domain toxin mass to the *nearest of many* active sources, which measures
source packing as much as diffusion: sparse-source snapshots own larger Voronoi
cells and report larger radii at identical diffusivity, and in the committed
Stage C results the active-source count explains most of the radius variation
(Spearman ≈ −0.96 against r50, ≈ −0.90 against r90) while source counts differed
across amplitude arms by up to eightfold at matched times.

Here the source geometry is controlled instead of corrected for: one release
event, at one position, identical across the amplitude arms of a seed.

## What is controlled

| Control | Mechanism |
|---|---|
| exactly one possible source | one cell carries ColE1; nothing else can release colicin |
| no second producer | every strain has `mu_max = 0`, so no division |
| bystanders never release | no bacteriocin locus ⇒ never SOS-induced (`FixBacteriocin::compute` skips agents without BI loci) |
| run does not end at lysis | three bystanders keep the global agent count above the population-stop threshold of 1 |
| no receptor sink | all `BtuB` expression is 0, so the profile is transport, not consumption |
| identical source position across amplitudes | initial placement is drawn before any amplitude-dependent state is consumed, so one seed gives one position; the analyzer asserts it to 1e-12 m rather than assuming it |
| radial support inside the domain | placement is confined to a 2 µm band at mid-depth, so every draw supports the predeclared 40 µm support |
| exact source time | the release window decays on a 300 s timescale, so the analyzer requires the provenance `event_time_s` of the lysis itself and refuses outputs that only carry provenance-export steps |

## Predeclared analysis

`assay_contract.json` is the contract, committed before any run exists, and
`analyze_assay.py` reads its definitions rather than restating them:

- **shell averaging** — fixed 2 µm spherical shells around the source; a shell
  value is the unweighted arithmetic mean of the concentrations of the voxel
  centres with radius in `[r_lo, r_hi)` (measure `dr`).
- **radial support** — 40 µm, and a run whose in-domain support is smaller is a
  blocker, never a silently trimmed profile.
- **boundary handling** — minimum image in periodic x and y, plain distance in
  non-periodic z, no image sources.
- **profile r50/r90** — shell means become a radial mass density `mean(r)·r²`
  accumulated over the support; `rXX` is where the normalized cumulative mass
  first reaches `XX/100`, linearly interpolated between shell centres.
- **analysis window** — `event_time_s < t ≤ event_time_s + 1500 s`; the snapshot
  coincident with the lysis step is excluded because the grid can be written
  before the burst's first chemistry application.
- **paired times** — amplitude arms are compared only at snapshot times all
  three arms observed; a missing paired time blocks the assay instead of
  changing each arm's support.

## Gate

`single_source_transport_gate` passes only if all 12 runs returned, every run
authenticates against the execution source SHA, `device_delivery`, one rank, its
seed and amplitude, `kd_corrinoid_btuB = 1e-4`, `b12_initial_conc = 1e-3`,
`burst_release_tau = 300 s`, the ColE1 pI/diffusion/burst-size identity and the
grid/provenance schedules; every producer run recorded exactly one strain-1
lysis with exact timing; the source position is shared across amplitudes; paired
times are complete; and at every paired time in every seed `r50` decreases
across amplitude 0 → 15 → 60 by at least 0.5 µm and `r90` by at least 0.1 µm,
while the 0–10 µm near-field mean increases.

The margins come from a single-seed pilot of this exact configuration (seed
20260913, 600 s horizon, run on a serial host build):

| amplitude | r50 (µm) | r90 (µm) | peak `bacteriocin_BtuB` |
|---|---|---|---|
| 0 | 28.72 | 37.25 | 1.15e-06 |
| 15 | 27.22 | 36.95 | 1.26e-05 |
| 60 | 25.44 | 36.51 | 4.65e-05 |

The pilot also shows what the assay is designed to give: one source at an
identical position and identical exact lysis time in all three arms, radii that
are stationary in time (the QSSA profile shape is set by amplitude, not by
snapshot age), and radii far enough inside the 40 µm support to be unclipped.
`r90` separates by only ~0.3 µm because the outer profile is dominated by the
shared support, hence its narrower predeclared margin.

`analyze_assay.py` exits 0 on pass, 2 when blocked (missing outputs, failed
authentication, missing paired times), 1 when the assay ran and the predeclared
ordering failed. The assay result validates the intrinsic transport law; it does **not** by itself calibrate or select a mucin-charge amplitude. A pass must be combined with the already-passed ecological Stage C population result (amplitudes 0/15/60) to support **amplitude 15 as a conditional D choice**. Record `C_transport_gate=true` only for the assay's exact execution SHA and image digest, then regenerate D with amplitude 15 on those same identities. Until that recorded pass, D remains blocked. Any code/image revision requires a new assay.

## Running it

Planning generation (no real SHA or digest, never submits anything):

```bash
cd experiments/single_source_transport_assay_v1
python3 prepare_assay.py --clean
```

Deployment generation requires the exact commit that is built:

```bash
EXECUTION_SOURCE_SHA=$(git rev-parse HEAD)
python3 prepare_assay.py --clean --deployment \
  --execution-source-sha "$EXECUTION_SOURCE_SHA" \
  --image-digest "<repo>@sha256:<64 hex>"
```

Each of the 12 jobs is one GPU, one MPI rank, a 100 µm cube on a 2 µm grid, 1 h
simulated time, with `bacteriocin_BtuB` grids every 120 s and provenance every
step. Return files as `generated/results/<array-index>/output.h5.gz`
(uncompressed `.h5` accepted), then:

```bash
python3 analyze_assay.py
```

Outputs land in `analysis/`: `assay_snapshots.csv` (every analyzed snapshot),
`assay_metrics.json` (authentication, source, profiles, ordering) and
`assay_gate.json` (the gate and its blockers).

## Requirement on the execution image

The provenance layer must carry `event_time_s`/`event_step`. Outputs produced
before that field existed are refused rather than reinterpreted, because the
provenance-export step is up to one export interval later than the event, and
under the 300 s release decay that error can inflate a source's release weight
several-fold.
